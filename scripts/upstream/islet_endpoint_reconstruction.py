#!/usr/bin/env python3
"""Bounded, source-only C018 islet endpoint reconstruction for v2.

All inputs are read-only. Generated files are confined to the v2 output tree.
Use ``extract`` after Comet searches and ``analyze`` after all desired runs finish.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import sys
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

BASE = Path(".")
V2 = Path(".")
MGF_DIR = V2 / "inputs/mgf_full"
PEP_DIR = V2 / "results/pepxml"
RESULT_DIR = V2 / "results"
EXACT_HLA_PATH = V2 / "inputs/stern_table_s1_hla.tsv"
AL_RE = re.compile(r"(?:HLA-)?([ABC])\*?(\d{1,2})[:*]?(\d{2})", re.I)
AA = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]{7,25}$")
DONOR_PASS = {"HP18101", "HP 20289", "R360", "R361", "R369"}
HISTORIC_RUN_DONOR = {
    "20180420_Stern_Islets_W632": "HP18101",
    "20200228_IsletsW632_OTOT_lumos": "R360",
    "20200302_IsletsW6322_OTOT_lumos": "R361",
    "20200724_Lumos_Stern_W632_OTOT_5": "R369",
    "20201202_Lumos_Stern_Islet_Set1_W632_5": "HP 20289",
    "20201202_Lumos_Stern_Islet_Set2_W632_10": "HP 20289",
    "20190405_Stern_Spleen_W632": "HP19026",
    "20200228_SpleenW632_OTOT_lumos": "R360",
    "20201202_Lumos_Stern_Spleen_W632_15": "HP 20289",
}
EXPECTED_MGF_BYTES = {
    "20180420_Stern_Islets_W632.mgf": 13185889,
    "20190405_Stern_Islet_Control_W632.mgf": 35322367,
    "20190405_Stern_Spleen_W632.mgf": 27623873,
    "20200228_IsletsW632_OTOT_lumos.mgf": 44988769,
    "20200228_SpleenW632_OTOT_lumos.mgf": 45664693,
    "20200302_IsletsW6322_OTOT_lumos.mgf": 41413600,
    "20200724_Lumos_Stern_W632_OTOT_5.mgf": 28246794,
    "20201202_Lumos_Stern_Islet_Set1_W632_5.mgf": 29769383,
    "20201202_Lumos_Stern_Islet_Set2_W632_10.mgf": 22595113,
    "20201202_Lumos_Stern_Spleen_W632_15.mgf": 102755608,
}


def norm_hla(value: str) -> set[str]:
    return {f"{g.upper()}*{int(a):02d}:{b}" for g, a, b in AL_RE.findall(str(value).upper().replace("HLA-", ""))}


def write_tsv(path: Path, rows: list[dict], columns: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if columns is None:
        columns = list(rows[0]) if rows else []
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig", errors="replace") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig", errors="replace") as f:
        return list(csv.DictReader(f))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def parse_fasta(path: Path) -> tuple[dict[str, str], dict[str, list[str]]]:
    seq_by_id: dict[str, str] = {}
    ids_by_seq: dict[str, list[str]] = defaultdict(list)
    header = None
    seq: list[str] = []
    def commit() -> None:
        nonlocal header, seq
        if header is not None:
            token = header.split()[0]
            value = "".join(seq).upper()
            seq_by_id[token] = value
            ids_by_seq[value].append(token)
    with path.open(encoding="ascii", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                commit()
                header, seq = line[1:], []
            else:
                seq.append(line)
        commit()
    return seq_by_id, ids_by_seq


def pepxml_rank1(path: Path) -> list[dict]:
    run = path.name[:-8] if path.name.endswith(".pep.xml") else path.stem
    rows: list[dict] = []
    for _, sq in ET.iterparse(path, events=("end",)):
        if sq.tag.rsplit("}", 1)[-1] != "spectrum_query":
            continue
        query = dict(sq.attrib)
        hits = [el for el in sq.iter() if el.tag.rsplit("}", 1)[-1] == "search_hit" and el.get("hit_rank") == "1"]
        if hits:
            hit = hits[0]
            protein = hit.get("protein", "")
            alt = [el.get("protein", "") for el in hit if el.tag.rsplit("}", 1)[-1] == "alternative_protein"]
            scores = {}
            mods: list[str] = []
            modified_peptide = ""
            for el in hit.iter():
                tag = el.tag.rsplit("}", 1)[-1]
                if tag == "search_score":
                    scores[el.get("name", "")] = el.get("value", "")
                elif tag == "modification_info":
                    modified_peptide = el.get("modified_peptide", "")
                    if el.get("mod_nterm_mass"):
                        mods.append(f"Nterm:{el.get('mod_nterm_mass')}")
                    if el.get("mod_cterm_mass"):
                        mods.append(f"Cterm:{el.get('mod_cterm_mass')}")
                elif tag == "mod_aminoacid_mass":
                    mods.append(f"{el.get('position','')}:{el.get('mass','')}")
            decoy = protein.startswith(("DECOY_", "REV_", "XXX_"))
            try:
                xcorr = float(scores.get("xcorr", "nan"))
            except ValueError:
                xcorr = math.nan
            scan = query.get("start_scan", "")
            rows.append({
                "run": run, "spectrum_id": query.get("spectrum", ""),
                "scan": scan, "charge": query.get("assumed_charge", ""),
                "peptide": hit.get("peptide", "").upper(),
                "modified_peptide": modified_peptide,
                "modifications": ";".join(mods),
                "rank": hit.get("hit_rank", "1"), "protein": protein,
                "alternative_proteins": ";".join(alt),
                "decoy": str(decoy).lower(), "xcorr": xcorr,
                "deltacn": scores.get("deltacn", ""),
                "sp": scores.get("sp", ""),
                "precursor_neutral_mass": query.get("precursor_neutral_mass", ""),
                "pep_xml": str(path), "legacy_accepted": "false",
                "qvalue": "", "tdc_fdr_at_score": "",
            })
        sq.clear()
    return rows


def assign_legacy(rows: list[dict]) -> None:
    by_run: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        by_run[r["run"]].append(r)
    for run_rows in by_run.values():
        run_rows.sort(key=lambda r: (-r["xcorr"] if math.isfinite(r["xcorr"]) else math.inf, r["scan"], r["spectrum_id"]))
        targets = decoys = 0
        for r in run_rows:
            if r["decoy"] == "true":
                decoys += 1
            else:
                targets += 1
            fdr = decoys / targets if targets else math.inf
            r["legacy_accepted"] = str(r["decoy"] == "false" and fdr <= 0.01).lower()


def assign_qvalues(rows: list[dict]) -> list[dict]:
    """Deduplicate Comet rank-1 charge queries to one winner per run/start_scan,
    then compute per-run concatenated target-decoy q-values with tied XCorr grouped.
    """
    by_scan: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for r in rows:
        by_scan[(r["run"], r["scan"])].append(r)
    winners = []
    for key, options in by_scan.items():
        options.sort(key=lambda r: (-r["xcorr"] if math.isfinite(r["xcorr"]) else math.inf,
                                   r["decoy"] == "true", r["peptide"], r["charge"], r["spectrum_id"]))
        winner = dict(options[0])
        winner["rank1_charge_rows_for_scan"] = len(options)
        winner["alternative_charge_psms"] = json.dumps([
            {"charge": x["charge"], "peptide": x["peptide"], "protein": x["protein"],
             "decoy": x["decoy"], "xcorr": x["xcorr"], "spectrum_id": x["spectrum_id"]}
            for x in options[1:]
        ], separators=(",", ":"))
        winners.append(winner)
    by_run: dict[str, list[dict]] = defaultdict(list)
    for r in winners:
        by_run[r["run"]].append(r)
    for run_rows in by_run.values():
        run_rows.sort(key=lambda r: (-r["xcorr"] if math.isfinite(r["xcorr"]) else math.inf))
        groups: list[list[dict]] = []
        for r in run_rows:
            if not groups or groups[-1][0]["xcorr"] != r["xcorr"]:
                groups.append([r])
            else:
                groups[-1].append(r)
        fdrs: list[float] = []
        targets = decoys = 0
        for group in groups:
            targets += sum(1 for r in group if r["decoy"] == "false")
            decoys += sum(1 for r in group if r["decoy"] == "true")
            fdrs.append(decoys / targets if targets else math.inf)
        q = math.inf
        for i in range(len(groups) - 1, -1, -1):
            q = min(q, fdrs[i])
            for r in groups[i]:
                r["tdc_fdr_at_score"] = fdrs[i]
                r["qvalue"] = q
    return winners


def extract(pep_dir: Path) -> None:
    pep_files = sorted(pep_dir.glob("*.pep.xml"))
    raw = []
    for p in pep_files:
        raw.extend(pepxml_rank1(p))
    assign_legacy(raw)
    winners = assign_qvalues(raw)
    fasta, aliases = parse_fasta(BASE / "work/msv/msv_db.fasta")
    for group in (raw, winners):
        for r in group:
            proteins = [r["protein"]] + [x for x in r.get("alternative_proteins", "").split(";") if x]
            target_aliases = [p.removeprefix("DECOY_") for p in proteins if not p.startswith("DECOY_")]
            r["target_aliases"] = ";".join(target_aliases)
            r["alias_sequences"] = ";".join(sorted({fasta.get(a, "") for a in target_aliases if fasta.get(a)}))
            r["exact_fasta_aliases_for_peptide"] = ";".join(aliases.get(r["peptide"], [])) if r["decoy"] == "false" else ""
            r["il_equivalent_candidate_match"] = ""
    cols = ["run", "spectrum_id", "scan", "charge", "peptide", "modified_peptide", "modifications",
            "rank", "protein", "alternative_proteins", "target_aliases", "alias_sequences",
            "exact_fasta_aliases_for_peptide", "decoy", "xcorr", "deltacn", "sp",
            "precursor_neutral_mass", "legacy_accepted", "tdc_fdr_at_score", "qvalue",
            "rank1_charge_rows_for_scan", "alternative_charge_psms", "pep_xml"]
    write_tsv(RESULT_DIR / "rank1_charge_rows.tsv", raw, cols)
    write_tsv(RESULT_DIR / "rank1_psms.tsv", winners, cols)
    write_tsv(RESULT_DIR / "pep_run_inventory.tsv", [
        {"run": p.name[:-8], "pep_xml": str(p), "rank1_charge_rows": sum(1 for r in raw if r["run"] == p.name[:-8]),
         "unique_scans": sum(1 for r in winners if r["run"] == p.name[:-8]),
         "multi_charge_scans": sum(1 for r in winners if r["run"] == p.name[:-8] and r["rank1_charge_rows_for_scan"] > 1)}
        for p in pep_files
    ])
    print(f"pepXML={len(pep_files)} rank1_charge_rows={len(raw)} unique_run_scans={len(winners)}")


def candidate_meta() -> tuple[dict[str, dict], set[str], dict[str, list[str]]]:
    rows = read_csv(BASE / "work/candidates/nc2025_all_raw.csv")
    meta: dict[str, dict] = {}
    all_sources: dict[str, list[str]] = defaultdict(list)
    for r in rows:
        p = r.get("peptide", "").strip().upper()
        if not AA.fullmatch(p):
            continue
        hla = norm_hla(r.get("HLA", ""))
        all_sources[p].append(json.dumps({"hla": sorted(hla), "ta_type": r.get("TA_type", ""),
                                          "cohort": r.get("cohort", ""), "gene": r.get("Gene_symbol", "")},
                                         ensure_ascii=False, separators=(",", ":")))
        if p not in meta:
            meta[p] = {"hla": hla, "ta": r.get("TA_type", ""), "cohort": r.get("cohort", ""),
                       "gene": r.get("Gene_symbol", ""), "n_sources": 0}
        meta[p]["n_sources"] += 1
    ae = {p for p, m in meta.items() if "aetsa" in m["ta"].lower()}
    return meta, ae, all_sources


def donor_inputs() -> tuple[dict[str, set[str]], dict[str, set[str]]]:
    ligands: dict[str, set[str]] = defaultdict(set)
    for r in read_csv(BASE / "work/msv/islet_ligands.csv"):
        d, p = r.get("Donor", "").strip(), r.get("Peptide", "").strip().upper()
        if d and p:
            ligands[d].add(p)
    hla: dict[str, set[str]] = defaultdict(set)
    for r in read_csv(BASE / "work/hla/cellrep_islet_classI.csv"):
        d = r.get("Donor", "").strip()
        if d:
            hla[d] |= norm_hla(r.get("HLA", ""))
    return ligands, hla


def exact_typed_donor_hla() -> tuple[dict[str, set[str]], list[dict[str, str]]]:
    rows = read_tsv(EXACT_HLA_PATH)
    exact = {}
    for r in rows:
        if r.get("typing_resolution", "").lower() == "two_field" and r.get("alleles_2field", "").strip():
            exact[r["donor"].strip()] = norm_hla(r["alleles_2field"])
    return exact, rows


def run_tissue(run: str) -> str:
    s = run.lower()
    if "spleen" in s:
        return "spleen"
    if "islet" in s:
        return "islet"
    return "unresolved"


def best_mapping(rates: dict[str, float]) -> tuple[str, float, float, bool]:
    ordered = sorted(rates.items(), key=lambda kv: (-kv[1], kv[0]))
    if not ordered:
        return "", 0.0, 0.0, False
    top, top_rate = ordered[0]
    second_rate = ordered[1][1] if len(ordered) > 1 else 0.0
    margin = top_rate - second_rate
    # Apply the historic mapping rule literally: unique highest recovery and >=70%.
    return top, top_rate, margin, (top_rate >= 0.70 and margin > 0.0)


def analyze() -> None:
    meta, ae, all_sources = candidate_meta()
    ligands, donor_hla = donor_inputs()
    canonical_rows = read_tsv(BASE / "out/step6/TABLE_S_evidence_grading.tsv")
    canonical = {r.get("peptide", "").upper(): r for r in canonical_rows}
    fasta, fasta_aliases = parse_fasta(BASE / "work/msv/msv_db.fasta")
    psm_path = RESULT_DIR / "rank1_psms.tsv"
    psms = read_tsv(psm_path)
    by_run: dict[str, list[dict]] = defaultdict(list)
    for r in psms:
        by_run[r["run"]].append(r)
    qaccepted = {run: {r["peptide"] for r in rows if r["decoy"] == "false" and float(r["qvalue"] or "inf") <= .01}
                 for run, rows in by_run.items()}
    legacy_accepted = {run: {r["peptide"] for r in rows if r["decoy"] == "false" and r["legacy_accepted"] == "true"}
                       for run, rows in by_run.items()}
    run_rows = []
    primary_map: dict[str, str] = {}
    for run in sorted(by_run):
        rates_legacy, rates_q = {}, {}
        for donor, controls in ligands.items():
            denom = len(controls)
            rates_legacy[donor] = len(controls & legacy_accepted.get(run, set())) / denom if denom else 0.0
            rates_q[donor] = len(controls & qaccepted.get(run, set())) / denom if denom else 0.0
        top_l, rate_l, margin_l, pass_l = best_mapping(rates_legacy)
        top_q, rate_q, margin_q, pass_q = best_mapping(rates_q)
        historical = HISTORIC_RUN_DONOR.get(run, "")
        tissue = run_tissue(run)
        mapping_status = "unresolved"
        chosen = ""
        # Use an all-donor competition and require the historical >=70% rule.
        if tissue == "islet" and pass_l and top_l == historical:
            mapping_status, chosen = "recomputed_unique_control_match", top_l
        elif tissue == "islet" and pass_l and not historical:
            mapping_status, chosen = "recomputed_unique_control_match_no_historical_label", top_l
        elif tissue == "islet" and historical:
            mapping_status = "historical_mapping_not_confirmed_by_current_controls"
        elif tissue == "spleen" and historical:
            mapping_status = "historical_spleen_pairing_nonjudgeable"
            chosen = historical
        if chosen and tissue == "islet" and chosen in DONOR_PASS and mapping_status.startswith("recomputed_unique"):
            primary_map[run] = chosen
        run_rows.append({
            "run": run, "donor": chosen, "eligible": str(chosen in DONOR_PASS and tissue == "islet" and mapping_status.startswith("recomputed_unique")).lower(),
            "tissue_from_filename": tissue, "historical_donor": historical,
            "legacy_competitor_top": top_l, "legacy_top_recovery": f"{rate_l:.6f}",
            "legacy_top_margin": f"{margin_l:.6f}", "legacy_unique_top_ge70": str(pass_l).lower(),
            "q01_competitor_top": top_q, "q01_top_recovery": f"{rate_q:.6f}",
            "q01_top_margin": f"{margin_q:.6f}", "q01_unique_top_ge70": str(pass_q).lower(),
            "mapped_donor_used": chosen, "mapping_status": mapping_status,
            "gate_judgeable": str(chosen in DONOR_PASS and tissue == "islet" and mapping_status.startswith("recomputed_unique")).lower(),
            "legacy_control_recovery_all_donors": json.dumps({d: round(v, 6) for d, v in sorted(rates_legacy.items())}, separators=(",", ":")),
            "q01_control_recovery_all_donors": json.dumps({d: round(v, 6) for d, v in sorted(rates_q.items())}, separators=(",", ":")),
        })
    write_tsv(RESULT_DIR / "run_donor_hla_evidence.tsv", run_rows)
    donor_gate_rows = []
    for donor in sorted(set(donor_hla) | set(ligands)):
        supporting = [r for r in run_rows if r["mapped_donor_used"] == donor and r["tissue_from_filename"] == "islet"]
        donor_gate_rows.append({
            "donor": donor, "v1_frozen_gate_status": "pass" if donor in DONOR_PASS else "fail_or_not_in_v1_pass_set",
            "historical_gate_scope": "5/6 pass; HP19026 was reported 5.1% and nonjudgeable",
            "recomputed_islet_runs": ";".join(r["run"] for r in supporting),
            "recomputed_legacy_control_recovery_max": max((r["legacy_top_recovery"] for r in supporting), default=""),
            "recomputed_q01_control_recovery_max": max((r["q01_top_recovery"] for r in supporting), default=""),
            "recomputed_gate_status": "pass" if supporting and any(float(r["legacy_top_recovery"]) >= 0.70 for r in supporting) else "not_pass_or_not_recomputed",
            "genotype_scope": "Table S5 ligand HLA union is not genotyping; see separate Table S1 exact-typed analysis",
        })
    write_tsv(RESULT_DIR / "donor_gate_recomputed.tsv", donor_gate_rows)

    # Donor-level HLA compatibility and two independent count rules.
    count_maps: dict[tuple[str, str, str], set[tuple[str, str]]] = defaultdict(set)
    for r in psms:
        if r["decoy"] != "false" or run_tissue(r["run"]) != "islet":
            continue
        p = r["peptide"].upper()
        if p not in ae:
            continue
        donor = primary_map.get(r["run"], "")
        if not donor:
            continue
        if r["legacy_accepted"] == "true":
            count_maps[(p, donor, "legacy")].add((r["run"], r["scan"]))
        if float(r["qvalue"] or "inf") <= .01:
            count_maps[(p, donor, "q01")].add((r["run"], r["scan"]))

    donor_rows = []
    per_candidate_rows = []
    psm_hits_by_method: dict[str, dict[str, dict[str, set[tuple[str, str]]]]] = {
        "legacy": defaultdict(lambda: defaultdict(set)), "q01": defaultdict(lambda: defaultdict(set))}
    # Candidate metadata compatibility is the first-row rule frozen by Step 5; all-source metadata retained.
    for p in sorted(ae):
        m = meta[p]
        for donor in sorted(donor_hla):
            compat = bool(m["hla"] & donor_hla[donor])
            mapped_runs = [run for run, d in primary_map.items() if d == donor]
            legacy_scans = set(count_maps.get((p, donor, "legacy"), set()))
            q_scans = set(count_maps.get((p, donor, "q01"), set()))
            # Keep the per-donor evidence and expose run+scan IDs.
            legacy_n, q_n = len(legacy_scans), len(q_scans)
            is_judged = donor in DONOR_PASS and bool(mapped_runs)
            status = ("not_HLA_compatible" if not compat else
                      "donor_not_judgeable_or_unmapped" if not is_judged else
                      "detected_ge2_spectra_legacy" if legacy_n >= 2 else
                      "single_spectrum_legacy" if legacy_n == 1 else
                      "detected_ge2_spectra_q01" if q_n >= 2 else
                      "single_spectrum_q01" if q_n == 1 else "no_qualifying_PSM")
            donor_rows.append({
                "peptide": p, "TA_type_first_row": m["ta"], "candidate_HLA_first_row": ",".join(sorted(m["hla"])),
                "candidate_HLA_all_sources": " || ".join(all_sources[p]), "donor": donor,
                "TableS5_ligand_HLA_union_not_genotype": ",".join(sorted(donor_hla[donor])),
                "HLA_compatible_under_v1_TableS5_union": str(compat).lower(),
                "v1_frozen_5donor_gate_set": str(donor in DONOR_PASS).lower(),
                "HLA_evidence_scope": "Table S5 per-ligand HLA union; not donor genotype; exact-typed branch separate",
                "mapped_islet_runs": ";".join(sorted(mapped_runs)),
                "legacy_distinct_spectra": legacy_n, "legacy_spectrum_ids": ";".join(f"{r}:{s}" for r, s in sorted(legacy_scans)),
                "q01_distinct_spectra": q_n, "q01_spectrum_ids": ";".join(f"{r}:{s}" for r, s in sorted(q_scans)),
                "evidence_status": status,
                "canonical_status_from_v1_table": canonical.get(p, {}).get("exact_canonical", "") +
                    ";IL=" + canonical.get(p, {}).get("il_canonical", ""),
            })
            if compat and is_judged:
                if legacy_n >= 2:
                    psm_hits_by_method["legacy"][p][donor] = legacy_scans
                if q_n >= 2:
                    psm_hits_by_method["q01"][p][donor] = q_scans
    write_tsv(RESULT_DIR / "aetsa_candidate_by_donor.tsv", donor_rows)

    # Preserve candidate-level PSMs with exact and I/L-equivalent alternatives.
    il_index: dict[str, set[str]] = defaultdict(set)
    for candidate in meta:
        il_index[candidate.replace("I", "J").replace("L", "J")].add(candidate)
    candidate_psm_rows, unresolved_psm_rows, spleen_psm_rows = [], [], []
    for r in psms:
        if r["decoy"] != "false":
            continue
        peptide = r["peptide"].upper()
        exact = peptide in meta
        il_alts = sorted(il_index.get(peptide.replace("I", "J").replace("L", "J"), set()) - {peptide})
        if not exact and not il_alts:
            continue
        run = r["run"]
        tissue = run_tissue(run)
        donor = primary_map.get(run, HISTORIC_RUN_DONOR.get(run, ""))
        canon = canonical.get(peptide, {})
        il_canon = [f"{p}:exact={canonical.get(p, {}).get('exact_canonical','')}:IL={canonical.get(p, {}).get('il_canonical','')}"
                    for p in il_alts]
        evidence = {
            "run": run, "tissue": tissue, "recomputed_donor": primary_map.get(run, ""),
            "historical_or_nonjudgeable_pairing": donor, "scan": r["scan"], "spectrum_id": r["spectrum_id"],
            "charge": r["charge"], "peptide": peptide, "exact_nc2025_candidate": str(exact).lower(),
            "exact_aeTSA_candidate": str(peptide in ae).lower(), "il_equivalent_candidate_alternatives": ";".join(il_alts),
            "protein_alias": r["protein"], "alternative_proteins": r["alternative_proteins"],
            "modified_peptide": r["modified_peptide"], "modifications": r["modifications"],
            "rank": r["rank"], "decoy": r["decoy"], "xcorr": r["xcorr"],
            "qvalue_per_run": r["qvalue"], "legacy_accepted": r["legacy_accepted"],
            "exact_canonical_in_existing_evidence": canon.get("exact_canonical", ""),
            "il_canonical_in_existing_evidence": canon.get("il_canonical", ""),
            "il_alternative_canonical_evidence": ";".join(il_canon),
            "canonical_scope": "Existing TABLE_S evidence; the v2 candidate FASTA is not a full human proteome search",
        }
        candidate_psm_rows.append(evidence)
        if tissue == "spleen":
            spleen_psm_rows.append(evidence)
        if tissue != "islet" or run not in primary_map or not exact or peptide not in ae:
            unresolved_psm_rows.append(evidence)
    candidate_psm_cols = ["run", "tissue", "recomputed_donor", "historical_or_nonjudgeable_pairing", "scan",
                          "spectrum_id", "charge", "peptide", "exact_nc2025_candidate", "exact_aeTSA_candidate",
                          "il_equivalent_candidate_alternatives", "protein_alias", "alternative_proteins",
                          "modified_peptide", "modifications", "rank", "decoy", "xcorr", "qvalue_per_run",
                          "legacy_accepted", "exact_canonical_in_existing_evidence", "il_canonical_in_existing_evidence",
                          "il_alternative_canonical_evidence", "canonical_scope"]
    write_tsv(RESULT_DIR / "nc2025_candidate_psm_evidence.tsv", candidate_psm_rows, candidate_psm_cols)
    write_tsv(RESULT_DIR / "unresolved_or_ambiguous_candidate_psms.tsv", unresolved_psm_rows, candidate_psm_cols)
    write_tsv(RESULT_DIR / "spleen_candidate_psms.tsv", spleen_psm_rows, candidate_psm_cols)

    # Audit each pepXML scan identifier against source MGF SCANS and preserve TITLE.
    scan_maps: dict[str, dict[str, list[dict]]] = {}
    scan_inventory = []
    for mgf in sorted(MGF_DIR.glob("*.mgf")):
        run = mgf.stem
        by_scan: dict[str, list[dict]] = defaultdict(list)
        blocks = 0
        in_block = False
        title = scan = charge = ""
        with mgf.open(encoding="utf-8", errors="replace") as f:
            for line in f:
                line = line.strip()
                if line == "BEGIN IONS":
                    blocks += 1
                    in_block, title, scan, charge = True, "", "", ""
                elif in_block and line.startswith("TITLE="):
                    title = line[6:]
                elif in_block and line.startswith("SCANS="):
                    scan = line[6:]
                elif in_block and line.startswith("CHARGE="):
                    charge = line[7:]
                elif line == "END IONS" and in_block:
                    if scan:
                        by_scan[scan].append({"block": blocks, "title": title, "charge": charge})
                    in_block = False
        scan_maps[run] = by_scan
        run_psms = [r for r in psms if r["run"] == run]
        matched = sum(1 for r in run_psms if r["scan"] in by_scan)
        scan_inventory.append({"run": run, "mgf_bytes": mgf.stat().st_size, "mgf_blocks": blocks,
                               "blocks_with_scans": sum(len(v) for v in by_scan.values()),
                               "unique_scan_ids": len(by_scan),
                               "duplicate_scan_ids": sum(1 for v in by_scan.values() if len(v) > 1),
                               "rank1_psms": len(run_psms), "rank1_scan_ids_matched": matched,
                               "all_rank1_scans_present": str(matched == len(run_psms)).lower()})
    write_tsv(RESULT_DIR / "mgf_scan_inventory.tsv", scan_inventory)
    scan_audit = []
    for r in psms:
        found = scan_maps.get(r["run"], {}).get(r["scan"], [])
        first = found[0] if found else {}
        audit = dict(r)
        audit.update({"mgf_scan_match": str(bool(found)).lower(), "mgf_block_number": first.get("block", ""),
                      "mgf_title": first.get("title", ""), "mgf_charge": first.get("charge", ""),
                      "mgf_scan_duplicate_occurrences": len(found)})
        scan_audit.append(audit)
    scan_cols = list(psms[0].keys()) + ["mgf_scan_match", "mgf_block_number", "mgf_title", "mgf_charge",
                                       "mgf_scan_duplicate_occurrences"] if psms else []
    write_tsv(RESULT_DIR / "rank1_psms_with_mgf_scan.tsv", scan_audit, scan_cols)

    endpoint_rows = []
    for method in ("legacy", "q01"):
        den_any = den_both = det_any = det_both = det_both_ge1 = 0
        for p in ae:
            # Frozen denominator is based on the five judgeable donor genotypes, not on
            # whether this replay happened to recover a run for each donor.
            compatible_judged = [d for d in DONOR_PASS if meta[p]["hla"] & donor_hla.get(d, set())]
            den_any += int(len(compatible_judged) >= 1)
            den_both += int(len(compatible_judged) >= 2)
            hit_donors = set(psm_hits_by_method[method].get(p, {}))
            det_any += int(bool(hit_donors))
            det_both += int(len(hit_donors) >= 2)
            donors_with_one_or_more = {d for d in DONOR_PASS
                                       if meta[p]["hla"] & donor_hla.get(d, set())
                                       and d in primary_map.values()
                                       and len(count_maps.get((p, d, method), set())) >= 1}
            det_both_ge1 += int(len(donors_with_one_or_more) >= 2)
        endpoint_rows.extend([
            {"method": method, "endpoint": "any", "aeTSA_unique_candidates": len(ae), "denominator": den_any, "detected": det_any,
             "detection_rule": "HLA compatible judgeable donor with >=2 distinct run+scan rank1 spectra"},
            {"method": method, "endpoint": "both", "aeTSA_unique_candidates": len(ae), "denominator": den_both, "detected": det_both,
             "detection_rule": ">=2 distinct run+scan rank1 spectra in each of >=2 HLA-compatible judgeable donors"},
            {"method": method, "endpoint": "both_ge1_spectrum_per_donor_sensitivity", "aeTSA_unique_candidates": len(ae),
             "denominator": den_both, "detected": det_both_ge1,
             "detection_rule": ">=1 distinct run+scan rank1 spectrum in each of >=2 HLA-compatible judgeable donors"},
        ])
    write_tsv(RESULT_DIR / "aetsa_endpoint_recomputed.tsv", endpoint_rows)

    frozen = {r["metric"]: r["value"] for r in read_tsv(BASE / "out/step8/FROZEN_NUMBERS.tsv")}
    legacy_rows = { (r["method"], r["endpoint"]): r for r in endpoint_rows }
    comparisons = []
    for endpoint, den_key, det_key in (("any", "islet_aetsa_den_any", "islet_aetsa_det_any"),
                                       ("both", "islet_aetsa_den_both", "islet_aetsa_det_both")):
        fden, fdet = int(frozen[den_key]), int(frozen[det_key])
        for method in ("legacy", "q01"):
            r = legacy_rows[(method, endpoint)]
            comparisons.append({"endpoint": endpoint, "method": method, "frozen_denominator": fden,
                                "recomputed_denominator": int(r["denominator"]),
                                "denominator_matches_frozen": str(int(r["denominator"]) == fden).lower(),
                                "frozen_detected": fdet, "recomputed_detected": int(r["detected"]),
                                "detection_matches_frozen": str(int(r["detected"]) == fdet).lower(),
                                "frozen_source": "out/step8/FROZEN_NUMBERS.tsv; v1 Step5 hard-coded islet det=0"})
        if endpoint == "both":
            for method in ("legacy", "q01"):
                r = legacy_rows[(method, "both_ge1_spectrum_per_donor_sensitivity")]
                comparisons.append({"endpoint": "both_ge1_spectrum_per_donor_sensitivity", "method": method,
                                    "frozen_denominator": fden, "recomputed_denominator": int(r["denominator"]),
                                    "denominator_matches_frozen": str(int(r["denominator"]) == fden).lower(),
                                    "frozen_detected": fdet, "recomputed_detected": int(r["detected"]),
                                    "detection_matches_frozen": str(int(r["detected"]) == fdet).lower(),
                                    "frozen_source": "out/step8/FROZEN_NUMBERS.tsv; sensitivity threshold stated explicitly"})
    write_tsv(RESULT_DIR / "frozen_vs_recomputed.tsv", comparisons)

    # Independent exact-typed Table S1 branch. Only donors with explicit two-field
    # alleles enter; serology/family-only donors remain unresolved in this branch.
    exact_hla, exact_hla_rows = exact_typed_donor_hla()
    write_tsv(RESULT_DIR / "exact_typed_donor_hla_source.tsv", exact_hla_rows)
    typed_candidate_rows = []
    typed_hit_donors: dict[str, dict[str, set[str]]] = {
        "legacy": defaultdict(set), "q01": defaultdict(set)}
    for p in sorted(ae):
        for donor, alleles in sorted(exact_hla.items()):
            compat = bool(meta[p]["hla"] & alleles)
            mapped_runs = [run for run, d in primary_map.items() if d == donor]
            legacy_scans = set(count_maps.get((p, donor, "legacy"), set()))
            q_scans = set(count_maps.get((p, donor, "q01"), set()))
            typed_candidate_rows.append({
                "peptide": p, "TA_type_first_row": meta[p]["ta"],
                "candidate_HLA_first_row": ",".join(sorted(meta[p]["hla"])), "donor": donor,
                "exact_two_field_HLA_from_TableS1": ";".join(sorted(alleles)),
                "HLA_compatible_exact_typed": str(compat).lower(),
                "mapped_islet_runs": ";".join(sorted(mapped_runs)),
                "legacy_distinct_spectra": len(legacy_scans),
                "legacy_spectrum_ids": ";".join(f"{run}:{scan}" for run, scan in sorted(legacy_scans)),
                "q01_distinct_spectra": len(q_scans),
                "q01_spectrum_ids": ";".join(f"{run}:{scan}" for run, scan in sorted(q_scans)),
                "source_text_line": next((r.get("source_text_line", "") for r in exact_hla_rows if r.get("donor") == donor), ""),
                "evidence_status": "not_HLA_compatible" if not compat else
                                   "mapped_donor_not_observed" if not mapped_runs else
                                   "ge2_legacy_spectra" if len(legacy_scans) >= 2 else
                                   "single_legacy_spectrum" if len(legacy_scans) == 1 else
                                   "ge2_q01_spectra" if len(q_scans) >= 2 else
                                   "single_q01_spectrum" if len(q_scans) == 1 else "no_qualifying_PSM",
            })
            if compat and mapped_runs:
                if len(legacy_scans) >= 2:
                    typed_hit_donors["legacy"][p].add(donor)
                if len(q_scans) >= 2:
                    typed_hit_donors["q01"][p].add(donor)
    write_tsv(RESULT_DIR / "aetsa_candidate_by_exact_typed_donor.tsv", typed_candidate_rows)
    typed_endpoint_rows = []
    typed_comparison_rows = []
    for method in ("legacy", "q01"):
        den_any = den_both = det_any = det_both = 0
        for p in ae:
            compatible = [d for d, hla in exact_hla.items() if meta[p]["hla"] & hla]
            den_any += int(len(compatible) >= 1)
            den_both += int(len(compatible) >= 2)
            hit_donors = typed_hit_donors[method].get(p, set())
            det_any += int(bool(hit_donors))
            det_both += int(len(hit_donors) >= 2)
        for endpoint, denominator, detected in (("any", den_any, det_any), ("both", den_both, det_both)):
            typed_endpoint_rows.append({"method": method, "endpoint": endpoint,
                                        "exact_two_field_donors": ";".join(sorted(exact_hla)),
                                        "denominator": denominator, "detected": detected,
                                        "detection_rule": "HLA-compatible exact-typed donor with >=2 distinct run+scan PSMs"})
            v1 = next(r for r in endpoint_rows if r["method"] == method and r["endpoint"] == endpoint)
            typed_comparison_rows.append({"method": method, "endpoint": endpoint,
                                          "v1_TableS5_union_denominator": v1["denominator"],
                                          "exact_typed_denominator": denominator,
                                          "v1_TableS5_union_detected": v1["detected"],
                                          "exact_typed_detected": detected,
                                          "note": "Not a frozen exact-typed endpoint; serology/family-only donors excluded"})
    write_tsv(RESULT_DIR / "aetsa_exact_typed_endpoint_recomputed.tsv", typed_endpoint_rows)
    write_tsv(RESULT_DIR / "aetsa_v1_union_vs_exact_typed.tsv", typed_comparison_rows)

    # Recompute the existing PDAC candidate entrapment endpoint from frozen score rows.
    pdac_scores = read_tsv(BASE / "out/pdac_entrap/scores.tsv")
    best: dict[tuple[str, str, str], float] = {}
    for r in pdac_scores:
        k = (r["sample"], r["peptide"].upper(), r["source"])
        best[k] = max(best.get(k, -math.inf), float(r["xcorr"]))
    per_sample: dict[str, dict[str, dict[str, float]]] = defaultdict(lambda: {"pdc": {}, "ent": {}})
    for (sample, peptide, source), score in best.items():
        if source in ("pdc", "ent"):
            per_sample[sample][source][peptide] = score
    compat_by_peptide: dict[str, set[str]] = defaultdict(set)
    pdac_alleles_by_peptide: dict[str, set[str]] = defaultdict(set)
    for r in read_tsv(BASE / "out/step4/STEP4_PDAC_HLA.tsv"):
        pdac_alleles_by_peptide[r["peptide"].upper()] |= norm_hla(r.get("allele", ""))
        if r["sample"] in DONOR_PASS and r["compat"].lower() == "true":
            compat_by_peptide[r["peptide"].upper()].add(r["sample"])
    exact_compat_by_peptide = {p: {d for d, hla in exact_hla.items() if hla & alleles}
                               for p, alleles in pdac_alleles_by_peptide.items()}
    pdac_rows = []
    for sample, group in per_sample.items():
        is_spleen = "spleen" in sample.lower()
        if not sample.startswith("I_"):
            continue
        stem = sample[2:]
        recomputed_donor = primary_map.get(stem, "")
        historical_donor = HISTORIC_RUN_DONOR.get(stem, "")
        donor = recomputed_donor or historical_donor
        pdc, ent = group["pdc"], group["ent"]
        for p, score in pdc.items():
            n_p = sum(1 for v in pdc.values() if v >= score)
            n_e = sum(1 for v in ent.values() if v >= score)
            fdr = (n_e + 1) / max(1, n_p)
            compat_donors = compat_by_peptide.get(p, set())
            exact_compat_donors = exact_compat_by_peptide.get(p, set())
            pdac_rows.append({"sample": sample, "run": stem, "tissue": "spleen" if is_spleen else "islet",
                              "donor_mapping": donor, "recomputed_donor": recomputed_donor,
                              "historical_donor": historical_donor, "peptide": p, "xcorr": score,
                              "n_pdc_ge": n_p, "n_ent_ge": n_e, "entrapment_fdr": fdr,
                              "entrap_q05": str(fdr <= .05).lower(), "compatible_judgeable_donors": ";".join(sorted(compat_donors)),
                              "mapped_donor_compatible": str(bool(recomputed_donor and recomputed_donor in compat_donors)).lower(),
                              "compatible_exact_typed_donors": ";".join(sorted(exact_compat_donors)),
                              "mapped_donor_exact_typed_compatible": str(bool(recomputed_donor and recomputed_donor in exact_compat_donors)).lower(),
                              "judgeable_islet_arm": str(not is_spleen and recomputed_donor in DONOR_PASS).lower()})
    write_tsv(RESULT_DIR / "pdac_existing_score_recompute.tsv", pdac_rows)
    ge1 = sum(1 for p, ds in compat_by_peptide.items() if ds)
    ge2 = sum(1 for p, ds in compat_by_peptide.items() if len(ds) >= 2)
    q05_pairs = [r for r in pdac_rows if r["entrap_q05"] == "true" and r["tissue"] == "islet"]
    q05_judged = [r for r in q05_pairs if r["judgeable_islet_arm"] == "true" and r["mapped_donor_compatible"] == "true"]
    q05_exact_typed = [r for r in q05_pairs if r["recomputed_donor"] in exact_hla and
                       r["mapped_donor_exact_typed_compatible"] == "true"]
    pdac_frozen = {r["metric"]: r["value"] for r in read_tsv(BASE / "out/step7/TABLE_2_pdac_islet.tsv")}
    write_tsv(RESULT_DIR / "pdac_endpoint_recomputed.tsv", [
        {"metric": "compatible_ge1", "recomputed": ge1, "frozen": pdac_frozen.get("compatible_ge1", ""), "matches": str(str(ge1) == pdac_frozen.get("compatible_ge1", "")).lower()},
        {"metric": "compatible_ge2", "recomputed": ge2, "frozen": pdac_frozen.get("compatible_ge2", ""), "matches": str(str(ge2) == pdac_frozen.get("compatible_ge2", "")).lower()},
        {"metric": "islet_q05_sample_peptide_pairs", "recomputed": len(q05_pairs), "frozen": pdac_frozen.get("q05_hits", ""), "matches": str(str(len(q05_pairs)) == pdac_frozen.get("q05_hits", "")).lower()},
        {"metric": "judgeable_compatible_q05_sample_peptide_pairs", "recomputed": len(q05_judged), "frozen": "not separately reported", "matches": ""},
        {"metric": "exact_typed_compatible_ge1", "recomputed": sum(1 for p, ds in exact_compat_by_peptide.items() if ds), "frozen": "not reported", "matches": ""},
        {"metric": "exact_typed_compatible_ge2", "recomputed": sum(1 for p, ds in exact_compat_by_peptide.items() if len(ds) >= 2), "frozen": "not reported", "matches": ""},
        {"metric": "exact_typed_mapped_q05_islet_sample_peptide_pairs", "recomputed": len(q05_exact_typed), "frozen": "not reported", "matches": ""},
    ])

    # Complete source manifest and per-run input inventory for reproducibility.
    source_paths = [BASE / "work/msv/comet.params", BASE / "work/msv/msv_db.fasta",
                    BASE / "work/msv/islet_ligands.csv", BASE / "work/candidates/nc2025_all_raw.csv",
                    BASE / "work/hla/cellrep_islet_classI.csv", BASE / "out/step8/FROZEN_NUMBERS.tsv",
                    BASE / "out/step5/STEP5_ENDPOINTS.tsv", BASE / "out/step6/TABLE_S_evidence_grading.tsv",
                    BASE / "out/pdac_entrap/scores.tsv", BASE / "out/step4/STEP4_PDAC_HLA.tsv",
                    BASE / "out/step7/TABLE_2_pdac_islet.tsv", EXACT_HLA_PATH]
    source_paths += sorted(MGF_DIR.glob("*.mgf")) + sorted(PEP_DIR.glob("*.pep.xml"))
    source_paths += [p for p in [V2 / "scripts/v2_islet_reconstruction.py", V2 / "scripts/run_comet_subset.sh",
                                 V2 / "README.md", V2 / "inputs/mgf_local_sha256.tsv"] if p.exists()]
    manifest_rows = []
    for p in source_paths:
        if p.exists():
            manifest_rows.append({"path": str(p), "bytes": p.stat().st_size, "sha256": sha256(p),
                                  "role": "spectrum_input" if p.suffix.lower() == ".mgf" else
                                         "Comet_search_output" if p.name.endswith(".pep.xml") else "source_or_frozen_input"})
    write_tsv(RESULT_DIR / "INPUT_SHA256.tsv", manifest_rows)
    local_hash_path = V2 / "inputs/mgf_local_sha256.tsv"
    local_hashes = {r.get("file", ""): r for r in read_tsv(local_hash_path)} if local_hash_path.exists() else {}
    transfer_rows = []
    for name, expected_bytes in sorted(EXPECTED_MGF_BYTES.items()):
        path = MGF_DIR / name
        present = path.exists()
        remote_bytes = path.stat().st_size if present else ""
        local_sha = local_hashes.get(name, {}).get("local_sha256", "")
        size_match = bool(present and remote_bytes == expected_bytes)
        remote_sha = sha256(path) if size_match else ""
        transfer_rows.append({
            "file": name,
            "source_local_path": name,
            "local_expected_bytes": expected_bytes, "server_bytes": remote_bytes,
            "server_file_present": str(present).lower(), "size_matches_local_source": str(size_match).lower(),
            "local_sha256": local_sha, "server_sha256": remote_sha,
            "sha256_matches_local_source": str(bool(local_sha and remote_sha and local_sha == remote_sha)).lower(),
            "status": "verified" if size_match and local_sha and local_sha == remote_sha else
                      "size_verified_hash_pending" if size_match else "missing_or_size_mismatch",
        })
    write_tsv(RESULT_DIR / "transfer_verification.tsv", transfer_rows)
    write_tsv(RESULT_DIR / "source_fasta_header_map.tsv", [
        {"header_alias": alias, "sequence": seq,
         "source_alias_class": "candidate" if alias.startswith("cand_") else
                               "islet_ligand" if alias.startswith("isl_") else "other"}
        for alias, seq in sorted(fasta.items())])
    candidate_source_rows = []
    for p in sorted(meta):
        candidate_source_rows.append({"peptide": p, "ta_type_first_row": meta[p]["ta"], "hla_first_row": ",".join(sorted(meta[p]["hla"])),
                                      "aetsa_first_row": str(p in ae).lower(), "source_row_count": meta[p]["n_sources"],
                                      "source_metadata_all_rows": " || ".join(all_sources[p])})
    write_tsv(RESULT_DIR / "nc2025_candidate_metadata.tsv", candidate_source_rows)
    print(f"aeTSA={len(ae)} mapped_judgeable_runs={len(primary_map)} donor_sets={len(DONOR_PASS)}")
    print(f"PDAC score rows={len(pdac_rows)} q05 islet sample-peptide pairs={len(q05_pairs)}; compatible={len(q05_judged)}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("action", choices=["extract", "analyze"])
    ap.add_argument("--base", type=Path, required=True,
                    help="Root containing work/, out/ and incoming/ source tables")
    ap.add_argument("--v2-root", type=Path, required=True,
                    help="Islet replay root containing inputs/ and results/")
    ap.add_argument("--pep-dir", type=Path,
                    help="Directory containing per-run Comet pepXML (defaults to <v2-root>/results/pepxml)")
    args = ap.parse_args()
    global BASE, V2, MGF_DIR, PEP_DIR, RESULT_DIR, EXACT_HLA_PATH
    BASE = args.base.resolve(strict=True)
    V2 = args.v2_root.resolve()
    MGF_DIR = V2 / "inputs/mgf_full"
    PEP_DIR = V2 / "results/pepxml"
    RESULT_DIR = V2 / "results"
    EXACT_HLA_PATH = V2 / "inputs/stern_table_s1_hla.tsv"
    PEP_DIR = args.pep_dir.resolve() if args.pep_dir else PEP_DIR
    if args.action == "extract":
        extract(PEP_DIR)
    else:
        analyze()


if __name__ == "__main__":
    main()

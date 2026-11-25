"""Report recovery of published human-islet HLA-I ligands in the C018 searches."""

import argparse
import csv
import gzip
import math
import re
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path


ALLELE = re.compile(r"([ABC])\*?(\d{2}):?(\d{2})(?!\d)")
PEPTIDE = re.compile(r"[ACDEFGHIKLMNPQRSTVWY]{8,12}")
Q_CUTOFF = 0.01


def read(path, delimiter=","):
    with path.open(encoding="utf-8-sig", newline="") as source:
        yield from csv.DictReader(source, delimiter=delimiter)


def write(path, rows):
    rows = list(rows)
    if rows:
        with path.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, rows[0], delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)


def alleles(value):
    return {f"{gene}*{family}:{second}" for gene, family, second in ALLELE.findall(value.upper())}


def donor_id(value):
    return re.sub(r"[^A-Z0-9]", "", value.upper())


def il_key(peptide):
    return peptide.replace("I", "L")


def pick(row, *keys):
    for key in keys:
        if key in row and row[key] not in (None, ""):
            return row[key]
    raise ValueError(f"None of {keys} available in {tuple(row)}")


def pdac_allele_union(base):
    by_sequence = defaultdict(set)
    for row in read(base / "out/step4/STEP4_PDAC_HLA.tsv", "\t"):
        by_sequence[row["peptide"].strip().upper()].update(alleles(row["allele"]))
    return by_sequence


def controls(base, exact_hla=None):
    all_candidates = set()
    restrictions = defaultdict(lambda: defaultdict(set))
    for row in read(base / "work/candidates/nc2025_all_raw.csv"):
        pep = row["peptide"].strip().upper()
        all_candidates.add(il_key(pep))
        if "aetsa" in row["TA_type"].lower():
            restrictions["NC2025_aeTSA"][pep].update(alleles(row["HLA"]))
    for pep, hla in pdac_allele_union(base).items():
        all_candidates.add(il_key(pep))
        restrictions["PDAC2025"][pep].update(hla)
    candidate_hla_by_length = {
        cohort: {
            length: set().union(*(hla for pep, hla in seqs.items() if len(pep) == length))
            for length in range(8, 13)
        } for cohort, seqs in restrictions.items()
    }

    donor_hla = defaultdict(set)
    if exact_hla:
        for row in read(exact_hla, "\t"):
            if row["typing_resolution"] == "two_field":
                donor_hla[donor_id(row["donor"])].update(alleles(row["alleles_2field"]))
    else:
        for row in read(base / "work/hla/cellrep_islet_classI.csv"):
            donor_hla[donor_id(row["Donor"])].update(alleles(row["HLA"]))

    reference = defaultdict(dict)
    for row in read(base / "work/msv/islet_ligands.csv"):
        pep, donor = row["Peptide"].strip().upper(), donor_id(row["Donor"])
        if not PEPTIDE.fullmatch(pep) or il_key(pep) in all_candidates:
            continue
        try:
            intensity = float(row["Intensity"])
        except ValueError:
            continue
        if not math.isfinite(intensity) or intensity <= 0:
            continue
        hla = alleles(row["HLA_norm"]) & donor_hla[donor]
        if not hla:
            continue
        previous = reference[donor].get(pep)
        if previous is None:
            reference[donor][pep] = (intensity, hla)
        else:
            reference[donor][pep] = (max(intensity, previous[0]), hla | previous[1])
    panel = []
    for donor, peptides in sorted(reference.items()):
        values = sorted(intensity for intensity, _ in peptides.values())
        quartiles = [values[math.ceil(len(values) * p / 4) - 1] for p in (1, 2, 3)]
        for pep, (intensity, hla) in sorted(peptides.items()):
            cohort_hla = {
                cohort: sorted(hla & by_length[len(pep)])
                for cohort, by_length in candidate_hla_by_length.items()
            }
            panel.append({
                "donor": donor, "peptide": pep, "length": len(pep),
                "author_intensity": intensity, "intensity_quartile": 1 + sum(intensity > x for x in quartiles),
                "hla": ";".join(sorted(hla)),
                "nc2025_matching_hla_and_length": ";".join(cohort_hla["NC2025_aeTSA"]),
                "pdac_matching_hla_and_length": ";".join(cohort_hla["PDAC2025"]),
            })
    return panel, restrictions, donor_hla


def load_psms(path):
    accepted = defaultdict(lambda: defaultdict(set))
    examined = set()
    for row in read(path, "\t"):
        run = pick(row, "run", "sample")
        examined.add(run)
        if "tdc_pseudocount" in row and (
            row["tdc_pseudocount"] != "1" or float(row["peptide_qvalue"]) > Q_CUTOFF
        ):
            continue
        if pick(row, "is_decoy", "decoy").strip().lower() in {"1", "true", "yes"}:
            continue
        qvalue = float(pick(row, "qvalue", "q_value", "q", "psm_qvalue"))
        if not math.isfinite(qvalue) or qvalue > Q_CUTOFF:
            continue
        pep = pick(row, "peptide", "sequence").strip().upper()
        scan = pick(row, "scan", "start_scan", "spectrum")
        accepted[run][pep].add(scan)
    return accepted, examined


def interval(k, n):
    if not n:
        return ("", "")
    z = 1.959963984540054
    p = k / n
    denominator = 1 + z * z / n
    center = (p + z * z / (2 * n)) / denominator
    half_width = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denominator
    return round(max(0, center - half_width), 6), round(min(1, center + half_width), 6)


def group_recovery(rows):
    groups = defaultdict(list)
    for row in rows:
        for panel in ("all_reference", "nc2025_matched", "pdac_matched"):
            if panel == "nc2025_matched" and not row["nc2025_matching_hla_and_length"]:
                continue
            if panel == "pdac_matched" and not row["pdac_matching_hla_and_length"]:
                continue
            for quartile in (0, row["intensity_quartile"]):
                groups[(panel, row["donor"], quartile)].append(row)
    for (panel, donor, quartile), values in sorted(groups.items()):
        for least in (1, 2):
            k = sum(row["distinct_scans"] >= least for row in values)
            n = len(values)
            low, high = interval(k, n)
            yield {
                "panel": panel, "donor": donor, "intensity_quartile": quartile or "all",
                "scan_threshold": least, "reference_n": n, "recovered_n": k,
                "recovery": round(k / n, 6), "descriptive_wilson_low": low,
                "descriptive_wilson_high": high,
            }


def probability(probabilities):
    none, exactly_one = 1.0, 0.0
    for p in probabilities:
        exactly_one = exactly_one * (1 - p) + none * p
        none *= 1 - p
    return min(1, max(0, 1 - none)), min(1, max(0, 1 - none - exactly_one))


def conditional_scenarios(recovery, restrictions, donor_hla, tested):
    rows = []
    positive_fraction = 0.20
    rates = {
        (r["panel"], r["donor"], r["intensity_quartile"], r["scan_threshold"]):
        (r["recovery"], r["descriptive_wilson_low"]) for r in recovery
    }
    for cohort, seqs in restrictions.items():
        compatible = {pep: sorted(d for d in tested if hla & donor_hla[d]) for pep, hla in seqs.items()}
        opportunities = [
            {"cohort": cohort, "peptide": pep, "n_compatible_donors": len(ds), "compatible_donors": ";".join(ds)}
            for pep, ds in sorted(compatible.items())
        ]
        yield "opportunities", opportunities
        panel = "nc2025_matched" if cohort == "NC2025_aeTSA" else "pdac_matched"
        for variant in ("JY_1to16_legacy", "JY_1to16_actual_HLA", "islet_Q1_recovery", "islet_Q1_Wilson_low"):
            for endpoint, min_donors in (("any", 1), ("both", 2)):
                estimates = []
                skipped = 0
                for ds in compatible.values():
                    if len(ds) < min_donors:
                        continue
                    if variant.startswith("JY"):
                        length = 2 if variant.endswith("legacy") else len(ds)
                        ps = [0.361] * length
                    else:
                        parts = [rates.get((panel, donor, 1, 2)) for donor in ds]
                        if any(value is None for value in parts):
                            skipped += 1
                            continue
                        ps = [value[int(variant.endswith("low"))] for value in parts]
                    independent = probability(ps)[min_donors - 1]
                    nested = sorted(ps, reverse=True)[min_donors - 1]
                    estimates.append((tuple(ds), independent, nested))
                probabilities = [p for _, p, _ in estimates]
                n = len(probabilities)
                if n:
                    p_zero = math.prod(1 - positive_fraction * p for p in probabilities)
                    grouped = defaultdict(float)
                    for profile, p, _ in estimates:
                        grouped[profile] = max(grouped[profile], p)
                    p_zero_grouped = math.prod(1 - positive_fraction * p for p in grouped.values())
                    nested_mean = sum(p for _, _, p in estimates) / n
                    rows.append({
                        "cohort": cohort, "scenario": variant, "endpoint": endpoint,
                        "included_candidates": n, "skipped_without_calibrant": skipped,
                        "assumed_positive_fraction": positive_fraction,
                        "conditional_detection_mean": round(sum(probabilities) / n, 6),
                        "expected_observed": round(positive_fraction * sum(probabilities), 6),
                        "conditional_detection_mean_if_nested_donor_recovery": round(nested_mean, 6),
                        "conditional_P_zero_independent_candidates": f"{p_zero:.8g}",
                        "n_distinct_HLA_opportunity_profiles": len(grouped),
                        "stress_P_zero_shared_outcome_within_profile": f"{p_zero_grouped:.8g}",
                    })
    yield "scenarios", rows


def mzid_reference(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rb") as source:
        root = ET.parse(source).getroot()
    peptides = {
        e.get("id"): (e.findtext("{*}PeptideSequence"), e.find("{*}Modification") is not None)
        for e in root.findall(".//{*}Peptide")
    }
    evidence = {e.get("id"): e.get("isDecoy") for e in root.findall(".//{*}PeptideEvidence")}
    spectra_data = {e.get("id"): e.get("location", "") for e in root.findall(".//{*}SpectraData")}
    for result in root.findall(".//{*}SpectrumIdentificationResult"):
        for item in result.findall("{*}SpectrumIdentificationItem"):
            if item.get("rank") != "1" or item.get("passThreshold", "").lower() != "true":
                continue
            ev = [evidence.get(e.get("peptideEvidence_ref")) for e in item.findall("{*}PeptideEvidenceRef")]
            if not ev or not any(v == "false" for v in ev) or any(v != "false" for v in ev):
                continue
            seq, modified = peptides.get(item.get("peptide_ref"), ("", False))
            if not PEPTIDE.fullmatch(seq or ""):
                continue
            yield {
                "author_mzid": path.name, "author_run": spectra_data.get(result.get("spectraData_ref"), ""),
                "author_spectrum": result.get("spectrumID", ""), "peptide": seq,
                "charge": item.get("chargeState", ""), "observed_mz": item.get("experimentalMassToCharge", ""),
                "calculated_mz": item.get("calculatedMassToCharge", ""),
                "modified": int(modified),
            }


def evaluate(base, output, psm_path, mapping, exact_hla):
    panel, restrictions, donor_hla = controls(base, exact_hla)
    write(output / "author_positive_panel.tsv", panel)
    if not psm_path:
        scenarios = conditional_scenarios([], restrictions, donor_hla,
                                          {"HP18101", "HP20289", "R360", "R361", "R369"})
        opportunities = []
        for kind, rows in scenarios:
            if kind == "opportunities":
                opportunities.extend(rows)
            else:
                write(output / "historical_donor_panel_jy_scenarios.tsv", rows)
        write(output / "historical_donor_panel_opportunities.tsv", opportunities)
        return
    accepted, examined = load_psms(psm_path)
    run_to_donor = {
        row["run"]: donor_id(row["donor"]) for row in read(mapping, "\t")
        if row["eligible"].lower() in {"1", "true", "yes"}
    }
    if exact_hla:
        run_to_donor = {run: donor for run, donor in run_to_donor.items() if donor_hla[donor]}
    if not run_to_donor or not set(run_to_donor) <= examined:
        raise ValueError(f"Missing matched runs in PSM export: {set(run_to_donor) - examined}")
    runs_by_donor = defaultdict(list)
    for run, donor in sorted(run_to_donor.items()):
        runs_by_donor[donor].append(run)
    donor_rows, run_rows = [], []
    for row in panel:
        donor, pep = row["donor"], row["peptide"]
        if donor not in runs_by_donor:
            continue
        ranked = [(run, accepted[run].get(pep, set())) for run in runs_by_donor[donor]]
        run_rows.extend({**row, "run": run, "distinct_scans": len(scans),
                         "scan_ids": ";".join(sorted(scans))} for run, scans in ranked)
        donor_rows.append({**row, "distinct_scans": sum(len(scans) for _, scans in ranked),
                           "n_runs_with_recovery": sum(bool(scans) for _, scans in ranked)})
    write(output / "reference_recovery_per_run.tsv", run_rows)
    write(output / "reference_recovery_per_donor.tsv", donor_rows)
    summaries = list(group_recovery(donor_rows))
    write(output / "tissue_recovery_summary.tsv", summaries)
    opportunities = []
    for kind, rows in conditional_scenarios(summaries, restrictions, donor_hla, set(runs_by_donor)):
        if kind == "opportunities":
            opportunities.extend(rows)
        else:
            write(output / "conditional_scenario_sensitivity.tsv", rows)
    write(output / "hla_opportunities.tsv", opportunities)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--psms", type=Path)
    parser.add_argument("--mapping", type=Path)
    parser.add_argument("--author-mzids", type=Path, nargs="*")
    parser.add_argument("--exact-hla", type=Path)
    args = parser.parse_args()
    if bool(args.psms) != bool(args.mapping):
        parser.error("--psms and --mapping are required together")
    args.out.mkdir(parents=True, exist_ok=True)
    evaluate(args.base, args.out, args.psms, args.mapping, args.exact_hla)
    if args.author_mzids:
        rows = [row for path in args.author_mzids for row in mzid_reference(path)]
        write(args.out / "source_author_positive_psms.tsv", rows)
        print("author positive identifications", len(rows), flush=True)


if __name__ == "__main__":
    main()

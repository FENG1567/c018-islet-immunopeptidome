"""Audit accepted matches from the fixed C018 positive/shuffle competition."""

import argparse
import math
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

from tissue_recovery import il_key, read, write


def psms(path):
    selected = {}
    for _, element in ET.iterparse(path, events=("end",)):
        if element.tag.split("}")[-1] != "spectrum_query":
            continue
        hit = element.find(".//{*}search_hit[@hit_rank='1']")
        if hit is not None:
            xcorr = hit.find("{*}search_score[@name='xcorr']")
            if xcorr is not None:
                peptide = hit.get("peptide")
                proteins = [hit.get("protein", "")]
                proteins.extend(e.get("protein", "") for e in hit.findall("{*}alternative_protein"))
                explicit_decoy = all(p.upper().startswith("DECOY_") for p in proteins)
                scan = element.get("start_scan")
                x = float(xcorr.get("value"))
                prior = selected.get(scan)
                if math.isfinite(x) and (prior is None or (x, explicit_decoy) > (prior["xcorr"], prior["is_decoy"])):
                    selected[scan] = {
                        "run": path.name.removesuffix(".pep.xml"),
                        "scan": scan, "peptide": peptide, "is_decoy": explicit_decoy,
                        "xcorr": x, "charge": element.get("assumed_charge"), "rank": 1,
                        "modified_peptide": (
                            hit.find("{*}modification_info").get("modified_peptide", "")
                            if hit.find("{*}modification_info") is not None else ""
                        ),
                    }
        element.clear()
    return list(selected.values())


def assign_q(rows, pseudocount=1):
    by_score = defaultdict(list)
    for row in rows:
        by_score[row["xcorr"]].append(row)
    decoys = targets = 0
    cumulative = {}
    for score in sorted(by_score, reverse=True):
        group = by_score[score]
        d = sum(row["is_decoy"] for row in group)
        decoys += d
        targets += len(group) - d
        cumulative[score] = min(1.0, (decoys + pseudocount) / max(1, targets))
    q = 1.0
    for score in sorted(by_score):
        q = min(q, cumulative[score])
        for row in by_score[score]:
            row["qvalue"] = q


def peptide_q(psm_rows, pseudocount=1):
    ranked_by_type = defaultdict(dict)
    for row in psm_rows:
        key = il_key(row["peptide"])
        group = ranked_by_type[row["is_decoy"]]
        if key not in group or row["xcorr"] > group[key]["xcorr"]:
            group[key] = row.copy()
    targets, decoys = ranked_by_type[False], ranked_by_type[True]
    mixed = targets.keys() & decoys.keys()
    ranked = [row for group in (targets, decoys) for key, row in group.items() if key not in mixed]
    assign_q(ranked, pseudocount)
    return {(il_key(r["peptide"]), r["is_decoy"]): r["qvalue"] for r in ranked}


def summarize(rows, path, panel, alpha, pseudocount):
    assign_q(rows, pseudocount)
    pq = peptide_q(rows, pseudocount)
    accepted = defaultdict(set)
    for row in rows:
        row["peptide_qvalue"] = pq.get((il_key(row["peptide"]), row["is_decoy"]), 1.0)
        if not row["is_decoy"] and row["qvalue"] <= alpha and row["peptide_qvalue"] <= alpha:
            accepted[il_key(row["peptide"])].add(row["scan"])
    pairs, aggregates = [], defaultdict(list)
    for row in panel:
        if not row["shuffle_id"]:
            continue
        target, shuffle = row["sequence"], row["shuffle_sequence"]
        pair = {"run": path.name.removesuffix(".pep.xml"), "tdc_pseudocount": pseudocount, **row,
                "target_scan_n": len(accepted[il_key(target)]), "shuffle_scan_n": len(accepted[il_key(shuffle)])}
        pairs.append(pair)
        labels = row["groups"].split(";")
        for label in labels:
            aggregates[label].append(pair)
    summary = []
    for group, evidence in sorted(aggregates.items()):
        for least in (1, 2):
            summary.append({
                "run": path.name.removesuffix(".pep.xml"), "tdc_pseudocount": pseudocount,
                "group": group, "scan_threshold": least,
                "panel_pairs": len(evidence), "accepted_targets": sum(r["target_scan_n"] >= least for r in evidence),
                "accepted_shuffles": sum(r["shuffle_scan_n"] >= least for r in evidence),
                "rank1_spectra": len(rows),
                "accepted_global_target_psms": sum(
                    not r["is_decoy"] and r["qvalue"] <= alpha and r["peptide_qvalue"] <= alpha for r in rows
                ),
            })
    return rows, pairs, summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--pepxml", type=Path, nargs="+", required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    panel = list(read(args.panel, "\t"))
    _psms, pairs, summaries = [], [], []
    for path in args.pepxml:
        raw_rows = psms(path)
        print(path.name, "spectra", len(raw_rows), "panel pairs", len(panel), flush=True)
        for pseudocount in (0, 1):
            psm_rows, pair_rows, summary_rows = summarize(
                [r.copy() for r in raw_rows], path, panel, 0.01, pseudocount)
            _psms.extend({**r, "tdc_pseudocount": pseudocount} for r in psm_rows)
            pairs.extend(pair_rows)
            summaries.extend(summary_rows)
    write(args.out / "matched_panel_psms.tsv", _psms)
    write(args.out / "matched_panel_pairs.tsv", pairs)
    write(args.out / "matched_panel_summary.tsv", summaries)


if __name__ == "__main__":
    main()

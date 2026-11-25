#!/usr/bin/env python3
"""Extract strict PRIDE-recipe rank hits for the four SIASAISMK spectra."""

import argparse
import csv
import xml.etree.ElementTree as ET
from pathlib import Path


TARGET = "SIASAISMK"
RUN_SCANS = {
    "SCP0004_IPP372_islet22_Slot2-1_1_355": {"57021"},
    "SCP0004_IPP372_islet22_Slot2-1_1_356": {"56976"},
    "SCP0004_IPP372_islet27_Slot2-2_1_358": {"58077"},
    "SCP0004_IPP372_islet27_Slot2-2_1_359": {"57601"},
}


def il_key(sequence):
    return sequence.replace("I", "J").replace("L", "J")


def classify(proteins):
    primary = proteins[0] if proteins else ""
    if primary.upper().startswith(("DECOY_", "REV_", "SHF_", "XXX_")):
        return "decoy_or_shuffle"
    if any(p.startswith("sp|") for p in proteins):
        return "reviewed_human_target"
    if any(p.startswith("REF_") for p in proteins):
        return "candidate_reference"
    return "other_target"


def write_tsv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) if rows else ["empty"]
    with path.open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--v2-base", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    root = args.v2_base / "tissue_calibration/pxd045211_pride_challenge"
    rows = []
    summaries = []
    for run, scans in RUN_SCANS.items():
        xml = root / run / "body" / f"{run}.candidate_index.pep.xml"
        tree = ET.parse(xml)
        for query in tree.findall(".//{*}spectrum_query"):
            scan = query.get("start_scan", "")
            if scan not in scans:
                continue
            hits = []
            for hit in query.findall(".//{*}search_hit"):
                score = hit.find("{*}search_score[@name='xcorr']")
                if score is None:
                    continue
                proteins = [hit.get("protein", "")]
                proteins.extend(node.get("protein", "") for node in hit.findall("{*}alternative_protein"))
                peptide = hit.get("peptide", "")
                row = {
                    "search": "strict_PRIDE_recipe",
                    "source_run": run,
                    "comet_scan": scan,
                    "hit_rank": int(hit.get("hit_rank", "0")),
                    "peptide": peptide,
                    "protein": hit.get("protein", ""),
                    "all_protein_mappings": ";".join(proteins),
                    "assignment_class": classify(proteins),
                    "xcorr": float(score.get("value", "nan")),
                    "I_L_equivalent_to_candidate": il_key(peptide) == il_key(TARGET),
                }
                hits.append(row)
            hits.sort(key=lambda r: r["hit_rank"])
            rows.extend(hits)
            candidate = next((r for r in hits if r["peptide"] == TARGET), None)
            human = [r for r in hits if r["assignment_class"] == "reviewed_human_target"
                     and il_key(r["peptide"]) != il_key(TARGET)]
            best_human = max(human, key=lambda r: r["xcorr"], default=None)
            summaries.append({
                "source_run": run,
                "comet_scan": scan,
                "top_output_lines": len(hits),
                "candidate_rank": candidate["hit_rank"] if candidate else "",
                "candidate_xcorr": candidate["xcorr"] if candidate else "",
                "best_reviewed_human_non_IL_peptide_within_top10": best_human["peptide"] if best_human else "",
                "best_reviewed_human_protein_within_top10": best_human["protein"] if best_human else "",
                "best_reviewed_human_rank_within_top10": best_human["hit_rank"] if best_human else "",
                "best_reviewed_human_xcorr_within_top10": best_human["xcorr"] if best_human else "",
                "candidate_minus_best_human_xcorr_within_top10": (
                    candidate["xcorr"] - best_human["xcorr"] if candidate and best_human else ""
                ),
                "scope_note": "rank competition only among stored top-output lines; not a candidate-subset FDR",
            })
        if scans:
            observed = {r["comet_scan"] for r in summaries if r["source_run"] == run}
            missing = scans - observed
            if missing:
                raise ValueError(f"{run}: missing SIASAISMK scans {sorted(missing)}")
        print(f"{run}: extracted {sum(r['source_run'] == run for r in summaries)} candidate spectrum", flush=True)
    write_tsv(args.out / "SIASAISMK_pride_top10_hits.tsv", rows)
    write_tsv(args.out / "SIASAISMK_pride_top10_human_competition.tsv", summaries)
    print(f"COMPLETE: {len(rows)} ranked hits; {len(summaries)} spectrum summaries", flush=True)


if __name__ == "__main__":
    main()

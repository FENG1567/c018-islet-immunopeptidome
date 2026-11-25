"""Stress-check repeat signals using the same matched null sequences across acquisitions."""

import argparse
import json
from collections import defaultdict
from pathlib import Path

from analyze_external_islets import candidate_map, eligible_inputs
from score_panel import assign_q, peptide_q, psms
from tissue_recovery import il_key, read, write


def accepted(path):
    rows = psms(path)
    assign_q(rows, 1)
    q_peptide = peptide_q(rows)
    counts = defaultdict(set)
    for row in rows:
        key = il_key(row["peptide"])
        if not row["is_decoy"] and not row["modified_peptide"] and row["qvalue"] <= 0.01:
            if q_peptide.get((key, False), 1.0) <= 0.01:
                counts[key].add(row["scan"])
    return {sequence: len(scans) for sequence, scans in counts.items()}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--base", type=Path, required=True)
    p.add_argument("--acquisition", type=Path, required=True)
    p.add_argument("--expanded", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    groups = eligible_inputs(args.acquisition)
    prep_by_run = {Path(name).stem: prep for prep, inputs in groups.items() for name, _ in inputs}
    by_run, pair_sets = {}, []
    for name in prep_by_run:
        run = args.expanded / name
        pair_sets.append({(r["group"], r["target"], r["null"])
                          for r in read(run / "expanded_entrapment_map.tsv", "\t")})
        cache = run / "unmodified_sequence_counts.json"
        by_run[name] = json.loads(cache.read_text()) if cache.is_file() else accepted(run / "body" / (name + ".pep.xml"))
        print(name, "passing target-like sequences", len(by_run[name]), flush=True)
    common_pairs = set.intersection(*pair_sets)
    restrictions = candidate_map(args.base)

    def observations(sequence):
        seq = il_key(sequence)
        by_prep = defaultdict(int)
        scans = defaultdict(int)
        for run, prep in prep_by_run.items():
            n = by_run[run].get(seq, 0)
            by_prep[prep] += int(n >= 1)
            scans[prep] += n
        return {
            "both_preps_at_least_one": int(len([n for n in by_prep.values() if n >= 1]) == 2),
            "both_preps_at_least_two_scans": int(len([n for n in scans.values() if n >= 2]) == 2),
            "all_four_technical_acquisitions": int(sum(by_prep.values()) == 4),
        }

    pair_rows, groups_rows = [], defaultdict(list)
    for group, target, null in sorted(common_pairs):
        key = ("NC2025_aeTSA" if group == "nc_aetsa" else "PDAC2025", target)
        families = {a.split(":")[0] for a in restrictions.get(key, set())}
        for sequence_type, sequence in (("target", target), ("matched_null", null)):
            row = {"list": key[0], "parent_target": target, "sequence_type": sequence_type,
                   "query_sequence": sequence, "family_possible": int(bool(families & {"A*02", "A*03"})),
                   **observations(sequence)}
            pair_rows.append(row)
            groups_rows[(key[0], "all", sequence_type)].append(row)
            if row["family_possible"]:
                groups_rows[(key[0], "A02_A03_family", sequence_type)].append(row)
    summary = []
    for (group, stratum, label), rows in sorted(groups_rows.items()):
        rows = {r["query_sequence"]: r for r in rows}
        summary.append({
            "list": group, "stratum": stratum, "type": label, "sequence_n": len(rows),
            **{field: sum(r[field] for r in rows.values()) for field in (
                "both_preps_at_least_one", "both_preps_at_least_two_scans", "all_four_technical_acquisitions")},
        })
    write(args.out / "joint_null_summary.tsv", summary)
    write(args.out / "joint_null_sequence_evidence.tsv", pair_rows)


if __name__ == "__main__":
    main()

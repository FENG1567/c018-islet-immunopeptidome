"""Join the two fixed search recipes by candidate, preparation and acquisition."""

import argparse
from collections import defaultdict
from pathlib import Path

from analyze_external_islets import candidate_map, eligible_inputs
from tissue_recovery import read, write


def collect(directory):
    by_candidate = defaultdict(list)
    seen_runs = set()
    for file in sorted(directory.glob("*/external_candidate_calls.tsv")):
        run = file.parent.name
        seen_runs.add(run)
        for row in read(file, "\t"):
            by_candidate[(row["group"], row["peptide"])].append(row)
    return by_candidate, seen_runs


def selected(rows, prep):
    passed = {r["source"]: int(r["threshold_pass_exact_unmodified_n"]) for r in rows
              if int(r["threshold_pass_exact_unmodified_n"]) > 0}
    counts = defaultdict(int)
    for source in passed:
        counts[prep[source]] += 1
    return passed, set(counts), {p for p, n in counts.items() if n == 2}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--base", type=Path, required=True)
    p.add_argument("--acquisition", type=Path, required=True)
    p.add_argument("--replay", type=Path, required=True)
    p.add_argument("--pride-recipe", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    prep = {Path(name).stem: preparation for preparation, files in eligible_inputs(args.acquisition).items()
            for name, _ in files}
    old, replay_runs = collect(args.replay)
    new, pride_runs = collect(args.pride_recipe)
    if set(prep) != replay_runs or set(prep) != pride_runs:
        raise ValueError(f"Runs missing: original {set(prep)-replay_runs}, PRIDE recipe {set(prep)-pride_runs}")
    calls = []
    for key, restrictions in sorted(candidate_map(args.base).items()):
        if len(old[key]) != len(prep) or len(new[key]) != len(prep):
            raise ValueError(f"Expected one candidate row per completed acquisition for {key}")
        runs_a, prep_a, pair_a = selected(old[key], prep)
        runs_b, prep_b, pair_b = selected(new[key], prep)
        compatible = sorted({a.split(":")[0] for a in restrictions} & {"A*02", "A*03"})
        calls.append({
            "list": key[0], "peptide": key[1], "candidate_hla_alleles": ";".join(sorted(restrictions)),
            "possible_A02_A03_family": ";".join(compatible), "exact_donor_HLA": "unknown",
            "legacy_recipe_exact_unmodified_psm_by_run": ";".join(f"{run}={k}" for run, k in sorted(runs_a.items())),
            "legacy_recipe_n_preparations": len(prep_a),
            "legacy_recipe_preparations_two_acquisitions": ";".join(sorted(pair_a)),
            "pride_recipe_exact_unmodified_psm_by_run": ";".join(f"{run}={k}" for run, k in sorted(runs_b.items())),
            "pride_recipe_n_preparations": len(prep_b),
            "pride_recipe_preparations_two_acquisitions": ";".join(sorted(pair_b)),
            "recurring_in_both_recipes_and_both_technical_acquisitions": ";".join(sorted(pair_a & pair_b)),
            "presentation_status": "unresolved candidate spectra" if runs_a or runs_b else "no qualifying candidate spectrum",
        })
    write(args.out / "external_candidate_all_summary.tsv", calls)
    write(args.out / "external_unresolved_candidate_signals.tsv",
          (r for r in calls if r["presentation_status"] == "unresolved candidate spectra"))
    methods = ("legacy_recipe", "pride_recipe")
    summaries = []
    for group in sorted({r["list"] for r in calls}):
        for stratum in ("all", "family_supported"):
            rows = [r for r in calls if r["list"] == group
                    and (stratum == "all" or r["possible_A02_A03_family"])]
            for method in methods:
                tested = len(rows)
                any_count = sum(r[f"{method}_n_preparations"] > 0 for r in rows)
                both_count = sum(r[f"{method}_n_preparations"] == 2 for r in rows)
                repeat_count = sum(bool(r[f"{method}_preparations_two_acquisitions"]) for r in rows)
                summaries.append({
                    "list": group, "stratum": stratum, "recipe": method, "candidate_sequences": tested,
                    "with_screen_signal_any_preparation": any_count,
                    "with_screen_signal_both_preparations": both_count,
                    "with_screen_signal_both_acquisitions_in_any_preparation": repeat_count,
                })
    write(args.out / "external_detection_summary.tsv", summaries)


if __name__ == "__main__":
    main()

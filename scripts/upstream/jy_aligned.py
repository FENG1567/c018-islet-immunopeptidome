"""Re-evaluate archived JY dilution searches at the same one- and two-scan thresholds."""

import argparse
import math
from collections import defaultdict
from pathlib import Path

from score_panel import assign_q, peptide_q, psms
from tissue_recovery import il_key, interval, probability, read, write


def source_jy(path):
    peptides = set()
    label = None
    parts = []
    for line in [*path.read_text().splitlines(), ">"]:
        if line.startswith(">"):
            if label and label.startswith("jy_"):
                peptides.add("".join(parts))
            label, parts = line[1:], []
        else:
            parts.append(line.strip())
    return peptides


def accepted(path, mode, jy):
    rows = psms(path)
    counts = defaultdict(set)
    if mode == "legacy_D_over_T_rank1":
        target, decoy = 0, 0
        for row in sorted(rows, key=lambda r: -r["xcorr"]):
            if row["is_decoy"]:
                decoy += 1
            else:
                target += 1
                if decoy / target <= 0.01 and row["peptide"] in jy:
                    counts[row["peptide"]].add(row["scan"])
    else:
        assign_q(rows, 1)
        per_peptide = peptide_q(rows)
        jy_il = {il_key(p) for p in jy}
        for row in rows:
            pep = il_key(row["peptide"])
            if (not row["is_decoy"] and pep in jy_il and row["qvalue"] <= 0.01
                    and per_peptide.get((pep, False), 1) <= 0.01):
                counts[pep].add(row["scan"])
    return counts, len(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--opportunities", type=Path, required=True)
    parser.add_argument("--original-fasta", type=Path, required=True)
    parser.add_argument("--pepxml", type=Path, nargs=3, required=True, metavar=("NEAT", "DILUTION4", "DILUTION16"))
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    recovery, candidate_sensitivity = [], []
    opportunities = defaultdict(list)
    for row in read(args.opportunities, "\t"):
        opportunities[row["cohort"]].append(row["compatible_donors"])
    jy = source_jy(args.original_fasta)
    if not jy:
        raise ValueError("No archived JY-only peptide labels found")
    for mode in ("legacy_D_over_T_rank1", "aligned_psm_and_peptide_q_plus1"):
        counts = [accepted(path, mode, jy) for path in args.pepxml]
        reference = counts[0][0]
        for dilution, file, (tested, n_rank1) in zip((1, 4, 16), args.pepxml, counts):
            for ref_need in (1, 2):
                reference_peptides = {pep for pep in reference if len(reference[pep]) >= ref_need}
                for test_need in (1, 2):
                    k = sum(len(tested.get(pep, ())) >= test_need for pep in reference_peptides)
                    n = len(reference_peptides)
                    low, high = interval(k, n)
                    recovery.append({"mode": mode, "dilution_factor": dilution, "file": file.name,
                                     "reference_scans_min": ref_need, "test_scans_min": test_need,
                                     "tested_rank1_spectra": n_rank1, "reference_n": n, "recovered_n": k,
                                     "recovery": round(k / n, 7), "descriptive_wilson_low": low,
                                     "descriptive_wilson_high": high})
                    if mode.startswith("aligned") and ref_need == 1 and test_need == 2:
                        for group, profiles in opportunities.items():
                            for end, min_donor in (("any", 1), ("both", 2)):
                                profile_rates = {}
                                values = []
                                for profile in profiles:
                                    m = len(profile.split(";")) if profile else 0
                                    if m >= min_donor:
                                        rate = probability([k / n] * m)[min_donor - 1]
                                        values.append(rate)
                                        profile_rates[profile] = rate
                                candidate_sensitivity.append({
                                    "cohort": group, "dilution_factor": dilution, "endpoint": end,
                                    "candidates": len(values), "empirical_JY_two_scan_recovery": round(k / n, 7),
                                    "conditional_probability_mean_under_JY_transfer": round(sum(values) / len(values), 7),
                                    "P_zero_at_20pct_if_candidates_independent": round(
                                        math.prod(1 - 0.20 * p for p in values), 10),
                                    "P_zero_at_20pct_shared_outcome_profile_stress": round(
                                        math.prod(1 - 0.20 * p for p in profile_rates.values()), 10),
                                })
    write(args.out / "jy_aligned_recovery.tsv", recovery)
    write(args.out / "jy_aligned_conditional_sensitivity.tsv", candidate_sensitivity)


if __name__ == "__main__":
    main()

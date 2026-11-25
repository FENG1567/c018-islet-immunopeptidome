"""Evaluate author-positive ligands against a different acquisition of one preparation."""

import argparse
import math
import re
from collections import defaultdict
from pathlib import Path

from tissue_recovery import interval, read, write


INDEX = re.compile(r"\bindex=(\d+)\b")


def signal_for_ids(mgf_path, requested):
    index, mz, charge, tic = None, None, None, 0.0
    signatures = {}
    with mgf_path.open("rt", encoding="ascii") as stream:
        for line in stream:
            if line.startswith("BEGIN IONS"):
                index, mz, charge, tic = None, None, None, 0.0
            elif line.startswith("TITLE="):
                match = INDEX.search(line)
                index = int(match.group(1)) if match else None
            elif index in requested:
                if line.startswith("PEPMASS="):
                    mz = float(line.split("=", 1)[1].split()[0])
                elif line.startswith("CHARGE="):
                    matched_charge = re.search(r"\d+", line)
                    charge = int(matched_charge.group()) if matched_charge else None
                elif line.startswith("END IONS"):
                    if index in signatures:
                        raise ValueError(f"Duplicate spectrum index {index} in {mgf_path}")
                    signatures[index] = (mz, charge, tic)
                elif line and line[0].isdigit():
                    _, intensity, *_ = line.split()
                    tic += float(intensity)
    return signatures


def source_names(row):
    return row["author_run"].replace("\\", "/").split("/")[-1]


def paired_benchmark(reference, ref_mgf, test_mgf, psms_path):
    by_source = defaultdict(list)
    for row in read(reference, "\t"):
        by_source[source_names(row)].append(row)
    source_rows = [row for row in by_source[ref_mgf.name]
                   if row["modified"] == "0" and INDEX.search(row["author_spectrum"])]
    if not source_rows or not by_source[test_mgf.name]:
        raise ValueError("Could not resolve both author-result data file references")
    checks = signal_for_ids(ref_mgf, {int(INDEX.search(row["author_spectrum"]).group(1)) for row in source_rows})
    strongest = {}
    verified_rows = 0
    for row in source_rows:
        index = int(INDEX.search(row["author_spectrum"]).group(1))
        observation = checks.get(index)
        if observation is None:
            continue
        mz, charge, tic = observation
        if not mz or charge != int(row["charge"]) or abs(float(row["observed_mz"]) - mz) > 50e-6 * mz:
            continue
        verified_rows += 1
        pep = row["peptide"]
        if pep not in strongest or tic > strongest[pep][0]:
            strongest[pep] = (tic, index, mz, charge)
    if not strongest or verified_rows < 0.95 * len(source_rows):
        raise ValueError(f"Index/precursor alignment below 95%: {verified_rows}/{len(source_rows)}")
    by_author_test = {row["peptide"] for row in by_source[test_mgf.name] if row["modified"] == "0"}
    accepted = defaultdict(set)
    for row in read(psms_path, "\t"):
        if row["run"] != test_mgf.stem:
            continue
        if (row["tdc_pseudocount"] == "1" and row["is_decoy"] == "False"
                and not row["modified_peptide"] and float(row["qvalue"]) <= 0.01
                and float(row["peptide_qvalue"]) <= 0.01):
            accepted[row["peptide"]].add(row["scan"])
    values = sorted(row[0] for row in strongest.values())
    edges = [values[math.ceil(len(values) * i / 4) - 1] for i in (1, 2, 3)]
    rows = []
    for pep, (tic, index, mz, charge) in sorted(strongest.items()):
        scans = accepted.get(pep, set())
        rows.append({
            "source_run": ref_mgf.stem, "test_run": test_mgf.stem, "peptide": pep,
            "reference_scan_index": index, "reference_mz": mz, "reference_charge": charge,
            "reference_ms2_fragment_tic": round(tic, 4), "source_tic_quartile": 1 + sum(tic > e for e in edges),
            "author_test_observation": int(pep in by_author_test), "c018_test_distinct_scans": len(scans),
            "c018_test_scan_ids": ";".join(sorted(scans)),
        })
    return rows, verified_rows, len(source_rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--ref-mgf", type=Path, required=True)
    parser.add_argument("--test-mgf", type=Path, required=True)
    parser.add_argument("--psms", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    rows, verified_rows, total_rows = paired_benchmark(
        args.reference, args.ref_mgf, args.test_mgf, args.psms)
    write(args.out / "paired_recovery_by_ligand.tsv", rows)
    summaries = []
    for quartile in (0, 1, 2, 3, 4):
        subset = [row for row in rows if not quartile or row["source_tic_quartile"] == quartile]
        n = len(subset)
        if not n:
            continue
        sources = {
            "author": sum(row["author_test_observation"] for row in subset),
            "c018_one_scan": sum(row["c018_test_distinct_scans"] >= 1 for row in subset),
            "c018_two_scans": sum(row["c018_test_distinct_scans"] >= 2 for row in subset),
        }
        for method, k in sources.items():
            low, high = interval(k, n)
            summaries.append({
                "source_run": args.ref_mgf.stem, "test_run": args.test_mgf.stem,
                "source_tic_quartile": quartile or "all", "method": method,
                "reference_positives": n, "recovered": k, "recovery": round(k / n, 6),
                "descriptive_wilson_low": low, "descriptive_wilson_high": high,
            })
    write(args.out / "paired_recovery_summary.tsv", summaries)
    print(f"author index / m/z / charge verification {verified_rows}/{total_rows}", flush=True)
    print(f"unique unmodified reference ligands {len(rows)}", flush=True)


if __name__ == "__main__":
    main()

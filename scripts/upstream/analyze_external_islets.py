"""Run the C018 external-islet search after each verified public file arrives."""

import argparse
import json
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

from matched_panel_search import create_panel, create_parameters, digest, search_one
from normalize_pasef_mgf import normalize
from tissue_recovery import alleles, il_key, pdac_allele_union, read, write


def eligible_inputs(acquisition):
    with (acquisition / "metadata/PXD045211.files.json").open() as source:
        content = json.load(source)
    if isinstance(content, dict):
        content = next(iter(content["_embedded"].values()))
    sizes = {row["fileName"]: int(row["fileSizeBytes"]) for row in content if row["fileName"].endswith(".mgf")}
    groups = defaultdict(list)
    metadata = acquisition / "metadata/PXD045211_community_annotated.sdrf.tsv"
    for row in read(metadata, "\t"):
        preparation = row["source name"]
        mgf = str(Path(row["comment[data file]"]).with_suffix(".mgf"))
        groups[preparation].append((mgf, sizes[mgf]))
    if set(sizes) != {name for group in groups.values() for name, _ in group}:
        raise ValueError("PRIDE MGF manifest and SDRF differ")
    return dict(groups)


def available(path, size):
    return (
        path.is_file() and path.stat().st_size == size
        and not path.with_name(path.name + ".aria2").is_file()
    )


def candidate_map(base):
    result = defaultdict(set)
    for row in read(base / "work/candidates/nc2025_all_raw.csv"):
        if "aetsa" in row["TA_type"].lower():
            result[("NC2025_aeTSA", row["peptide"].upper())].update(alleles(row["HLA"]))
    for peptide, restrictions in pdac_allele_union(base).items():
        result[("PDAC2025", peptide)].update(restrictions)
    return result


def candidate_calls(base, dataset, psms, crosswalk, output):
    mapping = {r["comet_scan"]: r for r in read(crosswalk, "\t")}
    rows = list(read(psms, "\t"))
    retained = defaultdict(list)
    for row in rows:
        if row["tdc_pseudocount"] != "1" or row["is_decoy"] != "False":
            continue
        retained[il_key(row["peptide"])].append(row)
    all_calls, evidence = [], []
    for (group, pep), restrictions in sorted(candidate_map(base).items()):
        compatibility = bool({s.split(":")[0] for s in restrictions} & {"A*02", "A*03"})
        matches = retained.get(il_key(pep), [])
        accepted = []
        for match in matches:
            key = match["scan"]
            native = mapping.get(key)
            if not native:
                raise ValueError(f"PSM scan {key} has no original spectrum key")
            q_psm = float(match["qvalue"])
            q_pep = float(match["peptide_qvalue"])
            passed = q_psm <= 0.01 and q_pep <= 0.01
            if passed:
                accepted.append((match, native))
            evidence.append({
                "source": dataset, "group": group, "peptide": pep,
                "psm_sequence": match["peptide"], "comet_scan": key, "author_index": native["author_index"],
                "native_scan": native["native_scans"], "alignment": "exact" if pep == match["peptide"] else "I/L-only",
                "xcorr": match["xcorr"], "psm_q": q_psm, "peptide_q": q_pep,
                "modified_peptide": match["modified_peptide"], "threshold_pass": int(passed),
            })
        all_calls.append({
            "source": dataset, "group": group, "peptide": pep,
            "hla": ";".join(sorted(restrictions)), "known_2_field_compatibility": "unknown",
            "possible_A02_A03_compatibility": int(compatibility),
            "rank1_psm_n": len(matches), "threshold_pass_psm_n": len(accepted),
            "threshold_pass_native_frame_n": len({m["native_scans"] for _, m in accepted}),
            "threshold_pass_exact_unmodified_n": sum(
                r["peptide"] == pep and not r["modified_peptide"] for r, _ in accepted),
        })
    write(output / "external_candidate_calls.tsv", all_calls)
    write(output / "external_candidate_psm_evidence.tsv", evidence)


def run_script(name, *args):
    program = Path(__file__).with_name(name)
    subprocess.run([sys.executable, "-u", str(program), *map(str, args)], check=True)


def analyze(base, sample, other, mzid, output, normalized_directory, pride_recipe, comet, params_template):
    outdir = output / sample.stem
    outdir.mkdir(parents=True, exist_ok=True)
    normalized = normalized_directory / sample.name
    if not normalized.is_file() or not normalized.with_suffix(".scan_map.tsv").is_file():
        normalize(sample, normalized)
    if not (outdir / "reference_targets_shuffles.fasta").is_file():
        create_panel(base, outdir, [mzid], other.name)
        create_parameters(params_template, outdir, pride_recipe)
    search_one(base, outdir, normalized, comet)
    run_script("score_panel.py", "--panel", outdir / "reference_targets_shuffles.tsv",
               "--out", outdir, "--pepxml", outdir / "body" / (sample.stem + ".pep.xml"))
    candidate_calls(base, sample.stem, outdir / "matched_panel_psms.tsv",
                    normalized.with_suffix(".scan_map.tsv"), outdir)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--acquisition", type=Path, required=True)
    parser.add_argument("--source-psms", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--normalized-from", type=Path)
    parser.add_argument("--comet", type=Path, required=True, help="Comet executable")
    parser.add_argument("--params-template", type=Path,
                        help="Archived raw-validation or PXD045211 PRIDE-recipe parameter snapshot")
    parser.add_argument("--pride-recipe", action="store_true")
    parser.add_argument("--watch-seconds", type=int, default=0)
    args = parser.parse_args()
    base, acquisition, output = (p.resolve(strict=True) for p in (args.base, args.acquisition, args.out.parent))
    output = output / args.out.name
    output.mkdir(parents=True, exist_ok=True)
    normalized_directory = (args.normalized_from or output / "normalized").resolve()
    if args.params_template:
        params_template = args.params_template.resolve(strict=True)
    elif args.pride_recipe:
        params_template = Path(__file__).resolve().parents[2] / "config/comet.pxd045211.pride_recipe.params"
    else:
        params_template = Path(__file__).resolve().parents[2] / "config/comet.pxd045211.raw_validation.params"
    grouped = eligible_inputs(acquisition)
    end_time = time.monotonic() + args.watch_seconds
    completed, paired, pair_failures = set(), set(), set()
    while True:
        for prep, inputs in grouped.items():
            if len(inputs) != 2:
                raise ValueError(f"{prep}: two acquisitions expected")
            ref = acquisition / "metadata" / f"peptides_1_1_0_{prep}.mzid.gz"
            for (name, size), (alternate, alt_size) in (inputs, inputs[::-1]):
                sample = acquisition / "raw" / name
                other = acquisition / "raw" / alternate
                if name not in completed and available(sample, size):
                    print("verified input bytes", size, "SHA256", digest(sample), name, flush=True)
                    if args.normalized_from and not (
                        (normalized_directory / name).is_file()
                        and (normalized_directory / name).with_suffix(".scan_map.tsv").is_file()
                    ):
                        continue
                    analyze(base, sample, other, ref, output, normalized_directory, args.pride_recipe,
                            args.comet.resolve(strict=True), params_template)
                    completed.add(name)
                if name in completed and name not in (paired | pair_failures) and available(other, alt_size):
                    outdir = output / sample.stem
                    try:
                        run_script("paired_recovery.py", "--reference", args.source_psms, "--ref-mgf", other,
                                   "--test-mgf", sample, "--psms", outdir / "matched_panel_psms.tsv", "--out", outdir)
                    except subprocess.CalledProcessError:
                        pair_failures.add(name)
                    else:
                        paired.add(name)
        if len(paired | pair_failures) == sum(len(g) for g in grouped.values()) or time.monotonic() >= end_time:
            break
        time.sleep(30)
    print("candidate_screen_completed", sorted(completed), "reciprocal_recovery_completed", sorted(paired),
          "reciprocal_recovery_failed", sorted(pair_failures), flush=True)


if __name__ == "__main__":
    main()

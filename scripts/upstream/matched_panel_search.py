"""Build and search a fixed positive/shuffled panel for the C018 sensitivity audit."""

import argparse
import hashlib
import random
import re
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

from tissue_recovery import PEPTIDE, donor_id, il_key, mzid_reference, read, write


SEED = "C018-v2-20261122"


def shuffled_sequence(peptide, label, blocked):
    generator = random.Random(SEED + ":" + label + ":" + peptide)
    middle = list(peptide[2:-1])
    for _ in range(200):
        generator.shuffle(middle)
        decoy = peptide[:2] + "".join(middle) + peptide[-1]
        if il_key(decoy) not in blocked:
            blocked.add(il_key(decoy))
            return decoy
    return ""


def sequence_groups(base, extra_results, author_run, excluded_control_donor=None):
    by_sequence = defaultdict(set)
    for row in read(base / "work/msv/islet_ligands.csv"):
        if excluded_control_donor and donor_id(row["Donor"]) == donor_id(excluded_control_donor):
            continue
        by_sequence[row["Peptide"].strip().upper()].add("stern_positive")
    for row in read(base / "work/candidates/nc2025_all_raw.csv"):
        by_sequence[row["peptide"].strip().upper()].add(
            "nc_aetsa" if "aetsa" in row["TA_type"].strip().lower() else "nc_other"
        )
    for row in read(base / "out/step6/TABLE_S_pdac_grading.tsv", "\t"):
        by_sequence[row["peptide"].strip().upper()].add("pdac_candidate")
    for path in extra_results:
        for row in mzid_reference(path):
            if author_run and not row["author_run"].endswith("/" + author_run):
                continue
            by_sequence[row["peptide"]].add("external_positive")
    if author_run and not any("external_positive" in tags for tags in by_sequence.values()):
        raise ValueError(f"No accepted source-author control for {author_run}")
    return {seq: labels for seq, labels in by_sequence.items() if PEPTIDE.fullmatch(seq)}


def create_panel(base, output, extra_results, author_run, excluded_control_donor=None):
    by_sequence = sequence_groups(base, extra_results, author_run, excluded_control_donor)
    blocked = {il_key(pep) for pep in by_sequence}
    by_il = defaultdict(set)
    for sequence in by_sequence:
        by_il[il_key(sequence)].add(sequence)
    rows = []
    with (output / "reference_targets_shuffles.fasta").open("w", encoding="ascii") as fasta:
        for idx, (sequence, labels) in enumerate(sorted(by_sequence.items()), 1):
            aliases = by_il[il_key(sequence)]
            family_labels = set().union(*(by_sequence[p] for p in aliases))
            groups = ";".join(sorted(family_labels))
            shuffled = shuffled_sequence(sequence, groups, blocked) if len(aliases) == 1 else ""
            target_id = f"REF_{idx:05d}"
            shuffled_id = f"SHF_{idx:05d}" if shuffled else ""
            fasta.write(f">{target_id}\n{sequence}\n")
            if shuffled:
                fasta.write(f">{shuffled_id}\n{shuffled}\n")
            rows.append({
                "target_id": target_id, "sequence": sequence, "groups": groups,
                "il_equivalent_count": len(aliases), "shuffle_id": shuffled_id,
                "shuffle_sequence": shuffled,
            })
    write(output / "reference_targets_shuffles.tsv", rows)
    print("panel counts", dict(Counter(g for row in rows for g in row["groups"].split(";"))),
          "target entries", len(rows), "shuffle entries", sum(bool(r["shuffle_id"]) for r in rows),
          flush=True)


def create_parameters(source, output, pride_recipe=False):
    raw = source.read_text(encoding="utf-8")
    settings = {
        "database_name": "../reference_targets_shuffles.fasta",
        "num_threads": "8",
        "output_pepxmlfile": "1",
        "output_txtfile": "0",
        "num_output_lines": "2",
    }
    if pride_recipe:
        settings.update({
            "peptide_mass_tolerance_upper": "10.0",
            "peptide_mass_tolerance_lower": "-10.0",
            "fragment_bin_tol": "0.03",
            "add_C_cysteine": "0.0",
            "variable_mod01": "15.994915 M 0 3 -1 0 0 0.0",
            "variable_mod02": "0.984016 NQ 0 2 -1 0 0 0.0",
            "variable_mod03": "42.010565 ^ 0 1 -1 0 0 0.0",
        })
    for parameter, value in settings.items():
        raw, count = re.subn(fr"^{parameter}[ \t]*=.*$", f"{parameter} = {value}", raw, flags=re.M)
        if count != 1:
            raise ValueError(f"{parameter} occurred {count} times")
    (output / "panel_search.params").write_text(raw, encoding="utf-8")


def digest(path):
    result = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            result.update(chunk)
    return result.hexdigest()


def search_one(base, output, source, comet_path):
    results = output / "body"
    results.mkdir(exist_ok=True)
    source = source.resolve(strict=True)
    link = results / source.name
    if not link.is_symlink():
        link.symlink_to(source)
    if link.resolve() != source:
        raise ValueError(f"existing input link points elsewhere: {link}")
    result = link.with_suffix(".pep.xml")
    completion = results / (source.stem + ".complete")
    config_hashes = digest(output / "panel_search.params") + " " + digest(output / "reference_targets_shuffles.fasta")
    if result.is_file() and result.stat().st_size and completion.is_file():
        stamp = completion.read_text().strip()
        if stamp.endswith(" " + config_hashes) and source.stat().st_mtime_ns <= completion.stat().st_mtime_ns:
            print("reuse verified input under same parameters", result, flush=True)
            return
    input_digest = digest(source)
    fingerprint = input_digest + " " + config_hashes
    if result.is_file() and result.stat().st_size and completion.is_file() and completion.read_text().strip() == fingerprint:
        print("reuse", result, flush=True)
        return
    print("start", source, source.stat().st_size, input_digest, flush=True)
    command = [str(comet_path), "-P" + str(output / "panel_search.params"), str(link)]
    with (results / (source.stem + ".comet.log")).open("w") as log:
        subprocess.run(command, check=True, cwd=results, stdout=log, stderr=subprocess.STDOUT)
    if not result.is_file() or not result.stat().st_size:
        raise RuntimeError(f"Comet did not produce {result}")
    completion.write_text(fingerprint + "\n", encoding="ascii")
    print("done", result.stat().st_size, result, flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--author-mzids", type=Path, nargs="*", default=[])
    parser.add_argument("--author-reference-run")
    parser.add_argument("--exclude-control-donor")
    parser.add_argument("--mgf", type=Path, nargs="*", default=[])
    parser.add_argument("--comet", type=Path, required=True, help="Comet executable used for the archived Comet 2026.03 search")
    parser.add_argument("--params-template", type=Path,
                        default=Path(__file__).resolve().parents[2] / "config/comet.local_competition.v2.params",
                        help="Archived parameter template; database_name is replaced for this output panel")
    parser.add_argument("--pride-recipe", action="store_true",
                        help="Apply the archived PXD045211 metadata adjustments to a base template")
    parser.add_argument("--reuse-panel", action="store_true")
    args = parser.parse_args()
    output, base = args.out.resolve(), args.base.resolve()
    output.mkdir(parents=True, exist_ok=True)
    if args.reuse_panel:
        if args.author_mzids or args.author_reference_run or args.exclude_control_donor:
            parser.error("--reuse-panel cannot change the reference")
        if not (output / "panel_search.params").is_file():
            parser.error("The requested panel does not exist")
    else:
        if (output / "body").exists():
            raise ValueError("Use a new output directory when changing the fixed search panel")
        create_panel(base, output, args.author_mzids, args.author_reference_run, args.exclude_control_donor)
        create_parameters(args.params_template.resolve(strict=True), output, pride_recipe=args.pride_recipe)
    print("panel FASTA SHA256", digest(output / "reference_targets_shuffles.fasta"), flush=True)
    for mgf in args.mgf:
        search_one(base, output, mgf, args.comet.resolve(strict=True))


if __name__ == "__main__":
    main()

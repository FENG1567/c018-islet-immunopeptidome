#!/usr/bin/env python3
"""Rebuild the supported C018 v4 figures from a verified supplement ZIP."""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import io
import os
from pathlib import Path, PurePosixPath
import sys
import tempfile
import zipfile


ROOT = Path(__file__).resolve().parents[1]
GENERATORS = ROOT / "vendor" / "figures"
REQUIRED = {
    "2": (
        "S2/archived_ligand_HLA_opportunities.tsv",
        "S2/aetsa_tissue_resolved_endpoints.tsv",
        "S2/aetsa_exact_typed_endpoint_recomputed.tsv",
        "S2/aetsa_endpoint_recomputed.tsv",
        "S2/lung_endpoints.tsv",
        "S2/ILMILQPQL_donor_HLA_S1_S2.tsv",
        "S2/ILMILQPQL_mapping_scenarios.tsv",
    ),
    "1": (
        "S1/local_run_inventory.tsv",
        "S1/lung_acquisition_eligibility.tsv",
        "S4/external_detection_summary.tsv",
        "S6/skin_reference_context.tsv",
    ),
    "3": (
        "S3/jy_aligned_recovery.tsv",
        "S3/HP20289_paired_recovery.tsv",
        "S3/pxd045211_pride_recipe_SCP0004_IPP372_islet22_Slot2-1_1_355.tsv",
        "S3/pxd045211_pride_recipe_SCP0004_IPP372_islet22_Slot2-1_1_356.tsv",
        "S3/pxd045211_pride_recipe_SCP0004_IPP372_islet27_Slot2-2_1_358.tsv",
        "S3/pxd045211_pride_recipe_SCP0004_IPP372_islet27_Slot2-2_1_359.tsv",
        "S3/pxd045211_raw_validation_SCP0004_IPP372_islet22_Slot2-1_1_355.tsv",
        "S3/pxd045211_raw_validation_SCP0004_IPP372_islet22_Slot2-1_1_356.tsv",
        "S3/pxd045211_raw_validation_SCP0004_IPP372_islet27_Slot2-2_1_358.tsv",
        "S3/pxd045211_raw_validation_SCP0004_IPP372_islet27_Slot2-2_1_359.tsv",
    ),
    "4": (
        "S2/archived_exact_panel_candidate_donor_endpoints.tsv",
        "S2/archived_exact_panel_candidate_run_matrix.tsv",
        "S2/archived_exact_panel_compatible_single_spectra.tsv",
    ),
    "5": (
        "S5/SIASAISMK_islet_spectrum_evidence.tsv",
        "S5/SIASAISMK_precursor_chemistry.tsv",
        "S5/SIASAISMK_tumor_reference_evidence.tsv",
    ),
    "S1": ("S3/jy_aligned_recovery.tsv",),
    "S2": ("S2/ILMILQPQL_donor_HLA_S1_S2.tsv",),
    "S3": (),
    "S4": (),
    "S5": (
        "S5/candidate_ion_peak_assignments.tsv",
        "S5/SIASAISMK_islet_spectrum_evidence.tsv",
        "S5/SIASAISMK_tumor_reference_evidence.tsv",
        "S5/tumour_hypothesis_by_ion.tsv",
    ),
    "S6": (),
}

AUTHOR_SUMMARIES = {
    "S3": "S3_auc.tsv",
    "S4": "S4_atlas_entries.csv",
    "S6": "S6_candidate_lists.csv",
}


def safe_member(name: str) -> str:
    path = PurePosixPath(name)
    if path.is_absolute() or not path.parts or any(part in ("", ".", "..") for part in path.parts):
        raise ValueError(f"Unsafe archive path: {name!r}")
    return path.as_posix()


def read_verified_archive(path: Path) -> dict[str, bytes]:
    try:
        archive = zipfile.ZipFile(path)
    except (OSError, zipfile.BadZipFile) as exc:
        raise SystemExit(f"Cannot open supplement ZIP: {exc}") from exc
    with archive:
        members: dict[str, zipfile.ZipInfo] = {}
        for info in archive.infolist():
            name = safe_member(info.filename)
            if name in members:
                raise SystemExit(f"Duplicate path in supplement ZIP: {name}")
            members[name] = info
        if "MANIFEST.tsv" not in members:
            raise SystemExit("Supplement ZIP has no MANIFEST.tsv")
        try:
            manifest_text = archive.read("MANIFEST.tsv").decode("utf-8-sig")
            rows = list(csv.DictReader(io.StringIO(manifest_text), delimiter="\t"))
        except (UnicodeDecodeError, csv.Error, KeyError) as exc:
            raise SystemExit(f"Cannot read supplement manifest: {exc}") from exc
        if not rows or not {"path", "sha256", "bytes"}.issubset(rows[0]):
            raise SystemExit("MANIFEST.tsv must have path, sha256, and bytes columns")

        verified: dict[str, bytes] = {}
        seen: set[str] = set()
        for row in rows:
            name = safe_member(row["path"])
            if name in seen:
                raise SystemExit(f"Duplicate path in supplement manifest: {name}")
            seen.add(name)
            if name not in members:
                raise SystemExit(f"Manifest file is absent from ZIP: {name}")
            data = archive.read(name)
            if len(data) != int(row["bytes"]):
                raise SystemExit(f"Size mismatch for {name}")
            digest = hashlib.sha256(data).hexdigest()
            if digest != row["sha256"].strip().lower():
                raise SystemExit(f"SHA-256 mismatch for {name}")
            verified[name] = data
        return verified


def load_generator(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise SystemExit(f"Cannot load bundled figure generator: {path.name}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--supplement", type=Path, required=True,
                        help="Path to the C018_Supplementary_Data_v4.zip file")
    parser.add_argument("--out", type=Path, required=True,
                        help="Directory for generated PDF figures and source-data CSVs")
    supported = ("1", "2", "3", "4", "5", "S1", "S2", "S3", "S4", "S5", "S6")
    parser.add_argument("--figures", nargs="*", choices=supported, default=supported,
                        help="Figures to render (default: all currently supported figures)")
    args = parser.parse_args()

    archive_data = read_verified_archive(args.supplement.expanduser().resolve())
    for number in args.figures:
        missing = [name for name in REQUIRED[number] if name not in archive_data]
        if missing:
            raise SystemExit(f"Supplement manifest lacks Figure {number} inputs: {', '.join(missing)}")

    args.out.expanduser().resolve().mkdir(parents=True, exist_ok=True)
    main_dir = args.out.expanduser().resolve() / "Figures"
    supp_dir = args.out.expanduser().resolve() / "Supplementary"
    with tempfile.TemporaryDirectory(prefix="c018_figure_inputs_") as temp:
        temp_root = Path(temp)
        os.environ.setdefault("MPLBACKEND", "Agg")
        os.environ.setdefault("MPLCONFIGDIR", str(temp_root / "mplconfig"))
        try:
            figure2 = load_generator("c018_v4_figure2", GENERATORS / "figure2_generator.py")
            figure4 = load_generator("c018_v4_figure4", GENERATORS / "figure4_generator.py")
        except ImportError as exc:
            raise SystemExit(f"Figure rendering requires matplotlib: {exc}") from exc

        for number in args.figures:
            for name in REQUIRED[number]:
                target = temp_root.joinpath(*PurePosixPath(name).parts)
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes(archive_data[name])
            if number in AUTHOR_SUMMARIES:
                summary_name = AUTHOR_SUMMARIES[number]
                summary_source = GENERATORS / "author_summaries" / summary_name
                summary_target = temp_root / "author_summaries" / summary_name
                summary_target.parent.mkdir(parents=True, exist_ok=True)
                summary_target.write_bytes(summary_source.read_bytes())
        figure_builders = {
            "1": figure2.figure_1,
            "2": figure2.figure_2,
            "3": figure2.figure_3,
            "4": figure4.figure_4,
            "5": figure4.figure_5,
            "S1": figure2.figure_s1,
            "S2": figure2.figure_s2,
            "S3": figure2.figure_s3,
            "S4": figure2.figure_s4,
            "S5": figure4.figure_s5,
            "S6": figure2.figure_s6,
        }
        for number in args.figures:
            figure_builders[number](temp_root, main_dir, supp_dir)

    print(f"Rendered C018 v4 Figure(s) {', '.join(args.figures)} to {args.out.expanduser().resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

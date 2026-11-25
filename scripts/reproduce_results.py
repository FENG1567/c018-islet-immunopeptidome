#!/usr/bin/env python3
"""Replay locked C018 results using only Python's standard library.

Usage (from the repository root)::

    python3 scripts/reproduce_results.py \
      --supplement /path/to/C018_Supplementary_Data_v4.zip \
      --out /path/to/replay-output

The supplement is read in place and is never extracted or modified. Every
manifest-listed file is checked for byte length and SHA-256 before any result
table is used. Derived results are compared with the author-owned
``data/expected_results.tsv`` values locked from v4R1. The output directory
receives ``reproduced_results.tsv``. No network access or third-party package
is used; the raw search and source data are not rerun.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import math
import re
import sys
import zipfile
from collections import defaultdict
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
EXPECTED_RESULTS = ROOT / "data" / "expected_results.tsv"
MANIFEST = "MANIFEST.tsv"


class ReplayError(Exception):
    """An input or result does not satisfy the locked replay contract."""


def read_tsv(archive: zipfile.ZipFile, member: str) -> list[dict[str, str]]:
    try:
        with archive.open(member, "r") as binary:
            stream = io.TextIOWrapper(binary, encoding="utf-8-sig", newline="")
            return list(csv.DictReader(stream, delimiter="\t"))
    except KeyError as exc:
        raise ReplayError(f"Required supplement member is missing: {member}") from exc
    except (UnicodeDecodeError, csv.Error) as exc:
        raise ReplayError(f"Could not parse {member} as UTF-8 TSV: {exc}") from exc


def verify_manifest(archive: zipfile.ZipFile) -> int:
    names = [info.filename for info in archive.infolist() if not info.is_dir()]
    if len(names) != len(set(names)):
        raise ReplayError("Supplement contains duplicate file names")
    if MANIFEST not in names:
        raise ReplayError(f"Supplement does not contain {MANIFEST}")

    manifest_rows = read_tsv(archive, MANIFEST)
    required_columns = {"path", "sha256", "bytes"}
    if not manifest_rows or not required_columns.issubset(manifest_rows[0]):
        raise ReplayError(f"{MANIFEST} must contain path, sha256, and bytes columns")

    expected_members: set[str] = set()
    for row in manifest_rows:
        member = row["path"]
        parts = Path(member).parts
        if (not member or "\\" in member or member.startswith("/") or
                any(part in {"", ".", ".."} for part in member.split("/"))):
            raise ReplayError(f"Unsafe path in {MANIFEST}: {member!r}")
        if member in expected_members or member == MANIFEST:
            raise ReplayError(f"Duplicate or self-referential path in {MANIFEST}: {member}")
        expected_members.add(member)
        try:
            expected_size = int(row["bytes"])
        except ValueError as exc:
            raise ReplayError(f"Invalid byte count in {MANIFEST}: {member}") from exc
        expected_hash = row["sha256"].lower()
        if not re.fullmatch(r"[0-9a-f]{64}", expected_hash):
            raise ReplayError(f"Invalid SHA-256 in {MANIFEST}: {member}")
        try:
            with archive.open(member, "r") as stream:
                digest = hashlib.sha256()
                size = 0
                while True:
                    block = stream.read(1024 * 1024)
                    if not block:
                        break
                    digest.update(block)
                    size += len(block)
        except KeyError as exc:
            raise ReplayError(f"Manifest member is missing from ZIP: {member}") from exc
        if size != expected_size:
            raise ReplayError(f"Byte-count mismatch for {member}: {size} != {expected_size}")
        if digest.hexdigest() != expected_hash:
            raise ReplayError(f"SHA-256 mismatch for {member}")

    actual_members = set(names) - {MANIFEST}
    if actual_members != expected_members:
        missing = sorted(expected_members - actual_members)
        unexpected = sorted(actual_members - expected_members)
        raise ReplayError(
            f"ZIP and manifest member sets differ (missing={missing[:5]}, "
            f"unexpected={unexpected[:5]})"
        )
    return len(manifest_rows)


def as_int(value: str, context: str) -> int:
    try:
        return int(value)
    except (TypeError, ValueError) as exc:
        raise ReplayError(f"Expected integer in {context}, got {value!r}") from exc


def as_float(value: str, context: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ReplayError(f"Expected number in {context}, got {value!r}") from exc
    if not math.isfinite(result):
        raise ReplayError(f"Non-finite number in {context}: {value!r}")
    return result


def flag(value: str, context: str) -> bool:
    normalized = value.strip().lower()
    if normalized in {"true", "1", "yes"}:
        return True
    if normalized in {"false", "0", "no", ""}:
        return False
    raise ReplayError(f"Expected Boolean in {context}, got {value!r}")


def fmt_number(value: float) -> str:
    return format(value, ".10g")


def add_result(
    results: dict[str, tuple[str, str, str]],
    metric: str,
    value: str | int | float,
    evidence_type: str,
    source: str,
) -> None:
    if metric in results:
        raise ReplayError(f"Internal replay produced duplicate metric: {metric}")
    if isinstance(value, float):
        rendered = fmt_number(value)
    else:
        rendered = str(value)
    results[metric] = (rendered, evidence_type, source)


def replay_pdac(archive: zipfile.ZipFile, results: dict[str, tuple[str, str, str]]) -> None:
    source = "S2/archived_exact_panel_candidate_run_matrix.tsv"
    rows = read_tsv(archive, source)
    by_peptide: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_peptide[row["peptide"]].append(row)

    if len(by_peptide) != 517 or len(rows) != 517 * 5:
        raise ReplayError(f"PDAC matrix is not 517 x 5: {len(by_peptide)} sequences, {len(rows)} rows")

    all_runs: set[str] = set()
    endpoints: dict[str, dict[str, object]] = {}
    for peptide, peptide_rows in by_peptide.items():
        if len(peptide_rows) != 5:
            raise ReplayError(f"PDAC candidate {peptide} has {len(peptide_rows)} run rows, expected 5")
        runs = [row["run"] for row in peptide_rows]
        if len(set(runs)) != 5:
            raise ReplayError(f"PDAC candidate {peptide} does not have five distinct runs")
        all_runs.update(runs)

        native_flags = {flag(row["native_unmodified_scope"], f"{source}:{peptide}") for row in peptide_rows}
        if len(native_flags) != 1:
            raise ReplayError(f"Native-scope status varies across runs for {peptide}")
        native_scope = native_flags.pop()

        compatible_by_donor: dict[str, set[bool]] = defaultdict(set)
        scans_by_donor: dict[str, int] = defaultdict(int)
        total_scans = 0
        compatible_scans = 0
        for row in peptide_rows:
            donor = row["donor"]
            compatible = flag(row["exact_typed_compatible"], f"{source}:{peptide}:{donor}")
            scan_n = as_int(row["accepted_exact_scan_n"], f"{source}:{peptide}:{row['run']}")
            if scan_n < 0:
                raise ReplayError(f"Negative scan count for {peptide} in {row['run']}")
            compatible_by_donor[donor].add(compatible)
            total_scans += scan_n
            if compatible:
                compatible_scans += scan_n
                scans_by_donor[donor] += scan_n

        if any(len(statuses) != 1 for statuses in compatible_by_donor.values()):
            raise ReplayError(f"HLA compatibility varies between runs for a candidate-donor pair: {peptide}")
        compatible_donors = {
            donor for donor, statuses in compatible_by_donor.items() if True in statuses
        }
        donor_hits = sum(count >= 2 for count in scans_by_donor.values())
        endpoints[peptide] = {
            "native": native_scope,
            "eligible_any": len(compatible_donors) >= 1,
            "eligible_both": len(compatible_donors) >= 2,
            "donors_ge2": donor_hits,
            "any_detected": donor_hits >= 1,
            "both_detected": donor_hits >= 2,
            "total_scans": total_scans,
            "compatible_scans": compatible_scans,
        }

    if len(all_runs) != 5:
        raise ReplayError(f"PDAC matrix has {len(all_runs)} distinct runs, expected 5")

    eligible_any = sum(bool(item["eligible_any"]) for item in endpoints.values())
    eligible_both = sum(bool(item["eligible_both"]) for item in endpoints.values())
    native_any = sum(bool(item["native"] and item["eligible_any"]) for item in endpoints.values())
    native_both = sum(bool(item["native"] and item["eligible_both"]) for item in endpoints.values())
    native_any_detected = sum(bool(item["native"] and item["any_detected"]) for item in endpoints.values())
    native_both_detected = sum(bool(item["native"] and item["both_detected"]) for item in endpoints.values())
    all_scans = sum(int(item["total_scans"]) for item in endpoints.values())
    compatible_scans = sum(int(item["compatible_scans"]) for item in endpoints.values())

    add_result(results, "pdac_matrix_rows", len(rows), "recomputed", source)
    add_result(results, "pdac_distinct_sequences", len(by_peptide), "recomputed", source)
    add_result(results, "pdac_distinct_runs", len(all_runs), "recomputed", source)
    add_result(results, "pdac_exact_hla_eligible_any", eligible_any, "recomputed", source)
    add_result(results, "pdac_exact_hla_eligible_both", eligible_both, "recomputed", source)
    add_result(results, "pdac_native_scope_exact_hla_eligible_any", native_any, "recomputed", source)
    add_result(results, "pdac_native_scope_exact_hla_eligible_both", native_both, "recomputed", source)
    add_result(results, "pdac_accepted_exact_scans_total", all_scans, "recomputed", source)
    add_result(results, "pdac_accepted_exact_scans_exact_hla_compatible", compatible_scans, "recomputed", source)
    add_result(results, "pdac_native_scope_any_endpoint_detected", native_any_detected, "recomputed", source)
    add_result(results, "pdac_native_scope_any_endpoint_denominator", native_any, "recomputed", source)
    add_result(results, "pdac_native_scope_both_endpoint_detected", native_both_detected, "recomputed", source)
    add_result(results, "pdac_native_scope_both_endpoint_denominator", native_both, "recomputed", source)

    # Cross-check row-wise aggregation against the separately supplied endpoint table.
    endpoint_source = "S2/archived_exact_panel_candidate_donor_endpoints.tsv"
    endpoint_rows = read_tsv(archive, endpoint_source)
    if {row["peptide"] for row in endpoint_rows} != set(endpoints):
        raise ReplayError("PDAC run matrix and donor-endpoint table contain different candidate sets")
    for row in endpoint_rows:
        item = endpoints[row["peptide"]]
        checks = {
            "native_unmodified_scope": bool(item["native"]),
            "any_endpoint": bool(item["any_detected"]),
            "two_endpoint": bool(item["both_detected"]),
        }
        for column, derived in checks.items():
            if flag(row[column], f"{endpoint_source}:{row['peptide']}:{column}") != derived:
                raise ReplayError(f"PDAC endpoint cross-check failed for {row['peptide']} ({column})")
        if as_int(row["compatible_donor_n"], f"{endpoint_source}:{row['peptide']}") != (
            2 if item["eligible_both"] else (1 if item["eligible_any"] else 0)
        ):
            raise ReplayError(f"PDAC compatible-donor count cross-check failed for {row['peptide']}")
        if as_int(row["donors_ge2_spectra_n"], f"{endpoint_source}:{row['peptide']}") != int(item["donors_ge2"]):
            raise ReplayError(f"PDAC donor scan-threshold cross-check failed for {row['peptide']}")
        # The donor endpoint table intentionally counts accepted scans only in
        # HLA-compatible donors. The run matrix also retains scans in untyped
        # or incompatible donors (for example FQAAPDRGF in R361); those are
        # included in the all-run diagnostic above, but not in this endpoint.
        if as_int(row["accepted_exact_scan_n"], f"{endpoint_source}:{row['peptide']}") != int(item["compatible_scans"]):
            raise ReplayError(f"PDAC compatible-donor scan-count cross-check failed for {row['peptide']}")


def replay_aetsa(archive: zipfile.ZipFile, results: dict[str, tuple[str, str, str]]) -> None:
    source = "S2/NC2025_exact_hla_tissue_resolved_candidates.tsv"
    rows = read_tsv(archive, source)
    endpoints = {
        "any": ("endpoint_any_eligible", "legacy_any_detected_ge2", "q01_any_detected_ge2"),
        "both": ("endpoint_both_eligible", "legacy_both_detected_ge2", "q01_both_detected_ge2"),
    }
    computed: dict[str, tuple[int, int, int]] = {}
    for name, (eligible_col, legacy_col, q01_col) in endpoints.items():
        eligible_rows = [row for row in rows if flag(row[eligible_col], f"{source}:{eligible_col}")]
        computed[name] = (
            len(eligible_rows),
            sum(as_int(row[legacy_col], f"{source}:{legacy_col}") > 0 for row in eligible_rows),
            sum(as_int(row[q01_col], f"{source}:{q01_col}") > 0 for row in eligible_rows),
        )

    summary_source = "S2/aetsa_tissue_resolved_endpoints.tsv"
    summary = {
        row["endpoint"]: row
        for row in read_tsv(archive, summary_source)
        if row["typing"] == "published_exact_HLA"
    }
    for name, (denominator, legacy_n, q01_n) in computed.items():
        if name not in summary:
            raise ReplayError(f"aeTSA summary is missing published_exact_HLA/{name}")
        if as_int(summary[name]["covered_NC2025"], f"{summary_source}:{name}") != denominator:
            raise ReplayError(f"aeTSA denominator cross-check failed for {name}")
        if as_int(summary[name]["eligible_donor_two_scan_detections"], f"{summary_source}:{name}") != q01_n:
            raise ReplayError(f"aeTSA detection cross-check failed for {name}")
        add_result(results, f"aetsa_exact_hla_{name}_eligible", denominator, "recomputed", source)
        add_result(results, f"aetsa_exact_hla_{name}_legacy_detected", legacy_n, "recomputed", source)
        add_result(results, f"aetsa_exact_hla_{name}_q01_detected", q01_n, "recomputed", source)


def recovery_and_models(archive: zipfile.ZipFile, results: dict[str, tuple[str, str, str]]) -> None:
    jy_source = "S3/jy_aligned_recovery.tsv"
    jy_rows = read_tsv(archive, jy_source)
    matches = [
        row for row in jy_rows
        if row["mode"] == "aligned_psm_and_peptide_q_plus1"
        and row["dilution_factor"] == "16"
        and row["reference_scans_min"] == "1"
        and row["test_scans_min"] == "2"
    ]
    if len(matches) != 1:
        raise ReplayError(f"Expected one 1:16 JY one-reference/two-test recovery row, found {len(matches)}")
    recovery_row = matches[0]
    recovered = as_int(recovery_row["recovered_n"], f"{jy_source}:recovered_n")
    reference_n = as_int(recovery_row["reference_n"], f"{jy_source}:reference_n")
    recovery = recovered / reference_n
    if not math.isclose(recovery, as_float(recovery_row["recovery"], f"{jy_source}:recovery"), abs_tol=5e-8):
        raise ReplayError("JY recovery fraction disagrees with recovered/reference counts")
    add_result(results, "jy_conditional_recovery_numerator", recovered, "carried-through", jy_source)
    add_result(results, "jy_conditional_recovery_denominator", reference_n, "carried-through", jy_source)
    add_result(results, "jy_conditional_recovery_fraction", recovery, "recomputed", jy_source)

    candidate_source = "S2/NC2025_exact_hla_tissue_resolved_candidates.tsv"
    candidates = read_tsv(archive, candidate_source)
    model_source = "S3/conditional_transfer_model.tsv"
    model_rows = [
        row for row in read_tsv(archive, model_source)
        if row["cohort_scope"] == "Tissue-resolved" and row["endpoint"] in {"any", "both"}
    ]
    model = {row["endpoint"]: row for row in model_rows}
    if set(model) != {"any", "both"}:
        raise ReplayError("Conditional-transfer table lacks tissue-resolved any/both rows")
    prevalences = {as_float(row["assumed_prevalence"], f"{model_source}:assumed_prevalence") for row in model.values()}
    if len(prevalences) != 1:
        raise ReplayError("Conditional-transfer endpoints use different prevalence assumptions")
    prevalence = prevalences.pop()
    add_result(results, "jy_model_assumed_prevalence", prevalence, "carried-through", model_source)

    any_rows = [row for row in candidates if flag(row["endpoint_any_eligible"], f"{candidate_source}:endpoint_any_eligible")]
    both_rows = [row for row in candidates if flag(row["endpoint_both_eligible"], f"{candidate_source}:endpoint_both_eligible")]

    def compatible_donor_count(row: dict[str, str]) -> int:
        value = row["compatible_tissue_resolved_donors"].strip()
        return len([name for name in value.split(";") if name.strip()]) if value else 0

    any_profile = defaultdict(int)
    for row in any_rows:
        donor_n = compatible_donor_count(row)
        if donor_n < 1:
            raise ReplayError(f"Eligible aeTSA candidate lacks a compatible donor: {row['peptide']}")
        any_profile[donor_n] += 1
    p0_any = 1.0
    for donor_n, candidate_n in any_profile.items():
        p_detect_if_presented = 1.0 - (1.0 - recovery) ** donor_n
        p0_any *= (1.0 - prevalence * p_detect_if_presented) ** candidate_n

    def probability_at_least_two_successes(n: int, probability: float) -> float:
        if n < 2:
            return 0.0
        return 1.0 - (1.0 - probability) ** n - n * probability * (1.0 - probability) ** (n - 1)

    p0_both = 1.0
    for row in both_rows:
        donor_n = compatible_donor_count(row)
        p_detect_if_presented = probability_at_least_two_successes(donor_n, recovery)
        p0_both *= 1.0 - prevalence * p_detect_if_presented

    for endpoint, eligible_n, p0 in (("any", len(any_rows), p0_any), ("both", len(both_rows), p0_both)):
        source_row = model[endpoint]
        if as_int(source_row["eligible_candidate_n"], f"{model_source}:{endpoint}") != eligible_n:
            raise ReplayError(f"Conditional-transfer candidate count cross-check failed for {endpoint}")
        if not math.isclose(p0, as_float(source_row["conditional_P_zero"], f"{model_source}:{endpoint}"), abs_tol=1e-8):
            raise ReplayError(f"Conditional-transfer model cross-check failed for {endpoint}: {p0:.10g}")
        add_result(results, f"jy_conditional_model_{endpoint}_eligible", eligible_n, "recomputed", candidate_source)
        add_result(results, f"jy_conditional_model_{endpoint}_P_zero", p0, "recomputed", candidate_source)


def replay_external(archive: zipfile.ZipFile, results: dict[str, tuple[str, str, str]]) -> None:
    candidate_source = "S4/external_candidate_all_summary.tsv"
    rows = [
        row for row in read_tsv(archive, candidate_source)
        if row["list"] == "PDAC2025" and row["possible_A02_A03_family"].strip()
    ]
    if len({row["peptide"] for row in rows}) != 122:
        raise ReplayError(f"Expected 122 family-supported PDAC external candidates, found {len(rows)}")
    if len(rows) != len({row["peptide"] for row in rows}):
        raise ReplayError("PDAC external candidate table contains duplicate candidate rows")

    counts = {
        "legacy_any_preparation": sum(as_int(row["legacy_recipe_n_preparations"], candidate_source) >= 1 for row in rows),
        "legacy_both_preparations": sum(as_int(row["legacy_recipe_n_preparations"], candidate_source) >= 2 for row in rows),
        "pride_any_preparation": sum(as_int(row["pride_recipe_n_preparations"], candidate_source) >= 1 for row in rows),
        "pride_both_preparations": sum(as_int(row["pride_recipe_n_preparations"], candidate_source) >= 2 for row in rows),
    }
    for metric, value in counts.items():
        add_result(results, f"pdac_external_{metric}", value, "recomputed", candidate_source)
    add_result(results, "pdac_external_family_supported_candidates", len(rows), "recomputed", candidate_source)

    evidence_source = "S4/joint_null_sequence_evidence.tsv"
    evidence_rows = [
        row for row in read_tsv(archive, evidence_source)
        if row["list"] == "PDAC2025" and flag(row["family_possible"], f"{evidence_source}:family_possible")
    ]
    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in evidence_rows:
        if row["sequence_type"] not in {"target", "matched_null"}:
            raise ReplayError(f"Unexpected sequence type in {evidence_source}: {row['sequence_type']}")
        grouped[(row["sequence_type"], row["parent_target"], row["query_sequence"])].append(row)
    target_pass = 0
    null_pass = 0
    target_n = 0
    null_n = 0
    for (kind, _parent, _query), group in grouped.items():
        passes = {flag(row["all_four_technical_acquisitions"], f"{evidence_source}:all_four_technical_acquisitions") for row in group}
        if len(passes) != 1:
            raise ReplayError(f"Inconsistent repeated joint-null rows for {kind} {_parent}/{_query}")
        if kind == "target":
            target_n += 1
            target_pass += int(passes.pop())
        else:
            null_n += 1
            null_pass += int(passes.pop())
    if target_n != 122 or null_n != 1220:
        raise ReplayError(f"Joint-null family-supported denominators differ: targets={target_n}, nulls={null_n}")
    add_result(results, "pdac_external_all_four_targets_passing", target_pass, "recomputed", evidence_source)
    add_result(results, "pdac_external_all_four_targets_denominator", target_n, "recomputed", evidence_source)
    add_result(results, "pdac_external_all_four_matched_nulls_passing", null_pass, "recomputed", evidence_source)
    add_result(results, "pdac_external_all_four_matched_nulls_denominator", null_n, "recomputed", evidence_source)


def replay_carried_summaries(archive: zipfile.ZipFile, results: dict[str, tuple[str, str, str]]) -> None:
    lung_source = "S2/lung_endpoints.tsv"
    lung_rows = read_tsv(archive, lung_source)
    keys = {
        (row["scenario"], row["endpoint"]): row for row in lung_rows
    }
    lung_expected = {
        ("S1_passing", "any"): ("primary", "any"),
        ("S1_passing", "both"): ("primary", "both"),
        ("S2_reversed", "any"): ("reversed", "any"),
        ("S2_reversed", "both"): ("reversed", "both"),
    }
    for key, (scenario_label, endpoint_label) in lung_expected.items():
        if key not in keys:
            raise ReplayError(f"Lung endpoint summary lacks {key[0]}/{key[1]}")
        row = keys[key]
        prefix = f"lung_{scenario_label}_{endpoint_label}"
        add_result(results, f"{prefix}_eligible_candidates", as_int(row["eligible_candidates"], lung_source), "carried-through", lung_source)
        add_result(results, f"{prefix}_signals_meeting_spectrum_criterion", as_int(row["signals_meeting_spectrum_criterion"], lung_source), "carried-through", lung_source)
        add_result(results, f"{prefix}_confirmed_presentation_calls", as_int(row["confirmed_presentation_calls"], lung_source), "carried-through", lung_source)

    auc_source = "S3/targeted_extraction_summary.tsv"
    auc_rows = read_tsv(archive, auc_source)
    auc_matches = [row for row in auc_rows if row["metric"] == "median_candidate_versus_shuffled_label_auc"]
    if len(auc_matches) != 1:
        raise ReplayError("Targeted-extraction table does not contain one median candidate-versus-shuffled AUC")
    add_result(results, "targeted_extraction_median_auc", as_float(auc_matches[0]["value"], auc_source), "carried-through", auc_source)


def read_expected() -> list[dict[str, str]]:
    try:
        with EXPECTED_RESULTS.open("r", encoding="utf-8-sig", newline="") as stream:
            rows = list(csv.DictReader(stream, delimiter="\t"))
    except OSError as exc:
        raise ReplayError(f"Could not read author-locked expected results {EXPECTED_RESULTS}: {exc}") from exc
    required = {"metric", "expected_value", "tolerance", "evidence_type", "source"}
    if not rows or not required.issubset(rows[0]):
        raise ReplayError(f"{EXPECTED_RESULTS} must contain {', '.join(sorted(required))}")
    if len({row["metric"] for row in rows}) != len(rows):
        raise ReplayError(f"Duplicate metric in {EXPECTED_RESULTS}")
    return rows


def compare_and_write(
    results: dict[str, tuple[str, str, str]],
    expected_rows: list[dict[str, str]],
    manifest_count: int,
    output_dir: Path,
) -> int:
    expected_by_metric = {row["metric"]: row for row in expected_rows}
    observed_metrics = set(results)
    missing = set(expected_by_metric) - observed_metrics
    extra = observed_metrics - set(expected_by_metric)
    if missing or extra:
        raise ReplayError(f"Expected/result metric sets differ (missing={sorted(missing)}, extra={sorted(extra)})")

    output_rows: list[dict[str, str]] = [{
        "metric": "supplement_manifest_entries_sha256_verified",
        "observed_value": str(manifest_count),
        "expected_value": "",
        "evidence_type": "recomputed",
        "source": MANIFEST,
        "status": "PASS",
    }]
    failures = 0
    for metric, expected in expected_by_metric.items():
        observed, evidence_type, source = results[metric]
        if evidence_type != expected["evidence_type"]:
            raise ReplayError(
                f"Evidence classification mismatch for {metric}: "
                f"script={evidence_type}, expected={expected['evidence_type']}"
            )
        tolerance = as_float(expected["tolerance"], f"{EXPECTED_RESULTS}:{metric}:tolerance")
        if tolerance == 0:
            matched = observed == expected["expected_value"]
        else:
            matched = math.isclose(
                as_float(observed, f"observed {metric}"),
                as_float(expected["expected_value"], f"expected {metric}"),
                rel_tol=0.0,
                abs_tol=tolerance,
            )
        if not matched:
            failures += 1
        output_rows.append({
            "metric": metric,
            "observed_value": observed,
            "expected_value": expected["expected_value"],
            "evidence_type": evidence_type,
            "source": source,
            "status": "PASS" if matched else "FAIL",
        })

    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "reproduced_results.tsv"
    with output_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=["metric", "observed_value", "expected_value", "evidence_type", "source", "status"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(output_rows)
    print(f"Verified {manifest_count} manifest-listed files; compared {len(expected_rows)} result metrics.")
    print(f"Replay table: {output_path}")
    if failures:
        print(f"FAILED: {failures} result metric(s) differ from the v4R1 lock.", file=sys.stderr)
    else:
        print("All result metrics match the v4R1 lock.")
    return failures


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Verify the locked C018 v4 supplement and replay result endpoints.")
    parser.add_argument("--supplement", required=True, type=Path, help="Path to the actual v4 final supplementary ZIP")
    parser.add_argument("--out", required=True, type=Path, help="Directory for reproduced_results.tsv")
    args = parser.parse_args(argv)

    try:
        if not args.supplement.is_file():
            raise ReplayError(f"Supplement ZIP does not exist: {args.supplement}")
        results: dict[str, tuple[str, str, str]] = {}
        with zipfile.ZipFile(args.supplement, "r") as archive:
            manifest_count = verify_manifest(archive)
            replay_pdac(archive, results)
            replay_aetsa(archive, results)
            recovery_and_models(archive, results)
            replay_external(archive, results)
            replay_carried_summaries(archive, results)
        expected_rows = read_expected()
        return 1 if compare_and_write(results, expected_rows, manifest_count, args.out) else 0
    except (ReplayError, OSError, zipfile.BadZipFile) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())

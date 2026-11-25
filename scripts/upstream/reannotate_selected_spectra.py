#!/usr/bin/env python3
"""Re-annotate selected C018 PXD045211 candidate MS/MS scans.

Reads the compact candidate-index MGF exports from the v2 server archive and
the two archived SIASAISMK tumour-reference PROXI peak lists. Emits observed
ion matches with mass errors and separately labelled tumour modification
hypotheses; it does not re-search the full dataset.
"""

import argparse
import csv
import json
import math
from pathlib import Path


RESIDUE = dict(zip("ACDEFGHIKLMNPQRSTVWY", [
    71.037114, 103.009185, 115.026943, 129.042593, 147.068414, 57.021464,
    137.058912, 113.084064, 128.094963, 113.084064, 131.040485,
    114.042927, 97.052764, 128.058578, 156.101111, 87.032028,
    101.047679, 99.068414, 186.079313, 163.063329,
]))
PROTON = 1.007276466
WATER = 18.010564684
TOL = 0.020
PAIRS = {
    "HLLQEELLL": "HLLQQKDIL",
    "SIASAISMK": "ALSSALQFA",
}
RUN_SCANS = {
    "SCP0004_IPP372_islet22_Slot2-1_1_355": {57021, 120696, 120813, 120865, 120950},
    "SCP0004_IPP372_islet22_Slot2-1_1_356": {56976, 119748, 119834},
    "SCP0004_IPP372_islet27_Slot2-2_1_358": {58077},
    "SCP0004_IPP372_islet27_Slot2-2_1_359": {57601},
}


def read_mgf(path: Path, requested: set[str]) -> dict[str, dict]:
    found = {}
    meta = None
    peaks = []
    with path.open("rt", encoding="ascii") as stream:
        for line in stream:
            line = line.strip()
            if line == "BEGIN IONS":
                meta, peaks = {}, []
            elif line == "END IONS":
                scan = meta.get("SCANS", "")
                if scan in requested:
                    meta["peaks"] = peaks
                    found[scan] = meta
                    if requested <= found.keys():
                        break
                meta = None
            elif meta is not None and "=" in line:
                key, value = line.split("=", 1)
                meta[key] = value
            elif meta is not None and line and line[0].isdigit():
                fields = line.split()
                if len(fields) >= 2 and float(fields[1]) > 0:
                    peaks.append((float(fields[0]), float(fields[1])))
    missing = requested - found.keys()
    if missing:
        raise ValueError(f"{path}: missing requested SCANS {sorted(missing)}")
    return found


def fragment_ions(sequence: str, precursor_charge: int, mods: dict[int, float] | None = None):
    mods = mods or {}
    residues = [RESIDUE[aa] + mods.get(index, 0.0) for index, aa in enumerate(sequence)]
    total = sum(residues)
    prefix = 0.0
    for cut, residue_mass in enumerate(residues[:-1], start=1):
        prefix += residue_mass
        suffix = total - prefix
        for charge in range(1, min(2, precursor_charge) + 1):
            yield f"b{cut}^{charge}", (prefix + charge * PROTON) / charge
            y_n = len(sequence) - cut
            yield f"y{y_n}^{charge}", (suffix + WATER + charge * PROTON) / charge


def best_peak(peaks, theoretical_mz):
    matches = [(abs(mz - theoretical_mz), -intensity, mz, intensity)
               for mz, intensity in peaks if abs(mz - theoretical_mz) <= TOL]
    if not matches:
        return None
    _, _, mz, intensity = min(matches)
    return mz, intensity, (mz - theoretical_mz) / theoretical_mz * 1e6


def load_tumour_peaks(path: Path):
    body = json.loads(path.read_text(encoding="utf-8"))
    spectrum = body[0] if isinstance(body, list) else body
    peaks = [(float(mz), float(intensity)) for mz, intensity in zip(spectrum["mzs"], spectrum["intensities"])
             if float(intensity) > 0]
    return spectrum.get("usi", ""), peaks


def write_tsv(path: Path, rows: list[dict]):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--v2-base", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    challenge = args.v2_base / "tissue_calibration/pxd045211_canonical_challenge"
    spectra = {}
    quality_path = args.v2_base / "tissue_calibration/pxd045211_fragment_evidence/external_fragment_psm_quality.tsv"
    native_scan_by_key = {}
    with quality_path.open("r", encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            if row.get("peptide") in PAIRS and row.get("native_scan"):
                native_scan_by_key[(row["source"], row["comet_scan"])] = row["native_scan"]
    for run, scans in RUN_SCANS.items():
        path = challenge / run / f"{run}.candidate_index.mgf"
        requested = {str(scan) for scan in scans}
        for scan, entry in read_mgf(path, requested).items():
            entry["run"] = run
            entry["mgf_title"] = entry.get("TITLE", "")
            key = (run, scan)
            if key not in native_scan_by_key:
                raise ValueError(f"No verified native_scan mapping in archived fragment evidence for {key}")
            entry["native_scan"] = native_scan_by_key[key]
            spectra[(run, scan)] = entry
        print(f"loaded {run}: {len(requested)} selected spectra from {path.name}", flush=True)

    assignment_rows = []
    summary_rows = []
    for (run, scan), spectrum in sorted(spectra.items()):
        precursor_charge = int(spectrum["CHARGE"].rstrip("+"))
        precursor_mz = float(spectrum["PEPMASS"].split()[0])
        rt = float(spectrum.get("RTINSECONDS", "nan"))
        for target, alternative in PAIRS.items():
            # Each selected scan belongs to exactly one peptide, as recorded in
            # the frozen candidate PSM evidence export.
            if (target == "HLLQEELLL" and int(scan) not in {120696, 120813, 120865, 120950, 119748, 119834}) or (
                target == "SIASAISMK" and int(scan) not in {57021, 56976, 58077, 57601}
            ):
                continue
            target_ions = list(fragment_ions(target, precursor_charge))
            alternative_ions = list(fragment_ions(alternative, precursor_charge))
            target_masses = [mz for _, mz in target_ions]
            alternative_masses = [mz for _, mz in alternative_ions]
            target_found, alternative_found = [], []
            for label, theoretical in target_ions:
                role = "shared" if min(abs(theoretical - mz) for mz in alternative_masses) <= 2 * TOL else "candidate_specific"
                peak = best_peak(spectrum["peaks"], theoretical)
                if peak and role == "candidate_specific":
                    target_found.append((label, peak[0], peak[1], peak[2]))
                assignment_rows.append({
                    "peptide": target, "competitor": alternative, "source_run": run,
                    "comet_scan": scan, "native_scan": spectrum["native_scan"],
                    "mgf_title": spectrum["mgf_title"],
                    "rt_seconds": rt, "precursor_mz": precursor_mz,
                    "precursor_charge": precursor_charge, "ion_sequence": "candidate",
                    "ion": label, "diagnostic_role": role,
                    "theoretical_mz": f"{theoretical:.6f}",
                    "observed_mz": f"{peak[0]:.6f}" if peak else "",
                    "observed_intensity": f"{peak[1]:.6f}" if peak else "",
                    "mass_error_ppm": f"{peak[2]:.3f}" if peak else "",
                    "fragment_tolerance_Da": TOL,
                })
            for label, theoretical in alternative_ions:
                role = "shared" if min(abs(theoretical - mz) for mz in target_masses) <= 2 * TOL else "alternative_specific"
                peak = best_peak(spectrum["peaks"], theoretical)
                if peak and role == "alternative_specific":
                    alternative_found.append((label, peak[0], peak[1], peak[2]))
                assignment_rows.append({
                    "peptide": target, "competitor": alternative, "source_run": run,
                    "comet_scan": scan, "native_scan": spectrum["native_scan"],
                    "mgf_title": spectrum["mgf_title"],
                    "rt_seconds": rt, "precursor_mz": precursor_mz,
                    "precursor_charge": precursor_charge, "ion_sequence": "competitor",
                    "ion": label, "diagnostic_role": role,
                    "theoretical_mz": f"{theoretical:.6f}",
                    "observed_mz": f"{peak[0]:.6f}" if peak else "",
                    "observed_intensity": f"{peak[1]:.6f}" if peak else "",
                    "mass_error_ppm": f"{peak[2]:.3f}" if peak else "",
                    "fragment_tolerance_Da": TOL,
                })
            summary_rows.append({
                "peptide": target, "competitor": alternative, "source_run": run,
                "comet_scan": scan, "native_scan": spectrum["native_scan"],
                "mgf_title": spectrum["mgf_title"], "precursor_mz": f"{precursor_mz:.6f}",
                "precursor_charge": precursor_charge,
                "candidate_specific_ions_observed": ";".join(x[0] for x in target_found),
                "candidate_specific_ion_n": len(target_found),
                "competitor_specific_ions_observed": ";".join(x[0] for x in alternative_found),
                "competitor_specific_ion_n": len(alternative_found),
                "fragment_tolerance_Da": TOL,
            })
    write_tsv(args.out / "candidate_ion_peak_assignments.tsv", assignment_rows)
    write_tsv(args.out / "candidate_ion_peak_summary.tsv", summary_rows)

    islet_sias = [(key, value) for key, value in sorted(spectra.items()) if int(key[1]) in {57021, 56976, 58077, 57601}]
    refs = [
        ("primary_C3L01031_scan16793", args.v2_base / "tissue_calibration/native_tumor_spectra/SIASAISMK_native_reference.json",
         2, [("unmodified_hypothesis", {}), ("M8_oxidized_hypothesis", {7: 15.994915})],
         "exact source PSM has uppercase sequence but no modification field"),
        ("PDO_PANFR0151_scan16495", args.v2_base / "tissue_calibration/native_tumor_spectra_alt/SIASAISMK_native_reference.json",
         2, [("M8_oxidized_source_PSM", {7: 15.994915})],
         "exact supplementary source PSM is SIASAISmK; M8 oxidation explicitly reported"),
    ]
    tumour_match_rows = []
    tumour_summary_rows = []
    for ref_name, ref_path, ref_charge, hypotheses, psm_note in refs:
        ref_usi, ref_peaks = load_tumour_peaks(ref_path)
        for hypothesis, mods in hypotheses:
            ref_ions = {label: (mz, best_peak(ref_peaks, mz))
                        for label, mz in fragment_ions("SIASAISMK", ref_charge, mods)}
            for (run, scan), islet in islet_sias:
                islet_charge = int(islet["CHARGE"].rstrip("+"))
                islet_ions = {label: (mz, best_peak(islet["peaks"], mz))
                              for label, mz in fragment_ions("SIASAISMK", islet_charge)}
                same_mass_observed = []
                shifted_labels = []
                all_shared_labels = sorted(
                    label for label in set(ref_ions) & set(islet_ions) if label.endswith("^1")
                )
                for label in all_shared_labels:
                    tumor_mz, tumor_peak = ref_ions[label]
                    islet_mz, islet_peak = islet_ions[label]
                    same_mz = abs(tumor_mz - islet_mz) <= TOL
                    if same_mz and tumor_peak and islet_peak:
                        same_mass_observed.append((label, islet_peak[1], tumor_peak[1]))
                    if not same_mz:
                        shifted_labels.append(label)
                    tumour_match_rows.append({
                        "tumour_reference": ref_name, "reference_PSM_hypothesis": hypothesis,
                        "reference_PSM_note": psm_note, "tumour_reference_USI": ref_usi,
                        "islet_run": run, "islet_scan": scan,
                        "islet_native_scan": islet["native_scan"],
                        "islet_mgf_title": islet["mgf_title"], "ion": label,
                        "islet_form": "unmodified; z1 fragments",
                        "tumour_form_hypothesis": hypothesis,
                        "islet_theoretical_mz": f"{islet_mz:.6f}",
                        "islet_observed_mz": f"{islet_peak[0]:.6f}" if islet_peak else "",
                        "islet_intensity": f"{islet_peak[1]:.6f}" if islet_peak else "",
                        "islet_error_ppm": f"{islet_peak[2]:.3f}" if islet_peak else "",
                        "tumour_theoretical_mz": f"{tumor_mz:.6f}",
                        "tumour_observed_mz": f"{tumor_peak[0]:.6f}" if tumor_peak else "",
                        "tumour_intensity": f"{tumor_peak[1]:.6f}" if tumor_peak else "",
                        "tumour_error_ppm": f"{tumor_peak[2]:.3f}" if tumor_peak else "",
                        "tumour_minus_islet_theoretical_mz_Da": f"{tumor_mz - islet_mz:.6f}",
                        "same_theoretical_mz": "yes" if same_mz else "no",
                        "observed_peak_at_each_hypothesized_mz": "yes" if tumor_peak and islet_peak else "no",
                        "fragment_tolerance_Da": TOL,
                    })
                conditional_similarity = ""
                # Only calculate a same-form similarity under the unmodified
                # primary-reference hypothesis. For M8-vs-unmodified comparisons
                # retain the per-ion shifted-mass table and do not pool a cosine.
                if ref_name == "primary_C3L01031_scan16793" and hypothesis == "unmodified_hypothesis":
                    numerator = sum(math.sqrt(i * t) for _, i, t in same_mass_observed)
                    denom = math.sqrt(
                        sum(v[1][1] for label, v in islet_ions.items() if label.endswith("^1") and v[1]) *
                        sum(v[1][1] for label, v in ref_ions.items() if label.endswith("^1") and v[1])
                    )
                    conditional_similarity = f"{numerator / denom:.6f}" if denom else ""
                tumour_summary_rows.append({
                    "tumour_reference": ref_name, "reference_PSM_hypothesis": hypothesis,
                    "reference_PSM_note": psm_note, "tumour_reference_USI": ref_usi,
                    "islet_run": run, "islet_scan": scan,
                    "islet_native_scan": islet["native_scan"],
                    "islet_mgf_title": islet["mgf_title"],
                    "by_label_singly_charged_ion_n": len(all_shared_labels),
                    "same_theoretical_mz_ion_n": sum(
                        abs(ref_ions[label][0] - islet_ions[label][0]) <= TOL for label in all_shared_labels
                    ),
                    "same_mz_ions_observed_in_both_n": len(same_mass_observed),
                    "same_mz_observed_ions": ";".join(x[0] for x in same_mass_observed),
                    "form_shifted_ions": ";".join(shifted_labels),
                    "conditional_same_form_sqrt_intensity_cosine": conditional_similarity,
                    "islet_charge": islet_charge, "tumour_charge_from_PSM": ref_charge,
                    "observed_tumour_precursor_mz": "not present in archived PROXI JSON",
                    "interpretation": "per-ion descriptive comparison; modification state of primary source scan is unresolved",
                })
    write_tsv(args.out / "tumour_hypothesis_by_ion.tsv", tumour_match_rows)
    write_tsv(args.out / "tumour_hypothesis_summary.tsv", tumour_summary_rows)
    print(f"COMPLETE: {len(assignment_rows)} candidate-ion rows; {len(summary_rows)} candidate spectra; "
          f"{len(tumour_match_rows)} tumour-ion rows; {len(tumour_summary_rows)} tumour comparisons", flush=True)


if __name__ == "__main__":
    main()

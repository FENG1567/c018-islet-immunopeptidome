#!/usr/bin/env python3
"""C018 v4 plotting implementations for the other main and supplement figures."""

import argparse
import csv
import math
import os
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch


NAVY = "#29465B"
BLUE = "#4E7794"
TEAL = "#3F8980"
GOLD = "#B27C34"
INK = "#24313A"
MID = "#697984"
LIGHT = "#EEF2F4"
PALE_BLUE = "#EDF3F7"
PALE_TEAL = "#EAF4F1"
PALE_GOLD = "#F8F2E8"
GRID = "#D8E0E4"
WHITE = "#FFFFFF"

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans", "sans-serif"],
    "font.size": 7,
    "axes.labelsize": 7,
    "axes.titlesize": 8,
    "axes.titleweight": "bold",
    "axes.spines.right": False,
    "axes.spines.top": False,
    "axes.linewidth": 0.75,
    "axes.edgecolor": MID,
    "xtick.labelsize": 6.2,
    "ytick.labelsize": 6.2,
    "xtick.color": INK,
    "ytick.color": INK,
    "text.color": INK,
    "axes.labelcolor": INK,
    "legend.frameon": False,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "savefig.dpi": 300,
    "figure.facecolor": WHITE,
    "savefig.facecolor": WHITE,
})


def read_tsv(root: Path, relative: str):
    path = root / relative
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(enumerate(csv.DictReader(handle, delimiter="\t"), start=2))


def read_csv(root: Path, relative: str):
    path = root / relative
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(enumerate(csv.DictReader(handle), start=2))


def first_match(rows, **criteria):
    for line, row in rows:
        if all(row.get(key) == str(value) for key, value in criteria.items()):
            return line, row
    raise KeyError(f"No matching source row for {criteria}")


def as_bool(value):
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def add_data(out, panel, mark, group, value, unit, source_path, selector, notes=""):
    out.append({
        "panel": panel,
        "mark": mark,
        "group": group,
        "value": value,
        "unit": unit,
        "source_path": source_path,
        "source_row_or_selector": selector,
        "notes": notes,
    })


def figure_title(fig, title, subtitle=None):
    fig.text(0.045, 0.965, title, fontsize=12, fontweight="bold", va="top")
    if subtitle:
        fig.text(0.047, 0.915, subtitle, fontsize=7.2, color=MID, va="top")


def box(fig, x, y, w, h, title, body, edge, fill, title_size=8.0, body_size=6.5):
    patch = FancyBboxPatch(
        (x, y), w, h, boxstyle="round,pad=0.008,rounding_size=0.012",
        linewidth=0.8, edgecolor=edge, facecolor=fill,
        transform=fig.transFigure, zorder=1,
    )
    fig.patches.append(patch)
    fig.text(x + 0.014, y + h - 0.035, title, fontsize=title_size,
             fontweight="bold", color=edge, va="top")
    fig.text(x + 0.014, y + h - 0.076, body, fontsize=body_size,
             color=INK, va="top", linespacing=1.32)


def save_figure(fig, number, main_dir, supp_dir, data_rows):
    is_supp = str(number).startswith("S")
    folder = supp_dir if is_supp else main_dir
    folder.mkdir(parents=True, exist_ok=True)
    stem = f"C018_Figure_{number}_v4"
    fig.savefig(folder / f"{stem}.pdf", bbox_inches="tight")
    # Optional alternate-backend renders; the release bundle uses PDF and CSV.
    if os.environ.get("C018_FIGURE_ALT_EXPORTS") == "1":
        fig.savefig(folder / f"{stem}.svg", bbox_inches="tight")
        fig.savefig(folder / f"{stem}.tiff", dpi=300, bbox_inches="tight",
                    pil_kwargs={"compression": "tiff_lzw"})
        fig.savefig(folder / f"{stem}.png", dpi=300, bbox_inches="tight")
    with (folder / f"C018_Figure_{number}_figure_data_v4.csv").open(
            "w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=[
            "panel", "mark", "group", "value", "unit", "source_path",
            "source_row_or_selector", "notes",
        ])
        writer.writeheader()
        writer.writerows(data_rows)
    plt.close(fig)
    print(f"Exported Figure {number}", flush=True)


def figure_1(data_root, main_dir, supp_dir):
    data = []
    fig = plt.figure(figsize=(7.15, 5.55))
    figure_title(fig, "Candidate lists and non-malignant comparison settings",
                 "List membership, biological sampling units and evidence layers are shown separately.")

    # Candidate sets
    box(fig, 0.055, 0.755, 0.425, 0.125, "NC2025",
        "589 source entries  |  220 aeTSA candidates used for donor endpoints",
        BLUE, PALE_BLUE, body_size=6.6)
    box(fig, 0.520, 0.755, 0.425, 0.125, "PDAC source list",
        "517 exact cancer-restricted sequences  |  33 assessed in source CTL experiments",
        TEAL, PALE_TEAL, body_size=6.6)

    # Cohort cards. Same-organ islets occupy the upper, visually central row.
    box(fig, 0.055, 0.492, 0.425, 0.208, "Archived pancreatic islets  ·  same organ",
        "7 archived acquisitions; 5 tissue-resolved, including 3 acquisitions\n"
        "from two two-field-typed donors\n"
        "Source study: cytokine-treated islets (IFNγ, IL-1β and TNFα, 48 h)\n"
        "W632/R369 acquisition: donor fingerprint supports R369; tissue link unresolved",
        TEAL, PALE_TEAL, body_size=6.25)
    box(fig, 0.520, 0.492, 0.425, 0.208, "Independent primary islets  ·  same organ",
        "2 untreated preparations × 2 technical acquisitions\n"
        "HLA-A*02 and HLA-A*03 family-level annotations\n"
        "Two preparations, not four biological replicates",
        TEAL, PALE_TEAL, body_size=6.25)
    box(fig, 0.055, 0.237, 0.425, 0.202, "Fibrotic-lung comparator",
        "12 acquisitions; 6 passed the reference gate, all from IPF samples\n"
        "Sample-to-donor/HLA mapping is inferred and was tested in both directions\n"
        "ILMILQPQL is retained as an unresolved signal",
        BLUE, PALE_BLUE, body_size=6.2)
    box(fig, 0.520, 0.237, 0.425, 0.202, "Skin reference  ·  background resource only",
        "25,518 source-table rows\n"
        "No skin endpoint or candidate overlap was assessed\n"
        "Row count is not a unique-peptide count",
        MID, LIGHT, body_size=6.4)

    # Evidence levels shown as a map rather than a serial pipeline.
    fig.text(0.055, 0.180, "EVIDENCE LAYERS", fontsize=7.0, color=MID, fontweight="bold")
    layer_labels = ["Comet spectrum search", "Shuffled / entrapment controls",
                    "HLA eligibility", "Targeted spectrum review"]
    layer_colors = [NAVY, BLUE, TEAL, GOLD]
    x0, gap, w = 0.055, 0.016, 0.211
    for i, (label, color) in enumerate(zip(layer_labels, layer_colors)):
        x = x0 + i * (w + gap)
        fig.patches.append(FancyBboxPatch(
            (x, 0.095), w, 0.060, boxstyle="round,pad=0.006,rounding_size=0.010",
            linewidth=0.65, edgecolor=color, facecolor=WHITE,
            transform=fig.transFigure,
        ))
        fig.text(x + w / 2, 0.125, label, fontsize=6.1, color=color,
                 ha="center", va="center", fontweight="bold")
    fig.text(0.055, 0.045,
             "Acquisitions are not independent biological units. Archived cytokine-treated islets and untreated external islets are distinct settings.",
             fontsize=6.15, color=MID)

    # Every displayed count has its own traceable row or source selector.
    ext = read_tsv(data_root, "S4/external_detection_summary.tsv")
    skin = read_tsv(data_root, "S6/skin_reference_context.tsv")
    local = read_tsv(data_root, "S1/local_run_inventory.tsv")
    lung = read_tsv(data_root, "S1/lung_acquisition_eligibility.tsv")
    add_data(data, "a", "source list size", "NC2025", 589, "candidate entries",
             "C018 v4R1 manuscript; author-supplied source-list aggregate",
             "all source rows; total list size recorded in manuscript source inventory",
             "The donor endpoint uses its 220-entry aeTSA subset.")
    for label, criteria, note in [
        ("NC2025 aeTSA", {"list": "NC2025_aeTSA", "stratum": "all", "recipe": "legacy_recipe"},
         "Endpoint list subset."),
        ("PDAC exact sequences", {"list": "PDAC2025", "stratum": "all", "recipe": "legacy_recipe"},
         "Source cancer-restricted list after exact-sequence definition."),
    ]:
        line, row = first_match(ext, **criteria)
        value = row["candidate_sequences"]
        add_data(data, "a", "candidate count", label, value, "sequences",
                 "S4/external_detection_summary.tsv",
                 line, note)
    add_data(data, "a", "CTL-tested subset", "PDAC", 33, "sequences",
             "C018 v4R1 manuscript; author-supplied source-study overlap aggregate", "all assessed sequences",
             "Exact overlap with the 517-sequence PDAC list.")
    add_data(data, "b", "archived acquisitions", "cytokine-treated islets", len(local), "acquisitions",
             "S1/local_run_inventory.tsv",
             "all run records", "Five acquisitions entered the tissue-resolved typed scenario.")
    add_data(data, "b", "tissue-resolved acquisitions", "archived islets", 5, "acquisitions",
             "C018 v4R1 manuscript; tissue-resolved cohort summary", "5 acquisitions",
             "One low-recovery control and the tissue-ambiguous W632 acquisition excluded.")
    add_data(data, "c", "biological preparations", "PXD045211", 2, "preparations",
             "S1/pxd045211_preparation_mapping.tsv", "islet22/islet27; files 355/356 and 358/359",
             "Untreated; two technical acquisitions per preparation.")
    add_data(data, "c", "lung cohort size", "fibrotic lung", sum(1 for _, r in lung if r["tissue"].lower() == "lung"), "acquisitions",
             "S1/lung_acquisition_eligibility.tsv", "all rows with tissue=lung",
             "All six gate-passing samples are from IPF; mapping is inferred.")
    add_data(data, "c", "gate-passing lung acquisitions", "IPF", sum(1 for _, r in lung if r["tissue"].lower() == "lung" and r["gate_status"].lower() == "pass"), "acquisitions",
             "S1/lung_acquisition_eligibility.tsv", "tissue=lung and gate_status=pass",
             "The four non-IPF controls did not pass the gate.")
    line, row = first_match(skin, source_description="Published skin HLA-I peptide table cited in the manuscript")
    add_data(data, "d", "source table rows", "skin context", row["source_rows"], "rows",
             "S6/skin_reference_context.tsv", line,
             "Background resource only; no skin endpoint or candidate overlap; row count is not a unique-peptide count.")
    for i, label in enumerate(layer_labels, start=1):
        add_data(data, "e", "evidence layer", label, "present", "analysis layer",
                 "C018 v4R1 manuscript Methods", f"four-layer evidence description {i}",
                 "These layers are not interpreted as a serial validation pipeline.")
    save_figure(fig, 1, main_dir, supp_dir, data)


def figure_2(data_root, main_dir, supp_dir):
    data = []
    fig = plt.figure(figsize=(7.15, 4.85))
    figure_title(fig, "HLA eligibility bounds the archived endpoint counts",
                 "Eligibility denominators are separate from qualifying calls; unresolved lung signals remain unclassified.")

    ax = fig.add_axes([0.16, 0.285, 0.39, 0.53])
    ax.set_title("a  Archived islet aeTSA candidates", loc="left", pad=9)
    scenarios = [
        ("Historical ligand-HLA union", 174, 49, "archive_data/S2/aetsa_endpoint_recomputed.tsv", 2, 3),
        ("Exact two-field; R369 retained", 93, 9, "archive_data/S2/aetsa_exact_typed_endpoint_recomputed.tsv", 2, 3),
        ("Tissue-resolved exact two-field", 49, 1, "local_conservative_verified_endpoints.tsv", 3, 4),
    ]
    y_centres = [2.45, 1.45, 0.45]
    colors = {"any donor": BLUE, "two donors": TEAL}
    for (label, any_n, both_n, source, any_line, both_line), y in zip(scenarios, y_centres):
        for sub_y, endpoint, value, source_line, color in [
            (y + 0.19, "any donor", any_n, any_line, colors["any donor"]),
            (y - 0.19, "two donors", both_n, both_line, colors["two donors"]),
        ]:
            ax.barh(sub_y, value, height=0.29, color=color, alpha=0.92, zorder=3)
            ax.text(value + (2.4 if value > 5 else 1.0), sub_y, str(value), va="center", fontsize=6.7, fontweight="bold")
            add_data(data, "a", "eligible candidates", f"{label}; {endpoint}", value, "of 220",
                     source, source_line,
                     "Zero candidates met the frozen ≥2-spectrum donor endpoint in this scenario.")
    ax.set_yticks(y_centres, [x[0] for x in scenarios])
    ax.set_xlim(0, 205)
    ax.set_ylim(-0.1, 2.9)
    ax.set_xlabel("Candidates with HLA-compatible donor(s)")
    ax.grid(axis="x", color=GRID, linewidth=0.45, zorder=0)
    ax.tick_params(axis="y", length=0)
    ax.text(0.5, -0.22, "0 qualifying aeTSA calls across all three HLA scenarios",
            transform=ax.transAxes, ha="center", va="top", fontsize=6.35,
            color=INK, fontweight="bold")
    ax.legend(handles=[
        plt.Line2D([0], [0], color=BLUE, lw=5, label="Any donor"),
        plt.Line2D([0], [0], color=TEAL, lw=5, label="At least two donors"),
    ], loc="lower right", bbox_to_anchor=(1.0, 1.01), ncol=2, fontsize=6.1)

    endpoint_rows = read_tsv(data_root, "frozen/STEP5_ENDPOINTS.tsv")
    lung_cells = [
        ("S1 main inferred mapping", "any", 71, 1, "unresolved", 2),
        ("S1 main inferred mapping", "both", 50, 1, "0 qualified", 3),
        ("S2 reversed inferred mapping", "any", 51, 1, "0 qualified", 4),
        ("S2 reversed inferred mapping", "both", 36, 0, "0 qualified", 5),
    ]
    fig.text(0.60, 0.815, "b  Lung mapping sensitivity", fontsize=8.0, fontweight="bold")
    tile_x, tile_y, tile_w, tile_h = 0.60, 0.645, 0.165, 0.137
    for i, (scenario, endpoint, denominator, code, audit_call, source_line) in enumerate(lung_cells):
        row_i, source_row = first_match(endpoint_rows, scenario="S1_passing" if scenario.startswith("S1") else "S2_reversed",
                                        cohort="lung", endpoint=endpoint)
        # Guard denominators/codes against the archived frozen endpoint table.
        assert int(source_row["denominator"]) == denominator and int(source_row["detected"]) == code
        col = i % 2
        row = i // 2
        x, y = tile_x + col * (tile_w + 0.02), tile_y - row * (tile_h + 0.025)
        fill = PALE_GOLD
        fig.patches.append(FancyBboxPatch(
            (x, y), tile_w, tile_h, boxstyle="round,pad=0.007,rounding_size=0.009",
            linewidth=0.65, edgecolor=GOLD, facecolor=fill,
            transform=fig.transFigure,
        ))
        fig.text(x + 0.009, y + tile_h - 0.028,
                 f"{scenario.split()[0]} · {endpoint}", fontsize=5.9, color=GOLD, fontweight="bold")
        fig.text(x + 0.009, y + 0.060,
                 f"eligible {denominator} · archive {code}", fontsize=5.5, color=INK)
        fig.text(x + 0.009, y + 0.028,
                 f"donor audit: {audit_call}", fontsize=5.5, color=INK)
        add_data(data, "b", "frozen lung endpoint row", f"{scenario}; {endpoint}",
                 f"{denominator}/{code}", "eligible/coded count",
                 "out/step5/STEP5_ENDPOINTS.tsv", row_i,
                 "The coded signal is ILMILQPQL, retained as unresolved; it is not classified positive or negative.")
        audit_source = "archive_data/S2/lung_correction_summary.tsv" if endpoint == "both" else "archive_data/S2/ILMILQPQL_S1_S2_endpoint_audit.tsv"
        audit_rows = read_tsv(data_root, audit_source)
        audit_idx = 2 if scenario.startswith("S1") else 3
        add_data(data, "b", "donor-mapped endpoint audit", f"{scenario}; {endpoint}", audit_call,
                 "qualified call status", audit_source, audit_idx,
                 "A zero qualified call is not a biological absence claim; ILMILQPQL stays unresolved.")
    fig.text(0.60, 0.305,
             "ILMILQPQL remains unresolved; archived codes are not positive/negative classifications.\n"
             "The S1 both-donor count is 0 because PF7 contributes one spectrum.\n"
             "S1/S2 sample-to-HLA aliases are inferred.",
             fontsize=6.35, color=INK, linespacing=1.45, va="top")
    add_data(data, "b", "mapping status", "S1 and S2", "inferred", "sample-to-HLA linkage",
             "archive_data/S2/ILMILQPQL_donor_HLA_S1_S2.tsv", "rows 2–7",
             "Under S1, PF5 has 2 spectra and PF7 has 1; under S2, PF5 maps to an untyped donor.")
    add_data(data, "b", "classification", "ILMILQPQL", "unresolved signal", "evidence status",
             "archive_data/S2/ILMILQPQL_S1_S2_endpoint_audit.tsv", "rows 2–3",
             "Excluded from positive and negative classifications.")
    save_figure(fig, 2, main_dir, supp_dir, data)


def figure_3(data_root, main_dir, supp_dir):
    data = []
    jy_rel = "S3/jy_aligned_recovery.tsv"
    jy_rows = read_tsv(data_root, jy_rel)
    groups = defaultdict(dict)

    # JY: one acquisition at 1:16, against source-identified reference peptides.
    for line, row in jy_rows:
        if (row["mode"] == "aligned_psm_and_peptide_q_plus1" and
                row["dilution_factor"] == "16" and row["reference_scans_min"] == "1" and
                row["test_scans_min"] in ("1", "2")):
            threshold = int(row["test_scans_min"])
            label = "JY eluate 1:16"
            groups[label][threshold] = {
                "n": int(row["recovered_n"]), "N": int(row["reference_n"]),
                "rate": float(row["recovery"]), "low": float(row["descriptive_wilson_low"]),
                "high": float(row["descriptive_wilson_high"]), "line": line,
                "source": "tissue_calibration/jy_aligned/jy_aligned_recovery.tsv",
                "notes": "One DDA acquisition per dilution; conditional on the source-reference set.",
            }

    # Archived cytokine-treated preparation, low-intensity reference quartile.
    stern_rel = "S3/HP20289_paired_recovery.tsv"
    stern_rows = read_tsv(data_root, stern_rel)
    for line, row in stern_rows:
        if row["panel"] != "searchable_low_intensity" or row["method"] not in ("c018_test_ge1_scan", "c018_test_ge2_scans"):
            continue
        threshold = 1 if row["method"] == "c018_test_ge1_scan" else 2
        src = "Set1" if "Set1" in row["source_run"] else "Set2"
        tst = "Set1" if "Set1" in row["test_run"] else "Set2"
        label = f"HP20289 low-intensity {src}→{tst}"
        groups[label][threshold] = {
            "n": int(row["recovered"]), "N": int(row["reference_n"]),
            "rate": float(row["recovery"]), "low": float(row["descriptive_wilson_low"]),
            "high": float(row["descriptive_wilson_high"]), "line": line,
            "source": stern_rel,
            "notes": "Two acquisitions from one HP20289 donor; source eluate identity is not verified.",
        }

    # Untreated primary islets: four paired directions, each under two search settings.
    pxd_files = sorted((data_root / "S3").glob("pxd045211_*.tsv"))
    for path in pxd_files:
        recipe = "archived settings" if "raw_validation" in path.name else "metadata-adjusted settings"
        protocol = "raw_validation" if "raw_validation" in path.name else "pride_recipe"
        for line, row in read_tsv(data_root, str(path.relative_to(data_root))):
            if row["source_tic_quartile"] != "1" or row["method"] not in ("c018_one_scan", "c018_two_scans"):
                continue
            threshold = 1 if row["method"] == "c018_one_scan" else 2
            prep = "islet22" if "islet22" in path.name else "islet27"
            src = row["source_run"].rsplit("_", 1)[-1]
            tst = row["test_run"].rsplit("_", 1)[-1]
            label = f"PXD045211 {prep} {src}→{tst} · {recipe}"
            groups[label][threshold] = {
                "n": int(row["recovered"]), "N": int(row["reference_positives"]),
                "rate": float(row["recovery"]), "low": float(row["descriptive_wilson_low"]),
                "high": float(row["descriptive_wilson_high"]), "line": line,
                "source": f"S3/{path.name}",
                "notes": "Lowest source fragment-TIC quartile; technical acquisition directions within two preparations.",
            }

    ordered = ["JY eluate 1:16"] + sorted(k for k in groups if k.startswith("HP20289")) + sorted(k for k in groups if k.startswith("PXD045211"))
    if any(set(groups[label]) != {1, 2} for label in ordered):
        raise ValueError("Each recovery comparison must have both one- and two-spectrum values")

    fig, ax = plt.subplots(figsize=(7.15, 6.25))
    figure_title(fig, "Conditional recovery of source-identified reference ligands",
                 "The same source-defined reference set is assessed under one- and two-spectrum rules; cohorts remain separate.")
    yvals = list(range(len(ordered)))[::-1]
    ax.set_position([0.335, 0.18, 0.54, 0.67])
    for y, label in zip(yvals, ordered):
        one, two = groups[label][1], groups[label][2]
        ax.plot([one["rate"] * 100, two["rate"] * 100], [y, y], color=GRID, lw=1.0, zorder=1)
        for threshold, item, color, offset in [
            (1, one, BLUE, -0.065), (2, two, TEAL, 0.065),
        ]:
            value = item["rate"] * 100
            lo = (item["rate"] - item["low"]) * 100
            hi = (item["high"] - item["rate"]) * 100
            ax.errorbar(value, y + offset, xerr=[[lo], [hi]], fmt="o", color=color,
                        ecolor=color, capsize=1.7, markersize=3.4, linewidth=0.85, zorder=3)
            add_data(data, "a", f"≥{threshold} distinct spectrum(s)", label,
                     f"{item['n']}/{item['N']} ({value:.3f})", "recovery fraction",
                     item["source"], item["line"], item["notes"])
    ax.set_yticks(yvals, ordered)
    ax.set_xlim(0, 82)
    ax.set_ylim(-0.8, len(ordered) - 0.15)
    ax.set_xlabel("Reference ligands recovered (%)")
    ax.grid(axis="x", color=GRID, linewidth=0.45, zorder=0)
    ax.tick_params(axis="y", length=0, pad=4)
    ax.legend(handles=[
        plt.Line2D([0], [0], marker="o", color="none", markerfacecolor=BLUE,
                   markeredgecolor=BLUE, markersize=4, label="≥1 spectrum"),
        plt.Line2D([0], [0], marker="o", color="none", markerfacecolor=TEAL,
                   markeredgecolor=TEAL, markersize=4, label="≥2 spectra"),
    ], loc="lower right", bbox_to_anchor=(1.0, 1.01), ncol=2, fontsize=6.2)
    fig.text(0.045, 0.103,
             "Error bars are descriptive Wilson 95% intervals. HP20289 and PXD045211 rows are within-donor/preparation acquisition pairs.\n"
             "JY 1:16 is eluate dilution in a single DDA acquisition and is not an islet candidate-specific detection limit.",
             fontsize=6.15, color=MID, linespacing=1.4)
    save_figure(fig, 3, main_dir, supp_dir, data)


def figure_4(data_root, main_dir, supp_dir):
    endpoint_rel = "S2/archived_exact_panel_candidate_donor_endpoints.tsv"
    run_rel = "S2/archived_exact_panel_candidate_run_matrix.tsv"
    single_rel = "S2/archived_exact_panel_compatible_single_spectra.tsv"
    endpoints = read_tsv(data_root, endpoint_rel)
    run_matrix = read_tsv(data_root, run_rel)
    single_rows = read_tsv(data_root, single_rel)

    # This table carries one row per exact PDAC-list sequence; the row count is
    # checked so that a truncated copy cannot silently change the denominator.
    if len(endpoints) != 517:
        raise ValueError(f"Expected 517 exact PDAC candidates, found {len(endpoints)}")
    compatible_any_all = [(line, row) for line, row in endpoints
                          if int(row["compatible_donor_n"]) >= 1]
    compatible_two_all = [(line, row) for line, row in endpoints
                          if int(row["compatible_donor_n"]) >= 2]
    native = [(line, row) for line, row in endpoints
              if as_bool(row["native_unmodified_scope"])]
    compatible_any = [(line, row) for line, row in native
                      if int(row["compatible_donor_n"]) >= 1]
    compatible_two = [(line, row) for line, row in native
                      if int(row["compatible_donor_n"]) >= 2]
    any_spectrum_endpoint = [(line, row) for line, row in compatible_any
                              if as_bool(row["any_endpoint"])]
    two_donor_endpoint = [(line, row) for line, row in compatible_any
                          if as_bool(row["two_endpoint"])]

    # Cross-reference the candidate/run matrix to the single-spectrum table.
    # The matrix keeps acquisition opportunities distinct from candidate counts.
    typed_runs = [(line, row) for line, row in run_matrix
                  if as_bool(row["native_unmodified_scope"])
                  and as_bool(row["exact_typed_compatible"])]
    typed_single_runs = [(line, row) for line, row in typed_runs
                         if int(row["accepted_exact_scan_n"] or 0) == 1]
    exported_single_runs = [(line, row) for line, row in single_rows
                            if as_bool(row["native_unmodified_scope"])
                            and as_bool(row["exact_typed_compatible"])
                            and int(row["accepted_exact_scan_n"] or 0) == 1]
    matrix_single_keys = {(row["peptide"], row["run"]) for _, row in typed_single_runs}
    exported_single_keys = {(row["peptide"], row["run"]) for _, row in exported_single_runs}
    if matrix_single_keys != exported_single_keys:
        raise ValueError("Single-spectrum audit rows disagree with the candidate/run matrix")

    expected = (len(endpoints), len(compatible_any_all), len(compatible_two_all),
                len(compatible_any), len(compatible_two), len(exported_single_runs),
                len(any_spectrum_endpoint), len(two_donor_endpoint))
    if expected != (517, 140, 19, 131, 19, 3, 0, 0):
        raise ValueError(f"Frozen same-organ audit counts changed: {expected}")

    data = []
    fig = plt.figure(figsize=(7.15, 4.45))
    figure_title(fig, "Same-organ PDAC–islet screen under exact donor and spectrum rules",
                 "The comparison uses non-malignant islets; no PDAC tumour specimens were analyzed in this arm.")
    ax = fig.add_axes([0.10, 0.30, 0.45, 0.53])
    stages = [
        ("Exact PDAC source list", len(endpoints), "all candidate rows"),
        ("≥1 compatible typed donor", len(compatible_any_all), "compatible_donor_n≥1"),
        ("≥2 compatible typed donors", len(compatible_two_all), "compatible_donor_n≥2"),
    ]
    yvals = list(range(len(stages)))[::-1]
    for y, (label, value, selector), color in zip(yvals, stages, [MID, BLUE, TEAL]):
        ax.barh(y, value, height=0.48, color=color, alpha=0.92, zorder=2)
        ax.text(value + 6, y, f"{value} / 517", va="center",
                fontsize=6.8, fontweight="bold", color=INK)
        add_data(data, "a", "HLA coverage", label, value, "sequences",
                 endpoint_rel, selector,
                 "Full fixed library sequence set; donor compatibility is an eligibility measure.")
    ax.set_title("a  HLA compatibility across the 517-sequence list", loc="left", pad=7)
    ax.set_yticks(yvals, [s[0] for s in stages])
    ax.set_xlim(0, 595)
    ax.set_xlabel("Exact PDAC-list sequences")
    ax.grid(axis="x", color=GRID, linewidth=0.45, zorder=0)
    ax.tick_params(axis="y", length=0, pad=4)

    tile_specs = [
        (0.61, 0.58, 0.17, len(compatible_any), "native unmodified; ≥1 typed donor", BLUE,
         endpoint_rel, "native_unmodified_scope=True; compatible_donor_n≥1", "sequences"),
        (0.80, 0.58, 0.17, len(compatible_two), "native unmodified; ≥2 typed donors", TEAL,
         endpoint_rel, "native_unmodified_scope=True; compatible_donor_n≥2", "sequences"),
        (0.61, 0.37, 0.17, len(exported_single_runs), "compatible single-spectrum records", GOLD,
         single_rel, "all exported rows; each has one accepted exact scan", "candidate–run records"),
        (0.80, 0.37, 0.17, len(any_spectrum_endpoint), "within-donor ≥2-spectrum endpoint", GOLD,
         endpoint_rel, "native_unmodified_scope=True; compatible_donor_n≥1; any_endpoint=True", "sequences"),
        (0.61, 0.16, 0.36, len(two_donor_endpoint), "two-donor ≥2-spectrum endpoint", GOLD,
         endpoint_rel, "native_unmodified_scope=True; compatible_donor_n≥1; two_endpoint=True", "sequences"),
    ]
    for x, y, width, value, label, color, source, selector, unit in tile_specs:
        box(fig, x, y, width, 0.15, str(value), label, color,
            PALE_GOLD if color == GOLD else (PALE_TEAL if color == TEAL else PALE_BLUE),
            title_size=9.0, body_size=5.0)
        add_data(data, "b", "same-organ screening result", label, value, unit,
                 source,
                 selector,
                 "Native-unmodified fixed-C scope is a search boundary, not a biological detection limit.")

    for line, row in exported_single_runs:
        add_data(data, "b", "single-spectrum candidate/run", row["peptide"],
                 row["accepted_exact_scan_n"], "accepted spectrum",
                 single_rel, line,
                 f"Compatible donor/run={row['donor']} / {row['run']}; not a two-spectrum endpoint.")
    add_data(data, "b", "eligible typed candidate/run opportunities", "native typed matrix",
             len(typed_runs), "candidate–run rows", run_rel,
             "native_unmodified_scope=True and exact_typed_compatible=True",
             "Acquisition rows are not independent biological replicates.")

    fig.text(0.10, 0.14,
             "A qualifying donor contributes at least two distinct spectra. The two-donor endpoint requires two qualifying donors.",
             fontsize=5.9, color=MID, va="top")
    fig.text(0.10, 0.075,
             "Candidate-level endpoint counts use exact sequence matching and the prespecified ≥2-distinct-spectrum rule; acquisition and donor counts are reported in Table S2.",
             fontsize=5.8, color=MID, va="top")
    save_figure(fig, 4, main_dir, supp_dir, data)


def figure_5(data_root, main_dir, supp_dir):
    evidence_rel = "S5/SIASAISMK_islet_spectrum_evidence.tsv"
    chemistry_rel = "S5/SIASAISMK_precursor_chemistry.tsv"
    tumour_rel = "S5/SIASAISMK_tumor_reference_evidence.tsv"
    evidence = read_tsv(data_root, evidence_rel)
    chemistry = read_tsv(data_root, chemistry_rel)
    tumour = read_tsv(data_root, tumour_rel)
    if len(evidence) != 4:
        raise ValueError(f"Expected four archived SIASAISMK islet spectra, found {len(evidence)}")

    evidence.sort(key=lambda item: (
        0 if item[1]["preparation"] == "islet22" else 1,
        item[1]["technical_acquisition"],
    ))
    chem_islet = [(line, row) for line, row in chemistry
                  if row["reference"] == "PXD045211 islet"]
    if not chem_islet:
        raise ValueError("No islet precursor chemistry row found")
    chem_line, chem_row = chem_islet[0]
    unmodified_mz = float(chem_row["theoretical_unmodified_mz_z1"])
    oxidized_mz = float(chem_row["theoretical_M8_oxidized_mz_z1"])
    primary_line, primary = first_match(tumour, reference_id="primary_C3L01031")
    pdo_line, pdo = first_match(tumour, reference_id="PDO_PANFR0151")

    data = []
    fig = plt.figure(figsize=(7.15, 5.15))
    figure_title(fig, "A recurrent islet signal remains unresolved at sequence and HLA level",
                 "SIASAISMK was assigned in both preparations. The comparison includes two tumour-source spectra.")

    # PSM score panel: two technical acquisitions are shown within each preparation.
    ax_score = fig.add_axes([0.12, 0.52, 0.37, 0.29])
    ax_score.set_title("a  Historical-recipe islet PSMs", loc="left", pad=7)
    yvals = list(range(len(evidence)))[::-1]
    colors = {"islet22": BLUE, "islet27": TEAL}
    for y, (line, row) in zip(yvals, evidence):
        score = float(row["xcorr"])
        ax_score.scatter([score], [y], color=colors[row["preparation"]], s=23, zorder=3)
        ax_score.text(score + 0.018, y, f"{score:.3f}  ·  q={float(row['peptide_q']):.4g}",
                      fontsize=5.8, va="center", color=INK)
        add_data(data, "a", "archived PSM XCorr", row["preparation"], score, "XCorr",
                 evidence_rel, line,
                 f"{row['technical_acquisition']}; peptide q={row['peptide_q']}; one acquisition within preparation.")
    ax_score.set_yticks(yvals, [f"{r['preparation']} · {r['technical_acquisition'].replace('technical acquisition ', 'T')}"
                                for _, r in evidence])
    ax_score.set_xlim(1.6, 2.65)
    ax_score.set_xlabel("Archived XCorr")
    ax_score.grid(axis="x", color=GRID, linewidth=0.45, zorder=0)
    ax_score.tick_params(axis="y", length=0, pad=4)

    # Precursor masses are plotted in the charge state actually observed in islets.
    ax_mass = fig.add_axes([0.59, 0.52, 0.34, 0.29])
    ax_mass.set_title("b  Islet precursor chemistry (z = 1)", loc="left", pad=7)
    obs_values = []
    for y, (line, row) in zip(yvals, evidence):
        obs = float(row["precursor_mz"])
        obs_values.append(obs)
        ax_mass.scatter([obs], [y], color=colors[row["preparation"]], s=23, zorder=3)
        add_data(data, "b", "observed islet precursor m/z", row["preparation"], obs, "m/z",
                 evidence_rel, line,
                 f"{row['technical_acquisition']}; z={row['charge']}; precursor remains unmodified in archived peptide field.")
    ax_mass.axvline(unmodified_mz, color=BLUE, lw=1.2, zorder=1)
    ax_mass.axvline(oxidized_mz, color=GOLD, lw=1.0, linestyle="--", zorder=1)
    ax_mass.text(unmodified_mz + 0.16, 3.72, f"unmodified {unmodified_mz:.4f}",
                 fontsize=5.55, color=BLUE, va="top")
    ax_mass.text(oxidized_mz - 0.18, 3.72, f"M8 {oxidized_mz:.4f}",
                 fontsize=5.55, color=GOLD, va="top", ha="right")
    ax_mass.set_xlim(min(obs_values) - 0.35, oxidized_mz + 0.55)
    ax_mass.set_ylim(-0.45, 3.95)
    ax_mass.set_yticks([])
    ax_mass.set_xlabel("Precursor m/z")
    ax_mass.grid(axis="x", color=GRID, linewidth=0.4, zorder=0)
    add_data(data, "b", "theoretical unmodified precursor", "islet z=1", unmodified_mz, "m/z",
             chemistry_rel, chem_line, "Theoretical [M+H]+ from the frozen precursor-chemistry audit.")
    add_data(data, "b", "theoretical M8-oxidized precursor", "islet z=1", oxidized_mz, "m/z",
             chemistry_rel, chem_line, "Theoretical [M+H]+ from the frozen precursor-chemistry audit.")

    # Keep tumour-source provenance separate: the primary scan's modification is unknown,
    # while the PDO source row explicitly records M8 oxidation.
    box(fig, 0.07, 0.20, 0.42, 0.205, "Primary PDAC reference · C3L-01031 scan 16793",
        f"Source PSM: {primary['source_PSM_sequence']}, z={primary['charge']}\n"
        "The archived row omits the modification field.\n"
        "Original precursor m/z and Met8 form are unresolved.", GOLD, PALE_GOLD, title_size=7.0, body_size=5.8)
    box(fig, 0.53, 0.20, 0.42, 0.205, "PDAC organoid reference · PANFR0151 scan 16495",
        f"Source PSM: {pdo['source_PSM_sequence']}, z={pdo['charge']}\n"
        "Met8 oxidation is explicit in the source row.\n"
        "The islet assignment uses an unmodified precursor.",
        BLUE, PALE_BLUE, title_size=7.0, body_size=5.8)
    add_data(data, "c", "tumour-reference provenance", "primary PDAC C3L-01031", "modification unresolved", "source status",
             tumour_rel, primary_line,
             f"Exact source status: {primary['modification_status_in_exact_source_row']}; observed precursor={primary['observed_precursor_mz_in_archive']}.")
    add_data(data, "c", "tumour-reference provenance", "PDO PANFR0151", "M8 oxidation explicit", "source status",
             tumour_rel, pdo_line,
             f"Source PSM={pdo['source_PSM_sequence']}; charge={pdo['charge']}; separate PDO spectrum.")

    fig.text(0.07, 0.105,
             "Four PSMs represent two technical acquisitions within each of two biological preparations.",
             fontsize=6.0, color=MID)
    fig.text(0.07, 0.062,
             "Sequence identity and HLA restriction remain unresolved.",
             fontsize=6.1, color=INK, fontweight="bold")
    save_figure(fig, 5, main_dir, supp_dir, data)


def figure_s1(data_root, main_dir, supp_dir):
    rel = "S3/jy_aligned_recovery.tsv"
    rows = read_tsv(data_root, rel)
    chosen = [(line, row) for line, row in rows if
              row["mode"] == "aligned_psm_and_peptide_q_plus1" and
              row["reference_scans_min"] == "1" and row["test_scans_min"] in ("1", "2")]
    chosen.sort(key=lambda item: (int(item[1]["dilution_factor"]), int(item[1]["test_scans_min"])))
    data = []
    fig, ax = plt.subplots(figsize=(6.7, 3.65))
    figure_title(fig, "JY source-reference recovery across eluate dilution",
                 "Recovery is conditional on peptides identified in the undiluted source reference set (n = 4,171).")
    colors = {1: BLUE, 2: TEAL}
    xloc = {"1": 0, "4": 1, "16": 2}
    for threshold in (1, 2):
        series = [(line, row) for line, row in chosen if int(row["test_scans_min"]) == threshold]
        xs = [xloc[row["dilution_factor"]] for _, row in series]
        ys = [float(row["recovery"]) for _, row in series]
        ax.plot(xs, ys, marker="o", color=colors[threshold], lw=1.3,
                markersize=4, label=f"≥{threshold} spectrum" + ("s" if threshold == 2 else ""), zorder=3)
        for x, (line, row), y in zip(xs, series, ys):
            lo = (y - float(row["descriptive_wilson_low"]))
            hi = (float(row["descriptive_wilson_high"]) - y)
            ax.errorbar(x, y, yerr=[[lo], [hi]], fmt="none", ecolor=colors[threshold],
                        capsize=2, linewidth=0.8, zorder=2)
            ax.annotate(f"{int(row['recovered_n'])}/{int(row['reference_n'])}\n{y:.1%}",
                        (x, y), xytext=(0, 6 if threshold == 1 else -25),
                        textcoords="offset points", ha="center", fontsize=5.8,
                        color=colors[threshold])
            add_data(data, "a", f"≥{threshold} spectrum(s)", f"1:{row['dilution_factor']} dilution",
                     f"{row['recovered_n']}/{row['reference_n']} ({row['recovery']})", "recovery fraction",
                     "S3/jy_aligned_recovery.tsv", line,
                     "One DDA acquisition per dilution; source-reference recovery, not candidate-specific sensitivity.")
    ax.set_xticks([0, 1, 2], ["Undiluted", "1:4", "1:16"])
    ax.set_ylim(0, 1.12)
    ax.set_xlim(-0.2, 2.2)
    ax.set_ylabel("Reference ligands recovered")
    ax.set_xlabel("JY eluate dilution")
    ax.grid(axis="y", color=GRID, linewidth=0.45)
    ax.legend(loc="lower left", fontsize=6.3)
    fig.text(0.05, 0.085,
             "JY HLA background and 400–650 m/z acquisition window differ from the islet cohorts; no tissue detection-limit transfer is made.",
             fontsize=6.1, color=MID)
    save_figure(fig, "S1", main_dir, supp_dir, data)


def figure_s2(data_root, main_dir, supp_dir):
    rel = "S2/ILMILQPQL_donor_HLA_S1_S2.tsv"
    rows = read_tsv(data_root, rel)
    data = []
    fig, axes = plt.subplots(1, 2, figsize=(7.15, 3.75), sharex=True)
    figure_title(fig, "ILMILQPQL remains an unresolved signal under both lung mappings",
                 "The mapping is inferred; bars are distinct MS/MS events within each sample.")
    status_color = {"True": TEAL, "False": MID, "unknown_unmapped": GOLD}
    for ax, scenario, title in zip(axes, ("S1_passing", "S2_reversed"),
                                   ("a  Main inferred mapping (S1)", "b  Reversed inferred mapping (S2)")):
        part = [(line, row) for line, row in rows if row["scenario"] == scenario]
        part.sort(key=lambda x: ["PF2", "PF5", "PF7"].index(x[1]["PF_label"]))
        y = list(range(len(part)))[::-1]
        ax.set_title(title, loc="left", pad=9)
        for yy, (line, row) in zip(y, part):
            n = int(row["spectra_count_frozen"])
            color = status_color[row["HLA_compatible"]]
            ax.barh(yy, n, height=0.38, color=color, alpha=0.92, zorder=3)
            ax.text(n + 0.08, yy, f"{n}  {row['mapped_donor']}", va="center", fontsize=6.0, color=INK)
            add_data(data, "a/b", "spectra per sample", f"{scenario}; {row['PF_label']}", n,
                     "distinct MS/MS events", rel,
                     line,
                     f"Mapped donor={row['mapped_donor']}; HLA={row['HLA_compatible']}; source alias verified={row['raw_sample_alias_verified']}.")
        ax.axvline(2, color=GOLD, linestyle=":", lw=0.9, zorder=1)
        ax.set_yticks(y, [f"{r['PF_label']}  ·  {r['HLA_compatible']}" for _, r in part])
        ax.set_xlim(0, 3.45)
        ax.set_xticks([0, 1, 2, 3])
        ax.set_xlabel("Distinct spectra in sample")
        ax.grid(axis="x", color=GRID, linewidth=0.4, zorder=0)
        ax.tick_params(axis="y", length=0)
    fig.subplots_adjust(left=0.13, right=0.98, top=0.78, bottom=0.24, wspace=0.35)
    fig.text(0.055, 0.103,
             "A 2-spectrum same-donor threshold is met by PF5 only in S1; S2 maps PF5 to an untyped donor.\n"
             "No confirmed positive or negative classification is assigned to ILMILQPQL.",
             fontsize=6.3, color=INK, linespacing=1.35)
    save_figure(fig, "S2", main_dir, supp_dir, data)


def figure_s3(data_root, main_dir, supp_dir):
    rows = read_tsv(data_root, "author_summaries/S3_auc.tsv")
    line, row = first_match(rows, metric="targeted_auc_median")
    auc = float(row["value"])
    data = []
    fig, ax = plt.subplots(figsize=(6.45, 3.15))
    figure_title(fig, "Archived candidate-label versus shuffled-label score ranking",
                 "The true presentation status of candidate labels was unknown.")
    ax.axvline(0.5, color=MID, linestyle="--", linewidth=0.9, label="Equal-ranking reference (AUC = 0.5)")
    ax.scatter([auc], [0], color=BLUE, s=38, zorder=3)
    ax.annotate(f"AUC {auc:.3f}", (auc, 0), xytext=(8, 7), textcoords="offset points",
                fontsize=7.1, fontweight="bold", color=BLUE)
    ax.set_xlim(0.35, 0.70)
    ax.set_ylim(-0.6, 0.6)
    ax.set_yticks([])
    ax.set_xlabel("AUC for distinguishing candidate labels from shuffled labels")
    ax.grid(axis="x", color=GRID, linewidth=0.45)
    ax.legend(loc="lower right", fontsize=6.1)
    fig.text(0.05, 0.09,
             "This AUC summarizes label-score separation in the archived extraction comparison; it is not presentation sensitivity or specificity.",
             fontsize=6.2, color=MID)
    add_data(data, "a", "candidate-versus-shuffle AUC", "median AUC", auc, "AUC",
             "author_summaries/S3_auc.tsv", line,
             "The candidate label is not a ground-truth presentation label.")
    add_data(data, "a", "equal-ranking reference", "AUC reference", 0.5, "AUC",
             "AUC definition", "equal-ranking reference", "Reference line only; no inferred confidence interval.")
    save_figure(fig, "S3", main_dir, supp_dir, data)


def figure_s4(data_root, main_dir, supp_dir):
    source = "author_summaries/S4_atlas_entries.csv"
    rows = read_csv(data_root, source)
    data = []
    fig = plt.figure(figsize=(6.6, 3.2))
    figure_title(fig, "Two atlas entries were not reproduced in PXD019643 RESULT files",
                 "Failure to reproduce in the queried RESULT files is not evidence of benign presentation.")
    for i, (line, row) in enumerate(rows):
        y = 0.62 - i * 0.30
        fig.patches.append(FancyBboxPatch(
            (0.08, y), 0.84, 0.20, boxstyle="round,pad=0.010,rounding_size=0.012",
            linewidth=0.7, edgecolor=GRID, facecolor=LIGHT,
            transform=fig.transFigure,
        ))
        fig.text(0.105, y + 0.145, row["mark"], fontsize=9.0, fontweight="bold", va="center")
        fig.text(0.37, y + 0.145, "not reproduced in PXD019643 RESULT files",
                 fontsize=6.9, color=GOLD, fontweight="bold", va="center")
        add_data(data, "a", "atlas sequence", row["mark"], row["value"], "reproduction status",
                 row["source_path"], row["source_line_or_selector"],
                 "Sequence-level audit only; no claim of benign tissue presentation.")
    save_figure(fig, "S4", main_dir, supp_dir, data)


def figure_s5(data_root, main_dir, supp_dir):
    ion_rel = "S5/candidate_ion_peak_assignments.tsv"
    evidence_rel = "S5/SIASAISMK_islet_spectrum_evidence.tsv"
    tumour_rel = "S5/SIASAISMK_tumor_reference_evidence.tsv"
    ion_rows = read_tsv(data_root, ion_rel)
    evidence = read_tsv(data_root, evidence_rel)
    tumour = read_tsv(data_root, tumour_rel)
    match_rel = "S5/tumour_hypothesis_by_ion.tsv"
    cross_matches = read_tsv(data_root, match_rel)
    match_counts = defaultdict(lambda: defaultdict(int))
    for _, row in cross_matches:
        hypothesis = row["reference_PSM_hypothesis"]
        match_counts[hypothesis][row["islet_scan"]] += row["observed_peak_at_each_hypothesized_mz"] == "yes"
    matched_ranges = {form: (min(items.values()), max(items.values())) for form, items in match_counts.items()}
    assert matched_ranges == {"unmodified_hypothesis": (1, 2), "M8_oxidized_hypothesis": (7, 8),
                              "M8_oxidized_source_PSM": (8, 9)}
    evidence.sort(key=lambda item: (
        0 if item[1]["preparation"] == "islet22" else 1,
        item[1]["technical_acquisition"],
    ))
    ion_order = [f"b{i}^1" for i in range(1, 9)] + [f"y{i}^1" for i in range(1, 9)]
    lookup = {}
    ion_source_rows = []
    for line, row in ion_rows:
        if row["peptide"] == "SIASAISMK" and row["ion_sequence"] == "candidate":
            key = (row["source_run"], row["ion"])
            if key in lookup:
                raise ValueError(f"Duplicate candidate-ion assignment for {key}")
            lookup[key] = (line, row)
            ion_source_rows.append((line, row))
    if len(evidence) != 4 or len(lookup) != 4 * len(ion_order):
        raise ValueError("Expected a complete 16-ion candidate map for four archived islet spectra")

    primary_line, primary = first_match(tumour, reference_id="primary_C3L01031")
    pdo_line, pdo = first_match(tumour, reference_id="PDO_PANFR0151")
    data = []
    fig = plt.figure(figsize=(7.15, 6.15))
    figure_title(fig, "SIASAISMK fragment assignments and tumour-reference provenance",
                 "Assigned islet fragments support a recurrent signal; sequence identity and matched tumour chemistry remain unresolved.")

    # Panel a: all theoretical b/y candidate ions are retained in the map.
    ax_map = fig.add_axes([0.12, 0.585, 0.82, 0.235])
    ax_map.set_title("a  Archived candidate-ion assignments across two preparations", loc="left", pad=7)
    prep_colors = {"islet22": BLUE, "islet27": TEAL}
    yvals = list(range(len(evidence)))[::-1]
    for y, (ev_line, ev_row) in zip(yvals, evidence):
        run = ev_row["source_run"]
        label = f"{ev_row['preparation']} · scan {ev_row['comet_scan']}"
        ax_map.set_yticks(yvals, [f"{r['preparation']} · scan {r['comet_scan']}" for _, r in evidence])
        for x, ion in enumerate(ion_order):
            line, row = lookup[(run, ion)]
            observed = bool(row["observed_mz"].strip())
            role = row["diagnostic_role"]
            ax_map.scatter([x], [y], marker="s", s=36, facecolor=WHITE,
                           edgecolor=GRID, linewidth=0.55, zorder=1)
            if observed:
                color = GOLD if role == "candidate_specific" else prep_colors[ev_row["preparation"]]
                ax_map.scatter([x], [y], marker="s", s=26, color=color,
                               edgecolor=WHITE, linewidth=0.25, zorder=2)
            status = "observed" if observed else "no assigned observed peak"
            add_data(data, "a", "candidate fragment assignment", f"{label}; {ion}", status,
                     "archived assignment status", ion_rel, line,
                     f"role={role} relative to {row['competitor']}; observed m/z={row['observed_mz'] or 'not retained'}.")
    ax_map.set_xticks(range(len(ion_order)), [x.replace("^1", "") for x in ion_order], fontsize=5.2)
    ax_map.axvline(7.5, color=GRID, linewidth=0.8)
    ax_map.set_xlim(-0.6, len(ion_order) - 0.4)
    ax_map.set_ylim(-0.6, len(evidence) - 0.4)
    ax_map.set_xlabel("Singly charged theoretical fragments")
    ax_map.tick_params(axis="y", length=0, pad=4)
    ax_map.tick_params(axis="x", length=0, pad=2)
    ax_map.legend(handles=[
        plt.Line2D([0], [0], marker="s", color="none", markerfacecolor=BLUE,
                   markeredgecolor=BLUE, markersize=4, label="Observed; shared fragment"),
        plt.Line2D([0], [0], marker="s", color="none", markerfacecolor=GOLD,
                   markeredgecolor=GOLD, markersize=4, label="Observed; candidate-specific vs ALSSALQFA"),
        plt.Line2D([0], [0], marker="s", color="none", markerfacecolor=WHITE,
                   markeredgecolor=GRID, markersize=4, label="No assigned observed peak"),
    ], loc="lower right", bbox_to_anchor=(1.0, 1.16), ncol=3, fontsize=5.2,
       handletextpad=0.25, columnspacing=0.75)

    # Panel b: only actually reported observed-mass errors are shown.
    ax_error = fig.add_axes([0.12, 0.345, 0.82, 0.17])
    ax_error.set_title("b  Mass error of observed candidate-specific fragments vs ALSSALQFA", loc="left", pad=6)
    run_positions = {row["source_run"]: i for i, (_, row) in enumerate(evidence)}
    ion_palette = {"b7^1": GOLD, "y2^1": BLUE, "b8^1": TEAL, "y8^1": NAVY,
                   "y7^1": "#8E6D9F", "b5^1": "#9B7852", "b6^1": "#657F61"}
    ion_offsets = {ion: (i - 3) * 0.035 for i, ion in enumerate(ion_palette)}
    plotted_error_n = 0
    for line, row in ion_source_rows:
        if row["diagnostic_role"] != "candidate_specific" or not row["observed_mz"].strip():
            continue
        ppm = float(row["mass_error_ppm"])
        x = run_positions[row["source_run"]] + ion_offsets.get(row["ion"], 0.0)
        ax_error.scatter([x], [ppm], s=20, color=ion_palette.get(row["ion"], MID),
                         edgecolor=WHITE, linewidth=0.35, zorder=3)
        add_data(data, "b", "candidate-specific fragment mass error", f"{row['source_run']}; {row['ion']}",
                 ppm, "ppm", ion_rel, line,
                 f"Observed m/z={row['observed_mz']}; role is relative to {row['competitor']}.")
        plotted_error_n += 1
    ax_error.axhline(0, color=MID, linewidth=0.65)
    ax_error.axhline(20, color=GRID, linestyle=":", linewidth=0.7)
    ax_error.axhline(-20, color=GRID, linestyle=":", linewidth=0.7)
    ax_error.set_xlim(-0.5, len(evidence) - 0.5)
    error_values = [float(row["mass_error_ppm"]) for _, row in ion_source_rows
                    if row["diagnostic_role"] == "candidate_specific"
                    and row["observed_mz"].strip()]
    error_limit = max(20.0, math.ceil(max(abs(x) for x in error_values) / 5.0) * 5.0 + 5.0)
    ax_error.set_ylim(-error_limit, error_limit)
    ax_error.set_xticks(range(len(evidence)),
                        [f"{r['preparation']} · {r['technical_acquisition'].replace('technical acquisition ', 'T')}"
                         for _, r in evidence])
    ax_error.set_ylabel("Mass error (ppm)")
    ax_error.tick_params(axis="x", length=0, labelsize=5.7)
    ax_error.grid(axis="y", color=GRID, linewidth=0.4, zorder=0)
    if plotted_error_n == 0:
        raise ValueError("No observed candidate-specific fragment errors to plot")

    box(fig, 0.07, 0.09, 0.42, 0.19, "Primary PDAC · C3L-01031 scan 16793",
        f"PSM {primary['source_PSM_sequence']}, z={primary['charge']}.\n"
        "Original modification state unresolved. Same-label matches:\n"
        "1–2 under unmodified, 7–8 under Met8 oxidation.",
        GOLD, PALE_GOLD, title_size=6.8, body_size=5.55)
    box(fig, 0.53, 0.09, 0.42, 0.19, "PDAC organoid · PANFR0151 scan 16495",
        f"PSM {pdo['source_PSM_sequence']}, z={pdo['charge']}.\n"
        "8–9 same-label matches to the unmodified islet spectra\n"
        "under the documented Met8-oxidised form.",
        BLUE, PALE_BLUE, title_size=6.8, body_size=5.55)
    add_data(data, "c", "tumour-reference evidence", "primary PDAC C3L-01031", "modification unresolved", "source status",
             tumour_rel, primary_line,
             f"{primary['modification_status_in_exact_source_row']}; archived precursor={primary['observed_precursor_mz_in_archive']}.")
    add_data(data, "c", "tumour-reference evidence", "PDO PANFR0151", "M8 oxidation explicit", "source status",
             tumour_rel, pdo_line,
             f"Source PSM={pdo['source_PSM_sequence']}; charge={pdo['charge']}.")
    for hypothesis, limits in matched_ranges.items():
        add_data(data, "c", "same-label matched fragments", hypothesis, f"{limits[0]}–{limits[1]}", "ions per islet comparison",
                 match_rel, f"reference_PSM_hypothesis={hypothesis}, group by islet_scan",
                 "Both observed peaks match their respective form-specific theoretical fragment mass.")
    fig.text(0.07, 0.045,
             "Fragment comparisons use form-specific theoretical masses. The precursor form of the primary PDAC reference remains unresolved.",
             fontsize=5.75, color=MID)
    save_figure(fig, "S5", main_dir, supp_dir, data)


def figure_s6(data_root, main_dir, supp_dir):
    source = "author_summaries/S6_candidate_lists.csv"
    rows = read_csv(data_root, source)
    vals = {row["mark"]: (line, row) for line, row in rows}
    data = []
    fig = plt.figure(figsize=(6.85, 3.65))
    figure_title(fig, "PDAC candidate-list construction and assay-tested subset",
                 "Source restriction predictions and islet donor typing remain separate evidence fields.")
    box(fig, 0.065, 0.52, 0.40, 0.24, "Source Table S6",
        "Cancer-restricted flag\nExact uppercase-sequence deduplication\n517 distinct candidate sequences",
        BLUE, PALE_BLUE, title_size=8.0, body_size=6.5)
    box(fig, 0.58, 0.52, 0.35, 0.24, "Source Table S10",
        "33 exact-overlap sequences\nassessed in source CTL experiments",
        TEAL, PALE_TEAL, title_size=8.0, body_size=6.5)
    fig.patches.append(FancyArrowPatch((0.47, 0.64), (0.57, 0.64), arrowstyle="-|>",
                                      mutation_scale=9, lw=0.9, color=MID,
                                      transform=fig.transFigure))
    fig.text(0.49, 0.68, "exact overlap", fontsize=5.8, color=MID, ha="center")
    fig.text(0.07, 0.36,
             "Candidate restrictions were predicted in the source study. Islet donor HLA typing defines a separate compatibility comparison.",
             fontsize=6.4, color=INK)
    fig.text(0.07, 0.25,
             "The 33 CTL-assessed sequences are a subset flag; the figure does not imply that all 517 candidates were functionally tested.",
             fontsize=6.25, color=MID)
    for key in ("CR set", "assay subset"):
        line, row = vals[key]
        add_data(data, "a", "candidate-list count", key, row["value"], "sequences",
                 row["source_path"], row["source_line_or_selector"], row["notes"])
    add_data(data, "a", "restriction type", "PDAC source list", "source-predicted", "HLA restriction",
             "C018_v2_evidence.zip::source Table S6", "candidate restriction annotations",
             "Predictions are distinct from measured recipient/donor HLA typing.")
    save_figure(fig, "S6", main_dir, supp_dir, data)


BUILDERS = {
    "1": figure_1,
    "2": figure_2,
    "3": figure_3,
    "4": figure_4,
    "5": figure_5,
    "S1": figure_s1,
    "S2": figure_s2,
    "S3": figure_s3,
    "S4": figure_s4,
    "S5": figure_s5,
    "S6": figure_s6,
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("figures", nargs="*", default=list(BUILDERS))
    args = parser.parse_args()
    main_dir = args.out / "Figures"
    supp_dir = args.out / "Supplementary"
    unknown = set(args.figures) - set(BUILDERS)
    if unknown:
        raise SystemExit(f"Not in this frozen pass: {sorted(unknown)}")
    for number in args.figures:
        BUILDERS[number](args.data_root, main_dir, supp_dir)
    print("Requested figure build completed.", flush=True)


if __name__ == "__main__":
    main()

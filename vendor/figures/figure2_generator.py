#!/usr/bin/env python3
"""C018 v4 Figure 2 plotting implementation used by the portable runner."""

import argparse
import csv
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
        "517 exact cancer-restricted sequences\n33 assessed in source CTL experiments",
        TEAL, PALE_TEAL, body_size=6.6)

    # Cohort cards. Same-organ islets occupy the upper, visually central row.
    box(fig, 0.055, 0.492, 0.425, 0.208, "Archived pancreatic islets  ·  same organ",
        "7 archived acquisitions. Five tissue-resolved.\n"
        "Three acquisitions from two donors with two-field typing.\n"
        "Source study: cytokine-treated islets (IFNγ, IL-1β and TNFα, 48 h)\n"
        "W632 fingerprint supports R369. Tissue link unresolved.",
        TEAL, PALE_TEAL, body_size=6.25)
    box(fig, 0.520, 0.492, 0.425, 0.208, "Independent primary islets  ·  same organ",
        "2 untreated preparations × 2 technical acquisitions\n"
        "HLA-A*02 and HLA-A*03 family-level annotations",
        TEAL, PALE_TEAL, body_size=6.25)
    box(fig, 0.055, 0.237, 0.425, 0.202, "Fibrotic-lung comparator",
        "12 acquisitions. Six passed the gate, all from IPF samples.\n"
        "Two inferred sample-to-donor mappings were tested.\n"
        "ILMILQPQL is retained as an unresolved signal",
        BLUE, PALE_BLUE, body_size=6.2)
    box(fig, 0.520, 0.237, 0.425, 0.202, "Skin reference  ·  background resource only",
        "25,518 source-table rows",
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
             "The two islet studies sampled distinct preparations under different cytokine conditions.",
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
             "all run records", "Five acquisitions were tissue-resolved, including three linked to donors with two-field HLA.")
    add_data(data, "b", "tissue-resolved acquisitions", "archived islets", 5, "acquisitions",
             "S1/archived_islet_sample_provenance.tsv", "tissue_resolved_for_islet_endpoint=yes",
             "One low-recovery control and the tissue-ambiguous W632 acquisition excluded.")
    add_data(data, "b", "two-field analysis acquisitions", "archived islets", 3, "acquisitions",
             "S1/archived_islet_sample_provenance.tsv", "included_in_exact_HLA_tissue_resolved_endpoint=yes",
             "Three eligible files from two typed donors.")
    add_data(data, "b", "two-field analysis donors", "archived islets", 2, "donors",
             "S1/archived_islet_sample_provenance.tsv", "distinct typed study donors with resolved tissue linkage",
             "R360 and HP20289.")
    add_data(data, "c", "biological preparations", "PXD045211", 2, "preparations",
             "PXD045211 SDRF and source-study Methods", "islet22/islet27; files 355/356 and 358/359",
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
    figure_title(fig, "HLA coverage and donor-level endpoints",
                 "Eligibility denominators are separate from qualifying calls; unresolved lung signals remain unclassified.")

    ax = fig.add_axes([0.178, 0.285, 0.365, 0.47])
    fig.text(0.178, 0.815, "a  Islet aeTSA candidates", fontsize=8.0, fontweight="bold")
    opportunities = [row for _, row in read_tsv(data_root, "S2/archived_ligand_HLA_opportunities.tsv")
                     if row["cohort"] == "NC2025_aeTSA"]
    ligand_any = sum(int(row["n_compatible_donors"]) >= 1 for row in opportunities)
    ligand_two = sum(int(row["n_compatible_donors"]) >= 2 for row in opportunities)
    assert (len(opportunities), ligand_any, ligand_two) == (220, 161, 14)
    tissue = read_tsv(data_root, "S2/aetsa_tissue_resolved_endpoints.tsv")
    exact_any_line, exact_any = first_match(tissue, typing="published_exact_HLA", endpoint="any")
    exact_two_line, exact_two = first_match(tissue, typing="published_exact_HLA", endpoint="both")
    typed = read_tsv(data_root, "S2/aetsa_exact_typed_endpoint_recomputed.tsv")
    typed_any_line, typed_any = first_match(typed, method="legacy", endpoint="any")
    typed_two_line, typed_two = first_match(typed, method="legacy", endpoint="both")
    historical = read_tsv(data_root, "S2/aetsa_endpoint_recomputed.tsv")
    hist_any_line, hist_any = first_match(historical, method="legacy", endpoint="any")
    hist_two_line, hist_two = first_match(historical, method="legacy", endpoint="both")
    scenarios = [
        ("Tissue-resolved exact two-field", int(exact_any["covered_NC2025"]), int(exact_two["covered_NC2025"]),
         "S2/aetsa_tissue_resolved_endpoints.tsv", exact_any_line, exact_two_line),
        ("Tissue-resolved ligand-HLA union", ligand_any, ligand_two, "S2/archived_ligand_HLA_opportunities.tsv",
         "NC2025_aeTSA, n_compatible_donors ≥1", "NC2025_aeTSA, n_compatible_donors ≥2"),
        ("Exact two-field, R369 retained", int(typed_any["denominator"]), int(typed_two["denominator"]),
         "S2/aetsa_exact_typed_endpoint_recomputed.tsv", typed_any_line, typed_two_line),
        ("Historical ligand-HLA union", int(hist_any["denominator"]), int(hist_two["denominator"]),
         "S2/aetsa_endpoint_recomputed.tsv", hist_any_line, hist_two_line),
    ]
    y_centres = [3.45, 2.45, 1.45, 0.45]
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
    ax.set_ylim(-0.1, 3.9)
    ax.set_xlabel("Candidates with HLA-compatible donor(s)")
    ax.grid(axis="x", color=GRID, linewidth=0.45, zorder=0)
    ax.tick_params(axis="y", length=0)
    ax.text(0.5, -0.28, "0 qualifying aeTSA calls under the within-donor two-spectrum rule",
            transform=ax.transAxes, ha="center", va="top", fontsize=6.35,
            color=INK, fontweight="bold")
    ax.legend(handles=[
        plt.Line2D([0], [0], color=BLUE, lw=5, label="Any donor"),
        plt.Line2D([0], [0], color=TEAL, lw=5, label="At least two donors"),
    ], loc="lower left", bbox_to_anchor=(0.0, 1.03), ncol=2, fontsize=6.1)

    endpoint_rows = read_tsv(data_root, "S2/lung_endpoints.tsv")
    lung_cells = [
        ("S1 main inferred mapping", "any", 71, 1, "Unresolved signal", 2),
        ("S1 main inferred mapping", "both", 50, 0, "0 qualifying calls", 3),
        ("S2 reversed inferred mapping", "any", 51, 0, "0 qualifying calls", 4),
        ("S2 reversed inferred mapping", "both", 36, 0, "0 qualifying calls", 5),
    ]
    fig.text(0.60, 0.815, "b  Lung mapping sensitivity", fontsize=8.0, fontweight="bold")
    tile_x, tile_y, tile_w, tile_h = 0.60, 0.645, 0.165, 0.137
    for i, (scenario, endpoint, denominator, code, audit_call, source_line) in enumerate(lung_cells):
        row_i, source_row = first_match(endpoint_rows, scenario="S1_passing" if scenario.startswith("S1") else "S2_reversed",
                                        endpoint=endpoint)
        # Guard denominators/codes against the archived frozen endpoint table.
        assert int(source_row["eligible_candidates"]) == denominator
        assert int(source_row["signals_meeting_spectrum_criterion"]) == code
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
                 f"Eligible candidates: {denominator}", fontsize=5.5, color=INK)
        fig.text(x + 0.009, y + 0.028,
                 f"{audit_call.capitalize()}", fontsize=5.5, color=INK)
        add_data(data, "b", "eligible lung candidates", f"{scenario}; {endpoint}",
                 denominator, "candidates", "S2/lung_endpoints.tsv", source_line,
                 "Candidate eligibility under the specified inferred linkage.")
        add_data(data, "b", "donor-level outcome", f"{scenario}; {endpoint}", audit_call,
                 "endpoint status", "S2/lung_endpoints.tsv", source_line,
                 "ILMILQPQL stays unresolved.")
    fig.text(0.60, 0.305,
             "ILMILQPQL remains an unresolved signal.\n"
             "PF5 contributes two spectra and PF7 contributes one.\n"
             "S1 and S2 are inferred sample-to-donor linkages.",
             fontsize=6.35, color=INK, linespacing=1.45, va="top")
    add_data(data, "b", "mapping status", "S1 and S2", "inferred", "sample-to-HLA linkage",
             "S2/ILMILQPQL_donor_HLA_S1_S2.tsv", "rows 2–7",
             "Under S1, PF5 has 2 spectra and PF7 has 1; under S2, PF5 maps to an untyped donor.")
    add_data(data, "b", "classification", "ILMILQPQL", "unresolved signal", "evidence status",
             "S2/ILMILQPQL_mapping_scenarios.tsv", "rows 2–3",
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
                "source": jy_rel,
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
    ax.set_position([0.335, 0.205, 0.54, 0.635])
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
    fig.text(0.045, 0.06,
             "Error bars are descriptive Wilson 95% intervals. HP20289 and PXD045211 rows are within-donor/preparation acquisition pairs.\n"
             "The JY point describes recovery after 1:16 eluate dilution in one selected DDA acquisition.",
             fontsize=6.15, color=MID, linespacing=1.4)
    save_figure(fig, 3, main_dir, supp_dir, data)


def figure_s1(data_root, main_dir, supp_dir):
    rel = "S3/jy_aligned_recovery.tsv"
    rows = read_tsv(data_root, rel)
    chosen = [(line, row) for line, row in rows if
              row["mode"] == "aligned_psm_and_peptide_q_plus1" and
              row["reference_scans_min"] == "1" and row["test_scans_min"] in ("1", "2")]
    chosen.sort(key=lambda item: (int(item[1]["dilution_factor"]), int(item[1]["test_scans_min"])))
    data = []
    fig, ax = plt.subplots(figsize=(6.7, 3.65))
    ax.set_position([0.13, 0.235, 0.83, 0.61])
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
    fig.text(0.05, 0.045,
             "The JY HLA background and 400–650 m/z acquisition window differ from those of the islet cohorts.",
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
    status_name = {"True": "compatible", "False": "incompatible", "unknown_unmapped": "HLA unknown"}
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
        ax.set_yticks(y, [f"{r['PF_label']}  ·  {status_name[r['HLA_compatible']]}" for _, r in part])
        ax.set_xlim(0, 3.45)
        ax.set_xticks([0, 1, 2, 3])
        ax.set_xlabel("Distinct spectra in sample")
        ax.grid(axis="x", color=GRID, linewidth=0.4, zorder=0)
        ax.tick_params(axis="y", length=0)
    fig.subplots_adjust(left=0.15, right=0.98, top=0.78, bottom=0.27, wspace=0.43)
    fig.text(0.055, 0.065,
             "PF5 reaches the two-spectrum rule in a compatible donor under S1. S2 maps PF5 to an untyped donor.\n"
             "The ILMILQPQL signal remains unresolved.",
             fontsize=6.3, color=INK, linespacing=1.35)
    save_figure(fig, "S2", main_dir, supp_dir, data)


def figure_s3(data_root, main_dir, supp_dir):
    rows = read_tsv(data_root, "author_summaries/S3_auc.tsv")
    line, row = first_match(rows, metric="targeted_auc_median")
    auc = float(row["value"])
    data = []
    fig, ax = plt.subplots(figsize=(6.45, 3.15))
    ax.set_position([0.12, 0.23, 0.83, 0.59])
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
    figure_title(fig, "Follow-up of two archived atlas entries")
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
             "Ely et al. Supplementary Table S6", "candidate restriction annotations",
             "Predictions are distinct from measured recipient/donor HLA typing.")
    save_figure(fig, "S6", main_dir, supp_dir, data)


BUILDERS = {
    "1": figure_1,
    "2": figure_2,
    "3": figure_3,
    "S1": figure_s1,
    "S2": figure_s2,
    "S3": figure_s3,
    "S4": figure_s4,
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
    print("Figure 4, Figure 5 and Figure S5 are built by their own generators.", flush=True)


if __name__ == "__main__":
    main()

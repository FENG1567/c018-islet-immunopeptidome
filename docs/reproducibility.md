# Figure reproducibility

The release runner regenerates the five main figures and six supplementary figures with the bundled C018 plotting code. It takes a locally obtained copy of `C018_Supplementary_Data_v4.zip`; the data archive is not vendored here. Before plotting, the runner checks the ZIP's full `MANIFEST.tsv` for file sizes and SHA-256 hashes.

```bash
python scripts/render_figures.py \
  --supplement /path/to/C018_Supplementary_Data_v4.zip \
  --out results/figures
```

To render a subset, for example Figures 2 and 4:

```bash
python scripts/render_figures.py \
  --supplement /path/to/C018_Supplementary_Data_v4.zip \
  --out results/figures \
  --figures 2 4
```

The output contains `Figures/C018_Figure_{1..5}_v4.pdf` and matching `_figure_data_v4.csv`, plus `Supplementary/C018_Figure_S{1..6}_v4.pdf` and matching CSVs. The runner requires Python 3.9 or later and Matplotlib 3.5 or later. PDF layout can vary slightly with Matplotlib and installed fonts; use the same Matplotlib/font versions when comparing renders pixel for pixel.

## Inputs by figure

All listed TSV inputs are read from the locally supplied ZIP and covered by its manifest. Files are extracted to a temporary directory only while rendering.

| Figure | Inputs and provenance |
|---|---|
| 1 | `S1/local_run_inventory.tsv`, `S1/lung_acquisition_eligibility.tsv`, `S4/external_detection_summary.tsv`, `S6/skin_reference_context.tsv`; the 589-entry NC2025 source-list size and 33-sequence CTL-tested overlap are author-supplied aggregates retained in the generator. |
| 2 | `S2/archived_ligand_HLA_opportunities.tsv`, `S2/aetsa_tissue_resolved_endpoints.tsv`, `S2/aetsa_exact_typed_endpoint_recomputed.tsv`, `S2/aetsa_endpoint_recomputed.tsv`, `S2/lung_endpoints.tsv`, `S2/ILMILQPQL_donor_HLA_S1_S2.tsv`, and `S2/ILMILQPQL_mapping_scenarios.tsv`. Counts and spectrum-criterion checks are read from these tables. |
| 3 | `S3/jy_aligned_recovery.tsv`, `S3/HP20289_paired_recovery.tsv`, and the eight `S3/pxd045211_{pride_recipe,raw_validation}_*.tsv` paired-recovery tables. |
| 4 | `S2/archived_exact_panel_candidate_donor_endpoints.tsv`, `S2/archived_exact_panel_candidate_run_matrix.tsv`, and `S2/archived_exact_panel_compatible_single_spectra.tsv`. The generator checks the 517-candidate denominator and frozen counts before saving. |
| 5 | `S5/SIASAISMK_islet_spectrum_evidence.tsv`, `S5/SIASAISMK_precursor_chemistry.tsv`, and `S5/SIASAISMK_tumor_reference_evidence.tsv`. |
| S1 | `S3/jy_aligned_recovery.tsv`. |
| S2 | `S2/ILMILQPQL_donor_HLA_S1_S2.tsv`. |
| S3 | The single author-derived median AUC value `0.527` in `vendor/figures/author_summaries/S3_auc.tsv`; this does not replay the archived score extraction. |
| S4 | Two author-summary lookup rows in `vendor/figures/author_summaries/S4_atlas_entries.csv`; the figure reports their lookup status and does not rerun the upstream atlas query. |
| S5 | `S5/candidate_ion_peak_assignments.tsv`, `S5/SIASAISMK_islet_spectrum_evidence.tsv`, `S5/SIASAISMK_tumor_reference_evidence.tsv`, and `S5/tumour_hypothesis_by_ion.tsv`. |
| S6 | Author-derived 517/33 list aggregates in `vendor/figures/author_summaries/S6_candidate_lists.csv`; the upstream cancer-restricted list construction and CTL overlap are not rerun. |

The author-summary files are limited plotting inputs for published aggregate results. They are not substitutes for the underlying source-list, atlas-query, or score-extraction analyses. Figure 1 likewise retains the 589 and 33 counts as author-supplied aggregates because the corresponding source lists are not included in the compact supplement.

## Manuscript table inputs

The displayed manuscript tables are authored summary JSON, not computed by the figure runner. The English and Chinese table specifications are maintained as separate JSON files, and the manuscript builder reads them together with the language-specific text, legends, and front/back matter. Those working-tree inputs are not part of this figure release.

## Public-sharing checks

The supplement includes acquisition labels, donor aliases, HLA typing, and inferred sample-to-donor linkage evidence. Those fields may remain sensitive genetic or health-related information even when names are absent. Review the original cohort consent, data-use terms, and repository licenses before publishing the supplement or rehosting its individual-level tables on GitHub. This release contains plotting code and compact author summaries only; it does not include the supplement ZIP or raw spectra.

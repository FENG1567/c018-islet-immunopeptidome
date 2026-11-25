# C018 immunopeptidomics reproducibility materials

This repository contains the analysis code and selected search configurations
for *HLA Coverage and Ligand Recovery in Human Pancreatic Islet Screens for
Pancreatic Cancer Antigen Candidates* (Feng et al.). The bundled tables,
figure-source data and replay checkpoints come from archived manuscript
revision **v4R1**. Figure numbering and layouts in the submitted manuscript
follow the later revision.

## Quick start

Use Python 3.12. The integrity and results checks use only the standard library.

```sh
python scripts/verify_release.py
```

The archived results and plotting commands below require a separate copy of the
historical **`C018_Supplementary_Data_v4.zip`** evidence archive. For this
optional replay, the archive's SHA-256 is
`d50e5dba07406057b5cf1223d1f39e9244814e13fb5476ffdd22dac88d9f392c`.
Use the intact archive as input, without extracting it first. The GitHub
package includes the summary tables and figure-source data. The current
manuscript supplement reports key supporting records in one document.

```sh
python scripts/reproduce_results.py --supplement input/C018_Supplementary_Data_v4.zip --out build/results
python -m pip install -r requirements-figures.txt
python scripts/render_figures.py --supplement input/C018_Supplementary_Data_v4.zip --out build/figures
```

## Contents and reproduction scope

| Path | Purpose |
| --- | --- |
| `data/tables/` | Machine-readable values and notes for all three main tables |
| `data/figure_data/` | Released mark-level evidence for Figures 1–5 and S1–S6 |
| `data/expected_results.tsv` | Numerical expectations fixed to manuscript v4R1 |
| `scripts/reproduce_results.py` | Reconstruct and check key results from archived evidence. Output identifies which values are recomputed and which are read from upstream summaries |
| `scripts/render_figures.py` | Run the available publication graphics generators from the submission supplement and packaged plotting inputs |
| `docs/reproducibility.md` | Per-figure input map and distinction between regenerated graphics and archived input summaries |
| `scripts/upstream/`, `config/` | Author search, scoring and recovery scripts with selected archived search recipes |
| `docs/source_datasets.tsv` | Accession numbers, source links, and input selectors |
| `docs/RAW_REPLAY.md` | Upstream input layout, search provenance, and the verified extent of raw-search replay |

The processed-data workflow starts with the candidate, acquisition, and
spectrum evidence in the historical archive. Raw datasets and source
publication supplements are obtained from the records listed in
[`docs/source_datasets.tsv`](docs/source_datasets.tsv).

Figure rendering reproduces the plotted content. PDF metadata and fonts can
vary by platform. Figures S3, S4 and S6 use the reported AUC, atlas lookup
and candidate-list summaries as inputs, respectively.

The raw-search components require the original spectral inputs, reference
lists, search databases, and Comet 2026.03. See the
[raw replay guide](docs/RAW_REPLAY.md) for branch-specific prerequisites and
validation limits.

## Licensing and citation

Author-written code is licensed under [MIT](LICENSE). The authors' original
summary results and figure-data annotations are licensed under
[CC BY 4.0](LICENSE-DATA.md). Third-party inputs retain their original terms.
Cite the project and the source studies used by the relevant workflow. See
[`CITATION.cff`](CITATION.cff) for repository metadata.

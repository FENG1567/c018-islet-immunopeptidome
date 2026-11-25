# Bundled figure generators

`figure2_generator.py` contains the v4-adjusted Figure 2 implementation derived from the submission build. `figure4_generator.py` contains the corresponding original plotting implementations for the other main and supplementary figures. `scripts/render_figures.py` calls the supported figure functions and supplies paths from the verified v4 supplementary archive.

The generators write v4 PDF and figure-data CSV names. Figure 2 reads the released S2 endpoint and HLA tables; Figure 4 reads the released exact candidate endpoint, candidate/run, and single-spectrum tables. For other figures, the runner maps the data tables listed in `docs/reproducibility.md` to the original plotting functions. It does not download upstream studies, raw spectra, or source candidate lists.

The three files in `author_summaries/` are compact values needed for Figures S3, S4, and S6. They document the exact author-supplied plotting inputs and do not claim to reproduce the upstream score extraction, atlas lookup, or candidate-list overlap calculation.

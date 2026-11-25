# C018 v4R1 source data

- `tables/Table_1.tsv` through `Table_3.tsv` are the editable, UTF-8 tabular
  values from the v4R1 English main tables. Their companion notes define the
  endpoints and units. They are author-produced summaries.
- `figure_data/` contains the released mark-level CSV for each of Figures 1–5
  and S1–S6. Columns identify the panel, plotted mark, group, reported value,
  unit, evidence table, and line or selector. Where the selector describes a
  logical subset it applies to all matching rows. These source selectors
  identify records in the historical v4 evidence archive used for v4R1.
- `expected_results.tsv` is the versioned numerical contract for the processed
  evidence replay. It identifies the values checked against C018 v4R1.

Counts refer to the unit reported by the corresponding table. A spectrum is
an acquisition-level observation. A donor endpoint is defined by the
candidate–donor eligibility and at least two qualifying spectra per donor.
Potential HLA compatibility, a screen-level match and sequence-confirmed
presentation are distinct claims. See the manuscript for sampling and
chemical restrictions.

[Licensing](../LICENSE-DATA.md) and [third-party sources](../docs/source_datasets.tsv).

# C018 raw-to-results replay guide

This guide accompanies the v4R1 release materials. The v4R1 manuscript reuses the existing analyses; its edits did not trigger a new spectrum search. The supported release path is therefore (1) audit the hashed S1–S6 tables and figure data already supplied, and (2) use the upstream modules below when the corresponding source spectra and historical analysis inputs are supplied.

The commands are modular and do not download data. They write only beneath the output directories passed to them. Use a fresh output directory for each rerun. The v4R1 values summarized from the archived PDAC candidate-by-acquisition matrix are handled by the separate top-level results runner; this guide does not recompute those aggregate values.

## Runtime and upstream files

- Search engine: Comet 2026.03 rev. 0 (`fa08489` in archived parameter headers). SHA-256 of the retained project Linux executable: `ad93b4cf60c2ed7afc2f41b8a3a05567c3938d999f17deb7b3676bc1fa91e7ed`. Obtain Comet separately from its [official releases](https://github.com/UWPR/Comet/releases).
- The six portable parameter snapshots are in `config/`. See `config/upstream_params_provenance.tsv` and `config/UPSTREAM_PARAMS.md` for hashes and path changes.
- `config/comet.archived_islet.primary.params` is the recovered primary-islet recipe. Its search FASTA `msv_db.fasta` has archived SHA-256 `907dccc82d81ef2135936f16df916f1e44b7aae30e8acd1d3a2d72d394d14dbf`.
- Python scripts use the standard library. The shell runners require Bash, GNU core utilities, and a compatible Comet executable.

Set these shell variables to user-supplied paths before running commands:

```bash
REPO=/path/to/C018-release/repo
C018_ROOT=/path/to/C018-source-tree
V2_ROOT=/path/to/v2_20261122/islet_reconstruction
COMET=/path/to/comet
```

`C018_ROOT` denotes the input tree with `work/`, `out/`, and source data tables. `V2_ROOT` denotes the islet replay tree with `inputs/mgf_full/`, `inputs/stern_table_s1_hla.tsv`, and writable `results/` directories. The internal v2 evidence ZIP is a compact code/evidence archive; it is not the raw input bundle.

## Archived-islet search and rank-one/q-value replay

The archived islet runner searches supplied MGF files, one at a time, with the supplied primary-search params, then the parser deduplicates rank-one charge queries by `(run, scan)`, computes per-run target–decoy q values, and writes source-scan records.

```bash
mkdir -p "$V2_ROOT/results/pepxml"
cp -n "$REPO/config/comet.archived_islet.primary.params" "$C018_ROOT/work/msv/comet.params"
"$REPO/scripts/upstream/run_islet_comet.sh" \
  "$COMET" "$C018_ROOT/work/msv/comet.params" \
  "$V2_ROOT/results/pepxml" "$V2_ROOT/inputs/mgf_full"/*.mgf

python "$REPO/scripts/upstream/islet_endpoint_reconstruction.py" extract \
  --base "$C018_ROOT" --v2-root "$V2_ROOT"
python "$REPO/scripts/upstream/islet_endpoint_reconstruction.py" analyze \
  --base "$C018_ROOT" --v2-root "$V2_ROOT"
```

The analysis step also reads the frozen source tables under `C018_ROOT/out/step5`, `out/step6`, `out/step7`, and `out/step8`, plus `work/candidates/nc2025_all_raw.csv`, `work/msv/islet_ligands.csv`, `work/msv/msv_db.fasta`, and `work/hla/cellrep_islet_classI.csv`. It cannot run from the released tables alone. The published v4 endpoint uses tissue-resolved sample linkage and source Table S1 two-field HLA typing; the older v2 replay also preserves its historical donor-union branch for comparison.

## Matched panel and Comet q values

`matched_panel_search.py` creates the positive/shuffled panel from supplied source ledgers, saves the panel FASTA and parameter file, and can search supplied MGF acquisitions. It expects the source layout used by the archived v2 code under `C018_ROOT` (islet ligand inventory, candidate list, and PDAC grading table). The panel seed is `C018-v2-20261122` in the source.

```bash
OUT=/path/to/new/local-panel-output
python "$REPO/scripts/upstream/matched_panel_search.py" \
  --base "$C018_ROOT" --out "$OUT" --comet "$COMET" \
  --params-template "$REPO/config/comet.local_competition.v2.params" \
  --mgf "$V2_ROOT/inputs/mgf_full/20200228_IsletsW632_OTOT_lumos.mgf" \
        "$V2_ROOT/inputs/mgf_full/20201202_Lumos_Stern_Islet_Set1_W632_5.mgf"

python "$REPO/scripts/upstream/score_panel.py" \
  --panel "$OUT/reference_targets_shuffles.tsv" --out "$OUT/summary" \
  --pepxml "$OUT/body/20200228_IsletsW632_OTOT_lumos.pep.xml" \
           "$OUT/body/20201202_Lumos_Stern_Islet_Set1_W632_5.pep.xml"
```

The external-islet branch below uses the opposite-acquisition source identifications when constructing each search panel. The metadata-adjusted recipe records 10 ppm precursor tolerance, a 0.03 Da fragment bin, no fixed cysteine modification, methionine oxidation, N/Q deamidation, and protein N-terminal acetylation. Candidate calls use whole-search PSM and I/L-grouped peptide q-value rules.

## External paired recovery and JY dilution

`analyze_external_islets.py` expects `--acquisition` to be a source directory containing `metadata/PXD045211.files.json`, `metadata/PXD045211_community_annotated.sdrf.tsv`, and the public MGF files under `raw/`. Its `--normalized-from`, when supplied, is a directory containing normalized MGF files and the matching `.scan_map.tsv` crosswalks. The script checks the file manifest against the SDRF before searching. `--source-psms` is the author PSM table used for opposite-acquisition recovery, and `--out` is a new output path whose parent already exists.

`paired_recovery.py` checks source MGF spectrum indices against author PSM precursor m/z and charge, then evaluates recovery in the opposite acquisition.

```bash
python "$REPO/scripts/upstream/analyze_external_islets.py" \
  --base "$C018_ROOT" --acquisition /path/to/PXD045211_acquisition_directory \
  --source-psms /path/to/source-author-psms.tsv --out "$OUT/external_pxd" \
  --normalized-from /path/to/normalized_mgf_directory --comet "$COMET" \
  --params-template "$REPO/config/comet.pxd045211.raw_validation.params"

# Add --pride-recipe and select comet.pxd045211.pride_recipe.params for the
# metadata-adjusted branch.

python "$REPO/scripts/upstream/paired_recovery.py" \
  --reference /path/to/source-author-psms.tsv \
  --ref-mgf /path/to/reference-acquisition.mgf \
  --test-mgf /path/to/test-acquisition.mgf \
  --psms /path/to/matched-panel-psms.tsv --out "$OUT/paired"
```

`jy_aligned.py` re-evaluates archived Comet PSMs from undiluted, 1:4, and 1:16 JY acquisitions. Supply the source opportunity table and original FASTA so the reference set is limited to original `jy_` entries observed in the undiluted run. The retained `jy_db.fasta` SHA-256 is `881dc3818bf78dff72a3fdfe38a026807cddaac6f86bc7a530b17b461de9f58f`.

```bash
python "$REPO/scripts/upstream/jy_aligned.py" \
  --out "$OUT/jy" --opportunities /path/to/jy_opportunities.tsv \
  --original-fasta /path/to/original_jy_database.fasta \
  --pepxml /path/to/undiluted.pep.xml /path/to/dilution_1to4.pep.xml \
           /path/to/dilution_1to16.pep.xml
```

## Selected-spectrum checks from v3

These scripts consume retained v2 source artifacts; they do not download or search spectra. The `v2-base` tree must contain the archived candidate-index MGF/pepXML files and the fragment-evidence mapping described in the script.

```bash
python "$REPO/scripts/upstream/reannotate_selected_spectra.py" \
  --v2-base /path/to/v2_20261122 --out "$OUT/selected_spectra"
python "$REPO/scripts/upstream/extract_pride_competition.py" \
  --v2-base /path/to/v2_20261122 --out "$OUT/pride_top10"
```

## PDAC results and search lineage

The exact v4R1 PDAC endpoint is replayed from the archived S2 candidate-by-acquisition matrix by the top-level results runner. This reproduces the manuscript-level summary without rerunning Comet:

```bash
python "$REPO/scripts/reproduce_results.py" \
  --supplement /path/to/C018_Supplementary_Data_v4.zip \
  --out "$OUT/result-replay"
```

The search-level method used the historical fixed-cysteine matched panel and is represented by `matched_panel_search.py`, `score_panel.py`, and `config/comet.local_competition.v2.params`. Replaying that search requires the original source tables, candidate/reference-panel inputs, five MGF acquisitions, Comet executable, and the archived Comet parameter context. Those full inputs and the per-run pepXML outputs are not in this release. The v3 no-static-cysteine sensitivity branch is excluded from the release code and is not the completed fixed-cysteine analysis used for the primary v4R1 endpoint; no completion claim is made for that sensitivity branch.

## What prevents a full raw-to-results rerun from this release alone

1. The release includes hashed result tables and selected PSM/peak evidence, not all source MGF/RAW files or the full per-run Comet pepXML outputs.
2. The six parameter files define recipes. Their referenced FASTA databases and complete source-input manifests are separate inputs.
3. The archived Comet binary fingerprint is listed above. The reader installs the executable separately.
4. The source human-reviewed protein FASTA, per-run full search outputs, and several source-side ledgers needed by the PDAC search are external inputs.
5. The PXD045211 MGF, mzIdentML/author PSMs, and candidate-index MGF/pepXML files needed for external recovery and selected-spectrum replay are not included in this repository. The evidence archive records checksums and source identifiers for some of these files, but not their contents.
6. The v4 S1–S6 package supports independent inspection of the cited row-level evidence, but it does not contain every intermediate table or an end-to-end workflow manifest linking each manuscript value to a script invocation.

Use public source identifiers and their repository deposit terms in `source_datasets.tsv` to obtain any missing source inputs. No private machine paths or prior run logs are required by the packaged commands.

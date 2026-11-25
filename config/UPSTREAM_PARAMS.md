# Archived Comet parameter provenance

Five parameter files in this directory are portable copies of parameter files listed in the internal `C018_v2_evidence.zip` archive. The archive SHA-256 is `b10939bbd705ba3d5df943e586c9f38803a7a723188b6bd7a4ada0cae8be5e7e`. The sixth copy comes from the retained primary-islet search input. The source and published-file SHA-256 values are recorded in `upstream_params_provenance.tsv`.

Every source file identifies **Comet 2026.03 rev. 0 (fa08489)** in its header. The archived v2 README and code identify the primary-islet template as the project’s `work/msv/comet.params`. Its copy was recovered from the retained project input tree. It is separate from the five branch-specific snapshots in the compact evidence archive.

## Files

- `comet.local_competition.v2.params`: fixed-cysteine local matched-panel search. This is the archived template used to seed the PDAC fixed-cysteine search branch.
- `comet.archived_islet.primary.params`: recovered primary archived-islet recipe. Its expected database is `msv_db.fasta`.
- `comet.pxd045211.raw_validation.params`: the archived raw-validation recipe for PXD045211.
- `comet.pxd045211.pride_recipe.params`: the PXD045211 metadata-adjusted recipe.
- `comet.stern_cross_acquisition.set1.params` and `comet.stern_cross_acquisition.set2.params`: opposite-acquisition reference-panel searches.

## Sanitization

Each copy uses a relative FASTA basename for `database_name`. Attribution comments were added and line endings normalized. All other parameter values are retained. At runtime, the search scripts may set `database_name` to the FASTA path for that run.

The Comet parameter template is © 2012–2026 Jimmy Eng. It is distributed
under Apache-2.0, with the [license text](COMET_APACHE_2.0.txt) included here.
See [the upstream license at fa08489](https://github.com/UWPR/Comet/blob/fa08489/LICENSE).
The project's original Python and shell code uses MIT separately.

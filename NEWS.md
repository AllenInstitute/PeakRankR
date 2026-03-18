# PeakRankR NEWS

## v1.1.0

### Bug fixes
- **Fixed hardcoded `"subclass"` column** (issue #7). The group column is now
  set via `group_by_column_name` (default `"cell_type"`), accepting any name.
- **Fixed incorrect dependencies** (issue #5). `DESCRIPTION` now accurately declares only `bedtoolsr` as a hard dependency; previously unused Bioconductor packages (`GenomicRanges`, `IRanges`, `rtracklayer`) were incorrectly listed as required.
- **Added input validation** (issue #6). `Peak_RankR()` checks all required
  columns, weight validity, and non-empty inputs with clear error messages.

### New features
- `check_bedtools()` — exported utility that verifies bedtools is accessible
  before running the pipeline, with install instructions if missing.
- Flexible column names via `chr_col`, `start_col`, `end_col`, `magnitude_col`.
- Replicate bigwig support — sensitivity averaged across all files per group.
- `bt.map` no-overlap `"."` placeholder handled correctly (converted to 0).
- Unit test suite via `testthat`.
- Vignette: `vignettes/getting-started.Rmd`.

### Removed
- `PeakRankR.R` (duplicate of `Peak_RankR.R` with hardcoded `cell.population`
  column and broken `dplyr` pipe references).
- `Peak_MACS2_rank.R`, `Peak_intersect_rank.R`, `Peak_coverage_rank.R`,
  `multiBigwig_summary_SS.R` — superseded by the unified `Peak_RankR.R`
  implementation. Their logic is incorporated into the main scoring pipeline.
- Root-level `test_file.tsv`, `bw_table`, `Astrocytes-1.bw`, `Astrocytes-2.bw`,
  `scaling.png`, `paper.md` — test/dev artefacts not part of the package.
- Stale `.tar.gz` build artefact.

## v1.0.0

- Initial release.

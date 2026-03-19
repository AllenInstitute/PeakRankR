# Changelog

## v0.1.0 (2026-02-26)

Initial release.

### Features
- `init` — initialize a deduplicated feature table from a BED/TSV peaks file
- `add-signal` — summarize BigWig signal over peaks (sum/mean/max)
- `add-gc` — compute GC fraction per peak from a reference FASTA
- `add-phylop` — compute mean PhyloP conservation scores; optional liftOver support
- `add-moments` — compute skewness, kurtosis, and bimodality per BigWig track
- CLI (`pypeakranker`) and Python API
- Unit test suite via `pytest`

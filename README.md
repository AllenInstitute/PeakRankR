# PyPeakRankR

<p align="center">
  <img src="PyPR.png" width="320" alt="PyPeakRankR logo"/>
</p>

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![Tests](https://github.com/AllenInstitute/PeakRankR/actions/workflows/tests.yml/badge.svg?branch=python-package)](https://github.com/AllenInstitute/PeakRankR/actions/workflows/tests.yml)

**PyPeakRankR** is a Python package that extracts quantitative features from
genomic peaks and assembles them into a reproducible **peak × feature matrix**.
It separates deterministic feature extraction from downstream ranking, so
prioritization can be performed reproducibly using any statistical or
machine-learning approach. PyPeakRankR does **not** perform peak calling.

## Installation

```bash
# From GitHub
pip install git+https://github.com/AllenInstitute/PeakRankR.git@python-package

# Or from source
git clone -b python-package https://github.com/AllenInstitute/PeakRankR
cd PeakRankR
pip install -e .
```

Requires Python >= 3.9. `pyBigWig` needs a C compiler and `zlib`.

## Quick start (CLI)

```bash
# Initialize a feature table from a BED file
pypeakranker init --peaks peaks.bed --out features.tsv

# Add BigWig signal summaries
pypeakranker add-signal --table features.tsv --bigwig-files sample1.bw sample2.bw \
  --stat mean --out features.tsv

# Add GC content
pypeakranker add-gc --table features.tsv --reference-fasta genome.fa --out features.tsv

# Add PhyloP conservation scores
pypeakranker add-phylop --table features.tsv --phylop-bw phyloP100way.bw --out features.tsv

# Add distribution moments (skewness, kurtosis, bimodality)
pypeakranker add-moments --table features.tsv --bigwig-files sample1.bw sample2.bw \
  --out features.tsv

# Rank peaks by cell-type specificity
pypeakranker rank-specificity --table features.tsv \
  --target-cols Astrocytes_mean \
  --background-cols Astrocytes_mean Neurons_mean Oligo_mean \
  --out ranked.tsv
```

## Quick start (Python API)

Every CLI subcommand has a matching Python function:

```python
from pypeakranker import (
    init_table, add_signal, add_gc, add_phylop, add_moments, rank_by_specificity
)

init_table("peaks.bed", "features.tsv")
add_signal("features.tsv", ["sample1.bw", "sample2.bw"], "features.tsv", stat="mean")
add_gc("features.tsv", "genome.fa", "features.tsv")
add_phylop("features.tsv", "phyloP.bw", "features.tsv")
add_moments("features.tsv", ["sample1.bw"], "features.tsv")
rank_by_specificity(
    "features.tsv",
    target_cols=["Astrocytes_mean"],
    background_cols=["Astrocytes_mean", "Neurons_mean", "Oligo_mean"],
    out_tsv="ranked.tsv",
)
```

## Subcommands and output columns

Every subcommand reads an existing TSV, appends columns, and writes out.

| Subcommand | Python function | Columns added |
|---|---|---|
| `init` | `init_table` | `chr`, `start`, `end` |
| `add-signal` | `add_signal` | `<sample>_summary` per BigWig |
| `add-gc` | `add_gc` | `GC_content` |
| `add-phylop` | `add_phylop` | `phyloP_mean` |
| `add-moments` | `add_moments` | `<sample>_skewness`, `_kurtosis`, `_bimodality` |
| `rank-specificity` | `rank_by_specificity` | `specificity_score`, `specificity_rank` |

`rank-specificity` scores each peak by the ratio of target-group signal to
mean background signal, then min-max normalises to [0, 1].

All functions accept `quiet=True` to suppress progress messages and
`allow_missing_chroms=True` (where applicable) to handle missing chromosomes.

## Cross-assembly scoring (liftOver)

`add-phylop` can score peaks from one assembly against a PhyloP BigWig in
another. Supply a UCSC chain file and PyPeakRankR calls `liftOver` automatically.

```bash
pypeakranker add-phylop \
  --table features.tsv \
  --phylop-bw hg38.phyloP100way.bw \
  --chain rheMac10ToHg38.over.chain.gz \
  --out features.tsv
```

<details>
<summary>Python API and additional options</summary>

```python
add_phylop(
    table_tsv="features.tsv",
    phylop_bw="hg38.phyloP100way.bw",
    out_tsv="features.tsv",
    chain_file="rheMac10ToHg38.over.chain.gz",
    liftover_exe="liftOver",       # path to binary if not on PATH
    max_len=5000,                  # intervals longer than this receive 0
    drop_lifted_coords=False,      # set True to omit chr_target columns
)
```

Peaks that fail to lift over receive a score of `0`. Pass
`--allow-missing-chroms` (CLI) or `allow_missing_chroms=True` (API) to
suppress errors for chromosomes absent from the PhyloP track.

**Requirements:** UCSC [`liftOver`](https://hgdownload.soe.ucsc.edu/admin/exe/)
binary on your `PATH` and a
[chain file](https://hgdownload.soe.ucsc.edu/goldenPath/rheMac10/liftOver/)
for your source → target assembly pair.

</details>

## Tests

```bash
pip install -e ".[test]"
pytest tests/
```

## Used in

- **Johansen et al. (2025)** — [Evaluating methods for the prediction of cell-type-specific enhancers in the mammalian cortex](https://doi.org/10.1016/j.xgen.2025.100879). *Cell Genomics.*
  PeakRankR ranked among the top 3 of 16 methods in the BICCN community challenge.

- **Wirthlin et al. (2026)** — [A Cross-Species Enhancer-AAV Toolkit for Cell Type-Specific Targeting Across the Basal Ganglia](https://doi.org/10.64898/2026.02.23.706695). *bioRxiv.*
  Used in the CERP pipeline for 514 candidate enhancers across basal ganglia cell types in mouse and macaque.

## Citation

> Somasundaram, S., Johansen, N.J., Bakken, T.E., and Miller, J.A. (2026).
> PyPeakRankR: Reproducible Peak-Level Feature Extraction for Regulatory Element Ranking.
> *Journal of Open Source Software* (submitted).
> https://github.com/AllenInstitute/PeakRankR/tree/python-package

## Authors

Saroja Somasundaram, Nelson J. Johansen, Trygve E. Bakken, Jeremy A. Miller
— Allen Institute for Brain Science

## License

MIT — see [LICENSE](LICENSE).

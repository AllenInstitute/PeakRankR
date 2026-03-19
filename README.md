# PyPeakRankR

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)

**PyPeakRankR** is a Python package for extracting quantitative features from
a predefined set of genomic peaks and assembling them into a reproducible,
analysis-ready feature table.

It generates a standardized **peak × feature matrix** enabling systematic
ranking and comparison of regulatory elements across cell types, conditions,
or species. PyPeakRankR does **not perform peak calling** — it standardizes
feature extraction so that downstream prioritization can be performed
reproducibly using any statistical or machine-learning approach.

---

## Features

| Module | Command | Output |
|---|---|---|
| Initialize table | `init` | Clean deduplicated peak table |
| BigWig signal | `add-signal` | Sum / mean / max per track |
| GC content | `add-gc` | `GC_content` column |
| Conservation | `add-phylop` | `phyloP_mean` column (optional liftOver) |
| Distribution moments | `add-moments` | Skewness, kurtosis, bimodality per track |
| ATAC specificity rank | `rank-specificity` | Specificity score and rank per peak |

---

## Installation

### Requirements

- Python >= 3.9
- `pyBigWig` requires a C build environment (`gcc`, `zlib`)

### Install from GitHub

```bash
pip install git+https://github.com/AllenInstitute/PeakRankR.git@python-package
```

### Install from source

```bash
git clone -b python-package https://github.com/AllenInstitute/PeakRankR
cd PeakRankR
pip install -e .
```

---

## Quick Start

```bash
# 1. Initialize a feature table from peaks
pypeakranker init \
  --peaks peaks.bed \
  --out features.tsv

# 2. Add BigWig signal summaries (mean across peak)
pypeakranker add-signal \
  --table features.tsv \
  --bigwig-files sample1.bw sample2.bw \
  --stat mean \
  --out features.tsv

# 3. Add GC content
pypeakranker add-gc \
  --table features.tsv \
  --reference-fasta genome.fa \
  --out features.tsv

# 4. Add PhyloP conservation scores
pypeakranker add-phylop \
  --table features.tsv \
  --phylop-bw phyloP100way.bw \
  --out features.tsv

# 5. Add distribution moments (skewness, kurtosis, bimodality)
pypeakranker add-moments \
  --table features.tsv \
  --bigwig-files sample1.bw sample2.bw \
  --out features.tsv
```

The resulting `features.tsv` contains:
- Peak coordinates (`chr`, `start`, `end`)
- One signal column per BigWig file
- `GC_content`
- `phyloP_mean`
- Per-track `_skewness`, `_kurtosis`, `_bimodality` columns

After `rank-specificity`, `ranked.tsv` additionally contains:
- `specificity_score` — normalised [0, 1], higher = more specific to target group
- `specificity_rank` — rank within the table (1 = best candidate to clone)

---

## Python API

```python
from pypeakranker import init_table, add_signal, add_gc, add_phylop, add_moments, rank_by_specificity

init_table("peaks.bed", "features.tsv")
add_signal("features.tsv", ["sample1.bw", "sample2.bw"], "features.tsv", stat="mean")
add_gc("features.tsv", "genome.fa", "features.tsv")
add_phylop("features.tsv", "phyloP.bw", "features.tsv")
add_moments("features.tsv", ["sample1.bw"], "features.tsv")
rank_by_specificity(
    "features.tsv",
    target_cols=["Astrocytes_mean"],
    background_cols=["Astrocytes_mean", "Neurons_mean", "Oligo_mean"],
    out_tsv="ranked.tsv"
)
```

---

## Cross-assembly scoring (liftOver)

To score peaks against a PhyloP track in a different genome assembly:

```bash
pypeakranker add-phylop \
  --table features.tsv \
  --phylop-bw hg38.phyloP100way.bw \
  --chain rheMac10ToHg38.over.chain.gz \
  --out features.tsv
```

Requires UCSC `liftOver` to be installed and on your `PATH`.

---

## CLI Reference

```
pypeakranker {init,add-signal,add-gc,add-phylop,add-moments,rank-specificity} --help
```

---

## Running Tests

```bash
pip install -e ".[test]"
pytest tests/
```

---

## Design Philosophy

PyPeakRankR separates **feature extraction** (deterministic, standardized)
from **peak ranking** (user-defined, flexible). This ensures ranking logic
remains transparent and adaptable to specific biological questions.

---


## Used in

PyPeakRankR was used in the following published studies:

- **Johansen et al. (2025)** — [Evaluating methods for the prediction of cell-type-specific enhancers in the mammalian cortex](https://doi.org/10.1016/j.xgen.2025.100879). *Cell Genomics.*
  PeakRankR ranked among the top 3 methods in the BICCN community challenge across 16 competing methods.

- **Wirthlin et al. (2026)** — [A Cross-Species Enhancer-AAV Toolkit for Cell Type-Specific Targeting Across the Basal Ganglia](https://doi.org/10.64898/2026.02.23.706695). *bioRxiv.*
  PyPeakRankR was used in the Cross-species Enhancer Ranking Pipeline (CERP) to compute ATAC specificity, PhyloP conservation, GC content, signal moments, and composite rankings for 514 candidate enhancers across basal ganglia cell types in mouse and macaque.

---

## Citation

If you use PyPeakRankR in your research, please cite:

> Somasundaram, S. and Johansen, N.J. (2026). PyPeakRankR: Reproducible Peak-Level
> Feature Extraction for Regulatory Element Ranking.
> *Journal of Open Source Software*.
> https://github.com/AllenInstitute/PeakRankR/tree/python-package

---

## Authors

Saroja Somasundaram — Allen Institute for Brain Science

Nelson J. Johansen — Allen Institute for Brain Science

## License

MIT License. See [LICENSE](LICENSE).

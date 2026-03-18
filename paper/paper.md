---
title: "PeakRankR: An R package for ranking enhancer peaks for cloning"
authors:
  - name: Saroja Somasundaram
    orcid: 0000-0002-3729-9849
    affiliation: Allen Institute
  - name: Nelson Johansen
    orcid: 0000-0002-4436-969X
    affiliation: Allen Institute
date: 2025
paperurl: https://github.com/AllenInstitute/PeakRankR
doi: 10.5281/zenodo.15238528
---

# PeakRankR: An R package for ranking enhancer peaks for cloning

**PeakRankR** is an R package developed to assist genomics researchers in
efficiently ranking genomic enhancer peaks for experimental cloning and
in-vivo targeting. Given single-cell ATAC-seq (scATAC-seq) peak data across
cell groups (e.g. cell types, subclasses), the package scores and ranks each
peak to identify the best candidates for targeting a specific group.

**PeakRankR** is quick and easy to run, making it accessible to researchers
with varying levels of R experience. It integrates directly into existing R
analysis pipelines and requires only a peak TSV file with group labels and a
table of bigwig file paths as input.

## Methods

PeakRankR computes a composite score for each peak per group using three
normalised components:

$$
\text{Score} = W_{\text{spec}} \times \text{Specificity}
             + W_{\text{sens}} \times \text{Sensitivity}
             + W_{\text{mag}}  \times \text{Magnitude}
$$

All three components are min-max normalised to [0, 1] before weighting.
Each weight defaults to 1 (equal importance) but is fully user-configurable.

- **Specificity**: Normalised ratio of the target group's mean bigwig signal
  to the mean signal across all background groups. Higher values indicate the
  peak is more active in the target group than in others.
- **Sensitivity**: Fraction of all groups that have NO signal at this peak.
  A peak active in only one group scores highest. A peak active in all groups
  scores 0. Computed as: 1 - (groups_with_signal / total_groups).
- **Magnitude**: Normalised mean bigwig signal level at the peak.

Bigwig signal extraction is performed via `rtracklayer::import()`, which reads
bigwig files directly without conversion. Final ranking uses either a rank-sum approach
(default) or direct composite score ordering.

## Features

- Flexible group column — works with any column name (cell type, subclass,
  cluster, etc.), not just hardcoded names.
- One bigwig file per group — each sample_id in the bw_table is treated as a distinct group.
- User-configurable weights for specificity, sensitivity, and magnitude.
- Input validation with clear, actionable error messages.
- `check_bedtools()` utility to verify the system dependency before running.

## Installation

```r
# Install bedtoolsr R interface
devtools::install_github("PhanstielLab/bedtoolsr")

# Install PeakRankR
devtools::install_github("AllenInstitute/PeakRankR", dependencies = TRUE)
```

## Usage

```r
library(PeakRankR)

tsv_file <- read.table("test_file.tsv", header = TRUE, sep = "\t")
bw_table  <- read.table("bw_table.txt",  header = TRUE, sep = "\t")

ranked <- Peak_RankR(
  tsv_file_df          = tsv_file,
  group_by_column_name = "cell_type",
  background_group     = unique(tsv_file$cell_type),
  bw_table             = bw_table,
  rank_sum             = TRUE,
  weights              = c(1, 1, 1)
)
```

## Citation

If you use PeakRankR in your research, please cite:

> Allen Institute. *PeakRankR: Package to rank enhancer peaks for cloning.*
> https://github.com/AllenInstitute/PeakRankR
> DOI: 10.5281/zenodo.15238528

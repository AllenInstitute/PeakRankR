---
title: "PeakRankR: A tool for ranking and visualizing genomic regions"
authors:
  - name: Saroja Somasundaram
    orcid: Your ORCID
    affiliation: Allen Institute
date: 2025-05-07
paperurl: https://github.com/AllenInstitute/PeakRankR
journal: Journal of Open Source Software
volume: XX
issue: YY
pages: ZZ-ZZ
doi: 10.21105/joss.XXXXX  # replace with the actual DOI when available
---

# PeakRankR: A tool for ranking genomic regions

**PeakRankR** is an R package developed to assist genomics researchers in efficiently ranking genomic regions, such as enhancers, based on metrics derived from sequencing data such as scATAC-seq (single-cell Assay for Transposase-Accessible Chromatin with sequencing). The package facilitates the analysis of peak data by offering a suite of feature extraction scripts that enable ranking and prioritization of enhancers for specific cell types or subclasses.

**PeakRankR** is quick and easy to run, making it accessible to researchers with varying levels of R experience. It derives a variety of informative features from genomic peak data conveniently using R, allowing for streamlined integration into existing analysis pipelines. At its core, the tool applies a linear model to score and rank peaks using user-defined weights applied to the extracted features. This flexible design lets users tailor the ranking process to their specific needs, although the default configuration has been shown to work well across several datasets.

The effectiveness of the default settings was validated in the study *"Evaluating Methods for the Prediction of Cell Type-Specific Enhancers in the Mammalian Cortex"*, where **PeakRankR** demonstrated strong performance in identifying relevant enhancers. By combining usability with analytical rigor, **PeakRankR** helps researchers identify candidate regulatory elements efficiently in complex, large-scale single-cell datasets.

## Methods

To be written 

## Features


### Features

- **MACS2 Rank**: Ranks genomic peaks based on the –log10(p-value) assigned by the peak caller MACS2, reflecting peak significance in the original data.
- **Intersection Rank**: Ranks peaks based on their specificity by calculating how often a peak overlaps with peaks from other cell types. Peaks unique to a given cell type are ranked higher.
- **Coverage Rank**: Ranks peaks by their normalized read coverage within a specific cell type relative to others. Peaks with disproportionately high coverage in one cell type receive higher ranks.

## Installation

PeakRankR is available from GitHub and can be installed using the `devtools` package in R:

```r
# Install the devtools package if not already installed
install.packages("devtools")

# Install PeakRankR from GitHub
devtools::install_github("AllenInstitute/PeakRankR")


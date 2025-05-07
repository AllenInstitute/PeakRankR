---
title: "PeakRankR: A tool for ranking and visualizing genomic regions"
authors:
  - name: Your Name
    orcid: Your ORCID
    affiliation: Your Affiliation
  - name: Co-Author Name (if applicable)
    orcid: Co-Author ORCID (if applicable)
    affiliation: Co-Author Affiliation (if applicable)
date: 2025-05-07
paperurl: https://github.com/AllenInstitute/PeakRankR
journal: Journal of Open Source Software
volume: XX
issue: YY
pages: ZZ-ZZ
doi: 10.21105/joss.XXXXX  # replace with the actual DOI when available
---

# PeakRankR: A tool for ranking and visualizing genomic regions

**PeakRankR** is an R package designed to help researchers in genomics efficiently rank genomic regions, such as enhancers or other DNA elements, based on user-defined metrics. The package supports the analysis of peak data and provides visualizations to aid interpretation.

## Summary

PeakRankR offers functionality for ranking genomic regions based on various scores and visualizing these rankings with customizable plots. This allows users to analyze and interpret data in the context of genome-wide experiments, particularly in epigenomics and gene regulation studies.

## Features

- **Ranking functionality**: Rank genomic regions by a given score, such as fold change, p-values, or other metrics.
- **Visualizations**: Create heatmaps, scatter plots, and other visualizations to explore the data.
- **Customizable output**: Support for flexible data input and output formats, with options for export.

## Installation

PeakRankR is available from GitHub and can be installed using the `devtools` package in R:

```r
# Install the devtools package if not already installed
install.packages("devtools")

# Install PeakRankR from GitHub
devtools::install_github("AllenInstitute/PeakRankR")


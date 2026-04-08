---
title: "PyPeakRankR: Reproducible Peak-Level Feature Extraction for Regulatory Element Ranking"
tags:
  - Python
  - ATAC-seq
  - genomics
  - regulatory elements
  - bioinformatics
  - BigWig
  - conservation
authors:
  - name: Saroja Somasundaram
    orcid: 0000-0002-3729-9849
    affiliation: 1
  - name: Nelson J. Johansen
    orcid: 0000-0002-4436-969X
    affiliation: 1
  - name: Jeremy A. Miller
    orcid: 0000-0003-4549-3747
    affiliation: 1
    corresponding: true
  - name: Trygve E. Bakken
    orcid: 0000-0003-1110-6431
    affiliation: 1
    corresponding: true
affiliations:
  - name: Allen Institute for Brain Science, Seattle, WA, USA
    index: 1
date: 2026-02-26
bibliography: paper.bib
doi: 10.5281/zenodo.15238527
---

# Summary

High-throughput chromatin accessibility assays such as ATAC-seq [@Buenrostro2013]
generate large sets of candidate regulatory elements called 'peaks'. Downstream analyses require
prioritizing peaks across cell types, experimental conditions, or species to
identify the most biologically relevant candidates for functional validation.
However, peak prioritization workflows are frequently implemented using ad hoc
scripts with inconsistent feature definitions and aggregation strategies, limiting
reproducibility and cross-study comparability.

**PyPeakRankR** is a Python package that standardizes quantitative feature
extraction for predefined genomic peaks. It produces a reproducible peak × feature
matrix — covering BigWig signal summaries, GC content, PhyloP conservation scores,
distribution moments, and ATAC specificity rankings — that can be used for
downstream statistical analysis, ranking, or machine learning applications. By
separating deterministic feature generation from user defined ranking logic,
PyPeakRankR promotes transparency, modularity, and reproducibility. It runs
in minutes on thousands of peaks, making it a practical first step before downstream modelling.

# Statement of need

Prioritizing genomic peaks across cell types requires combining multiple
features: signal intensity, sequence properties such as GC content,
evolutionary conservation (PhyloP [@Pollard2010]), and higher-order signal
statistics. These features are typically computed using custom per-project
scripts that vary across laboratories, complicating benchmarking and
cross-study integration. PyPeakRankR addresses this gap for computational
biologists working with single-cell ATAC-seq (scATAC-seq) or bulk ATAC-seq
data who need to systematically prioritize peaks for experimental follow-up,
particularly for enhancer discovery or adeno-associated virus (AAV) tool design.

# State of the field

Existing genomics tools each address part of the problem but none assembles
a unified, portable feature matrix. Peak callers such as MACS2 [@Zhang2008]
identify open chromatin regions but rank peaks only by fold change or p-value,
reflecting signal strength rather than cell-type specificity — a peak with
high fold change may be active across many cell types and therefore a poor
candidate for cell-type targeted AAV tools. Differential accessibility tools
such as ArchR [@Granja2021] test for cell-type enrichment but operate within
their own data model and do not produce portable, tool-agnostic feature tables.
Annotation tools such as GREAT [@McLean2010] link peaks to genes but do not
score chromatin features. At the library level, `pyBigWig` [@Ramirez2020pyBigWig]
provides low-level BigWig access without peak-level aggregation, `deepTools`
[@Ramírez2016] computes matrix summaries oriented toward visualization, and
`pyfaidx` [@Shirley2015] enables FASTA access without a genomics feature
pipeline.

Table 1 summarises feature coverage across these tools.

| Tool | Peak-level signal | Cell-type specificity | GC content | PhyloP conservation | Signal moments | Portable TSV output | CLI + Python API | Cross-assembly (liftOver) |
|------|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| PyPeakRankR | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| ArchR | ✓ | ✓ | – | – | – | partial | – | – |
| MACS2 | ✓ | – | – | – | – | ✓ | ✓ | – |
| deepTools | ✓ | partial | – | – | – | ✓ | ✓ | – |
| GREAT | – | – | – | – | – | ✓ | web only | – |
| pyBigWig | low-level | – | – | – | – | – | ✓ | – |
| pyfaidx | – | – | via FASTA | – | – | – | ✓ | – |

: Feature coverage across genomics tools. ✓ = supported natively;
partial = limited or indirect support; – = not supported. ArchR portable
output is partial because outputs are tied to the ArchR project object. {#tbl:tools}

PyPeakRankR fills this gap by combining `pyBigWig`, `pyfaidx`, and
`scipy` [@Virtanen2020] into a flexible CLI pipeline that assembles
heterogeneous features into a single reproducible TSV table.

# Software design

PyPeakRankR is built around three core design decisions:

**1. Table-first, flexible pipeline.** Every subcommand reads an existing tab-separated values (TSV) file
and appends columns without modifying the peak coordinates. Users can run any
subset of steps or add custom columns; the table remains valid throughout.
The table-first design enables incremental extension.

**2. Separation of feature extraction from ranking.** Feature extraction is
deterministic: given the same inputs and the same BigWig files, the same table
is always produced. Ranking is deliberately left to the user or to the `rank-specificity` subcommand, which implements one well-defined ranking formula
but is not the only option. This separation means benchmarking studies can
compare ranking strategies using the same upstream feature matrix, which is
exactly how PyPeakRankR was used in the Brain Initiative Cell Census Network (BICCN) challenge [@Johansen2025].

**3. Command-line interface (CLI) + Python API parity.** Every subcommand wraps a
public Python function (`init_table`, `add_signal`, `add_gc`, `add_phylop`,
`add_moments`, `rank_by_specificity`), so the tool works equally in shell
pipelines and Python notebooks without reimplementing logic.

The specificity ranking formula computes the ratio of target group signal to
mean background signal, then min-max normalises to [0, 1]. This matches the
CERP pipeline [@Wirthlin2026] and the ATAC-specificity metric validated in
the BICCN challenge [@Johansen2025].

Each feature has a distinct biological rationale. GC content is lower in active
enhancers than in promoters or bulk genomic DNA, reflecting differences in
nucleosome occupancy. PhyloP conservation [@Pollard2010] identifies peaks under
cross-species purifying selection. Signal distribution moments — kurtosis
(sharpness), skewness (asymmetry), and bimodality (Sarle's coefficient) — are
motivated by Lu et al. [@Lu2015], who showed these shape features distinguish
enhancers from promoters in ChIP-seq data more reliably than signal intensity
alone. The table-first design is directly extensible: future columns could
include sequence model importance scores (e.g., from Borzoi or Enformer) or
spatially resolved scores from multiplexed error-robust fluorescence in situ
hybridization (MERFISH), integrating epigenomic and spatial context in one
reproducible matrix.

![Features collected by PyPeakRankR for each candidate peak: GC content,
PhyloP conservation, ATAC specificity, and signal distribution moments
(kurtosis, skewness, bimodality). Figure adapted from Wirthlin et al.
(2026) [@Wirthlin2026].](figure6_panelA.png){ width=100% }

# Benchmarking against the BICCN challenge

PyPeakRankR's specificity-ratio approach was validated in the Brain
Initiative Cell Census Network (BICCN) Community Challenge
[@Johansen2025], where PeakRankR ranked third among 16 methods for
predicting functional cell-type-specific enhancers across 19 cortical
cell types. Figure 2 illustrates the approach using three validated L5
extratelencephalic (ET) enhancers from that challenge.

The specificity metric correctly distinguishes peaks whose signal is
concentrated in the target cell type from peaks with comparable total
signal spread across multiple cell types. In vivo AAV screening
confirmed that the high-specificity peaks drove on-target expression
while the low-specificity peak did not, demonstrating that the
specificity ratio captures biologically meaningful enrichment that
signal-magnitude ranking alone misses [@Johansen2025].

![ATAC-seq signal tracks across five cortical cell types for three
validated L5 ET enhancers from the BICCN Community Challenge
[@Johansen2025]. Each column shows per-cell-type read-pileup profiles;
dashed vertical lines mark the MACS2 summit. AiE0456m and AiE0463m
show signal concentrated in L5 ET neurons, while AiE0460m shows
comparable signal across all five cell types. In vivo validation from
adeno-associated virus (AAV) screening in mouse motor cortex confirmed
on-target expression for the high-specificity peaks
[@Johansen2025].](biccn_three_peaks.png)

# Research impact statement

PyPeakRankR extends the R package PeakRankR, which used a minimal set of
three features (ATAC specificity, signal magnitude, and peak coverage) in
the BICCN Community Challenge [@Johansen2025]. PyPeakRankR re-implements
and expands this approach in Python to integrate directly with sequence
models and modern genomics workflows. In the BICCN challenge, PeakRankR
ranked among the top three methods out of 16 competing approaches for
predicting functional cell-type specific enhancers across 19 cortical cell
types, with direct experimental validation via in vivo AAV screening.

In a recent basal ganglia study from our group [@Wirthlin2026], PyPeakRankR was used
within the Cross-species Enhancer Ranking Pipeline (CERP) across multiple
basal ganglia (BG) cell types in mouse and macaque. The composite feature rankings
outperformed conventional fold-change approaches, and the resulting
enhancer-AAV tools achieved >70% on-target specificity across cell types,
with exemplary enhancers exceeding 90%.

These results establish direct experimental utility.

# Implementation

PyPeakRankR is implemented in Python (>=3.9) with the following dependencies:
`pandas` [@Reback2020] for tabular data handling, `numpy` [@Harris2020] for
numerical computation, `pyBigWig` [@Ramirez2020pyBigWig] for BigWig signal extraction,
`pyfaidx` [@Shirley2015] for FASTA sequence access, and `scipy` [@Virtanen2020]
for statistical distribution metrics. Installable via pip from GitHub; includes a `pypeakranker` CLI, unit tests,
and example data (`tests/test.bed`). Source:
<https://github.com/AllenInstitute/PeakRankR/tree/python-package> (MIT).

# AI usage disclosure

Generative AI tools (Claude, Anthropic) assisted with code scaffolding,
paper drafting, and test scaffolding. All outputs were reviewed and validated
by the authors, who take full responsibility for all submitted materials.
All core design decisions were made by the human authors.

# Acknowledgements

The authors thank the bioinformatics and enhancer adeno-associated virus (AAV) teams at the Allen Institute for Brain Science for feedback on feature definitions and pipeline design.

This research was supported by the Allen Institute, founded by Jody Allen, chair and co‑founder of Allen Family Philanthropies, and the late Paul G. Allen, investor, philanthropist, and co‑founder of Microsoft. We gratefully acknowledge their vision and generosity, which make this work possible. This research was also supported by U.S. National Institutes of Health (NIH) BRAIN Initiative Human and Mammalian Brain Atlas (HMBA) BICAN grant UM1MH130981. The content of this study is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.

# References

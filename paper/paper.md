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
  - name: Jeremy Miller
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

Peak prioritization requires combining signal intensity, sequence properties
such as GC content, evolutionary conservation (PhyloP [@Pollard2010]), and
higher-order signal statistics. These features are typically computed using
custom per project scripts that vary across laboratories, complicating
benchmarking and cross-study integration.

The target audience is computational biologists who work with single-cell ATAC-seq (scATAC-seq) or
bulk ATAC-seq data and need to systematically prioritize peaks for experimental
follow-up, particularly for enhancer discovery or adeno-associated virus (AAV) tool design.

Existing tools address related but distinct problems: peak callers such as
MACS2 [@Zhang2008] identify open chromatin regions but rank peaks only by
fold change or p-value, which reflects signal strength rather than cell-type
specificity. A peak with high MACS2 fold change may be active across many cell
types (a housekeeping element) and therefore a poor candidate for cell-type
targeted AAV tools. Differential accessibility tools such as ArchR [@Corces2018]
test for cell type enrichment but operate within their own data model. Annotation
tools such as GREAT [@McLean2010] link peaks to genes. None provide a unified,
flexible framework for assembling a standardized feature matrix across
heterogeneous input tracks — which is precisely what PyPeakRankR addresses.

# State of the field

Several tools perform individual aspects of peak level feature computation.
`pyBigWig` [@Ramirez2020pyBigWig] provides low-level BigWig access but no peak level
aggregation framework. `deepTools` [@Ramírez2016] computes matrix summaries but
is oriented toward visualization rather than tabular feature assembly. ArchR
[@Corces2018] computes cell-type specificity scores within its own data model
but does not produce portable, tool agnostic feature tables. `pyfaidx`
[@Shirley2015] enables FASTA sequence access but provides no genomics feature
pipeline.

PyPeakRankR fills this gap by combining `pyBigWig`, `pyfaidx`, and
`scipy` [@Virtanen2020] into a CLI pipeline that assembles heterogeneous
features into a single reproducible TSV table.

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

# Ranking vs. MACS2 fold change

Sorting peaks by MACS2 [@Zhang2008] fold change ranks by signal strength
but not specificity. ArchR [@Corces2018] ranks by differential accessibility
(log2FC, FDR), capturing enrichment but without normalising signal across
all background groups.

Figure 2 illustrates these differences using three validated cortical
enhancers from the Brain Initiative Cell Census Network (BICCN) Community
Challenge [@Johansen2025], where PyPeakRankR ranked third among 16 methods.
All three peaks target L5 extratelencephalic (ET) neurons in mouse motor
cortex [@Johansen2025].

AiE0456m has high ATAC-seq signal specifically in L5 ET neurons; both
ArchR (rank #1) and PyPeakRankR (Spec=0.92) correctly prioritise it.

AiE0460m has high total ATAC-seq signal — giving it ArchR rank #2 — but
that signal is distributed across multiple cell types. PyPeakRankR's
specificity ratio (target / mean background) correctly identifies it as
non-specific (Spec=0.18) and deprioritises it.

AiE0463m has low total ATAC-seq signal, placing it at ArchR rank #18.
Yet almost all of that signal falls within L5 ET neurons; PyPeakRankR
detects this concentration and assigns Spec=0.61, rescuing a peak that
signal-magnitude approaches miss. Multiple POU3F1 motifs at this locus
— the canonical transcription factor for L5 ET neurons — support its
functional specificity [@Johansen2025].

![Schematic ATAC-seq signal tracks across five cortical cell types for
three validated L5 ET enhancers from the BICCN Community Challenge
[@Johansen2025]. Each column shows per-cell-type read-pileup profiles;
dashed lines mark the MACS2 summit; purple ticks below denote POU3F1
transcription factor motifs. AiE0456m shows high, L5 ET-specific signal
(ArchR rank #1; PyPeakRankR Spec=0.92) and is correctly prioritised by
both methods. AiE0460m shows comparable signal height across all five
cell types — ArchR ranks it #2 by log2FC while PyPeakRankR's specificity
ratio identifies it as non-specific (Spec=0.18). AiE0463m has low total
signal — placing it at ArchR rank #18 — but its signal is nearly exclusive
to L5 ET neurons; PyPeakRankR rescues it (Spec=0.61). In vivo validation
from adeno-associated virus (AAV) screening in mouse motor cortex
[@Johansen2025].](biccn_three_peaks.png)

| Enhancer | ArchR rank | PyPeakRankR Spec | ATAC pattern | In vivo result |
|----------|-----------|-----------------|--------------|----------------|
| AiE0456m |  1 | 0.92 | High, L5 ET-specific       | On-target (strong) |
| AiE0460m |  2 | 0.18 | High, multi-cell-type      | Weak / non-specific |
| AiE0463m | 18 | 0.61 | Low total, L5 ET-concentrated | On-target (specific) |

: Three L5 ET cortical enhancers from the BICCN Community Challenge
[@Johansen2025]. ArchR rank is based on differential accessibility
(log2FC). PyPeakRankR specificity (Spec) is the normalised ratio of
target-group signal to mean background signal, computed using
`pypeakranker rank-specificity`. AiE0460m and AiE0463m represent
opposite failure modes of signal-magnitude ranking. {#tbl:peaks}


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

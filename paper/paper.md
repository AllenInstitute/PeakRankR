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
affiliations:
  - name: Allen Institute for Brain Science, Seattle, WA, USA
    index: 1
date: 2026-02-26
bibliography: paper.bib
doi: 10.5281/zenodo.15238527
---

# Summary

High-throughput chromatin accessibility assays such as ATAC-seq [@Buenrostro2013]
generate large sets of candidate regulatory elements. Downstream analyses require
prioritizing peaks across cell types, experimental conditions, or species to
identify the most biologically relevant candidates for functional validation.
However, peak prioritization workflows are frequently implemented using ad hoc
scripts with inconsistent feature definitions and aggregation strategies, limiting
reproducibility and cross-study comparability.

**PyPeakRankR** is a Python package that standardizes quantitative feature
extraction for predefined genomic peaks. It produces a reproducible peak × feature
matrix — covering BigWig signal summaries, GC content, PhyloP conservation scores,
distribution moments, and ATAC specificity rankings — that can be used for
downstream statistical analysis, ranking, or machine-learning applications. By
separating deterministic feature generation from user-defined ranking logic,
PyPeakRankR promotes transparency, modularity, and reproducibility in regulatory
element analysis.

# Statement of need

Peak prioritization for experimental validation (e.g., cloning enhancers into
AAV vectors for in vivo targeting) requires combining signal intensity from
BigWig tracks, sequence-based properties such as GC content, evolutionary
conservation metrics such as PhyloP [@Pollard2010], and higher-order signal
distribution statistics. These features are typically computed using custom
per-project scripts that vary across laboratories. Such variability complicates
benchmarking, replication, and cross-study integration.

The target audience is computational biologists and genomics researchers who
work with scATAC-seq or bulk ATAC-seq data and need to systematically prioritize
peaks for experimental follow-up — particularly for enhancer discovery, AAV tool
design, or regulatory element ranking across cell types or species.

Existing tools address related but distinct problems: peak callers such as
MACS2 [@Zhang2008] identify open chromatin regions but rank peaks only by
fold-change or p-value, which reflects signal strength rather than cell-type
specificity. A peak with high MACS2 fold-change may be active across many cell
types (a housekeeping element) and therefore a poor candidate for cell-type
targeted AAV tools. Differential accessibility tools such as ArchR [@Corces2018]
test for cell-type enrichment but operate within their own data model. Annotation
tools such as GREAT [@McLean2010] link peaks to genes. None provide a unified,
composable framework for assembling a standardized feature matrix across
heterogeneous input tracks — which is precisely what PyPeakRankR addresses.

# State of the field

Several tools perform individual aspects of peak-level feature computation.
`pyBigWig` [@Ramirez2020pyBigWig] provides low-level BigWig access but no peak-level
aggregation framework. `deepTools` [@Ramírez2016] computes matrix summaries but
is oriented toward visualization rather than tabular feature assembly. ArchR
[@Corces2018] computes cell-type specificity scores within its own data model
but does not produce portable, tool-agnostic feature tables. `pyfaidx`
[@Shirley2015] enables FASTA sequence access but provides no genomics feature
pipeline.

PyPeakRankR fills this gap by composing these libraries — `pyBigWig` for
signal extraction, `pyfaidx` for sequence access, `scipy` [@Virtanen2020]
for distribution metrics — into a composable CLI pipeline that assembles
heterogeneous features into a single reproducible TSV table.

# Software design

PyPeakRankR is built around three core design decisions:

**1. Table-first, composable pipeline.** Every subcommand reads an existing TSV
and appends columns — it never modifies the peak coordinates. This means a user
can run any subset of steps, add custom columns from other tools, and the table
remains valid throughout. The alternative — a single monolithic script that
computes all features at once — requires users to rerun everything when adding
a new feature type. The table-first design enables incremental extension.

**2. Separation of feature extraction from ranking.** Feature extraction is
deterministic: given the same inputs and the same BigWig files, the same table
is always produced. Ranking is deliberately left to the user — or to the
`rank-specificity` subcommand, which implements one well-defined ranking formula
but is not the only option. This separation means benchmarking studies can
compare ranking strategies using the same upstream feature matrix, which is
exactly how PyPeakRankR was used in the BICCN challenge [@Johansen2025].

**3. CLI + Python API parity.** Every subcommand is a thin wrapper over a
public Python function (`init_table`, `add_signal`, `add_gc`, `add_phylop`,
`add_moments`, `rank_by_specificity`). This means the tool works equally well
in shell pipelines and in Python notebooks or workflows, without reimplementing
logic.

The specificity ranking formula is: for each peak, compute the ratio of the
target group's signal to the mean signal across all background groups, then
min-max normalise to [0, 1] across all peaks. This definition matches the CERP
pipeline used in [@Wirthlin2026] and is consistent with the ATAC-specificity
metric validated in the BICCN challenge [@Johansen2025].

Figure 1 illustrates the features collected by PyPeakRankR for each candidate
peak: sequence conservation (PhyloP), GC content, ATAC-seq specificity across
cell groups, and signal distribution shape metrics (kurtosis, skewness, and
bimodality).

![Features collected by PyPeakRankR for each candidate enhancer peak. PhyloP
quantifies cross-species sequence conservation. GC content is computed from the
reference genome. ATAC specificity captures differential chromatin accessibility
between cell groups. Kurtosis, skewness, and bimodality describe the shape of
the per-base signal profile within the peak. Figure adapted from Wirthlin et al.
(2026) [@Wirthlin2026].](figure6_panelA.png)

# Ranking vs. MACS2 fold-change

Sorting peaks by MACS2 [@Zhang2008] fold-change ranks by signal strength
but not specificity — a broadly accessible peak is a poor candidate for
cell-type targeted experiments.

Figure 2 illustrates this using ten MACS2 narrowPeak calls (Table 1)
scored against four ENCODE human tissue ATAC-seq BigWig tracks
(GRCh38): colonic mucosa [@ENCODE2020] (ENCFF557AZH, ENCSR970UNF),
liver (ENCFF160VHY, ENCSR802GEV), heart left ventricle
(ENCFF455AFI, ENCSR117PYB), and lung (ENCFF210HIS, ENCSR647AOY).
Signal was extracted using PyPeakRankR's `add-signal` command.
Specificity scores (colon / mean across tissues, normalised [0, 1])
diverge substantially from MACS2 fold-change ranks. P6 (chrX,
FC=14.7, rank #7) achieves the highest specificity (1.000) with colon
mean=320 versus liver=7.6 and heart=7.1. P4 (chr2, FC=16.6, rank #3)
scores lowest (0.000) with signal spread across all tissues. P1
(chr4) is ambiguous: high colon signal but substantial heart
signal (mean=83.8). Table 1 lists all ten peaks.

![ATAC-seq signal tracks across four human tissues (colonic mucosa, liver,
heart left ventricle, lung) for selected peaks from Table 1. Left panel:
top peaks by MACS2 fold-change, showing broadly accessible signal across
all tissues. Right panel: top peaks by PyPeakRankR specificity score,
showing concentrated signal in the target tissue. Signal values are MACS2
p-value BigWig tracks from ENCODE (GRCh38); dashed lines mark MACS2
summits. Specificity scores are the ratio of target to mean background
signal across all four tissues, min-max normalised to [0,
1].](browser_tracks_comparison.png)


| Peak | Coordinates (GRCh38) | FC | FC rank | Colon | Liver | Heart | Lung | Spec | Spec rank |
|------|----------------------|----|---------|-------|-------|-------|------|------|-----------|
| P1  | chr4:49,149,084–49,149,919   | 11.66 | 10 | 444.8 |  4.4 | 83.8 | 49.4 | 0.105 | 8  |
| P2  | chr18:12,947,298–12,948,641  | 15.87 |  5 | 327.3 | 14.8 | 19.5 | 55.2 | 0.383 | 4  |
| P3  | chr1:244,451,043–244,452,405 | 17.03 |  2 | 339.2 |  9.7 | 21.4 | 50.3 | 0.650 | 3  |
| P4  | chr2:178,450,461–178,452,049 | 16.59 |  3 | 242.1 | 18.2 | 16.4 | 43.8 | 0.000 | 10 |
| P5  | chr5:119,267,962–119,269,594 | 16.40 |  4 | 212.0 | 14.3 | 18.3 | 35.7 | 0.015 | 9  |
| P6  | chrX:1,391,993–1,393,130     | 14.75 |  7 | 319.8 |  7.6 |  7.1 | 49.0 | 1.000 | 1  |
| P7  | chr2:197,498,901–197,500,709 | 13.91 |  9 | 241.8 | 13.2 |  8.7 | 39.8 | 0.526 | 5  |
| P8  | chr8:144,826,667–144,827,949 | 15.51 |  6 | 296.6 |  6.1 | 18.6 | 47.3 | 0.628 | 2  |
| P9  | chr10:132,536,759–132,537,934| 14.17 |  8 | 357.5 | 11.6 | 27.6 | 70.4 | 0.127 | 7  |
| P10 | chr10:69,123,341–69,124,698  | 17.35 |  1 | 248.7 | 10.0 | 17.9 | 35.3 | 0.533 | 6  |

: Ten MACS2 narrowPeak calls used in Figure 2 (GRCh38, colonic mucosa
ENCSR970UNF [@ENCODE2020]). MACS2 fold-change (FC) and PyPeakRankR
specificity score both range 0–1 within this set; ranks frequently
diverge, confirming that signal strength alone does not predict
tissue-specific accessibility. Specificity computed using `pypeakranker rank-specificity` with
four ENCODE tissue BigWig tracks. {#tbl:peaks}



# Research impact statement

PyPeakRankR has demonstrated realized research impact in two published studies.
In the BICCN Community Challenge [@Johansen2025], PeakRankR — which uses the
same specificity, magnitude, and coverage features as PyPeakRankR — ranked
among the top three methods out of 16 competing approaches for predicting
functional cell-type-specific enhancers across 19 cortical cell types. The
challenge provided direct experimental validation of ranked candidates via
in vivo AAV screening, establishing that ATAC-seq specificity-based feature
ranking identifies functional enhancers at a level competitive with sequence
deep learning models.

In the companion basal ganglia study [@Wirthlin2026], PyPeakRankR was used
within the Cross-species Enhancer Ranking Pipeline (CERP) across multiple
BG cell types in mouse and macaque. The composite feature rankings
outperformed conventional fold-change approaches, and the resulting
enhancer-AAV tools achieved >70% on-target specificity across cell types,
with exemplary enhancers exceeding 90%.

These applications demonstrate that PyPeakRankR's feature extraction produces
rankings with direct experimental utility across species and brain regions.

# Implementation

PyPeakRankR is implemented in Python (>=3.9) with the following dependencies:
`pandas` [@Reback2020] for tabular data handling, `numpy` [@Harris2020] for
numerical computation, `pyBigWig` [@Ramirez2020pyBigWig] for BigWig signal extraction,
`pyfaidx` [@Shirley2015] for FASTA sequence access, and `scipy` [@Virtanen2020]
for statistical distribution metrics. The package is installable via pip from GitHub, includes a `pypeakranker`
CLI, unit tests, and example data (`tests/test.bed`). Source code:
<https://github.com/AllenInstitute/PeakRankR/tree/python-package> (MIT).

# AI usage disclosure

Generative AI tools (Claude, Anthropic) assisted with code scaffolding,
paper drafting, and test scaffolding. All outputs were reviewed and validated
by the authors, who take full responsibility for all submitted materials.
All core design decisions were made by the human authors.

# Acknowledgements

Development was supported by the Allen Institute for Brain Science.
The authors thank the bioinformatics and enhancer AAV teams for feedback
on feature definitions and pipeline design.

# References

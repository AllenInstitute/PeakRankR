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

PyPeakRankR fills the gap between these low-level libraries and higher-level
analysis frameworks. Rather than competing with or duplicating existing tools,
it composes them: `pyBigWig` for signal extraction, `pyfaidx` for sequence
access, `scipy` [@Virtanen2020] for distribution metrics. The key design
contribution is the composable CLI pipeline that assembles heterogeneous features
into a single, reproducible TSV table that any downstream tool or statistical
model can consume.

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

A common approach to peak prioritization is to sort peaks by MACS2
[@Zhang2008] fold-change or significance score. While this ranks peaks by
signal strength, it does not capture cell-type specificity. A peak with high
fold-change may be broadly accessible across many cell types — a housekeeping
element — and therefore a poor candidate for cell-type targeted experiments.

Figure 2 illustrates this distinction using ten MACS2 peaks from a real
ATAC-seq experiment. When ranked by fold-change (left), peak P1 (chr4) is
placed last because it has the lowest fold-change (11.7) of the set. When
ranked by PyPeakRankR specificity score (right), P1 is ranked first because
its signal is concentrated in the target cell type relative to background.
Conversely, peaks P3 (chr1) and P10 (chr10) rank at the top by fold-change
but near the bottom by specificity, consistent with broad chromatin
accessibility across cell types. This divergence illustrates why fold-change
alone is insufficient for selecting cell-type specific regulatory elements.

![Comparison of MACS2 fold-change ranking (left) versus PyPeakRankR
specificity ranking (right) for ten MACS2 narrowPeak calls from a real
ATAC-seq experiment (test.bed). Peaks with the highest MACS2 fold-change
are not necessarily the most cell-type specific. P1 (chr4, green border)
ranks last by fold-change (FC = 11.7) but first by specificity. P3 (chr1)
and P10 (chr10) rank first and second by fold-change (FC = 17.0 and 17.4)
but near the bottom by specificity, consistent with broad chromatin
accessibility across cell types. Specificity scores are the ratio of target
to mean background ATAC signal, min-max normalised to [0, 1]; rank 1 is
the peak most exclusively active in the target
group.](browser_tracks_comparison.png)

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
rankings with direct experimental utility, spanning multiple species and brain
regions. The software is openly available and documented for community reuse.

# Implementation

PyPeakRankR is implemented in Python (>=3.9) with the following dependencies:
`pandas` [@Reback2020] for tabular data handling, `numpy` [@Harris2020] for
numerical computation, `pyBigWig` [@Ramirez2020pyBigWig] for BigWig signal extraction,
`pyfaidx` [@Shirley2015] for FASTA sequence access, and `scipy` [@Virtanen2020]
for statistical distribution metrics. The package is installable via pip from
GitHub, provides a `pypeakranker` CLI entry point, and includes unit tests
covering all core functions. Source code is available at
<https://github.com/AllenInstitute/PeakRankR/tree/python-package> under the
MIT license.

# AI usage disclosure

Generative AI tools (Claude, Anthropic) were used to assist with: code
scaffolding and refactoring of module structure, drafting sections of this
paper and the README, and test scaffolding. All AI-assisted outputs were
reviewed, edited, and validated by the authors. All core design decisions —
the table-first pipeline architecture, the composable CLI structure, the
specificity ranking formula, and the feature set — were made by the human
authors. The authors take full responsibility for the accuracy and content of
all submitted materials.

# Acknowledgements

Development was supported by the Allen Institute for Brain Science and informed
by regulatory genomics workflows developed in the Human Cell Types program.
The authors thank the Allen Institute bioinformatics and enhancer AAV teams for
feedback on feature definitions and pipeline design.

# References

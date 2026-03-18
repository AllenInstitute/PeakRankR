# PeakRankR

> **Rank enhancer peaks for cloning** — an R package from the [Allen Institute](https://alleninstitute.org)

[![License](https://img.shields.io/badge/License-Allen_Institute-blue.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15238528.svg)](https://doi.org/10.5281/zenodo.15238528)

PeakRankR scores and ranks genomic enhancer peaks across user-defined cell groups (e.g. cell types, subclasses, clusters) to prioritise candidates for experimental cloning and in-vivo targeting.

---

## How it works

Each peak receives a composite score from three normalised components:

```
PeakRankR_score = W(specificity) × Specificity
              + W(sensitivity) × Sensitivity
              + W(magnitude)   × Magnitude
```

| Component | Definition |
|---|---|
| **Specificity** | Normalised ratio of target-group signal to mean background signal |
| **Sensitivity** | Fraction of groups with NO signal at this peak — fewer groups sharing a peak = higher score |
| **Magnitude** | Normalised mean signal level |

All three components are min-max normalised to [0, 1] before weighting. By default each weight is 1 (equal importance).

---

## Installation

### 1. Install system dependency: bedtools

```bash
# macOS (Homebrew)
brew install bedtools

# Ubuntu / Debian
sudo apt-get install bedtools

# Conda
conda install -c bioconda bedtools
```

Verify installation:
```bash
which bedtools   # macOS / Linux
where bedtools   # Windows
```

### 2. Install R dependencies

```r
install.packages("devtools")
devtools::install_github("PhanstielLab/bedtoolsr")
```

### 3. Install PeakRankR

```r
devtools::install_github("AllenInstitute/PeakRankR", dependencies = TRUE)
```

---

## Quick Start

```r
library(PeakRankR)

# Optional: set bedtools path if not on system PATH
# options(bedtools.path = "/usr/local/bin/")

# Verify bedtools is accessible
check_bedtools()

# Load the built-in example files
tsv_file <- read.table(
  system.file("extdata", "test_file.tsv", package = "PeakRankR"),
  header = TRUE, sep = "\t"
)

# Build bw_table using the bundled bigwig files
extdata_path <- system.file("extdata", package = "PeakRankR")
bw_table <- read.table(
  file.path(extdata_path, "bw_table.txt"),
  header = TRUE, sep = "\t"
)
# Resolve filenames to full paths
bw_table$file_path <- file.path(extdata_path, bw_table$file_path)

# Run ranking
# Note: only Astrocyte has bundled bigwig files — other groups will
# rank by magnitude only. For full scoring, supply your own bigwig files.
ranked <- Peak_RankR(
  tsv_file_df          = tsv_file,
  group_by_column_name = "cell_type",
  background_group     = unique(tsv_file$cell_type),
  bw_table             = bw_table,
  rank_sum             = TRUE,
  weights              = c(1, 1, 1)
)

head(ranked[order(ranked$PeakRankR_rank), ])
```

---

## Input file formats

### Peak TSV (`tsv_file_df`)

Tab-separated. Required columns (names are configurable via function arguments):

| Column | Default arg | Default name | Description |
|---|---|---|---|
| Chromosome | `chr_col` | `chr` | e.g. `chr1` |
| Start | `start_col` | `start` | 0-based start coordinate |
| End | `end_col` | `end` | End coordinate |
| Magnitude | `magnitude_col` | *(optional)* | Signal strength — if omitted, auto-computed from bigwig |
| Group | `group_by_column_name` | `cell_type` | Cell type or group label |

The bundled example (`inst/extdata/test_file.tsv`) uses `cell.population` as the group column and has no magnitude column — magnitude is computed automatically from the bigwig signal:
```
chr	start	end	cell.population
chr6	31131000	31131500	Astrocytes-1
chr6	31126700	31127200	Astrocytes-1
```

If your data has different column names or includes a pre-computed magnitude column:
```r
Peak_RankR(
  tsv_file_df          = my_peaks,
  group_by_column_name = "subclass",     # your group column
  magnitude_col        = "peak.score"    # optional: omit to auto-compute from bw
)
```

### Bigwig table (`bw_table`)

Two-column tab-separated file. `sample_id` must match values in your group column.

```
file_path	sample_id
/data/bw/Excitatory_rep1.bw	Excitatory
/data/bw/Excitatory_rep2.bw	Excitatory
/data/bw/Inhibitory_rep1.bw	Inhibitory
```

One bigwig file per group. Each `sample_id` in your bw_table is treated as a distinct group.

---

## Function reference

### `Peak_RankR()`

| Argument | Default | Description |
|---|---|---|
| `tsv_file_df` | — | Peak data frame |
| `group_by_column_name` | `"cell_type"` | Group column name |
| `background_group` | all groups | Specificity background |
| `bw_table` | — | Bigwig path table |
| `rank_sum` | `TRUE` | Rank by rank-sum vs composite score |
| `weights` | `c(1,1,1)` | Weights: specificity, sensitivity, magnitude |
| `chr_col` | `"chr"` | Chromosome column |
| `start_col` | `"start"` | Start column |
| `end_col` | `"end"` | End column |
| `magnitude_col` | `"magnitude"` | Magnitude column |

### `check_bedtools()`

Verifies the `bedtools` binary is accessible before running the pipeline.

---

## Citation

> Allen Institute. *PeakRankR: Package to rank enhancer peaks for cloning.*  
> https://github.com/AllenInstitute/PeakRankR  
> DOI: [10.5281/zenodo.15238528](https://doi.org/10.5281/zenodo.15238528)

---

## License

BSD 2-Clause. See [LICENSE](LICENSE).

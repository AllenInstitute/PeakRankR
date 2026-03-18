#' PeakRankR: Rank Enhancer Peaks for Cloning
#'
#' PeakRankR prioritises a list of enhancer peaks from different cell groups
#' (e.g. cell types, subclasses, clusters) for experimental cloning and
#' targeting. It scores each peak using three complementary measures:
#' specificity, sensitivity, and magnitude, then returns a ranked list
#' per group.
#'
#' ## Quick start
#'
#' ```r
#' library(PeakRankR)
#'
#' # check bedtools is on your PATH
#' check_bedtools()
#'
#' # load your data
#' tsv_file <- read.table("your_peaks.tsv", header = TRUE, sep = "\t")
#' bw_table  <- read.table("your_bw_table.txt", header = TRUE, sep = "\t")
#'
#' ranked <- Peak_RankR(
#'   tsv_file_df          = tsv_file,
#'   group_by_column_name = "cell_type",  # or "subclass", "cluster", etc.
#'   background_group     = unique(tsv_file$cell_type),
#'   bw_table             = bw_table,
#'   rank_sum             = TRUE,
#'   weights              = c(1, 1, 1)
#' )
#'
#' head(ranked[order(ranked$PeakRankR_rank), ])
#' ```
#'
#' ## Scoring formula
#'
#' \deqn{Score = W_{spec} \cdot Specificity + W_{sens} \cdot Sensitivity + W_{mag} \cdot Magnitude}
#'
#' All three components are min-max normalised to [0, 1] before weighting.
#'
#' ## System requirement
#'
#' This package requires the `bedtools` binary on your system PATH:
#' - Install: <https://bedtools.readthedocs.io>
#' - R interface: `devtools::install_github("PhanstielLab/bedtoolsr")`
#' - Set a custom path: `options(bedtools.path = "/path/to/bedtools/bin/")`
#'
#' @keywords internal
"_PACKAGE"

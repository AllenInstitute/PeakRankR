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
#' ## Dependencies
#' - Bigwig signal is read via `rtracklayer` (Bioconductor).
#' - Genomic ranges are handled by `GenomicRanges`, `IRanges`, and `S4Vectors`.
#' - The `bedtools` system binary is used by `check_bedtools()` to verify the
#'   genomic toolchain. Install from <https://bedtools.readthedocs.io>.
#'
#' @keywords internal
"_PACKAGE"

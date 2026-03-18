#' Rank Enhancer Peaks Across Cell Groups
#'
#' Scores and ranks a set of genomic peaks by specificity, sensitivity, and
#' magnitude of signal across user-defined cell groups (e.g. cell types,
#' subclasses, clusters). Returns the input data frame with additional ranking
#' columns appended.
#'
#' The composite score is:
#' \deqn{PeakRankR_score = W_{spec} \cdot Specificity + W_{sens} \cdot Sensitivity + W_{mag} \cdot Magnitude}
#'
#' All three components are min-max normalised to [0, 1] before weighting.
#' Each weight defaults to 1 (equal importance).
#'
#' @param tsv_file_df A data frame of peaks. Must contain columns for
#'   chromosome, start, end, magnitude, and a group column. Column names are
#'   controlled by the chr_col, start_col, end_col, magnitude_col, and
#'   group_by_column_name arguments.
#' @param group_by_column_name Character. Name of the column in tsv_file_df
#'   that identifies the cell group (e.g. "cell_type", "subclass", "cluster").
#'   Default: "cell_type".
#' @param background_group Character vector. Groups to use as background when
#'   computing specificity. If NULL (default), all unique values of
#'   group_by_column_name are used.
#' @param bw_table A data frame with two columns: file_path (path to bigwig
#'   file) and sample_id (must match values in the group column).
#' @param rank_sum Logical. If TRUE (default), final rank is based on the
#'   sum of per-metric ranks. If FALSE, rank is based on the composite score.
#' @param weights Numeric vector of length 3: weights for
#'   c(specificity, sensitivity, magnitude). Default: c(1, 1, 1).
#' @param chr_col Character. Name of the chromosome column. Default: "chr".
#' @param start_col Character. Name of the start position column. Default: "start".
#' @param end_col Character. Name of the end position column. Default: "end".
#' @param magnitude_col Character. Name of the magnitude/signal column in
#'   \code{tsv_file_df}. If \code{NULL} (default), magnitude is computed
#'   automatically from the mean bigwig signal of the target group at each peak.
#' @param verbose Logical. Print progress messages. Default: TRUE.
#'
#' @return The input data frame with five additional columns per group:
#'   \describe{
#'     \item{specificity_score}{Normalised ratio of target group signal to background.}
#'     \item{sensitivity_score}{Fraction of groups with NO signal at this peak (1 = peak active in only one group, 0 = peak active in all groups).}
#'     \item{magnitude_score}{Normalised mean signal magnitude.}
#'     \item{PeakRankR_score}{Weighted composite score.}
#'     \item{PeakRankR_rank}{Rank within the group (1 = best candidate).}
#'   }
#'
#' @examples
#' \dontrun{
#' tsv_file <- read.table(
#'   system.file("extdata", "test_file.tsv", package = "PeakRankR"),
#'   header = TRUE, sep = "\t"
#' )
#' bw_table <- read.table(
#'   system.file("extdata", "bw_table.txt", package = "PeakRankR"),
#'   header = TRUE, sep = "\t"
#' )
#' ranked <- Peak_RankR(
#'   tsv_file_df          = tsv_file,
#'   group_by_column_name = "cell.population",
#'   background_group     = unique(tsv_file$cell.population),
#'   bw_table             = bw_table,
#'   rank_sum             = TRUE,
#'   weights              = c(1, 1, 1)
#' )
#' head(ranked[order(ranked$PeakRankR_rank), ])
#' }
#'
#' @export
Peak_RankR <- function(tsv_file_df,
                        group_by_column_name = "cell_type",
                        background_group     = NULL,
                        bw_table,
                        rank_sum             = TRUE,
                        weights              = c(1, 1, 1),
                        chr_col              = "chr",
                        start_col            = "start",
                        end_col              = "end",
                        magnitude_col        = NULL,
                        verbose              = TRUE) {

  .validate_inputs(
    tsv_file_df          = tsv_file_df,
    group_by_column_name = group_by_column_name,
    bw_table             = bw_table,
    weights              = weights,
    chr_col              = chr_col,
    start_col            = start_col,
    end_col              = end_col,
    magnitude_col        = magnitude_col
  )

  check_bedtools(quiet = !verbose)

  if (is.null(background_group)) {
    background_group <- unique(tsv_file_df[[group_by_column_name]])
  }

  df <- tsv_file_df
  names(df)[names(df) == chr_col]              <- "chr"
  names(df)[names(df) == start_col]            <- "start"
  names(df)[names(df) == end_col]              <- "end"
  if (!is.null(magnitude_col))
    names(df)[names(df) == magnitude_col]      <- "magnitude"
  names(df)[names(df) == group_by_column_name] <- "group_col"

  # Warn once if background_group contains values not present in the data
  unmatched_bg <- setdiff(background_group, unique(df$group_col))
  if (length(unmatched_bg) > 0)
    warning("background_group contains values not found in the data: ",
            paste(unmatched_bg, collapse = ", "),
            ". These will be ignored.", call. = FALSE)

  groups <- unique(df$group_col)
  if (verbose) message("PeakRankR: scoring ", length(groups), " group(s)...")

  results_list <- lapply(groups, function(grp) {
    if (verbose) message("  -> ", grp)
    .score_group(df, grp, background_group, bw_table, weights, rank_sum)
  })

  ranked_df <- do.call(rbind, results_list)

  names(ranked_df)[names(ranked_df) == "chr"]       <- chr_col
  names(ranked_df)[names(ranked_df) == "start"]     <- start_col
  names(ranked_df)[names(ranked_df) == "end"]       <- end_col
  if (!is.null(magnitude_col))
    names(ranked_df)[names(ranked_df) == "magnitude"] <- magnitude_col
  names(ranked_df)[names(ranked_df) == "group_col"] <- group_by_column_name

  if (verbose) message("PeakRankR: done.")
  ranked_df
}

#' @noRd
.score_group <- function(df, target_group, background_group,
                         bw_table, weights, rank_sum) {
  target_df     <- df[df$group_col == target_group, ]
  background_df <- df[df$group_col %in% background_group, ]
  target_bw     <- bw_table[bw_table$sample_id == target_group, "file_path"]

  if (length(target_bw) == 0) {
    warning("No bigwig files found for group '", target_group,
            "'. Ranking by magnitude only.", call. = FALSE)
    target_df$specificity_score <- NA_real_
    target_df$sensitivity_score <- NA_real_
    mag_raw_fb <- if ("magnitude" %in% names(target_df)) target_df$magnitude else rep(0, nrow(target_df))
    target_df$magnitude_score   <- .normalise(mag_raw_fb)
    target_df$PeakRankR_score   <- weights[3] * target_df$magnitude_score
    target_df$PeakRankR_rank    <- rank(-target_df$PeakRankR_score,
                                        ties.method = "min")
    return(target_df)
  }

  spec_scores <- .compute_specificity(target_df, background_df,
                                       bw_table, target_group)
  # Sensitivity: peaks active in fewer groups rank higher
  all_groups  <- unique(df$group_col)
  sens_scores <- .compute_sensitivity(target_df, bw_table, all_groups)

  # Magnitude: use provided column if available, otherwise compute from bw
  if ("magnitude" %in% names(target_df)) {
    mag_raw <- target_df$magnitude
  } else {
    mag_raw <- .get_mean_signal(target_df, target_bw)
  }
  mag_scores  <- .normalise(mag_raw)
  composite   <- weights[1]*spec_scores + weights[2]*sens_scores +
                 weights[3]*mag_scores

  target_df$specificity_score <- spec_scores
  target_df$sensitivity_score <- sens_scores
  target_df$magnitude_score   <- mag_scores
  target_df$PeakRankR_score   <- composite

  if (rank_sum) {
    target_df$PeakRankR_rank <- rank(
      rank(-spec_scores, ties.method = "min") +
      rank(-sens_scores, ties.method = "min") +
      rank(-mag_scores,  ties.method = "min"),
      ties.method = "min"
    )
  } else {
    target_df$PeakRankR_rank <- rank(-composite, ties.method = "min")
  }
  target_df
}

#' @noRd
.compute_specificity <- function(target_df, background_df,
                                  bw_table, target_group) {
  # target_signal: vector of length n_peaks — A1 signal at each peak
  target_signal <- .get_mean_signal(
    target_df,
    bw_table[bw_table$sample_id == target_group, "file_path"]
  )

  bg_groups <- unique(background_df$group_col)

  # Per-peak background: for each peak, get signal from every background group
  # Result: matrix n_peaks x n_bg_groups
  bg_signal_matrix <- do.call(cbind, lapply(bg_groups, function(grp) {
    bw_paths <- bw_table[bw_table$sample_id == grp, "file_path"]
    if (length(bw_paths) == 0) return(rep(0, nrow(target_df)))
    .get_mean_signal(target_df, bw_paths)
  }))

  if (!is.matrix(bg_signal_matrix))
    bg_signal_matrix <- matrix(bg_signal_matrix, nrow = nrow(target_df))

  # Per-peak average background across all groups
  avg_bg_per_peak <- rowMeans(bg_signal_matrix, na.rm = TRUE)

  # Avoid division by zero
  avg_bg_per_peak[avg_bg_per_peak == 0] <- NA

  # Ratio: how much stronger is target signal vs average background at each peak
  ratio <- target_signal / avg_bg_per_peak
  ratio[is.na(ratio)] <- 0

  .normalise(ratio)
}

#' @noRd
#' Sensitivity = fraction of ALL groups that have signal > 0 at this peak,
#' INVERTED so that peaks active in FEWER groups score HIGHER.
#' A peak active in only 1 group gets score 1.0 (most sensitive/specific).
#' A peak active in all groups gets score close to 0.
.compute_sensitivity <- function(target_df, bw_table, all_groups) {
  n_peaks <- nrow(target_df)

  # For each group, check if any of its bw files has signal > 0 at each peak
  group_has_signal <- do.call(cbind, lapply(all_groups, function(grp) {
    bw_paths <- bw_table[bw_table$sample_id == grp, "file_path"]
    if (length(bw_paths) == 0) return(rep(0, n_peaks))
    # Mean signal across bw files for this group
    grp_signal <- .get_mean_signal(target_df, bw_paths)
    as.numeric(grp_signal > 0)
  }))

  if (!is.matrix(group_has_signal))
    group_has_signal <- matrix(group_has_signal, nrow = n_peaks)

  # Count how many groups have signal at each peak
  n_groups_with_signal <- rowSums(group_has_signal, na.rm = TRUE)
  n_total_groups       <- length(all_groups)

  # Invert: fewer groups with signal = higher score
  # Score = 1 - (n_groups_with_signal / n_total_groups)
  # Peak in only 1 group: 1 - 1/N (high)
  # Peak in all groups:   1 - N/N = 0 (low)
  1 - (n_groups_with_signal / n_total_groups)
}

#' @noRd
.get_mean_signal <- function(peak_df, bw_paths) {
  # Write peaks to a temp BED3 file
  tmp_bed <- tempfile(fileext = ".bed")
  on.exit(unlink(tmp_bed), add = TRUE)

  write.table(
    data.frame(chr   = peak_df$chr,
               start = as.integer(peak_df$start),
               end   = as.integer(peak_df$end)),
    tmp_bed, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
  )

  n_peaks      <- nrow(peak_df)
  signals_list <- lapply(bw_paths, function(bw_path) {
    if (!file.exists(bw_path)) {
      warning("Bigwig file not found: ", bw_path, call. = FALSE)
      return(rep(NA_real_, n_peaks))
    }
    tryCatch({
      # bedtools cannot read .bw directly — use rtracklayer to read signal
      if (!requireNamespace("rtracklayer", quietly = TRUE))
        stop("Please install rtracklayer: BiocManager::install('rtracklayer')",
             call. = FALSE)

      # Build GRanges for all peaks (1-based, as required by Bioc)
      sel <- GenomicRanges::GRanges(
        seqnames = peak_df$chr,
        ranges   = IRanges::IRanges(
          start = as.integer(peak_df$start) + 1L,
          end   = as.integer(peak_df$end)
        )
      )

      # Import signal overlapping our peaks from the bigwig
      hits <- rtracklayer::import(bw_path, which = sel, format = "BigWig")

      if (length(hits) == 0) return(rep(0, n_peaks))

      # For each peak, mean signal across overlapping bigwig intervals
      # Use findOverlaps on all peaks at once against the imported signal
      ols <- GenomicRanges::findOverlaps(sel, hits)

      scores <- vapply(seq_len(n_peaks), function(i) {
        idx <- S4Vectors::subjectHits(ols[S4Vectors::queryHits(ols) == i])
        if (length(idx) == 0) return(0)
        mean(hits$score[idx], na.rm = TRUE)
      }, numeric(1))

      return(scores)

    }, error = function(e) {
      warning("Error reading '", bw_path, "': ", conditionMessage(e),
              call. = FALSE)
      rep(NA_real_, n_peaks)
    })
  })

  signal_mat <- do.call(cbind, signals_list)
  if (!is.matrix(signal_mat)) signal_mat <- matrix(signal_mat, nrow = n_peaks)
  rowMeans(signal_mat, na.rm = TRUE)
}

#' @noRd
.normalise <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) return(rep(0, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

#' @noRd
.validate_inputs <- function(tsv_file_df, group_by_column_name, bw_table,
                              weights, chr_col, start_col, end_col,
                              magnitude_col) {
  if (!is.data.frame(tsv_file_df))
    stop("`tsv_file_df` must be a data frame.", call. = FALSE)

  required_cols <- c(chr_col, start_col, end_col, group_by_column_name)
  if (!is.null(magnitude_col)) required_cols <- c(required_cols, magnitude_col)
  miss <- setdiff(required_cols, names(tsv_file_df))
  if (length(miss) > 0)
    stop("Required columns missing from `tsv_file_df`: ",
         paste(miss, collapse = ", "),
         "\nUse chr_col/start_col/end_col/magnitude_col/group_by_column_name ",
         "to specify your actual column names.", call. = FALSE)

  if (!is.data.frame(bw_table))
    stop("`bw_table` must be a data frame.", call. = FALSE)

  miss_bw <- setdiff(c("file_path", "sample_id"), names(bw_table))
  if (length(miss_bw) > 0)
    stop("Required columns missing from `bw_table`: ",
         paste(miss_bw, collapse = ", "),
         "\n`bw_table` needs columns 'file_path' and 'sample_id'.",
         call. = FALSE)

  if (!is.numeric(weights) || length(weights) != 3)
    stop("`weights` must be a numeric vector of length 3 (specificity, sensitivity, magnitude).", call. = FALSE)

  if (any(weights < 0))
    stop("All `weights` must be non-negative.", call. = FALSE)

  if (nrow(tsv_file_df) == 0)
    stop("`tsv_file_df` has no rows.", call. = FALSE)

  if (nrow(bw_table) == 0)
    stop("`bw_table` has no rows.", call. = FALSE)

  # Warn if the data has a column literally named "group_col", "chr", "start",
  # "end", or "magnitude" that would clash with internal rename targets
  internal_names <- c("chr", "start", "end", "magnitude", "group_col")
  user_col_args  <- c(chr_col, start_col, end_col, magnitude_col,
                      group_by_column_name)
  clashing <- setdiff(
    intersect(names(tsv_file_df), internal_names),
    user_col_args
  )
  if (length(clashing) > 0)
    warning("Column(s) '", paste(clashing, collapse = "', '"),
            "' in tsv_file_df share names with PeakRankR internal variables. ",
            "Consider renaming them to avoid unexpected behaviour.",
            call. = FALSE)

  # Check no two column arguments point to the same column
  col_args <- c(chr_col, start_col, end_col, magnitude_col, group_by_column_name)
  dupes <- col_args[duplicated(col_args)]
  if (length(dupes) > 0)
    stop("Column name conflict: '", paste(dupes, collapse = "', '"),
         "' is used by more than one column argument. ",
         "Each of chr_col, start_col, end_col, magnitude_col, and ",
         "group_by_column_name must refer to a different column.",
         call. = FALSE)

  invisible(TRUE)
}

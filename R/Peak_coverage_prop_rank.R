#' Peak coverage proportion and ranking function
#'
#' @description For a set of genomic regions, this function:
#' - Takes precomputed per-sample signal from BigWig files (multiBigwig summary)
#' - Normalizes each region's coverage into proportions across samples
#' - Scales those proportions per sample (0–1)
#' - Ranks the normalized proportions per region per sample (descending)
#'
#' @param multiBigwig_summary_df A data.frame from multiBigwig_summary_SS() with coverage per sample
#' @param bw_table A data.frame with columns `bw_path` and `sample_id`
#'
#' @return A data.frame of the same regions with coverage ranks per sample
#' @export

Peak_coverage_rank <- function(multiBigwig_summary_df, bw_table) {
  message("Running peak coverage proportion and rank...")
  
  stopifnot(is.data.frame(multiBigwig_summary_df),
            is.data.frame(bw_table),
            all(bw_table$sample_id %in% colnames(multiBigwig_summary_df)))
  
  sample_cols <- as.character(bw_table$sample_id)
  
  ## Step 1: Normalize to per-peak proportions
  coverage_matrix <- multiBigwig_summary_df[, sample_cols, drop = FALSE]
  row_sums <- rowSums(coverage_matrix, na.rm = TRUE)
  prop_matrix <- sweep(coverage_matrix, 1, row_sums, FUN = "/")
  
  cov_prop_matrix <- multiBigwig_summary_df
  cov_prop_matrix[, sample_cols] <- prop_matrix
  
  ## Step 2: Scale each column to [0, 1]
  range01 <- function(x) {
    rng <- range(x, na.rm = TRUE)
    if (diff(rng) == 0) return(rep(0, length(x)))
    (x - rng[1]) / (rng[2] - rng[1])
  }
  
  scaled_matrix <- as.data.frame(lapply(cov_prop_matrix[, sample_cols, drop = FALSE], range01))
  cov_df_norm <- cov_prop_matrix
  cov_df_norm[, sample_cols] <- scaled_matrix
  
  ## Step 3: Rank each column (higher value gets lower rank)
  ranked_matrix <- as.data.frame(lapply(cov_df_norm[, sample_cols, drop = FALSE],
                                        function(x) rank(-x, ties.method = "max")))
  cov_df_ranked <- cov_df_norm
  cov_df_ranked[, sample_cols] <- ranked_matrix
  
  return(cov_df_ranked)
}

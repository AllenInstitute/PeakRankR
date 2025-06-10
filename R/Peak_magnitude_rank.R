#' Peak signal magnitude scaling and ranking
#'
#' @description For a set of genomic regions, this function:
#' - Takes already summarized signal magnitudes from BigWig files
#' - Scales those signals per sample to [0–1]
#' - Ranks the scaled magnitudes per region per sample (descending)
#'
#' @param multiBigwig_summary_df A data.frame from multiBigwig_summary_SS()
#' @param bw_table A data.frame with columns `bw_path` and `sample_id`
#'
#' @return A data.frame of the same regions with magnitude ranks per sample
#' @export

Peak_magnitude_rank <- function(multiBigwig_summary_df, bw_table) {
  message("Running peak magnitude scaling and rank...")
  
  stopifnot(is.data.frame(multiBigwig_summary_df),
            is.data.frame(bw_table),
            all(bw_table$sample_id %in% colnames(multiBigwig_summary_df)))
  
  sample_cols <- as.character(bw_table$sample_id)
  
  ## Step 1: Scale each column to [0, 1]
  range01 <- function(x) {
    rng <- range(x, na.rm = TRUE)
    if (diff(rng) == 0) return(rep(0, length(x)))
    (x - rng[1]) / (rng[2] - rng[1])
  }
  
  scaled_matrix <- as.data.frame(lapply(multiBigwig_summary_df[, sample_cols, drop = FALSE], range01))
  mag_df_norm <- multiBigwig_summary_df
  mag_df_norm[, sample_cols] <- scaled_matrix
  
  ## Step 2: Rank scaled values per column (higher signal = rank 1)
  ranked_matrix <- as.data.frame(lapply(mag_df_norm[, sample_cols, drop = FALSE],
                                        function(x) rank(-x, ties.method = "max")))
  mag_df_norm_ranked <- mag_df_norm
  mag_df_norm_ranked[, sample_cols] <- ranked_matrix
  
  return(mag_df_norm_ranked)
}

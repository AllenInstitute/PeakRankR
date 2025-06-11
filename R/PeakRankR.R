#' Process BigWig Signal Pipeline to Rank Peaks and Extract Top Regions
#'
#' This function runs a multi-step pipeline to summarize BigWig signal over genomic regions,
#' calculate peak coverage proportion and magnitude ranks, combine these ranks (optionally weighted),
#' and then extract the top N ranked peaks per sample. Optionally, it adds UCSC Genome Browser URLs for visualization.
#'
#' @param bw_table A data.frame or tibble containing BigWig file paths and associated metadata.
#' @param regions_df A data.frame of genomic regions to summarize (must include chr, start, end, Peak_ID).
#' @param summary_type Character string specifying the summary statistic to use for BigWig signal aggregation (default: "mean").
#' @param n_top_values Integer specifying the number of top peaks to extract per sample (default: 500).
#' @param make_ucsc_url Logical indicating whether to add UCSC Genome Browser URL links to the output (default: TRUE).
#' @param base_url Character string of the base URL for UCSC Genome Browser links (no default; required if `make_ucsc_url = TRUE`).
#' @param weights Numeric vector of length 2 to weight coverage and magnitude ranks respectively (default: c(1, 1)).
#'
#' @return A data.frame with top ranked peaks per sample, including columns: chr, start, end, Peak_ID, subclass (sample),
#' and optionally UCSC_URL if `make_ucsc_url = TRUE`.
#'
#' @import dplyr
#' @import tidyr
#' @export
process_bigwig_pipeline <- function(bw_table,
                                    regions_df,
                                    summary_type = "mean",
                                    n_top_values = 500,
                                    make_ucsc_url = TRUE,
                                    base_url = NULL,
                                    weights = c(1, 1)) {
  if (make_ucsc_url && (is.null(base_url) || base_url == "")) {
    stop("Argument 'base_url' must be provided if 'make_ucsc_url' is TRUE.")
  }
  
  if (!is.numeric(weights) || length(weights) != 2) {
    stop("Argument 'weights' must be a numeric vector of length 2 (e.g., c(1, 1))")
  }
  
  message("Running multiBigwig_summary_SS...")
  multi_summary <- multiBigwig_summary_SS(
    data_table = bw_table,
    data_frame_of_regions = regions_df,
    summary_type = summary_type,
    parallel = TRUE
  )
  
  ## Running Peak_coverage_rank
  cov_df_norm_ranked <- Peak_coverage_rank(
    multiBigwig_summary_df = multi_summary,
    bw_table = bw_table
  )
  
  ## Running Peak_magnitude_rank
  mag_df_norm_ranked <- Peak_magnitude_rank(
    multiBigwig_summary_df = multi_summary,
    bw_table = bw_table
  )
  
  message("Combining ranks with weights")
  coord_cols <- c("chr", "start", "end")
  sample_cols <- setdiff(colnames(cov_df_norm_ranked), coord_cols)
  
  ## Apply weights to coverage and magnitude ranks
  sum_mat <- weights[1] * cov_df_norm_ranked[, sample_cols] +
    weights[2] * mag_df_norm_ranked[, sample_cols]
  
  ## Combine with coordinates
  sum_meta <- cbind(cov_df_norm_ranked[, coord_cols], sum_mat)
  
  ## Function to extract top n rows per column
  top_n_rows <- function(column, n = n_top_values) {
    sorted_indices <- order(column)
    sum_meta[sorted_indices[1:min(n, length(sorted_indices))], coord_cols, drop = FALSE]
  }
  
  ## Apply top_n_rows to each sample
  sorted_list <- lapply(sample_cols, function(col) {
    top_n_rows(sum_meta[[col]], n_top_values)
  })
  names(sorted_list) <- sample_cols
  
  ## Add subclass info
  list_with_subclass <- Map(function(df, group) {
    cbind(df, subclass = group)
  }, sorted_list, names(sorted_list))
  
  ## Combine and add rank
  df <- do.call(rbind, lapply(list_with_subclass, function(df_group) {
    df_group$rank <- seq_len(nrow(df_group))
    df_group
  }))
  rownames(df) <- NULL
  
  ## Convert types
  df$start <- as.integer(as.character(df$start))
  df$end <- as.integer(as.character(df$end))
  
  ## Optional UCSC URL
  if (make_ucsc_url) {
    df$UCSC_URL <- paste0(
      base_url,
      df$chr, "%3A", df$start, "%2D", df$end
    )
  }
  
  return(df)
}

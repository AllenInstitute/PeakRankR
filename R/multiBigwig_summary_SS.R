#' Summarize signal from BigWig files over genomic regions
#'
#' @description Computes summary statistics (mean, median, min, or max)
#' across genomic regions from multiple BigWig files.
#' Inspired by deeptools' multiBigwigSummary, bedtoolsr (Phanstiel lab)
#'
#' @param data_table A data.frame with columns `bw_path` and `sample_id`
#' @param data_frame_of_regions A data.frame with columns `chr`, `start`, `end`
#' @param summary_type Summary type: "mean", "median", "min", or "max" (default = "mean")
#' @param parallel Logical. Use parallel processing? Default is TRUE.
#' @return A data.frame with signal summaries per region per sample
#' @export

multiBigwig_summary_SS <- function(data_table,
                                   data_frame_of_regions,
                                   summary_type = "mean",
                                   parallel = TRUE) {
  
  ## Validate input
  stopifnot(is.data.frame(data_table),
            all(c("bw_path", "sample_id") %in% colnames(data_table)))
  
  ## Convert input regions to GRanges
  peaks_gr <- GenomicRanges::makeGRangesFromDataFrame(
    df = data_frame_of_regions,
    seqnames.field = "chr",
    start.field = "start",
    end.field = "end",
    keep.extra.columns = TRUE
  )
  
  ## Prepare BigWig files
  bw_files <- setNames(as.character(data_table$bw_path),
                       data_table$sample_id)
  bw_list <- rtracklayer::BigWigFileList(bw_files)
  
  ## Function to summarize one BigWig over all regions
  summarize_one <- function(sample_id) {
    summary_df <- as.data.frame(rtracklayer::summary(bw_list[[sample_id]], peaks_gr, type = summary_type))
    ## Use :: as separator to avoid splitting seqnames
    summary_df$region <- paste(summary_df$seqnames, summary_df$start, summary_df$end, sep = "::")
    rownames(summary_df) <- summary_df$region
    summary_df <- summary_df["score", drop = FALSE]
    colnames(summary_df) <- sample_id
    return(summary_df)
  }
  
  
  
  ## Apply summarization across all samples
  sample_ids <- names(bw_list)
  
  summary_list <- if (parallel) {
    BiocParallel::bplapply(sample_ids, summarize_one)
  } else {
    lapply(sample_ids, summarize_one)
  }
  
  ## Combine into matrix
  count_mat <- do.call(cbind, summary_list)
  count_mat <- tibble::rownames_to_column(count_mat, "region")
  count_mat <- tidyr::separate(count_mat, region, into = c("chr", "start", "end"), sep = "::", convert = TRUE)
  return(count_mat)
}

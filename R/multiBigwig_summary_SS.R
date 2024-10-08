#' Enrichments at genomics regions (Adapted from ALPS)
#'
#' @description \code{multiBigwig_summary} is a function to calculate
#' enrichments from a set of given bigwig files
#' and a set of genomics regions.
#' This function is similar to
#' \code{deeptools multiBigwigSummary} python package.
#' @param data_table a dataframe that contains \code{bw_path}, \code{sample_id}.
#' \code{sample_id} ids will be used in the final result matrix
#' @param summary_type whether to calculate mean, median, min or max
#' for each genomic region in within consensus peak-set
#' from the bigwig files. Default is \code{mean}
#' Added by SS
#' @param data_frame_of_regions df containing chr, start, end for group of interest
#' @param parallel logical. Whether to parallelize the calculation process,
#' default is \code{TRUE}
#' @return data.frame of enrichments within given genomic regions
#' @export
multiBigwig_summary_SS <- function(data_table,
                                summary_type = "mean",
				# Added by SS
				data_frame_of_regions,
                                parallel = TRUE) {

    assertthat::assert_that(is.data.frame(data_table), msg = "Please provide `data_table`")
    assertthat::assert_that(assertthat::has_name(data_table, "bw_path"))
    assertthat::assert_that(assertthat::has_name(data_table, "sample_id"))
    # Added by SS
    #assertthat::assert_that(assertthat::has_name(data_table, "bed_path"))

    ## create a consensus peak-set from all
    ## files
    # Added by SS
    #all_beds <- data_table$bed_path %>% as.character()
    # Added by SS
    #peaks_gr <- merge_GR(x = all_beds)
    # Added by SS
    peaks_gr <- GenomicRanges::makeGRangesFromDataFrame(df = data_frame_of_regions, seqnames.field = "chr", start.field = "start", end.field = "end", keep.extra.columns = TRUE)

    ## bw list
    all_bw_files <- data_table$bw_path %>% as.character()
    names(all_bw_files) <- data_table$sample_id %>% as.character()

    bwL <- rtracklayer::BigWigFileList(all_bw_files)

    if (parallel) {

        .par_fun <- function(x) {
          
          # Added by SS
          # to give the string as column name to dplyr
           #x <- sym(x)
          x_bw <- bwL[[x]]
            
            x_res <- rtracklayer::summary(x_bw,
                peaks_gr, type = summary_type) %>%
                as.data.frame() %>% dplyr::mutate(region = paste(seqnames, start, end, sep = "_")) %>%
                tibble::column_to_rownames(var = "region") %>%
                dplyr::select(score) 
            #%>%
             #   dplyr::rename(`:=`(!!x, score))
            # Added by SS
              colnames(x_res) <- x

            return(x_res)
        }

        all_pid <- bwL %>% names

        all_pid_reslst <- BiocParallel::bplapply(all_pid, .par_fun)

        count_mat <- all_pid_reslst %>% as.data.frame()
          # Added by SS
          colnames(count_mat) <- all_pid
        
          count_mat <- count_mat %>% tibble::rownames_to_column(var = "region") %>%
            tidyr::separate(region, into = c("chr", "start", "end"))

    } else {

        count_mat <- peaks_gr %>% as.data.frame() %>%
            dplyr::mutate(region = paste(seqnames, start, end, sep = "_")) %>%
            dplyr::select(region)

        for (i in 1:length(bwL)) {

            sample_id <- bwL[i] %>% names
            sample_path <- bwL[[i]]

            sample_res <- rtracklayer::summary(sample_path, peaks_gr, type = summary_type) %>%
              as.data.frame() %>%
              dplyr::mutate(region = paste(seqnames, start, end, sep = "_")) %>%
              dplyr::select(region, score) %>%
              dplyr::rename(`:=`(!!sample_id, score))

            suppressMessages(count_mat <- dplyr::left_join(count_mat, sample_res, by = "region"))
        }
        count_mat <- count_mat %>%
          tidyr::separate(region, into = c("chr", "start", "end"))
    }
    return(count_mat)
}


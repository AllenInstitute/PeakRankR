#' A Peak_annotation_rank function 
#'
#' This function allows you to extract the peaks for a group by category along with its MACS2 rank from the CERP tsv file and places the output in a new directory called peak_subset_files with two files per group, 
#' one as is and one sorted using bedtools in two directories 
#' 
#' @param peak_table.list
#' @param group_by_column_name
#'
#' @keywords MACS2
#'
#' @export
peak_annotation_rank = function(peak_table, group_by_column_name, rank_columns=c("rank")){
    ##
    print("Computing peak annotation / MACS2 ranks.")
    ##
    peak_table.list = list()
    ## Loop through all groups and build grange objects to be used for next steps
    for(group in unique(peak_table[[group_by_column_name]])){
            ## Gather Rank info for group of interest
            peak_table.group = peak_table %>% filter(.data[[group_by_column_name]] == group) %>% select(union(c("chr", "start", "end"), rank_columns))
            peak_table.group = peak_table.group %>% dplyr::rename_at(rank_columns, ~ paste0(gsub("\\.","_", rank_columns), "_rank_up"))
            ##
            peak_table.list[[group]] = makeGRangesFromDataFrame(peak_table.group, keep.extra.columns=TRUE)
    }
    ##
    return(peak_table.list)
}
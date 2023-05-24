#' A Peak_intersect_rank function 
#'
#' This Peak_intersect() function intersects each group peaks by background_group peak files using bedtools and 
#' Peak_intersect_rank() ranks for specificity (group peaks with low intersection). The results are stored in peak_intersect_rank folder.
#' 
#' @param peak_table.list
#'
#' @keywords bedtools intersect
#'
#' @export
peak_intersect = function(peak_table.list){
    ##
    print("Computing peak intersect rank.")
    ##
    peak_table.overlaps = list()
    ## Intersect each group against all others
    for(group in names(peak_table.list)){
        ## Build a list of granges for all background groups
        background_groups = GRangesList(peak_table.list[-which(names(peak_table.list) == group)])
        ## Intersect group of interest against all background groups and count the number of hits
        peak_table.list[[group]]$peak_intersect_rank_down = countOverlaps(peak_table.list[[group]], background_groups)
    }
    return(peak_table.list)
}
#' A Peak_coverage_rank function 
#'
#' Peak_coverage_rank() function ranks (descending) the coverage proportion per peak of a group across groups 
#' Needs the bigwig files  
#' 
#' @param tsv_file_df CERP tsv_file_df REQUIRED
#' @param group_by_column column name in the tsv_file that has the groups 
#' @param group unique group in the tsv_file that needs to be ranked
#' @param background_group list of groups from tsv_file with which the group peaks should be intersected 
#' @param bw_table table of sample name and bigwig files
#' @export
Peak_coverage_rank <- function(tsv_file_df, group_by_column,  group, background_group, bw_table){
  
  print("running peak coverage rank...")
  tsv_file_df <- as.data.frame(tsv_file_df)
  # df of coordinates for group
  bed_regions <- tsv_file_df[tsv_file_df[[group_by_column]] == group, c("chr","start","end",group_by_column)]
  
  # coverage matrix for coordinates
  final_coverage <- multiBigwig_summary_SS(data_table = bw_table, summary_type = "mean", data_frame_of_regions = bed_regions, parallel = TRUE  )

  
   # prop calculation
  # first three are chr, start, end
  final_coverage$coverage_sum <- rowSums(final_coverage[,c(4:dim(final_coverage)[2])])
  final_coverage$coverage_prop <- final_coverage[,colnames(final_coverage) == group]/final_coverage$coverage_sum
  setDT(final_coverage)[order(-coverage_prop), peak_cov_rank := rleid(coverage_prop)]
 

  # sorted output
  final_coverage$cell.population <- group
  final_coverage_df <- as.data.frame(final_coverage)
  cov_df <- final_coverage_df[,c("chr","start","end",group_by_column,"peak_cov_rank")]
  cov_df_sorted <- bt.sort(cov_df)
  colnames(cov_df_sorted) <- colnames(cov_df)

 
  return(cov_df_sorted)
  
  
}
  

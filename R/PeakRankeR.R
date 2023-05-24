#' A Peak_RankeR function 
#'
#' This function calculates the final rank for each peak in tsv file given the columns of *ranks produced using PeakRankeR package 
#' 
#' @param rank_df data_frame which has the peaks for a specific group and newly added rank columns generated using PeakRankeR package
#'
#' @keywords PeakRankeR 
#'
#' @export
PeakRankeR = function(peak_table, group_by_column_name, bigwig.dir){

  ##
  peakRankeR.gr = peak_annotation_rank(peak_table, "cell.population")
  peakRankeR.gr = peak_intersect(peakRankeR.gr)
  peakRankeR.gr = peak_coverage(peakRankeR.gr, peak_table, bigwig.dir)

  ## Convert the list of granges into a data.frame
  peakranker.df = grList_to_df(peakRankeR.gr)

  ## Normalizing rank
  final_df_norm = peakranker.df %>% 
                    group_by(cell.population) %>% 
                    mutate_at(vars(grep("rank_down", colnames(peakranker.df))), function(x) {x=rank(x); range01(x)}) %>%
                    mutate_at(vars(grep("rank_up", colnames(peakranker.df))), function(x) {x=rank(-x); range01(x)})

  # ## Adding extra weight maybe a good option in second release???
  final_df_norm$rank_sum = rowSums(final_df_norm[,grep("rank", colnames(peakranker.df)),drop=F])

  # ##
  # final_df_norm = final_df_norm %>%
  #   group_by(cell.population) %>% 
  #   arrange(rank_sum, .by_group = TRUE) %>% as.data.frame()

  return(final_df_norm)
}

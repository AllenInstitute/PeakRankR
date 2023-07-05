#' A Peak_intersect_rank function 
#'
#' Peak_intersect_rank() ranks for specificity (peaks of a group having low intersection with other groups). 
#' 
#' @param tsv_file_df CERP tsv_file_df REQUIRED
#' @param group_by_column column name in the tsv_file that has the groups 
#' @param group unique group in the tsv_file for that needs to be ranked 
#' @param background_group list of groups from tsv_file with which the group peaks should be intersected 
#' @export

Peak_intersect_rank <- function(tsv_file_df, group_by_column_name, group, background_group){
  
  print("running peak intersect rank...")
  
  # making sure we work on a df
  df <- as.data.frame(tsv_file_df)
  
  # group of interest bed frame
  a <- df[df[[group_by_column_name]] == group,c("chr","start","end",group_by_column_name) ]
  a_sort <- bt.sort(a)
  colnames(a_sort) <- colnames(a)
  
  # background bedframe
  b <- df[df[[group_by_column_name]] %in% background_group, c("chr","start","end") ]
  b_sort <- bt.sort(b)
  colnames(b_sort) <- colnames(b)
  
  # intersection
  df_intersected <- bt.intersect(a_sort,b_sort,wao = TRUE)
  
  # count and rank
  agg_df <- df_intersected  %>% 
    group_by(V1,V2,V3,V4) %>%
    count() %>% 
    as.data.frame()
  colnames(agg_df) <-  c("chr","start","end",group_by_column_name,"peak_intersect_rank")
  
  # return the df per group with the intersect rank
  return(agg_df)
  
}

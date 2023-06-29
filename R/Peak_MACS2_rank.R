#' A Peak_MACS2_rank function 
#'
#' Add rank to the tsv file based on peak magnitude given by MACS2
#' 
#' @param tsv_file_df CERP tsv_file_df REQUIRED
#' @param group_by_column column name in the tsv_file that has the groups 
#' @keywords MACS2
#' @export
#' @examples 
#' Peak_MACS2_rank()
#' 

#########################################################################################################################################################


Peak_MACS2_rank <- function(tsv_file_df, group_by_column_name){
  
  print("Running peak rank...")
  # to give the string as column name to dplyr
  group_by_column_name <- sym(group_by_column_name)
  
  # making sure we work on a df
  df <- as.data.frame(tsv_file_df)
  
  # assign ranks to EVERY group based on magnitude
  df_final <- df %>%
    group_by(!!group_by_column_name) %>%
    mutate(MACS2_rank = as.numeric(as.factor(rank(rank)))) %>%
      #-peak.magnitude)))) %>%
    as.data.frame()
  
  # Returning the entire table as is with the new column
  return(df_final)
  
}
 

#' range function 
#'
#' This function scales a vector ranging from 0 - 1
#'
#' @param x vector that needs to be scaled
#' @export
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

#' Peak_RankeR function
#'
#' This function calculates the final rank for group of interest against background group 
#' 
#' @param tsv_file_df CERP tsv_file_df REQUIRED
#' @param group_by_column column name in the tsv_file that has the groups 
#' @param background_group list of groups from tsv_file with which the group peaks should be intersected 
#' @param bw_table table of sample name and bigwig files
#' @param rank_sum if the rank sum (0 - 1) should be provided in the output
#' @param weights vector of the features MACS2, intersect, coverage in that order
#' @export 
Peak_RankeR <- function(tsv_file_df, group_by_column_name, background_group, bw_table, rank_sum, weights){
   #print(group)
  weights_vector <- weights
  
   # MACS2 rank
   MACS2_df <- Peak_MACS2_rank(tsv_file_df, group_by_column_name)
   MACS2_df_sub <- MACS2_df[,c("chr","start","end",group_by_column_name,"MACS2_rank")]
  
   
   # intersect rank
   # intersected_df <- Peak_intersect_rank(tsv_file,group_by_column_name,group, background_group)
   intersected_df_list <- list()
   intersected_df_list <- lapply(unique(tsv_file_df[[group_by_column_name]]), function(x) Peak_intersect_rank(tsv_file_df,group_by_column_name,x, background_group))
   intersected_df <-  bind_rows(intersected_df_list, .id = "column_label")

   
    # cov rank
    # cov_df <- Peak_coverage_rank(tsv_file, group_by_column_name,group ,background_group, bw_table)
    cov_df_list <- list()
    cov_df_list <- lapply(unique(tsv_file_df[[group_by_column_name]]), function(x) Peak_coverage_rank(tsv_file_df, group_by_column_name,x ,background_group, bw_table))
    cov_df <-  bind_rows(cov_df_list, .id = "column_label")
    

   df1 <- MACS2_df %>%
      left_join(intersected_df) 
    #%>%
     #select("chr","start","end",group_by_column_name,"MACS2_rank","peak_intersect_rank")
   # interchanged MACS2 and intersect df
    PP_df <- df1 %>%
     left_join(cov_df)

   
    # From here
    # Commenting from here for test
     PP_df_norm <- PP_df %>% group_by(cell.population) %>%
     mutate_at(c("MACS2_rank","peak_intersect_rank","peak_cov_rank"), function(x) range01(x)) %>%
      mutate_all(~ifelse(is.nan(.), 0, .))
     
  
     
    # multiplying weights
     sum_of_weights <- sum(weights)
     # Divide by sum_of_weights
     PP_df_norm$MACS2_rank <- (weights_vector[1]/sum_of_weights) * PP_df_norm$MACS2_rank
     PP_df_norm$peak_intersect_rank <- (weights_vector[2]/sum_of_weights) * PP_df_norm$peak_intersect_rank
     PP_df_norm$peak_cov_rank <- (weights_vector[3]/sum_of_weights) * PP_df_norm$peak_cov_rank

     # multiplying weights
     #PP_df_norm$MACS2_rank <- (weights_vector[1]) * PP_df_norm$MACS2_rank
     #PP_df_norm$peak_intersect_rank <- (weights_vector[2]) * PP_df_norm$peak_intersect_rank
     #PP_df_norm$peak_cov_rank <- (weights_vector[3]) * PP_df_norm$peak_cov_rank
     

    PP_df_norm$rank_sum  <- rowSums(PP_df_norm[c("MACS2_rank","peak_intersect_rank","peak_cov_rank")])
    PP_df_norm <- as.data.frame(PP_df_norm)
    
    PP_df_final <- PP_df_norm %>% 
      group_by(cell.population) %>%
      mutate(PeakRankeR_rank = as.numeric(as.factor(frank(rank_sum, ties.method = "min"))))
   
    PP_df_final <- as.data.frame(PP_df_final)
    if(rank_sum) {
      PP_df_norm_sub <- PP_df_final[,!(names(PP_df_final) %in% c("MACS2_rank","peak_intersect_rank","peak_cov_rank","column_label"))]
      return(PP_df_norm_sub)
    }else{
      PP_df_norm_sub <- PP_df_final[,!(names(PP_df_final) %in% c("MACS2_rank","peak_intersect_rank","peak_cov_rank","column_label","rank_sum"))]
      return(PP_df_norm_sub)
     }

}


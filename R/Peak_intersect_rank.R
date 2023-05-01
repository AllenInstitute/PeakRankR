#' A Peak_intersect_rank function 
#'
#' This Peak_intersect() function intersects each group peaks by background_group peak files using bedtools and 
#' Peak_intersect_rank() ranks for specificity (group peaks with low intersection). The results are stored in peak_intersect_rank folder.
#' 
#' @param group unique group in the tsv_file for which the ranked peak intersect file needs to be created
#' @param background_group list of groups from tsv_file with which the group peaks should be intersected 
#' @param group_intersected_file_df the intersected file for group of interest produced by bedtools 
#' @keywords bedtools intersect
#' @export
#' @examples 
#' Peak_intersect()
#' Peak_intersect_rank()

#########################################################################################################################################################

# setwd("/allen/programs/celltypes/workgroups/rnaseqanalysis/sarojaS/230418_PeakRanker_package/PRP/R")
# library(dplyr)

#########################################################################################################################################################


# Function 1
Peak_intersect <- function(group, background_group){
  
  print(group)
 
  intersecting_groups <-  paste0("./peak_subset_ranks/",background_group,".peak.rank.tsv")
  intersecting_groups_string <- paste( unlist(intersecting_groups), collapse=' ')
  # intersecting sorted peak_files
  cmd2 <- paste0("bedtools intersect -a ./peak_subset_ranks/", group, ".peak.rank.tsv -b ", intersecting_groups_string," -wao -filenames > ./peak_intersect_files/", group, ".intersect.tsv.tmp", collapse = "")
  system(cmd2)
}

# Function 2
# Input of the intersected file as a df object
Peak_intersect_rank <- function(group_intersected_file_df){
 
  # Taking only the first 4 columns, chr, start, end, MACS2 rank
  df_intersect_sub <- group_intersected_file_df[,c(1,2,3,4)]
  
  agg_tbl <- df_intersect_sub %>% group_by(V4) %>% 
    summarise(total_count=n(),
              .groups = 'drop')
  df2 <- agg_tbl %>% as.data.frame()
  colnames(df2) <- c("V4","peak_intersect_count")
  
  df_intersect_final <- df_intersect_sub %>%
    left_join(df2)
  peak_intersect_count_df <- unique(df_intersect_final)
  colnames(peak_intersect_count_df) <- c("chr","start","end","peak_rank","peak_intersect_rank")
  
  return(peak_intersect_count_df)
}

#########################################################################################################################################################

# testing the function

# Peak_intersect("Trhr",c("Astrocytes-2","Astrocytes-1","Trhr"))
# test_df <- read.csv("./peak_intersect_files/Trhr.intersect.tsv.tmp", sep = "\t",header = F)
# test_result <- Peak_intersect_rank(test_df)
# write.table(test_result, "./peak_intersect_rank/Trhr.MACS2.intersect.rank.tsv", sep = "\t", col.names = T, row.names = F, quote = F)

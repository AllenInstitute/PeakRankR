#' A Peak_MACS2_rank function 
#'
#' This function allows you to extract the peaks for a group by category along with its MACS2 rank from the CERP tsv file and places the output in a new directory called peak_subset_files with two files per group, 
#' one as is and one sorted using bedtools in two directories 
#' 
#' @param tsv_file_df CERP tsv_file_df REQUIRED
#' @param group_by_column column name in the tsv_file that has the groups 
#' @param group unique group in the tsv_file for which the peak rank file needs to be created, for eg. Pvalb
#' @keywords MACS2
#' @export
#' @examples 
#' Peak_MACS2_rank()
#' 

#########################################################################################################################################################

# setwd("/allen/programs/celltypes/workgroups/rnaseqanalysis/sarojaS/230418_PeakRanker_package/PRP/R")

# To do:
# change the bedtools to granges and don't write out files

#########################################################################################################################################################



Peak_MACS2_rank <- function(tsv_file_df, group_by_column_name, group){
  
  print(group)
  
  # Check for the column!
  # change it to to a df explicitly
  df <- tsv_file_df[tsv_file_df[[group_by_column_name]] == group, c("chr","start","end","rank")]
  
  # writing peak rank files
  write.table(df, file = paste0("./peak_subset_files/", group, ".peak.rank.tsv.tmp", collapse = ""), sep = "\t", row.names = F, col.names = F, quote = F)
  
  cmd1 <- paste0("bedtools sort -i ./peak_subset_files/", group, ".peak.rank.tsv.tmp > ./peak_subset_ranks/", group, ".peak.rank.tsv", collapse = "")
  system(cmd1)
  
  # error if no rows (warning/rows function)
}
 
#########################################################################################################################################################

# testing the function

# test_file <- read.csv("./../../PeakRankeRPackage/subclass_annotated_markerPeaks.tsv", sep = "\t", header = T)
# Peak_MACS2_rank(test_file, "cell.population", "Trhr")
#########################################################################################################################################################

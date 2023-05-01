#' A Peak_coverage_rank function 
#'
#' The Peak_coverage() function calculates coverage for a groups' peak  across background groups list provided using multiBigWigSummary..
#' and Peak_coverage_rank() function ranks (ascending) the coverage proportion per peak of a group across groups 
#' Needs the bigwig files  
#' 
#' @param group unique group  in the tsv_file for which the coverage across groups needs to be extracted 
#' @param background_group list of groups from tsv_file from which coverage needs to be extracted for the gorup peaks
#' @param group_coverage_df result of multiBigWigSummary for group of interest
#' @param group_bigwig_file_name group of interest bigwig file name 
#' @param group_macs2_intersect_rank_df df of ranks calculated for other features
#' @keywords multiBigWigSummary 
#' @export
#' @examples 
#' Peak_coverage()
#' Peak_coverage_rank()

#########################################################################################################################################################

# setwd("/allen/programs/celltypes/workgroups/rnaseqanalysis/sarojaS/230418_PeakRanker_package/PRP/R")
# library(data.table)
# library(dplyr)

#########################################################################################################################################################



Peak_coverage <- function(group, background_group){
  
  peak_file <- paste0("./peak_subset_ranks/",group,".peak.rank.tsv")
  
  # Name shud be consistent with tsv file for bigwigs
  bwfiles <- paste0("./bwfiles/",background_group,"*.bw")
  bwfiles_string <- paste( unlist(bwfiles), collapse=' ')
  
  print(group)
  
  cmd3 <- paste0(" multiBigwigSummary BED-file  -b ", bwfiles_string, " --outFileName ./peak_coverage_files/", group, ".npz --outRawCounts ./peak_coverage_files/", group, ".scores_per_bin.tab  --BED ", peak_file , collapse = "")
  system(cmd3)
  
}
  
Peak_coverage_rank <- function(group_coverage_df,group_bigwig_file_name, group_macs2_intersect_rank_df){
  
  # for R column name match
   subclass_edit1 <- gsub("-",".",group_bigwig_file_name)
   subclass_edit2 <-gsub("_",".",subclass_edit1)
   subclass_final <- gsub("$",".",subclass_edit2)
  

  cov_df$coverage_sum <- rowSums(group_coverage_df[,c(4:dim(group_coverage_df)[2])])
  print(colnames(cov_df)[grepl(subclass_final, colnames(cov_df), fixed=TRUE)])
  cov_df$coverage_prop <- cov_df[,grepl(subclass_final, colnames(cov_df), fixed=TRUE)]/cov_df$coverage_sum
  setDT(cov_df)[order(-coverage_prop), peak_cov_rank := rleid(coverage_prop)]
  colnames(cov_df)[c(1:3)] <- c("chr","start","end")
  
  # joining the rank
  final_df <- group_macs2_intersect_rank_df %>%
    left_join(cov_df)
  final_df_sub <- final_df[,c("chr","start","end","peak_rank","peak_intersect_rank","peak_cov_rank")]
  return(final_df_sub)
  
  
}


#########################################################################################################################################################

# testing the function

# Peak_coverage("Trhr", c("Astrocytes-1","Trhr", "Astrocytes-2"))

# notice header is TRUE!
# test_df <- read.csv("./peak_coverage_files/Trhr.scores_per_bin.tab", sep = "\t", header = T)
# test_rank_df <- read.csv("./peak_intersect_rank/Trhr.MACS2.intersect.rank.tsv", sep = "\t", header = T)

# result_df <- Peak_coverage_rank(test_df, "Trhr-TileSize-25-normMethod-ReadsInTSS-ArchR.bw", test_rank_df)
# write.table(result_df, "./peak_coverage_rank/Trhr.MACS2.intersect.coverage.rank.tsv", sep = "\t", col.names = T, row.names = F, quote = F)

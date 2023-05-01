# Take tsv file
# Add a column
# populate with NA
# For group input, get the rank files and normalize the rank and rank


# Improvements:

# bed files as input rather than file names (named list of files)
# how to generalize it for second version for features
# Give me the df object for any number of features to calculate the rank,
# take the df and add the columns to it or make a df additive

# For version 1
# bedtools to granges
# only bigwig files needed give the ArchR project and point to bigwig file..
# Test data: /allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/NHP

####################################################################################################################################

#' A Peak_RankeR function 
#'
#' This function calculates the final rank for each peak in tsv file given the columns of *ranks produced using PeakRankeR package 
#' 
#' @param rank_df data_frame which has the peaks for a specific group and newly added rank columns generated using PeakRankeR package
#' @keywords PeakRankeR 
#' @export
#' @examples 
#' Peak_RankeR()
#' range01()
#' 

# rank range
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

Peak_RankeR <- function(group_rank_df){
 
   # colnames(rank_df) <- c("chr","start","end","cov_rank","intersect_rank","peak_rank")
 

  # Normalizing rank
  final_df_norm <- group_rank_df %>%
    mutate_at(vars(colnames(group_rank_df[c(4:dim(group_rank_df)[2])])), function(x) range01(x))
  
  # Adding extra weight maybe a good option in second release???
  
  final_df_norm$rank_sum <- rowSums(final_df_norm[c(4:dim(final_df_norm)[2])])
  final_df_norm_sort <- final_df_norm[order(final_df_norm$rank_sum),]
  
  return(final_df_norm_sort[,c(1,2,3)])

}

####################################################################################################################################

# Testing the function

# PRP <- Peak_RankeR(result_df)


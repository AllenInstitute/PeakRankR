# PeakRankeR

The package consists of functions to rank enhancer peaks in a group of cells given a tsv file input from CERP pipeline. 

Install:

   Clone this github and source the files

The workflow is as follows:

COMING SOON ....

Steps to run the script to get a ranked list of peaks in a group for cloning 

1. Install the package from here
2. Set the working directory containing the bwfiles under ./bwfiles/*.bw
3. Create directories peak_subset_files, peak_subset_ranks, peak_coverage_files
4. Run the following lines of code:
   
   Note: subclass_annotated_markerPeaks.tsv is available here for reference, We are ranking peaks for "Trhr" (subclass) in the file

        # Get the MACS2 ranks for Trhr
        
        test_file <- read.csv("./subclass_annotated_markerPeaks.tsv", sep = "\t", header = T)
        Peak_MACS2_rank(test_file, "cell.population", "Trhr")
        
        # Get the intersect ranks for Trhr
        # To get the intersect rank for a group e.g.Trhr against a background group apply Peak_MACS2_rank on the background groups
        
        Peak_intersect("Trhr",c("Astrocytes-2","Astrocytes-1","Trhr"))
        test_intersect_df <- read.csv("./peak_intersect_files/Trhr.intersect.tsv.tmp", sep = "\t",header = F)
        test_intersect_result <- Peak_intersect_rank(test_intersect_df)
        write.table(test_intersect_result, "./peak_intersect_rank/Trhr.MACS2.intersect.rank.tsv", sep = "\t", col.names = T, row.names = F, quote = F)

        # Get the coverage ranks for Trhr
        
        Peak_coverage("Trhr", c("Astrocytes-1","Trhr", "Astrocytes-2"))
        test_coverage_df <- read.csv("./peak_coverage_files/Trhr.scores_per_bin.tab", sep = "\t", header = T)
        test_coverage_rank_df <- read.csv("./peak_intersect_rank/Trhr.MACS2.intersect.rank.tsv", sep = "\t", header = T)
        result_df <- Peak_coverage_rank(test_coverage_df, "Trhr-TileSize-25-normMethod-ReadsInTSS-ArchR.bw", test_coverage_rank_df)
        
        # Calculate final ranks and produce ranked peak coordinates
        
        PRP <- Peak_RankeR(result_df)
        

        
        

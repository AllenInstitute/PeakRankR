# PeakRankR

The R Package can be used to prioritize a list of enahncers/peaks from different groups (e.g. celltypes, subclasses) for cloning and targeting the group of interest. It takes in tsv file with coordinates and group (at bare minimum) and two column file (refer below table for sample) listing the bigwig file path and sample id as input and returns the same tsv file with peak ranks calculated per group as output

## Installation

Clone this github repo or install via: 

install.packages("/allen/programs/celltypes/workgroups/rnaseqanalysis/sarojaS/230418_PeakRanker_package/PeakRankeR_2023-06-08/PeakRankeR/PeakRankR_0.0.0.9000.tar.gz", repos = NULL) 

## Workflow

The workflow is as follows:

![workflow](https://github.com/AllenInstitute/PeakRankeR/blob/master/workflow.png)

## Pictorial representation on scaling for rank sum


![scaling](https://github.com/AllenInstitute/PeakRankeR/blob/master/scaling.png)

# Run

```
tsv_file <- "input peaks file with coordinates and group name (cell.population) columns"
bw_table <- "path to bigwig table" (example given below)
Ranked_peaks_file <- Peak_RankeR(tsv_file_df         = ArchR_tsv_file,
				group_by_column_name = "cell.population",
				background_group     = unique(ArchR_tsv_file$"group_by_column"),
				bw_table             = example_bw_table, 
				rank_sum             = TRUE,
				weights              = c(1,1,1)))
```

example_bw_table:

bw_path                                    | sample_id
-------------------------------------------| -------------
/allen/programs/celltypes/Astrocytes-1.bw  | Astrocytes-1
/allen/programs/celltypes/Astrocytes-2.bw  | Astrocytes-2


       
 

        

        

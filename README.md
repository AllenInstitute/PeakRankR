# PeakRankR

The R Package can be used to prioritize a list of enahncers/peaks from different groups (e.g. celltypes, subclasses) for cloning and targeting the group of interest. It takes in tsv file with coordinates only/and group (at bare minimum) and two column file (refer below table for sample) listing the bigwig file path and sample id as input and returns the same tsv file with peak ranks calculated per group as output

## Installation

git clone git@github.com:AllenInstitute/PeakRankR.git


install.packages("devtools")


library(devtools)


build("~/put/the/package/path/here")

## Algorithm

PeakRankRscore  =  W(specificity)SpecificityPeak +
      	           W(sensitivity)SensitivityPeak +	 
      		   W(magnitude)MagnitudePeak


where W stands for the weight of each feature. By default each weight variable is set to 1 indicating equal importance for all three features

## Run

```
tsv_file <- "input peaks file with coordinates only/and group name (cell.population) columns"
bw_table <- "path to bigwig table" (example given below)


If group name is given:

Ranked_peaks_file <- Peak_RankeR(tsv_file_df         = ArchR_tsv_file,
				group_by_column_name = "cell.population",
				background_group     = unique(ArchR_tsv_file$"group_by_column"),
				bw_table             = example_bw_table, 
				rank_sum             = TRUE,
				weights              = c(1,1,1))
```

example_bw_table:

bw_path                                    | sample_id
-------------------------------------------| -------------
/allen/programs/celltypes/Astrocytes-1.bw  | Astrocytes-1
/allen/programs/celltypes/Astrocytes-2.bw  | Astrocytes-2

Note: 1. The sample_id column should match the group_by_column_name values
      2. All arguments to the function are mandatory

       
 

        

        

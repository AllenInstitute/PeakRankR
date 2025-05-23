[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15238528.svg)](https://doi.org/10.5281/zenodo.15238528)

# PeakRankR

The R Package can be used to prioritize a list of enahncers/peaks from different groups (e.g. celltypes, subclasses) for cloning and targeting the group of interest. It takes in tsv file with coordinates, group and magnitude (at bare minimum) and two column file (refer below table for sample) listing the bigwig file path and sample id as input and returns the same tsv file with peak ranks calculated per group as output

## Required tools to install and run PeakRankR

```
# Identiy location of bedtools using:
system("which bedtools")
# If not installed please install bedtools using https://bedtools.readthedocs.io/en/latest/content/installation.html#
options(bedtools.path = "/path/to/")

# bedtoolsr installation 
install.packages("devtools")
library(devtools)
devtools::install_github("PhanstielLab/bedtoolsr")

```
## Installation of PeakRanR

```
devtools::install_github("AllenInstitute/PeakRankR", dependencies = T)
library(PeakRankR)

```


## PeakRankR Algorithm

PeakRankRscore  =  W(specificity)SpecificityPeak +
      	           W(sensitivity)SensitivityPeak +	 
      		   W(magnitude)MagnitudePeak


where W stands for the weight of each feature. By default each weight variable is set to 1 indicating equal importance for all three features

## Running PeakRankR

```
tsv_file <- read.table("test_file.tsv",header=TRUE)  # input peaks file with coordinates only/and group name (cell.population) columns (example: test_file.tsv)
bw_table <- read.table("bw_table.txt",header=TRUE) # path to bigwig table (example: bw_table)


# If group name is given:

Ranked_peaks_file <- Peak_RankeR(tsv_file_df         = tsv_file,
				group_by_column_name = "cell.population",
				background_group     = unique(tsv_file$"cell.population"),
				bw_table             = bw_table, 
				rank_sum             = TRUE,
				weights              = c(1,1,1))

```


### Note: 

1. In the bw_table file, the sample_id column should match the group_by_column_name values

2. All arguments to the function are mandatory

       
## License
The license for this package is available on Github at: https://github.com/AllenInstitute/PeakRankR/blob/master/LICENSE

## Level of Support
We are planning on occasionally updating this repo with no fixed schedule, but likely several times per year. Community involvement is encouraged through both issues and pull requests. 

        

        

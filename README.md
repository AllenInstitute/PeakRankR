[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15238528.svg)](https://doi.org/10.5281/zenodo.15238528)

# PeakRankR

The R Package can be used to prioritize a list of enahncers/peaks from different groups (e.g. celltypes, subclasses) for cloning and targeting the group of interest. It takes in 3 coloumn bed file with coordinates (header required) and two column file listing the bigwig file path and sample id as input and returns the tsv file with n_top (user provided) prioritized peaks for each group

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
git clone --branch saroja-updated-branch https://github.com/AllenInstitute/PeakRankR.git
devtools::install(".") # From the directory with the description file

```


## PeakRankR Algorithm

PeakRankRscore  =  W(specificity)SpecificityPeak + W(magnitude)MagnitudePeak


where W stands for the weight of each feature. By default each weight variable is set to 1 indicating equal importance for both features

## Running PeakRankR

```
tsv_file <- read.table("test_file.tsv",header=TRUE)  # input 3 column bed file with coordinates and header. If no header, then add colnames c("chr", "start, "end")
bw_table <- read.table("bw_table.txt",header=TRUE) #  bigwig and sample table example (bw_table)




PeakRankR::process_bigwig_pipeline(bw_table = bw_table,
                                   regions_df = tsv_file,
                                   summary_type = "mean",
                                   n_top_values = 5,
                                   make_ucsc_url = FALSE,
                                   weights = c(1,1))

```


       
## License
The license for this package is available on Github at: https://github.com/AllenInstitute/PeakRankR/blob/master/LICENSE

## Level of Support
We are planning on occasionally updating this repo with no fixed schedule, but likely several times per year. Community involvement is encouraged through both issues and pull requests. 

        

        

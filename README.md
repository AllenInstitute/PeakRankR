# PeakRankeR

The R Package can be used to prioritize a list of enahncers/peaks from different groups (e.g. celltypes, subclasses) for cloning and targeting the group of interest. It takes in tsv file with coordinates and group (at bare minimum) and two column file (refer test directory for sample) listing the bigwig file path and sample id as input and returns the same tsv file with peak ranks calculated per group as output

## Installation

Clone this github repo and source the files or 

--COMING SOON--

## Workflow

The workflow is as follows:

![workflow](https://github.com/AllenInstitute/PeakRankeR/blob/master/workflow.png)

## Pictorial representation on scaling for rank sum


![scaling](https://github.com/AllenInstitute/PeakRankeR/blob/master/scaling.png)

# Run

```
tsv_file <- "input peaks file with coordinates and group name columns
bw_table <- "path to bigwig table"
Ranked_peaks_file <- Peak_RankeR(tsv_file,
				"group by column",
				unique(tsv_file$"group_by_column"),
				bw_table, 
				TRUE,
				c(1,1,1)))
```
Replace the string in "" with the right values
       
 

        

        

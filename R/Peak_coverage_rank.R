#' A peak_coverage_rank function 
#'
#' The Peak_coverage() function calculates coverage for a groups' peak  across background groups list provided using multiBigWigSummary..
#' and Peak_coverage_rank() function ranks (ascending) the coverage proportion per peak of a group across groups 
#' Needs the bigwig files  
#' 
#' @param peak_table.list
#' @param peak_table
#' @param bigwig.dir
#'
#' @keywords multiBigWigSummary 
#'
#' @export
peak_coverage = function(peak_table.list, peak_table, bigwig.dir){

  ##
  print("Computing peak coverage rank.")
  
  ## Create BED file for deeptools
  peak.file = file.path(bigwig.dir, "peak.bed")
  ## Remove duplicates and keep only essentials for coverage calculation (chr, start, end)
  write.table(peak_table %>% filter(!duplicated(coordinates)) %>% select(chr, start, end), sep="\t", row.names=F, col.names=F, quote=F, file=peak.file)
  
  ## Name should be consistent with tsv file for bigwigs (some big assumptions here!)
  bwfiles = list.files(bigwig.dir, pattern="*.bw")
  bwfiles = paste( unlist(bwfiles), collapse=' ')
  bwfiles = file.path(bigwig.dir, bwfiles)
  
  ## Run multiBigwigSummary, use labels to specify correct subclass
  cmd = paste0(" multiBigwigSummary BED-file  -b ", bwfiles, 
                    " --outFileName ", file.path(bigwig.dir, "peakCoverages.npz"), 
                    " --outRawCounts ", file.path(bigwig.dir, "peakCoverages.scores_per_bin.tab"), 
                    " --BED ", peak.file , 
                    collapse = "")
  system(cmd)

  ## Load in results
  coverage.df = read_tsv(file.path(bigwig.dir, "peakCoverages.scores_per_bin.tab")) %>% as.data.frame()
  colnames(coverage.df) = gsub("'|#", "", colnames(coverage.df))

  ## Assumptions in file names again... We should have the user pass a list of named bigwigs???
  colnames(coverage.df) = unlist(lapply(strsplit(colnames(coverage.df), "-"), "[[", 1))

  ## Convert the multiBigwigSummary into a granges
  coverage.df = makeGRangesFromDataFrame(coverage.df, keep.extra.columns=TRUE)

  ## Map coverage results onto group peak tables
  for(group in colnames(mcols(coverage.df))){
    ##
    print(group)
    group.bw = gsub("\\.", "-", group)
    ## Get the coverage results and per-group peak tables in the same order
    overlaps = findOverlaps(peak_table.list[[group.bw]], coverage.df, type="within")
    coverage.df.group = coverage.df[subjectHits(overlaps),]

    ## Quick sanity check
    if(!all(coverage.df.group == peak_table.list[[group.bw]])){ stop("Something has gone wrong with ordering of granges, could not proceed.") }

    ## Extract coverages from grange metadata fields.
    coverage.values = Reduce(cbind, mcols(coverage.df.group)@listData); colnames(coverage.values) = colnames(mcols(coverage.df.group))

    ## Compute coverage proportion w.r.t. within-group accessibility.
    coverage.summary = rowSums(coverage.values[,setdiff(colnames(mcols(coverage.df.group)), group)])
    peak_table.list[[group.bw]]$peak_coverage_rank_up = coverage.values[,group] / coverage.summary
  }
  return(peak_table.list)
}
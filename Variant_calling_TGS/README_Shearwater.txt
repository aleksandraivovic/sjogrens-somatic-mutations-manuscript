Variants in targeted deep sequencing were called using Shearwater algorithm.

“shearwater_pipe_WithPredefinedBAMList.R” is the script for calling variants. It takes as input the tumor and normal bam files, and a bed file with bait capture regions (SJOGRENS_BAITS.bed). 
The tumor samples are all lymphocytes sorted from salivary gland, found “2094_metadata_simple.csv”. The normal samples are all the fibroblast samples. 

“wrapper_shearwaterML_multiplebams.R” is a wrapper for the main Shearwater script when running across many samples. 

“annotate_filter_shearwater.R” takes the input from the variant caller to annotate with statistics and filter the variants that pass QC thresholds.

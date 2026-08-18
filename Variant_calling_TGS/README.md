Variants in targeted deep sequencing were called using Shearwater algorithm, as described in the original publication https://academic.oup.com/bioinformatics/article/30/9/1198/234115

Shearwater is based on the deepSNV R package: https://github.com/im3sanger/deepSNV

The scripts provided were used to run variant calling on targeted deep sequencing from sorted lymphocytes from salivary gland biopsies.

“shearwater_pipe_WithPredefinedBAMList.R” is the script for calling variants. It takes as input the tumor and normal bam files, and a bed file with bait capture regions (provided as SJOGRENS_BAITS.bed for this study). 
The tumor samples are all lymphocytes sorted from salivary gland, found “2094_metadata_simple.csv”. The normal samples are all the fibroblast samples. 

“wrapper_shearwaterML_multiplebams.R” is a wrapper for the main Shearwater script when running across many samples. 

“annotate_filter_shearwater.R” takes the input from the variant caller to annotate with statistics and filter the variants that pass QC thresholds.

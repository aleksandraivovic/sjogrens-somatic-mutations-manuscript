# HDP no priors - extract components
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(hdp)
})


#chlist <- vector("list", 7)
chlist <- list()

for (i in c(1:8)){
  filename <- paste0("output/out_noprior/hdpout_bySample_noprior_10jun24_",i,".rds")
  chlist <- append(chlist, readRDS(filename))
}

mut_example_multi <- hdp_multi_chain(chlist)

# extract signatures
mut_example_multi <- hdp_extract_components(mut_example_multi)

saveRDS(mut_example_multi, "./output/out_noprior/hdpout_bySample_noprior_8chains_object_150524.rds")
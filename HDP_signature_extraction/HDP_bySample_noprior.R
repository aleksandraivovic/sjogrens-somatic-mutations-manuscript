##### HDP for signature extraction on WGS samples from LCM
## Run with min 200 mutations per sample for sig extraction and NO priors

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(reshape2)
  library(RColorBrewer)
  library(stringr)
  library(hdp)
  library(sigfit)
  # library(ggtree)
  library(lsa)
  library(lattice)
})


## load COSMIC signatures
data("cosmic_signatures_v3")

cosmic_signatures_v3[1:4,1:4]

##------------------------------------------------------------
## load counts matrix
matrices = read.delim("HDPinpT_all.tsv",  sep = "\t", header = T, check.names = F) # 459
matrices = t(matrices)

matrices[1:6,1:4]

cosmic_signatures_v3[1:6,1:4]

##------------------------------------------------------------

# filter mut counts matrix to keep only clones with minimum 200 mutations

matrices_min200 = matrices[,colSums(matrices)>=200]
matrices_min200[1:6,1:4] # 389 samples

# make rownames match cosmic - order and naming
colnames(cosmic_signatures_v3)

# change naming
rownames(matrices_min200) = paste0(substr(rownames(matrices_min200), 5, 5), substr(rownames(matrices_min200), 1, 1), substr(rownames(matrices_min200), 7, 7), ">",
                                   substr(rownames(matrices_min200), 5, 5), substr(rownames(matrices_min200), 3, 3), substr(rownames(matrices_min200), 7, 7))

matrices_min200 = as.data.frame(matrices_min200)

matrices_min200$MutationType <- rownames(matrices_min200)

matrices_min200[1:6,1:4]

# change order
matrices_min200$MutationType==colnames(cosmic_signatures_v3)

matrices_min200$MutationType = factor(matrices_min200$MutationType, levels = colnames(cosmic_signatures_v3))

matrices_min200 = matrices_min200[order(matrices_min200$MutationType),]

matrices_min200 = matrices_min200[,c(390, 1:389)]

# check
matrices_min200$MutationType==colnames(cosmic_signatures_v3)
matrices_min200[1:6,1:4]

##------------------------------------------------------------
# load metadata

sample_metadata = read.csv("2107_sample_metadata.csv", header = T)
head(sample_metadata)

sample_metadata %>% select(Diagnosis, Cell_type) %>% table
# Diagnosis  Acinus B cells Duct T cells
# Normal       46       1   66       2
# Sjogrens     96      46  198       1


# order mutcounts matrix to be in the same sample order as metadata
length(sample_metadata$Sample_ID) # 456
length(colnames(matrices_min200)) # 390 (incl. MutationType)

sample_metadata_filtered = sample_metadata[sample_metadata$Sample_ID %in% colnames(matrices_min200),] # dim(sample_metadata_filtered) 389 10

sample_names <- sapply(sample_metadata_filtered$Sample_ID, as.character)
matrices_ordered = matrices_min200[, c("MutationType", sample_names)]
matrices_ordered[1:6,1:4]

colnames(matrices_ordered)==c("MutationType", sample_names)

##-------------------------------------------------------------------------------


rownames(matrices_ordered) = matrices_ordered$MutationType
matrices_ordered = matrices_ordered %>% select(!MutationType)

write.table(matrices_ordered, "matriced_ordered_byPt_min200.tsv", sep = "\t", row.names = T, col.names = T, quote = F)

# transpose counts matrix: samples (clones) in rows
matrices_ordered_t = t(matrices_ordered)

# sort metadata by patient
sample_metadata_filtered <- sample_metadata_filtered %>% arrange(Patient) 

##-------------------------------------------------------------------------------
# prepare hierarchy by setting up dp node counts

# number of prior sigs=0
nps <- 0
Nsample <- nrow(matrices_ordered_t) # 389
parent = 21 # number of parent nodes (does not include grandparent or prior nodes) 1+number of patients

idx = table(sample_metadata_filtered$Patient)[sample_metadata_filtered$Patient %>% unique] %>% cumsum() %>% data.frame()
colnames(idx) = 'start'
idx$end = idx$start
idx <- idx %>% mutate_at(c("start"), tibble::lst("start"=lag), n = 1 )
idx['PD42759','start'] <- 0
idx$start <- idx$start+1
idx$dp_start = idx$start+nps+parent+1 # priors + parents + grandparent dps
idx$dp_end = idx$end+nps+parent+1
for (i in 1:nrow(idx)) {
  idx$length[i] = length(idx$start[i]:idx$end[i])
}

# manually assign node values for sample nodes
idx$cp = c(3:23)
idx$pp = c(2:22)

idx

#         start end dp_start dp_end length cp pp
# PD42759     1  17       23     39     17  3  2
# PD42760    18  32       40     54     15  4  3
# PD42761    33  49       55     71     17  5  4
# PD42763    50  68       72     90     19  6  5
# PD42764    69  74       91     96      6  7  6
# PD42765    75  88       97    110     14  8  7
# PD42766    89 106      111    128     18  9  8
# PD42767   107 114      129    136      8 10  9
# PD42768   115 129      137    151     15 11 10
# PD42769   130 140      152    162     11 12 11
# PD42770   141 171      163    193     31 13 12
# PD42771   172 194      194    216     23 14 13
# PD42772   195 206      217    228     12 15 14
# PD42773   207 233      229    255     27 16 15
# PD42774   234 255      256    277     22 17 16
# PD45527   256 286      278    308     31 18 17
# PD45528   287 301      309    323     15 19 18
# PD45529   302 318      324    340     17 20 19
# PD45530   319 342      341    364     24 21 20
# PD45531   343 347      365    369      5 22 21
# PD45532   348 389      370    411     42 23 22


cp_index = c()
for (i in 1:nrow(idx)) {
  a = rep(idx$cp[i], idx$length[i])
  cp_index = c(cp_index, a)
}


pp_index = c()
for (i in 1:nrow(idx)) {
  a = rep(idx$pp[i], idx$length[i])
  pp_index = c(pp_index, a)
}



#-------------------------------------------------------------------------------------------------------------
# initialize cosmic prior. REQUIRES R 3.6
# set up hdp no prior per Yichen's script

### 1. Run HDP chains
out = "hdp_noprior_"

i <-as.numeric(commandArgs(T)[1])

#for (i in 1:4){

print(paste0("Starting job ",i))


hdp_mut <- hdp_init(ppindex = c(0, rep(1,21), pp_index), # index of all nodes
                    cpindex = c(1, rep(2,21), cp_index), # index of the CP to use
                    hh = rep(1, 96), # prior is uniform over 96 categories
                    alphaa = c(rep(1,max(cp_index))), # shape hyperparams for new CPs
                    alphab = c(rep(1,max(cp_index))))  # rate hyperparameters for 2 CPs

hdp_mut <- hdp_setdata(hdp_mut,  
                       dpindex = (1+nps+parent)+1:Nsample, ## dp indices of sample nodes
                       matrices_ordered_t) # samples should be in rows

hdp_activated <- dp_activate(hdp_mut, 1:numdp(hdp_mut), initcc=10, seed=i*200)

chlist <- hdp_posterior(hdp_activated,
                        burnin=20000,
                        n=100,
                        space=1000,
                        cpiter=3,
                        seed=i*1e3)




# save data.
output_dir = "out_noprior"

saveRDS(chlist, paste0(output_dir, "/hdpout_bySample_noprior_10jun24_",i,".rds"))








###### HDP no prior postprocess - include HDP component 4 as a novel sig

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

setwd("out_noprior/")


#----------------------------------------------------------
# load cosmic refs
data("cosmic_signatures_v3")
ref = cosmic_signatures_v3

ref = t(ref)
features = rownames(ref)

# load HDP components 
hdp_sigs = read.delim("noprior_output_hdp_sigs_topN_8chains_130624.txt", sep = " ", header = T, check.names = F)
hdp_sigs[1:6,1:4]

# check that hdp_sigs and cosmic ref tables are in the same order
rownames(ref)==rownames(hdp_sigs)

# make new ref matrix combining HDP Sig 4 and just the desired cosmic sigs
gdsigs = c("SBS1", "SBS2", "SBS5", "SBS8", "SBS9", "SBS13", "SBS17b", "SBS18", "SBS40")

ref_wComp4 = ref[,gdsigs]

ref_wComp4 = cbind(ref_wComp4, "Sig4" = hdp_sigs[,"Signature4"])

### Refit reference signatures with Sigfit
##-------------------------------------------------------------------------------
library(sigfit)
library(RColorBrewer)

## Get ordered mutation counts for all samples with >100 muts
## load counts matrix

matrices_ordered = read.delim("matriced_ordered_byPt_min100.tsv", sep = "\t", header = T, check.names = F)
matrices_ordered_t = t(matrices_ordered)

matrices_ordered_t[1:6,1:4]
cosmic_signatures_v3[1:6,1:4]
ref_wComp4[1:6,1:4] # transpose?

colnames(matrices_ordered_t)==colnames(cosmic_signatures_v3)

# read in original hdp2ref but manually change Component 4 to match Component 4
hdp2ref=readRDS("./initial_decomposed_sigs_130624.rds") # hdp2ref (or sigs_deconv_R2) is a list of which cosmic refs each hdp sig breaks down into
hdp2ref$Signature4 <- "Sig4"


##-------------------------------------------------------------------------------
#### now run signature filtering and Sigfit 

# samples to run Sigfit on
samples_min100 = read.delim("matriced_ordered_byPt_min100.tsv", header = ) %>% colnames()
# fix samples names to replace "." with "-"
samples_min100 = str_replace(string = samples_min100, pattern = "\\.", replacement = "-")


# hdp_exposures is filtered to keep only HDP output that is present in at least 5% of 5 samples 
hdp_exposures=read.table("./noprior_output_hdp_exposures_8chains_130624.txt", header = T, check.names = F) # hdp_exposures is the output of hdp on samples with >=200 muts - sample vs hdp component contribution
# nrow 389

# fix samples names to replace "." with "-"
rownames(hdp_exposures) = str_replace(string = rownames(hdp_exposures), pattern = "\\.", replacement = "-")

hdp_counts = matrices_ordered_t # counts for samples with >100 muts
samples_min100 = rownames(hdp_counts)
# nrow 448

#rownames(hdp_exposures)==rownames(hdp_counts) #hdp_exposures has samples with >200 muts while counts has >100

colnames(hdp_exposures)=paste0("Signature",colnames(hdp_exposures)) # rename hdp components to "Signature0:9" - all the original hdp exposures incl 0

my_sigs = c(unlist(hdp2ref) %>% unique)
sig_order = 1:length(my_sigs)
names(sig_order) = my_sigs

# subset reference sigs for just sigs of interest - already did this previously - ref_wComp4
# final_sigs = cosmic_signatures_v3[my_sigs,]  ## "SBS1"   "SBS2"   "SBS5"   "SBS8"   "SBS9"   "SBS13"  "SBS17b" "SBS18"  "SBS40"  "Sig4"  
final_sigs = t(ref_wComp4)

colourCount = nrow(final_sigs)
getPalette = colorRampPalette(brewer.pal(10, "Set3"))
all_cols=getPalette(8)
all_cols=c(all_cols,"firebrick","magenta")
final_sigs = final_sigs[names(sort(sig_order[rownames(final_sigs)])),]
names(all_cols)=rownames(final_sigs)
all_cols
patients=unique(substr(rownames(hdp_counts),1,7))

###--Sigfit refitting for each sample--

# filtering hdp output (hdp_exposures) if desired to remove comp 0
hdp_exposures_tmp = hdp_exposures[,2:ncol(hdp_exposures)] 

### This was already done in the postprocessing:
# common_signatures = colSums(hdp_exposures>0.05)>5 # Keep only the signatures that make up at least 5% of at least 5 samples
# n_signatures = sum(common_signatures)
# hdp_exposures_tmp = hdp_exposures[,common_signatures]

exposure_names = colnames(hdp_exposures_tmp)
#exposure_names = paste0("Signature",colnames(hdp_exposures_tmp)) # rename hdp components to "Signature0:9" - all the original hdp exposures incl 0
hdp_names = names(hdp2ref)

# Check to see if the two sets of signature names are the same. We start with exposure_names[2] to ignore "Signature0".
same_names = setequal(exposure_names, hdp_names)
if (!same_names){
  stop('Signature names within `hdp_names` do not match `exposure_names`. Were they filtered differently?')
}

bad_samples <- list() # keep track of the samples which don't have enough signatures that pass the thresholds, below


#fit_samples settings
iter = 20000
warmup = 10000
chains = 5

for (patient in patients[1:length(patients)]){
  # Check to see if the two sets of signature names are the same. We start with exposure_names[2] to ignore "Signature0".
  # Duplicating this from above in case we just re-run the `for` loop
  same_names = setequal(exposure_names, hdp_names)
  if (!same_names){
    stop('Signature names within `hdp_names` do not match `exposure_names`. Were they filtered differently?')
  }
  # hdp_exposures_tmp_pt=hdp_exposures_tmp[rowSums(hdp_counts)>100,] # in Yichen's: the final signatures in each tree make up at least 5% of mutations in at least one branch, with branch length >200 -> changed this to at least 10% of mutations in branches of >100 muts
  hdp_exposures_tmp_pt=hdp_exposures_tmp
  
  # which hdp sigs to be used
  hdp_sigs_final=colnames(hdp_exposures_tmp_pt)[colSums(hdp_exposures_tmp_pt[grepl(paste0(patient,"a_"),rownames(hdp_exposures_tmp_pt)),]>0.05)>0] #0.05 or 0.1
  # hdp_sigs_final is sigs in a given patient that occur in at least one branch >5%
  
  for (sample in grep(patient,samples_min100,value=T)){
    # hdp_sigs_present=colnames(hdp_exposures_tmp)[hdp_exposures_tmp[sample,]>0.05] # sigs present at >5% in a given sample
    hdp_sigs_present = hdp_sigs_final  # just use sigs present in the patient, not in each individual sample
    if (length(hdp_sigs_present) == 0){
      print(paste0('None of the current signatures make up at least 5% of Sample ', sample))
      bad_samples <- append(bad_samples, sample)
      next
    } else if (sum(hdp_sigs_present %in% hdp_sigs_final) == 0){
      print(paste0('None of the signatures making up at least 5% of Sample ', sample, ' are present in the final list of signatures'))
      bad_samples <- append(bad_samples, sample)
      next
    }
    
    hdp_sigs_present=hdp_sigs_present[hdp_sigs_present %in% hdp_sigs_final]
    ref_sigs_present=unique(c(unlist(hdp2ref[hdp_sigs_present])))
    if (length(ref_sigs_present) == 0){
      stop(paste('Signature(s)', hdp_sigs_present, 'not found in hdp2ref. Did you filter hdp2ref differently than hdp_exposures?'))
    }
    
    ref_sigs_present = ref_sigs_present[ref_sigs_present %in% rownames(final_sigs)]
    counts_patient=as.data.frame(hdp_counts)[sample,]
    # colnames(counts_patient)=colnames(cosmic_signatures_v3) # sets colnames to mutation classes
    ref_sigs_present=names(sort(sig_order[ref_sigs_present])) # just getting the name of each signature? already is...
    
    if (length(ref_sigs_present) == 1){
      exposure_matrix <- matrix(1.0)
      colnames(exposure_matrix) <- c(sample)
      rownames(exposure_matrix) <- ref_sigs_present
    } else {
      fit=fit_signatures(counts=counts_patient, 
                         signatures = final_sigs[ref_sigs_present,,drop=FALSE],  # "drop" necessary in case only one signature present; returns dataframe instead of just the row
                         iter = iter, # 20000
                         warmup = warmup, #  10000
                         model="poisson",
                         chains = chains) # 5
      pars <- retrieve_pars(fit,par = "exposures",hpd_prob = 0.95)
      exposure_matrix<-t(pars$mean[,ref_sigs_present])
    }
    ref_sigs_present_final=rownames(exposure_matrix)[exposure_matrix[,]>0.05]
    
    if (sum(exposure_matrix[rownames(exposure_matrix) %in% c('SBS2','SBS13'),])>0.05){
      ref_sigs_present_final=unique(c(ref_sigs_present_final,"SBS2","SBS13"))
    }
    
    # ref_sigs_present_final=names(sort(sig_order[ref_sigs_present_final])) # just getting the name of each signature? already is...
    if (length(ref_sigs_present_final) == 1){
      exposure_matrix <- matrix(1.0)
      colnames(exposure_matrix) <- c(sample)
      rownames(exposure_matrix) <- ref_sigs_present_final
    } else if (length(ref_sigs_present_final) != length(ref_sigs_present)){
      fit=fit_signatures(counts=counts_patient,
                         signatures = final_sigs[ref_sigs_present_final,,drop=FALSE],  # "drop" necessary in case only one signature present; returns dataframe instead of just the row
                         iter = iter,  # iter = 20000,warmup = 10000,model="poisson",chains = 5
                         warmup = warmup,
                         model="poisson",
                         chains = chains) 
      pars <- retrieve_pars(fit,par = "exposures",hpd_prob = 0.95)
      exposure_matrix<-t(pars$mean[,ref_sigs_present_final])
    }
    write.table(exposure_matrix,paste0("./sigfit_output_wComp4_oct2024/",sample,"_contrib.txt"),quote = F,sep='\t')
    #  if (length(ref_sigs_present_final) > 1){
    #    plot_all(mcmc_samples = fit,
    #             out_path = "./sigfit_output", 
    #             prefix = sample)
    #  }
  }
}

n_failed <- length(bad_samples)
if (n_failed > 0){
  message(paste(n_failed, 'sample(s) did not have enough signatures that passed the thresholds above:'))
  for (sample in bad_samples) {
    message(sample)
  }
} else {
  print('All samples successfully passed the provided thresholds.')
}





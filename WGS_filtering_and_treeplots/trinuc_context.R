## Generate trinucleotide context matrix for signature extraction

#setwd("~/Mounts/ai5_network/scratch117/2107/caveman_jul22/")
setwd("/nfs/users/nfs_a/ai5/scratch117/2107/caveman_jul22")

options(stringsAsFactors=FALSE)

library("GenomicRanges")
library("Rsamtools")
library("MASS")
library(dplyr)
library(tidyr)
library(stringr)
#library(lme4)
#library(lmerTest)
#library(sjPlot)
library(reshape2)
library(RColorBrewer)

genomeFile = "~/scratch117/hg19/genome.fa"

n <- 30
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
cols=c("grey80","peachpuff","forestgreen","firebrick","steelblue","pink2",
       "turquoise1",'orange2',"chartreuse","mediumorchid3","grey20",
       "yellow2","mediumaquamarine","tomato4")

metadata = read.table("lcm_all_metadata.tsv", sep = "\t", header = T, stringsAsFactors = F)
#meta = metadata %>% filter(Tissue=="Epithelium")

samples = metadata$Sample_ID # 456

#all_vafs <- data.frame()

#for (sample in samples) {
#  patient=substr(sample, 1, 7)  
#  NV <- read.table(paste0(patient,"/bb_filt/","snp_NV_filtered_all.txt"), sep = " ", header = T,  stringsAsFactors = F)
#  colnames(NV) <- str_replace(colnames(NV), "[.]", "-")
#  NV$mut_id <- rownames(NV)
#  NV <- NV[order(rownames(NV)),]
#  
#  NR <- read.table(paste0(patient,"/bb_filt/","snp_NR_filtered_all.txt"), sep = " ", header = T,  stringsAsFactors = F)
#  colnames(NR) <- str_replace(colnames(NR), "[.]", "-")
#  NR$mut_id <- rownames(NR)
#  NR <- NR[order(rownames(NR)),]
#  
#  VAF <- data.frame("Sample"=sample, "mut_id"=NV$mut_id, "NV" = NV[[sample]], "Depth" = NR[[sample]], "vaf"=NV[[sample]]/NR[[sample]])
#  VAF <- VAF %>% filter(vaf > 0)
#  VAF <- VAF %>% drop_na()
#  all_vafs <- rbind(all_vafs, VAF)
#}

#write.table(all_vafs, "all_samples_vafs.tsv", row.names = F, col.names = T, quote = F, sep = "\t")

all_vafs = read.delim("../all_samples_vafs.tsv", header = T)

#for (sample in samples) {
mutations <- data.frame(separate(all_vafs, col = mut_id, into = c("chr","pos","ref","mut"), sep = "_"))
mutations$pos=as.numeric(mutations$pos)
mutations = mutations[(mutations$ref %in% c("A","C","G","T")) & (mutations$mut %in% c("A","C","G","T")) & mutations$chr %in% c(1:22,"X","Y"),]
mutations$trinuc_ref = as.vector(scanFa(genomeFile, GRanges(mutations$chr, IRanges(as.numeric(mutations$pos)-1, as.numeric(mutations$pos)+1))))    

# 2. Annotating the mutation from the pyrimidine base
ntcomp = c(T="A",G="C",C="G",A="T")
mutations$sub = paste(mutations$ref,mutations$mut,sep=">")
mutations$trinuc_ref_py = mutations$trinuc_ref
for (j in 1:nrow(mutations)) {
  if (mutations$ref[j] %in% c("A","G")) { # Purine base
    mutations$sub[j] = paste(ntcomp[mutations$ref[j]], ntcomp[mutations$mut[j]], sep=">")
    mutations$trinuc_ref_py[j] = paste(ntcomp[rev(strsplit(mutations$trinuc_ref[j],split="")[[1]])],collapse="")
  }
}  

#  # 3. Counting subs
#  freqs = table(paste(mutations$sub,paste(substr(mutations$trinuc_ref_py,1,1),substr(mutations$trinuc_ref_py,3,3),sep="-"#),sep=","))
#  sub_vec = c("C>A","C>G","C>T","T>A","T>C","T>G")
#  ctx_vec = paste(rep(c("A","C","G","T"),each=4),rep(c("A","C","G","T"),times=4),sep="-")
#  full_vec = paste(rep(sub_vec,each=16),rep(ctx_vec,times=6),sep=",")
#  freqs_full = freqs[full_vec];
#  freqs_full[is.na(freqs_full)] = 0;
#  names(freqs_full) = full_vec
#  write.table(freqs_full, paste0(sample, "_hdp_inp.tsv"), sep = "\t", col.names = T, quote = F, row.names = F)
##}

for (sample in samples) {
  muts <- mutations[mutations$Sample==sample,]
  freqs = table(paste(muts$sub,paste(substr(muts$trinuc_ref_py,1,1),substr(muts$trinuc_ref_py,3,3),sep="-"),sep=","))
  sub_vec = c("C>A","C>G","C>T","T>A","T>C","T>G")
  ctx_vec = paste(rep(c("A","C","G","T"),each=4),rep(c("A","C","G","T"),times=4),sep="-")
  full_vec = paste(rep(sub_vec,each=16),rep(ctx_vec,times=6),sep=",")
  freqs_full = freqs[full_vec];
  freqs_full[is.na(freqs_full)] = 0;
  names(freqs_full) = full_vec
  write.table(freqs_full, paste0(sample, "_hdp_inp.tsv"), sep = "\t", col.names = T, quote = F, row.names = F)
}

write.table(mutations, "mut_counts_all.tsv", sep = "\t", col.names = T, quote = F, row.names = F)

### get transposed matrix of counts for HDP input
HDPinp <- as.data.frame(freqs_full)
colnames(HDPinp) <- c("trinuc", "Freq")
HDPinp <- HDPinp %>% dplyr::select(trinuc)

for (sample in c(samples)) {
  hdpinp <- read.table(paste0(sample, "_hdp_inp.tsv"), sep = '\t', header = T)
  HDPinp[,sample] <- hdpinp$Freq
}

rownames(HDPinp) <- HDPinp$trinuc
HDPinpT <- HDPinp %>% dplyr::select(-trinuc) %>% t()

write.table(HDPinpT, "HDPinpT_all.tsv", sep = '\t', row.names = T, col.names = T, quote = F)

###############################

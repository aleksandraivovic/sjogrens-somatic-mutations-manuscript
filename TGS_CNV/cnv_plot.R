#setwd("~/Desktop/output-all/")
## on farm:
setwd("/nfs/users/nfs_a/ai5/scratch126/scratch117_migration/ai5/2094/cnvkit/output-all/")
options(stringsAsFactors=FALSE)
library(dplyr)
library(reshape2)
library(ggplot2)

tumors <- readLines("../tumor_IDs.txt") # all biopsy samples
pbmc <- readLines("../pbmc_IDs.txt") # all pbmc samples
all <- append(tumors, pbmc)

all = all[all != "PD30967m_hum"]

monosomy <- data.frame()

# Get mean depth across X and autosomes for all samples
for (sample in all) {
  cnr <- read.table(paste0(sample, ".sample.dupmarked.cnr"), sep = "\t", header = T)
  #  cnr <- cnr %>% filter(gene != "Antitarget")
  cnr <- cnr %>% filter(chromosome!='Y')
  Xchr <- cnr %>% filter(chromosome == 'X')
  autosome <- cnr %>% filter(chromosome != 'X')
  
  samp <- data.frame("sample" = sample,
                     "Xmean" = mean(Xchr$log2),
                     "Xsd" = sd(Xchr$log2),
                     "autoMean" = mean(autosome$log2),
                     "autoSd" = sd(autosome$log2),
                     "MannWhitney_pval" = wilcox.test(Xchr$log2, autosome$log2, alternative = "less")$p.value,
                     "t_test_pval"=t.test(Xchr$log2, autosome$log2, alternative = "less")$p.value
  )
  
  monosomy <- rbind(monosomy, samp)
}

metadata <- read.table("../../2094_sample_metadata.tsv", header = T, sep = "\t")

# Add cell type and diagnosis meta data
for (i in 1:nrow(monosomy)){
  for (j in 1:nrow(metadata)) {
    if(as.character(monosomy$sample)[i]==as.character(metadata$Sample)[j]){
      monosomy$cell_type[i]=as.character(metadata$CellType[j]);
      monosomy$origin[i]=as.character(metadata$Origin[j]);
      monosomy$focus_group[i]=as.character(metadata$FocusGroup[j]);
      monosomy$diagnosis[i]=as.character(metadata$Diagnosis[j]);
    }
  }
}

# Remove male samples
monosomy$patient <- substr(monosomy$sample, 1, 7) # dim 312
males <- c("PD30975", "PD30978", "PD42054", "PD42069", "PD42072", "PD42081", "PD42082", "PD42086", "PD42087")
monosomy <- monosomy %>% filter(!patient %in% males) %>% filter(!sample=="PD42084j_hum") #dim 261

# FDR adjusted q-values for Mann-Whitney and T-test
monosomy$qvalMW <- p.adjust(p = monosomy$MannWhitney_pval, method = "fdr", n = cnr$gene[!cnr$gene=='Antitarget'] %>% length())
monosomy$qvalTT <- p.adjust(p = monosomy$t_test_pval, method = "fdr", n = cnr$gene[!cnr$gene=='Antitarget'] %>% length())


# Effect size of difference between autosomal and X chr copy ratio
monosomy$effect <- (monosomy$autoMean-monosomy$Xmean)/sqrt(((monosomy$autoSd)^2 + (monosomy$Xsd)^2)/2)

# plot q-value and effect size cut-offs
ggplot(monosomy, aes(x = effect, y = qvalMW)) + geom_rect(mapping = aes(xmin=0.1, xmax=1.25, ymin=-0.02, ymax=0.01, fill='Signif'), alpha=0.5, color='red') + geom_hline(yintercept = 0.01, col = 'red') + geom_vline(xintercept = 0.1, col='red') + geom_point(alpha=0.6) + xlim(-0.6,1.25) + ylim(-0.02,1) + theme_light() + labs(x="Effect size", y="q-value", size=12) + theme(legend.position = 'none', axis.text = element_text(size = 12))

## CN and error bar calculations
monosomy$Xratio <- 2^monosomy$Xmean
monosomy$Xratio_ymin <- 2^(monosomy$Xmean-monosomy$Xsd)
monosomy$Xratio_ymax <- 2^(monosomy$Xmean+monosomy$Xsd)
monosomy$autoRatio <- 2^monosomy$autoMean
monosomy$autoRatio_ymin <- 2^(monosomy$autoMean-monosomy$autoSd)
monosomy$autoRatio_ymax <- 2^(monosomy$autoMean+monosomy$autoSd)

# filter significant samples by q-value and effect size
monosomy$sample <- monosomy$sample %>% factor(ordered = T, levels = monosomy$sample[order(monosomy$Xmean, decreasing = F)])
signif <- monosomy %>% filter(qvalTT < 0.01 & effect > 0.1)

X <- monosomy[c("sample","Xratio","Xratio_ymin","Xratio_ymax")]
colnames(X) <- c("sample","av","ymin","ymax")
X$chr <- "X"
A <- monosomy[c("sample","autoRatio","autoRatio_ymin","autoRatio_ymax")]
colnames(A) <- c("sample","av","ymin","ymax")
A$chr <- "A"

ratios <- rbind(X,A)
ratios$chr <- factor(ratios$chr,levels=c("X","A"))

sigs <- data.frame("sample"=monosomy$sample,"y"=2.75)
sigs$y[monosomy$effect < 0.2 | monosomy$qvalTT > 0.01] <- NA


### plot CN of autosomes vs X chr
#ggplot(ratios, aes(color=chr, y=2*av, x=sample)) +
#  geom_point(position=position_dodge(.5), stat="identity") +
#  geom_errorbar(aes(ymin=2*ymin, ymax=2*ymax), width=0,position=position_dodge(.5)) +
#  geom_text(data=sigs, aes(y = y), label="*", color='black', size = 6) +
#  theme_light() +
#  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  + labs(x="Sample", y= #"Copy number") + scale_fill_discrete(labels = c("X Chromosome","Autosomes")) + theme(legend#.title = element_blank())

## plot CN of just CD8 samples
monosomyCD8 <- monosomy %>% filter(cell_type=="Cytotoxic T Cells")
monosomyCD8$diff <- monosomyCD8$autoMean-monosomyCD8$Xmean
monosomyCD8$sample <- factor(monosomyCD8$sample, levels = monosomyCD8$sample[order(monosomyCD8$diff, decreasing = T)])

X <- monosomyCD8[c("sample","Xratio","Xratio_ymin","Xratio_ymax")]
colnames(X) <- c("sample","av","ymin","ymax")
X$chr <- "X"
A <- monosomyCD8[c("sample","autoRatio","autoRatio_ymin","autoRatio_ymax")]
colnames(A) <- c("sample","av","ymin","ymax")
A$chr <- "A"

ratios <- rbind(X,A)
ratios$chr <- factor(ratios$chr,levels=c("X","A"))

sigs <- data.frame("sample"=monosomyCD8$sample,"y"=2.75)
sigs$y[monosomyCD8$effect < 0.1 | monosomyCD8$qvalTT > 0.01] <- NA

ggplot(ratios, aes(color=chr, y=2*av, x=sample)) +
  geom_point(size=5, position=position_dodge(.5), stat="identity") +
  scale_color_brewer(palette = "Set2", labels = c("X Chromosome","Autosomes")) +
  geom_errorbar(aes(ymin=2*ymin, ymax=2*ymax), width=0.3,position=position_dodge(.5)) +
  geom_text(data=sigs, aes(y = y), label="*", color='black', size = 12) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=16), axis.title = element_text(size = 20)) +
  theme(axis.text.y = element_text(size=20), axis.title = element_text(size = 20)) +
  labs(x="Sample", y= "Copy number") +
  theme(legend.title = element_blank(), legend.text = element_text(size=20), legend.position = "top")

monosomyCD8$clonesize <- (monosomyCD8$autoRatio-monosomyCD8$Xratio)*2 ## add clone size of X loss in CD8

## plot CN of just CD4 samples
remove(X, A, sigs, ratios)
monosomyCD4 <- monosomy %>% filter(cell_type=="Helper T cells")
monosomyCD4$diff <- monosomyCD4$autoMean-monosomyCD4$Xmean
monosomyCD4$sample <- factor(monosomyCD4$sample, levels = monosomyCD4$sample[order(monosomyCD4$diff, decreasing = T)])

X <- monosomyCD4[c("sample","Xratio","Xratio_ymin","Xratio_ymax")]
colnames(X) <- c("sample","av","ymin","ymax")
X$chr <- "X"
A <- monosomyCD4[c("sample","autoRatio","autoRatio_ymin","autoRatio_ymax")]
colnames(A) <- c("sample","av","ymin","ymax")
A$chr <- "A"

ratios <- rbind(X,A)
ratios$chr <- factor(ratios$chr,levels=c("X","A"))

sigs <- data.frame("sample"=monosomyCD4$sample,"y"=2.75)
sigs$y[monosomyCD4$effect < 0.1 | monosomyCD4$qvalTT > 0.01] <- NA

ggplot(ratios, aes(color=chr, y=2*av, x=sample)) +
  geom_point(size=5, position=position_dodge(.5), stat="identity") +
  scale_color_brewer(palette = "Set2", labels = c("X Chromosome","Autosomes")) +
  geom_errorbar(aes(ymin=2*ymin, ymax=2*ymax), width=0.3,position=position_dodge(.5)) +
  geom_text(data=sigs, aes(y = y), label="*", color='black', size = 12) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=16), axis.title = element_text(size = 20)) +
  theme(axis.text.y = element_text(size=20), axis.title = element_text(size = 20)) +
  labs(x="Sample", y= "Copy number") +
  theme(legend.title = element_blank(), legend.text = element_text(size=20), legend.position = "top")


## plot CN of B/plasma/plasmablast samples
remove(X, A, sigs, ratios)
monosomyB <- monosomy %>% filter(cell_type %in% c("B Cells", "Plasmablasts", "Plasma Cells"))
monosomyB$diff <- monosomyB$autoMean-monosomyB$Xmean
monosomyB$sample <- factor(monosomyB$sample, levels = monosomyB$sample[order(monosomyB$diff, decreasing = T)])

X <- monosomyB[c("sample","Xratio","Xratio_ymin","Xratio_ymax")]
colnames(X) <- c("sample","av","ymin","ymax")
X$chr <- "X"
A <- monosomyB[c("sample","autoRatio","autoRatio_ymin","autoRatio_ymax")]
colnames(A) <- c("sample","av","ymin","ymax")
A$chr <- "A"

ratios <- rbind(X,A)
ratios$chr <- factor(ratios$chr,levels=c("X","A"))

sigs <- data.frame("sample"=monosomyB$sample,"y"=2.75)
sigs$y[monosomyB$effect < 0.1 | monosomyB$qvalTT > 0.01] <- NA

ggplot(ratios, aes(color=chr, y=2*av, x=sample)) +
  geom_point(size=5, position=position_dodge(.5), stat="identity") +
  scale_color_brewer(palette = "Set2", labels = c("X Chromosome","Autosomes")) +
  geom_errorbar(aes(ymin=2*ymin, ymax=2*ymax), width=0.3,position=position_dodge(.5)) +
  geom_text(data=sigs, aes(y = y), label="*", color='black', size = 12) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=16), axis.title = element_text(size = 20)) +
  theme(axis.text.y = element_text(size=20), axis.title = element_text(size = 20)) +
  labs(x="Sample", y= "Copy number") +
  theme(legend.title = element_blank(), legend.text = element_text(size=20), legend.position = "top")


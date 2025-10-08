### LME model from adjusted mutation burden and Sigit mutational signatures 
# adapted Sept 2025

#----------------------------------------------------------

library("GenomicRanges")
library("Rsamtools")
library("MASS")
library(dplyr)
library(tidyr)
library(stringr)
library(lme4)
library(lmerTest)
library(sjPlot)
library(reshape2)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)

### Format Sigfit output per sample into one matrix for LME and other analysis ###

# list output files from Sigfit with Comp4

setwd("~/Desktop/Code_SjD_paper/Mutation_burden_plots/sigfit_output_wComp4_oct2024_old/")
#length(list.files())

# pull original mutcounts file for sample IDs
matrices_ordered = read.delim("~/Desktop/Code_SjD_paper/Mutation_burden_plots/matriced_ordered_byPt_min100.tsv", sep = "\t", header = T, check.names = F)
hdp_counts = t(matrices_ordered) # should be 448 samples with >100 mutations

# samples with >=100 mutations
samples = rownames(hdp_counts)
patients=unique(substr(samples,1,7))

# all reference signatures used for Sigfit 
my_sigs = c("SBS1", "SBS2", "SBS5", "SBS8", "SBS9", "SBS13", "SBS17b", "SBS18", "SBS40", "Sig4")

# merge Sigfit contributions into one data frame
output_all = data.frame()

filenames = list.files(pattern='*_contrib.txt')

for (filename in filenames) {
  sample = gsub("_contrib.txt","",filename)
  output_single = read.delim(filename, header = T)
  output_single = t(output_single) %>% as.data.frame()
  output_all = bind_rows(output_all, output_single)
}

output_all <- output_all %>% dplyr::select(all_of(my_sigs))
output_all[is.na(output_all)] <- 0

# fix names
rownames(output_all) = str_replace(string = rownames(output_all), pattern = "\\.", replacement = "-")

output_all$Sample_ID = rownames(output_all)

#----------------------------------------------------------
## Load metadata file with VAF, Depth, and number of mutations
sample_metadata = read.delim("~/Desktop/Code_SjD_paper/Mutation_burden_plots/lcm_all_metadata_wVAF_wSensitivity.tsv", header = T) # nrow 456

# fix names
sample_metadata$Sample_ID = str_replace(string = sample_metadata$Sample_ID, pattern = "\\.", replacement = "-")

# filter metadata to only include samples with >100 mutations as analyzed by Sigfit
sample_metadata_filt = sample_metadata %>% filter(Sample_ID %in% rownames(output_all)) # nrow

# merge metadata with signatures
sample_metadata_filt = full_join(sample_metadata_filt, output_all, by = "Sample_ID")



#----------------------------------------------------------
#######
### Prepare LME model analysis


# mutation burden adjusted by sensitivity
sample_metadata_filt$adj_MutBurden <- as.integer(sample_metadata_filt$Mutations/sample_metadata_filt$sensitivity)

# rename Comp4 (Sig4) as SBS_A
colnames(sample_metadata_filt)[25] <- "SBS_A"

# add number of muts per signature
my_sigs = c("SBS1", "SBS2", "SBS5", "SBS8", "SBS9", "SBS13", "SBS17b", "SBS18", "SBS40", "SBS_A")

sample_metadata_filt$Muts_SBS1 = sample_metadata_filt$SBS1*sample_metadata_filt$adj_MutBurden

for (sig in my_sigs) {
  sample_metadata_filt[paste0("Muts_", sig)] <- sample_metadata_filt[sig]*sample_metadata_filt$adj_MutBurden
}

# aggregate certain mutation types
sample_metadata_filt$nonAPOBEC_muts <- rowSums(sample_metadata_filt[,c("Muts_SBS1", "Muts_SBS5", "Muts_SBS8", "Muts_SBS9", "Muts_SBS17b", "Muts_SBS18", "Muts_SBS40", "Muts_SBS_A")])
sample_metadata_filt$APOBEC_muts <- rowSums(sample_metadata_filt[,c("Muts_SBS2", "Muts_SBS13")])
sample_metadata_filt$Clocklike <- rowSums(sample_metadata_filt[,c("Muts_SBS1", "Muts_SBS5", "Muts_SBS40")])

#write.table(sample_metadata_filt, "2107_all_meta_wSigs_Oct2024.tsv", sep = "\t", col.names = T, quote = F, row.names = F)

###
# filtering samples to use for LME by VAF*depth or by sensitivity
sample_metadata_filt %>% filter(VAFxDepth>=3) %>% dim #364 out of 448
sample_metadata_filt %>% filter(VAFxDepth>=3.5) %>% dim #303 out of 448
sample_metadata_filt %>% filter(sensitivity>0.4) %>% dim #335 out of 448

# filter to analyze only glandular tissue or only lymphocytes
sample_metadata_filt_epith = sample_metadata_filt %>% filter(Tissue=="Epithelium")
sample_metadata_filt_lympth = sample_metadata_filt %>% filter(Tissue=="Lymphocyte" & Cell_type != "T cells") # exclude the 1 T cell sample

sample_metadata_filt_epith = sample_metadata_filt_epith %>% filter(VAFxDepth>=3.5) #293
#sample_metadata_filt_epith = sample_metadata_filt_epith %>% filter(sensitivity>0.4) #335


#----------------------------------------------------------
### Run LME analysis

# total adjusted mutation burden LME
lmer_model <- lmer(adj_MutBurden ~ Age + (Age - 1|Patient) + Diagnosis + Cell_type, data = sample_metadata_filt_epith, REML = F)
lme4:::drop1.merMod(lmer_model, test = "Chisq")
## model with updated diagnosis:
# Single term deletions
# 
# Model:
#   adj_MutBurden ~ Age + (Age - 1 | Patient) + Diagnosis + Cell_type
# npar    AIC     LRT   Pr(Chi)    
# <none>         4336.6                      
# Age          1 4363.1 28.4893 9.422e-08 ***
#   Diagnosis    1 4334.6  0.0174  0.895063    
# Cell_type    1 4343.0  8.4288  0.003693 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##########
#OLD Model:
##  adj_MutBurden ~ Age + (Age - 1 | Patient) + Diagnosis + Cell_type
##npar    AIC    LRT   Pr(Chi)    
##<none>         4333.8                     
##  Age          1 4365.1 33.295 7.918e-09 ***
##  Diagnosis    1 4334.6  2.799  0.094331 .  
##  Cell_type    1 4341.2  9.337  0.002246 ** 

print(summary(lmer_model))

# Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: adj_MutBurden ~ Age + (Age - 1 | Patient) + Diagnosis + Cell_type
#    Data: sample_metadata_filt_epith
# 
#      AIC      BIC   logLik deviance df.resid 
#   4336.6   4358.7  -2162.3   4324.6      287 
# 
# Scaled residuals: 
#     Min      1Q  Median      3Q     Max 
# -1.9823 -0.5413 -0.1328  0.3311  6.3256 
# 
# Random effects:
#  Groups   Name Variance  Std.Dev.
#  Patient  Age  5.853e+00   2.419 
#  Residual      1.398e+05 373.953 
# Number of obs: 293, groups:  Patient, 21
# 
# Fixed effects:
#                   Estimate Std. Error       df t value Pr(>|t|)    
# (Intercept)         -3.597    159.150   39.030  -0.023  0.98208    
# Age                 17.918      2.519   37.034   7.113    2e-08 ***
# DiagnosisSjogrens  -10.252     77.523   20.906  -0.132  0.89606    
# Cell_typeDuct     -140.333     47.653  289.873  -2.945  0.00349 ** 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#             (Intr) Age    DgnssS
# Age         -0.915              
# DgnssSjgrns -0.464  0.225       
# Cell_typDct -0.244  0.073 -0.055

##############
### OLD
##  Age                 19.086      2.361   35.295   8.084 1.52e-09 ***
##  DiagnosisSjogrens  136.740     75.811   17.793   1.804  0.08824 .  
##  Cell_typeDuct     -148.747     47.486  291.751  -3.132  0.00191 ** 

##---------------------------------------------
####
## compare with previous analysis on this dataset - it's the same
# mutburden_epi_old = read.delim("~/Desktop/Somatic_paper/LCM_WGS/caveman_jul22/all_info_epi_wSigs.tsv", header = T)
# mutburden_epi_old_filt = mutburden_epi_old %>% filter(VAFxDepth>=3.5) # 293 samples
# 
# lmer_model <- lmer(adj_MutBurden ~ Age + (Age - 1|Patient) + Diagnosis + Cell_type, data = mutburden_epi_old_filt, REML = F)
# lme4:::drop1.merMod(lmer_model, test = "Chisq")
# print(summary(lmer_model))


##---------------------------------------------
set.seed(0)
n <- 30
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
cols=c("grey80","peachpuff","forestgreen","firebrick","steelblue","pink2",
       "turquoise1",'orange2',"chartreuse","mediumorchid3","grey20",
       "yellow2","mediumaquamarine","tomato4")


plot_model(lmer_model, type = "pred", terms = "Age", title = "Salivary epithelium: Mutation Burden", show.data = T)

pdf(file = "~/Desktop/Code_SjD_paper/Mutation_burden_plots/Mut_burden_Age_2107_Epith.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 5.5)

plot_model(lmer_model, type = "pred", terms = "Age", show.data = F, colors ='black') + 
  geom_point(data = sample_metadata_filt_epith, aes(x = Age, y = adj_MutBurden, colour = Diagnosis), alpha=0.75, size=2.5) + 
  scale_color_manual(values = c("#377EB8", "#E41A1C"), labels = c("Biopsy-negative control", "Sjögren disease")) + 
  labs(y = "Adjusted Number of Substitutions", colour = "Diagnosis", title = "Salivary Epithelium: Mutation Burden") + theme_light() # + stat_cor() 

dev.off() 

#ggsave("Epith_wgs_age_burden_scatter.pdf", width = 10, height = 7)

pdf(file = "~/Desktop/Code_SjD_paper/Mutation_burden_plots/Mut_burden_Age_2107_Epith_byPt.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 5.5)

plot_model(lmer_model, type = "pred", terms = "Age", show.data = F) + 
  geom_point(data = sample_metadata_filt_epith, aes(x = Age, y = adj_MutBurden, colour = Patient), alpha=0.75, size=2.5) + 
  scale_color_manual(values = col_vector) + 
  labs(y = "Adjusted Number of Substitutions", colour = "Patient", title = "Salivary Epithelium: Mutation Burden") + theme_light()

dev.off()


#ggsave("Epith_wgs_age_burden_scatter_bypatient.png", width = 10.5, height = 7)

#---------------------------------------------------------------------------
#### adjusted mutation burden of only APOBEC muts
lmer_model_3 <- lmer(APOBEC_muts ~ Age + (Age - 1|Patient) + Diagnosis + Cell_type, data = sample_metadata_filt_epith, REML = F)

lme4:::drop1.merMod(lmer_model_3, test = "Chisq")

lmer_model_3 <- lmer(APOBEC_muts ~ Age + (Age - 1|Patient) + Cell_type, data = sample_metadata_filt_epith, REML = F)
print(summary(lmer_model_3))

# Fixed effects:
#   Estimate Std. Error       df t value Pr(>|t|)    
# (Intercept)   -134.789     90.343   78.913  -1.492 0.139696    
# Age              2.147      1.525   45.151   1.408 0.165952    
# Cell_typeDuct  118.199     33.173  292.858   3.563 0.000428 ***

# APOBEC-only mutations by cell type
pdf(file = "~/Desktop/Code_SjD_paper/Mutation_burden_plots/Apobec_Mut_burden_Age_2107_Epith_byPt.pdf",   # The directory you want to save the file in
    width = 4, # The width of the plot in inches
    height = 5)

ggplot(sample_metadata_filt_epith, aes(x = Cell_type, y = APOBEC_muts)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(aes(color=Cell_type), cex=2, position = position_jitter(seed=1)) + 
  geom_point(color='black', cex=2, position = position_jitter(seed=1), shape=1) + 
  scale_color_manual(values = col_vector[5:6]) + labs(x = "cell type", y = "Number APOBEC muts") + 
  theme_light() + theme(legend.position = "bottom", axis.text = element_text(size=12)) + 
  stat_compare_means(method = "wilcox.test", label.x = 1.2, label = "p.format")

dev.off()

# APOBEC fraction by cell type
pdf(file = "~/Desktop/Code_SjD_paper/Mutation_burden_plots/Apobec_fraction_2107_Epith_byPt.pdf",   # The directory you want to save the file in
    width = 4, # The width of the plot in inches
    height = 5)

ggplot(sample_metadata_filt_epith, aes(x = Cell_type, y = SBS13+SBS2)) + 
  # geom_boxplot(outlier.shape = NA) + 
  geom_point(aes(color=Cell_type), cex=2, position = position_jitter(seed=1)) + 
  geom_point(color='black', cex=2, position = position_jitter(seed=1), shape=1) + 
  scale_color_manual(values = col_vector[5:6]) + labs(x = "cell type", y = "apobec sig fraction") + 
  theme_light() + theme(legend.position = "bottom", axis.text = element_text(size=12)) + 
  stat_compare_means(method = "wilcox.test", label.x = 1.2, label = "p.format")

dev.off()

#---------------------------------------------------------------------------
melt_sigs <- melt(sample_metadata_filt_epith, id.vars = c("Sample_ID", "Diagnosis", "Cell_type"), measure.vars = c("SBS5", "SBS1", "SBS40", "SBS9", "SBS8", "SBS18", "SBS17b", "SBS_A", "SBS2", "SBS13"), variable.name = "Signature", value.name = "Exposure")

#write.table(melt_sigs, "epithel_info_melted.tsv", sep = "\t", col.names = T, row.names = F, quote = F)

melt_sigs$Signature <- factor(melt_sigs$Signature, levels = c("SBS5", "SBS1", "SBS40", "SBS9", "SBS8", "SBS18", "SBS17b", "SBS_A", "SBS2", "SBS13"))
melt_sigs <- melt_sigs[order(melt_sigs$Signature), ]

# order by increasing APOBEC contribution
sig13 = melt_sigs %>% filter(Signature=="SBS13")
melt_sigs$Sample_ID <- factor(melt_sigs$Sample_ID, levels = sig13[order(sig13$Exposure),]$Sample_ID)

# barplot faceted by diagnosis and cell type
pdf(file = "~/Desktop/Code_SjD_paper/Mutation_burden_plots/Sigs_barplot_ordered_byApobec.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 6)

ggplot(melt_sigs, aes(x = Sample_ID, y = Exposure, fill = Signature)) + 
  geom_bar(stat="identity", position = position_fill(reverse=T)) + 
  theme(axis.text.x=element_blank()) +
  scale_fill_manual(values = col_vector, labels = c("SBS5", "SBS1", "SBS40", "SBS9", "SBS8", "SBS18", "SBS17b", "SBS_A", "SBS2", "SBS13")) + 
  facet_wrap(Diagnosis ~ Cell_type, scales = "free_x")

dev.off()

#----------------------------------------------------------
# WGS Epi depth

meta2 = sample_metadata_filt_epith
meta2$Sample_ID = factor(meta2$Sample_ID, levels = meta2[order(meta2$Depth),]$Sample_ID)

# coverage
ggplot(meta2, aes(x = reorder(Sample_ID, Depth), y = Depth, color = Cell_type)) + geom_point() + theme_light() + labs(x = "Samples", y = "Depth") + 
  theme_minimal() +
  theme(axis.title = element_text(size = 14), axis.text.y = element_text(size = 12), axis.text.x = element_blank()) + scale_color_manual(values = col_vector[5:6])

#ggsave("Epith_WGS_Depth_colbyCelltype.png", width = 8, height = 7)

# median VAF
ggplot(meta2, aes(x = reorder(Sample_ID, median_VAF), y = median_VAF, color = Cell_type)) + geom_point() + theme_light() + labs(x = "Samples", y = "median_VAF") + theme(axis.title = element_text(size = 14), axis.text.y = element_text(size = 12), axis.text.x = element_blank()) + theme_minimal() + scale_color_manual(values = col_vector[5:6])


# boxplot coverage vs cell type
ggplot(meta2, aes(x = Cell_type, y = Depth)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(aes(color=Cell_type), cex=2, position = position_jitter(seed=1))  +
  theme_light() + scale_color_manual(values = col_vector[5:6])


# boxplot VAF vs cell type
ggplot(meta2, aes(x = Cell_type, y = median_VAF)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(aes(color=Cell_type), cex=2, position = position_jitter(seed=1))  +
  scale_color_manual(values = col_vector[5:6]) +
  theme_light()




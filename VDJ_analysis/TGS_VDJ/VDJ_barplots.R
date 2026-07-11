library(ggplot2)
library(dplyr)
library(paletteer)

## Mixcr for 2094 TGS

setwd("2094_Mixcr_output/")

metadata = read.csv("2094_metadata_simple.csv", header = T, check.names = F)
metadata$patient = substr(metadata$`Sanger ID`, 1, 7)

patients = metadata$patient %>% unique()
samples = metadata$`Sanger ID`

groups = c("B cells", "Plasmablasts", "Plasma cells", "Helper T cells", "Cytotoxic T cells", "bulk B cells", "bulk T cells")

# PD30967b = read.delim("all_IGH/PD30967b_hum_IGH_clones.txt", check.names = F, header = T)
# PD30967b = PD30967b %>% select("cloneId", "cloneCount", "cloneFraction", "allVHitsWithScore", "allDHitsWithScore", "allJHitsWithScore", "allCHitsWithScore", "aaSeqCDR3")
# 
# PD30967b$sample = "PD30967b"

# ggplot(PD30967b, aes(x=sample, y=cloneCount, color = cloneId)) + geom_bar( position = "stack", stat="identity", color = cols_blue[1:nrow(PD30967b)]) + 
 # coord_flip() + scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, 1)) 


###############
## B CELL SAMPLES
bcell_samples = metadata$`Sanger ID`[metadata$`Cell type`=="B cells"]
bcell_samples_short = list()

# remove samples with empty files
for (sample in  bcell_samples) {
  if (file.info(paste0("all_IGH/", sample, "_hum_IGH_clones.txt"))$size > 417) {
    bcell_samples_short <- c(bcell_samples_short,sample)
  }
}

all_igh = data.frame()

for (sample in bcell_samples_short) {
  sample_igh = read.delim(paste0("all_IGH/", sample, "_hum_IGH_clones.txt"), check.names = F, header = T)
  sample_igh = sample_igh %>% select("cloneId", "cloneCount", "cloneFraction", "allVHitsWithScore", "allDHitsWithScore", "allJHitsWithScore", "allCHitsWithScore", "aaSeqCDR3")
  sample_igh$sample = sample
  sample_igh = sample_igh %>% head(20)
  all_igh = rbind(all_igh, sample_igh)
}



samples_ordered_by_max_cloneCount <-
  all_igh %>% group_by(sample) %>% slice(which.max(cloneCount)) %>% arrange(cloneCount) %>% pull(sample)

all_igh_by_cloneCount <-
  all_igh %>% arrange(sample) %>% mutate(
    sample = factor(sample, levels =
                      samples_ordered_by_max_cloneCount),
    cloneCount = factor(cloneCount)
  )

pdf("IGH_all_barplot.pdf", width = 10, height = 10)

ggplot(all_igh_by_cloneCount,
       aes(x = sample, y = cloneCount, fill = cloneCount)) + geom_bar(position = "stack",
                                                                      stat = "identity",
                                                                      color = 'grey')  + scale_fill_brewer(palette = "Blues", name = "Clone size \n(# transcripts)") + 
  theme_light() +
  scale_y_discrete(name = "Clones", breaks = seq(1, 20, by = 1)) +
  ggtitle("IGH clones for all B cell samples") +
  coord_flip()

dev.off()


###############
## PLASMABLAST SAMPLES
plasmablasts_samples = metadata$`Sanger ID`[metadata$`Cell type`=="Plasmablasts"]
plasmablasts_samples_short = list()

# remove samples with empty files
for (sample in  plasmablasts_samples) {
  if (file.info(paste0("all_IGH/", sample, "_hum_IGH_clones.txt"))$size > 417) {
    plasmablasts_samples_short <- c(plasmablasts_samples_short,sample)
  }
}

all_igh = data.frame()

for (sample in plasmablasts_samples_short) {
  sample_igh = read.delim(paste0("all_IGH/", sample, "_hum_IGH_clones.txt"), check.names = F, header = T)
  sample_igh = sample_igh %>% select("cloneId", "cloneCount", "cloneFraction", "allVHitsWithScore", "allDHitsWithScore", "allJHitsWithScore", "allCHitsWithScore", "aaSeqCDR3")
  sample_igh$sample = sample
  sample_igh = sample_igh %>% head(20)
  all_igh = rbind(all_igh, sample_igh)
}



samples_ordered_by_max_cloneCount <-
  all_igh %>% group_by(sample) %>% slice(which.max(cloneCount)) %>% arrange(cloneCount) %>% pull(sample)

all_igh_by_cloneCount <-
  all_igh %>% arrange(sample) %>% mutate(
    sample = factor(sample, levels =
                      samples_ordered_by_max_cloneCount),
    cloneCount = factor(cloneCount)
  )

pdf("plasmablast_IGH_all_barplot.pdf", width = 10, height = 10)

ggplot(all_igh_by_cloneCount,
       aes(x = sample, y = cloneCount, fill = cloneCount)) + geom_bar(position = "stack",
                                                                      stat = "identity",
                                                                      color = 'grey')  + scale_fill_brewer(palette = "Blues", name = "Clone size \n(# transcripts)") + 
  theme_light() +
  scale_y_discrete(name = "Clones", breaks = seq(1, 20, by = 1)) +
  ggtitle("IGH clones for all plasmablast samples") +
  coord_flip()

dev.off()


###############
## PLASMA CELL SAMPLES
plasmacell_samples = metadata$`Sanger ID`[metadata$`CellType_simple`=="Plasma Cells"]
plasmacell_samples_short = list()

# remove samples with empty files
for (sample in  plasmacell_samples) {
  if (file.info(paste0("all_IGH/", sample, "_hum_IGH_clones.txt"))$size > 417) {
    plasmacell_samples_short <- c(plasmacell_samples_short,sample)
  }
}

all_igh = data.frame()

for (sample in plasmacell_samples_short) {
  sample_igh = read.delim(paste0("all_IGH/", sample, "_hum_IGH_clones.txt"), check.names = F, header = T)
  sample_igh = sample_igh %>% select("cloneId", "cloneCount", "cloneFraction", "allVHitsWithScore", "allDHitsWithScore", "allJHitsWithScore", "allCHitsWithScore", "aaSeqCDR3")
  sample_igh$sample = sample
  sample_igh = sample_igh %>% head(20)
  all_igh = rbind(all_igh, sample_igh)
}



samples_ordered_by_max_cloneCount <-
  all_igh %>% group_by(sample) %>% slice(which.max(cloneCount)) %>% arrange(cloneCount) %>% pull(sample)

all_igh_by_cloneCount <-
  all_igh %>% arrange(sample) %>% mutate(
    sample = factor(sample, levels =
                      samples_ordered_by_max_cloneCount),
    cloneCount = factor(cloneCount)
  )

pdf("plasma_IGH_all_barplot.pdf", width = 10, height = 10)

ggplot(all_igh_by_cloneCount, aes(x = sample, y = cloneCount, fill = cloneCount)) + 
  geom_bar(position = "stack", stat = "identity", color = 'grey')  + 
  # scale_fill_brewer(palette = "Blues", name = "Clone size \n(# transcripts)") + 
  scale_fill_manual(values = rev(paletteer_c("grDevices::Blues 3", 5)), name = "Clone size \n(# transcripts)") +
  theme_light() +
  scale_y_discrete(name = "Clones", breaks = seq(1, 20, by = 1)) +
  ggtitle("IGH clones for all plasma cell samples") +
  coord_flip()

dev.off()


###############
## CD4 T CELL SAMPLES
CD4_samples = metadata$`Sanger ID`[metadata$`CellType_simple`=="Helper T cells"]
CD4_samples_short = list()

# remove samples with empty files
for (sample in  CD4_samples) {
  if (file.info(paste0("all_TRB/", sample, "_hum_TRB_clones.txt"))$size > 417) {
    CD4_samples_short <- c(CD4_samples_short,sample)
  }
}

all_trb = data.frame()

for (sample in CD4_samples_short) {
  sample_trb = read.delim(paste0("all_TRB/", sample, "_hum_TRB_clones.txt"), check.names = F, header = T)
  sample_trb = sample_trb %>% select("cloneId", "cloneCount", "cloneFraction", "allVHitsWithScore", "allDHitsWithScore", "allJHitsWithScore", "allCHitsWithScore", "aaSeqCDR3")
  sample_trb$sample = sample
  sample_trb = sample_trb %>% head(20)
  all_trb = rbind(all_trb, sample_trb)
}



samples_ordered_by_max_cloneCount <-
  all_trb %>% group_by(sample) %>% slice(which.max(cloneCount)) %>% arrange(cloneCount) %>% pull(sample)

all_trb_by_cloneCount <-
  all_trb %>% arrange(sample) %>% mutate(
    sample = factor(sample, levels =
                      samples_ordered_by_max_cloneCount),
    cloneCount = factor(cloneCount)
  )

pdf("CD4_TRB_all_barplot.pdf", width = 10, height = 10)

ggplot(all_trb_by_cloneCount,
       aes(x = sample, y = cloneCount, fill = cloneCount)) + geom_bar(position = "stack",
                                                                      stat = "identity",
                                                                      color = 'grey')  + scale_fill_manual(values = rev(paletteer_c("grDevices::Greens 3",))  , name = "Clone size \n(# transcripts)")  + 
  theme_light() +
  scale_y_discrete(name = "Clones", breaks = seq(1, 20, by = 1)) +
  ggtitle("TRB clones for all CD4 T cell samples") +
  coord_flip()

dev.off()


###############
## CD8 T CELL SAMPLES
CD8_samples = metadata$`Sanger ID`[metadata$`CellType_simple`=="Cytotoxic T cells"]
CD8_samples_short = list()

# remove samples with empty files
for (sample in  CD8_samples) {
  if (file.info(paste0("all_TRB/", sample, "_hum_TRB_clones.txt"))$size > 417) {
    CD8_samples_short <- c(CD8_samples_short,sample)
  }
}

all_trb = data.frame()

for (sample in CD8_samples_short) {
  sample_trb = read.delim(paste0("all_TRB/", sample, "_hum_TRB_clones.txt"), check.names = F, header = T)
  sample_trb = sample_trb %>% select("cloneId", "cloneCount", "cloneFraction", "allVHitsWithScore", "allDHitsWithScore", "allJHitsWithScore", "allCHitsWithScore", "aaSeqCDR3")
  sample_trb$sample = sample
  sample_trb = sample_trb %>% head(20)
  all_trb = rbind(all_trb, sample_trb)
}



samples_ordered_by_max_cloneCount <-
  all_trb %>% group_by(sample) %>% slice(which.max(cloneCount)) %>% arrange(cloneCount) %>% pull(sample)

all_trb_by_cloneCount <-
  all_trb %>% arrange(sample) %>% mutate(
    sample = factor(sample, levels =
                      samples_ordered_by_max_cloneCount),
    cloneCount = factor(cloneCount)
  )

pdf("CD8_TRB_all_barplot.pdf", width = 10, height = 10)

ggplot(all_trb_by_cloneCount,
       aes(x = sample, y = cloneCount, fill = cloneCount)) + geom_bar(position = "stack",
                                                                      stat = "identity",
                                                                      color = 'grey')  + scale_fill_manual(values = rev(paletteer_c("grDevices::Greens 3", 12)), name = "Clone size \n(# transcripts)") + 
  theme_light() +
  scale_y_discrete(name = "Clones", breaks = seq(1, 20, by = 1)) +
  ggtitle("TRB clones for all CD8 T cell samples") +
  coord_flip()

dev.off()








### Immunarch analysis for diversity comparison - June 2025library(immunarch)
## Basic analysis: https://immunarch.com

library(immunarch)
library(stringr)
library(rlist)
library(gridExtra)
options(stringsAsFactors = F)

setwd("~/Desktop/Somatic_paper/Figures_pdf/TGS_bulk_lymphocytes/VDJ_calls/2094_Mixcr_output/")


######## ####  IGH analysis ####  #### 

all_clns_igh <- repLoad("all_IGH/")

# filtering all_clns_igh object to exclude samples with less than 5 total clones
new_clns <- list()
for (name in names(all_clns_igh$data)) {
  d <- all_clns_igh$data[[name]]
  if (0 %in% dim(d)) {
    next
  }
  if ((d %>% select(Clones) %>% sum) >= 5) {
    new_clns$data[[name]] <- d
  }
}

new_clns$meta <- all_clns_igh$meta %>% filter(Sample %in% names(new_clns$data))
all_clns_igh = new_clns

# some changes to metadata
all_clns_igh$meta$Diagnosis = str_replace(all_clns_igh$meta$Diagnosis, "Amb", "PSS")
all_clns_igh$meta$Diagnosis = str_replace(all_clns_igh$meta$Diagnosis, "Yes", "PSS")
all_clns_igh$meta$Diagnosis = str_replace(all_clns_igh$meta$Diagnosis, "No", "nonPSS")
all_clns_igh$meta$Diagnosis = str_replace(all_clns_igh$meta$Diagnosis, "N", "nonPSS")


#######
## subset dataset to exclude PBMC samples
biopsy_samples = all_clns_igh$meta$Sample[all_clns_igh$meta$Origin=="biopsy" & all_clns_igh$meta$CellType != "Bulk_B"]

all_clns_igh_biopsy = all_clns_igh$data[biopsy_samples]
all_clns_igh_biopsy_meta = all_clns_igh$meta %>% filter(Origin=="biopsy" & CellType != "Bulk_B")

# Number of clones
exp_vol <- repExplore(all_clns_igh_biopsy, .method = "clones", .coding = T )
p1 = vis(exp_vol, .by = c("CellType"), .meta = all_clns_igh_biopsy_meta) + scale_fill_manual(values = paletteer_d("colorBlindness::Blue2Gray8Steps") ) # #0099CCFF #66E5FFFF #99FFFFFF
p2 = vis(exp_vol, .by = c("Diagnosis"), .meta = all_clns_igh_biopsy_meta) + scale_fill_manual(values = c("#377EB8", "#E41A1C"))
  

# Percentrage of repertoire occupied by clones of varying size
repClonality(all_clns_igh_biopsy, .method = "homeo", .clone.types = c(Small = .05, Medium = .15, Large = .25, Hyperexpanded = .5)) %>% vis()

# Inverse Simpson diversity index
div_inv = repDiversity(all_clns_igh_biopsy, .method = "inv.simp")
p3 = vis(div_inv, .by = c("CellType"), .meta = all_clns_igh_biopsy_meta) + scale_fill_manual(values = paletteer_d("colorBlindness::Blue2Gray8Steps") ) # #0099CCFF #66E5FFFF #99FFFFFF
p4 = vis(div_inv, .by = c("Diagnosis"), .meta = all_clns_igh_biopsy_meta) + scale_fill_manual(values = c("#377EB8", "#E41A1C"))

pdf("~/Desktop/Somatic_paper/Figures_pdf/TGS_bulk_lymphocytes/VDJ_calls/BCR_IGH_byCellType.pdf", height = 11, width = 6)
grid.arrange(p1, p3, ncol = 1, nrow = 2) # Arrange plots side-by-side
dev.off()

pdf("~/Desktop/Somatic_paper/Figures_pdf/TGS_bulk_lymphocytes/VDJ_calls/BCR_IGH_byDiagnosis.pdf", height = 11, width = 5)
grid.arrange(p2, p4, ncol = 1, nrow = 2) # Arrange plots side-by-side
dev.off()


################################################################################

######## ####  TRB analysis ####  #### 

all_clns_trb <- repLoad("all_TRB/")
                        
# filtering all_clns object to exclude samples with less than 5 clones
 new_clns <- list()
 for (name in names(all_clns$data)) {
   d <- all_clns$data[[name]]
   if (0 %in% dim(d)) {
     next
   }
   if ((d %>% select(Clones) %>% sum) >= 5) {
     new_clns$data[[name]] <- d
   }
 }

new_clns$meta <- all_clns_trb$meta %>% filter(Sample %in% names(new_clns$data))
all_clns_trb = new_clns
# some changes to metadata
all_clns_trb$meta$Diagnosis = str_replace(all_clns_trb$meta$Diagnosis, "Amb", "PSS")
all_clns_trb$meta$Diagnosis = str_replace(all_clns_trb$meta$Diagnosis, "Yes", "PSS")
all_clns_trb$meta$Diagnosis = str_replace(all_clns_trb$meta$Diagnosis, "No", "nonPSS")
all_clns_trb$meta$Diagnosis = str_replace(all_clns_trb$meta$Diagnosis, "N", "nonPSS")


## subset dataset to exclude PBMC samples
biopsy_samples = all_clns_trb$meta$Sample[all_clns_trb$meta$Origin=="biopsy"]

all_clns_trb_biopsy = all_clns_trb$data[biopsy_samples]
all_clns_trb_biopsy_meta = all_clns_trb$meta %>% filter(Origin=="biopsy")

# Number of clones
exp_vol <- repExplore(all_clns_trb_biopsy, .method = "clones", .coding = T )
p5 = vis(exp_vol, .by = c("CellType"), .meta = all_clns_trb_biopsy_meta) + scale_fill_manual(values = paletteer_d("LaCroixColoR::Lime")[c(1,3)] ) # #2CB11BFF #BDDE9BFF 
p6 = vis(exp_vol, .by = c("Diagnosis"), .meta = all_clns_trb_biopsy_meta) + scale_fill_manual(values = c("#377EB8", "#E41A1C"))

# Percentrage of repertoire occupied by clones of varying size
repClonality(all_clns_trb_biopsy, .method = "homeo", .clone.types = c(Small = .05, Medium = .15, Large = .25, Hyperexpanded = .5)) %>% vis()

# Inverse Simpson diversity index
div_inv = repDiversity(all_clns_trb_biopsy, .method = "inv.simp")
p7 = vis(div_inv, .by = c("CellType"), .meta = all_clns_trb_biopsy_meta) + scale_fill_manual(values = paletteer_d("LaCroixColoR::Lime")[c(1,3)] ) # #2CB11BFF #BDDE9BFF 
p8 = vis(div_inv, .by = c("Diagnosis"), .meta = all_clns_trb_biopsy_meta) + scale_fill_manual(values = c("#377EB8", "#E41A1C"))


pdf("~/Desktop/Somatic_paper/Figures_pdf/TGS_bulk_lymphocytes/VDJ_calls/TCR_TRB_byCellType.pdf", height = 11, width = 6)
grid.arrange(p5, p7, ncol = 1, nrow = 2) # Arrange plots side-by-side
dev.off()

pdf("~/Desktop/Somatic_paper/Figures_pdf/TGS_bulk_lymphocytes/VDJ_calls/TCR_TRB_byDiagnosis.pdf", height = 11, width = 5)
grid.arrange(p6, p8, ncol = 1, nrow = 2) # Arrange plots side-by-side
dev.off()


################################################################################

######## ####  Combined IGH and TRB analysis ####  #### 

all_clns_combined = c(all_clns_trb_biopsy, all_clns_igh_biopsy)

all_clns_combined_meta = rbind(all_clns_trb_biopsy_meta, all_clns_igh_biopsy_meta)

# all_clns_combined_meta$CellType = factor(all_clns_combined_meta$CellType, levels = c("B Cells", "Plasmablasts", "Plasma" ,"CD4", "CD8"))
# 
# all_clns_combined_meta = all_clns_combined_meta[order(all_clns_combined_meta$CellType),]

order_CellType <- c("B Cells", "Plasmablasts", "Plasma" ,"CD4", "CD8")
order_Diagnosis <- c()

# Number of clones
exp_vol <- repExplore(all_clns_combined, .method = "clones", .coding = T )

p9 <- vis(exp_vol, .by = c("CellType"), .meta = all_clns_combined_meta) + 
  scale_x_discrete(limits = order_CellType) +
  scale_fill_manual(values = c("#0099CCFF", "#2CB11BFF", "#BDDE9BFF", "#66E5FFFF", "#99FFFFFF"))
  
  
p10 = vis(exp_vol, .by = c("Diagnosis"), .meta = all_clns_combined_meta) + scale_fill_manual(values = c("#377EB8", "#E41A1C"))

# Percentrage of repertoire occupied by clones of varying size
repClonality(all_clns_combined, .method = "homeo", .clone.types = c(Small = .05, Medium = .15, Large = .25, Hyperexpanded = .5)) %>% vis()

# Inverse Simpson diversity index
div_inv = repDiversity(all_clns_combined, .method = "inv.simp")

p11 = vis(div_inv, .by = c("CellType"), .meta = all_clns_combined_meta) + 
  scale_x_discrete(limits = order_CellType) +
  scale_fill_manual(values = c("#0099CCFF", "#2CB11BFF", "#BDDE9BFF", "#66E5FFFF", "#99FFFFFF"))

p12 = vis(div_inv, .by = c("Diagnosis"), .meta = all_clns_combined_meta) + scale_fill_manual(values = c("#377EB8", "#E41A1C"))

pdf("~/Desktop/Somatic_paper/Figures_pdf/TGS_bulk_lymphocytes/VDJ_calls/All_clones_IGH_TRB_byCellType.pdf", height = 11, width = 6)
grid.arrange(p9, p11, ncol = 1, nrow = 2) # Arrange plots side-by-side
dev.off()

pdf("~/Desktop/Somatic_paper/Figures_pdf/TGS_bulk_lymphocytes/VDJ_calls/All_clones_IGH_TRB_byDiagnosis.pdf", height = 11, width = 5)
grid.arrange(p10, p12, ncol = 1, nrow = 2) # Arrange plots side-by-side
dev.off()


################################################################################
### Notes: 
#  If data is grouped, then statistical tests for comparing means of groups will be performed, unless .test = FALSE is supplied. 
#  In case there are only two groups, the Wilcoxon rank sum test (https://en.wikipedia.org/wiki/Wilcoxon_signed-rank_test) is performed 
#  (R function wilcox.test with an argument exact = FALSE) for testing if there is a difference in mean rank values between two groups. 
#  In case there more than two groups, the Kruskal-Wallis test (https://en.wikipedia.org/wiki/Kruskal A significant Kruskal-Wallis test 
#   indicates that at least one sample stochastically dominates one other sample. 
#  Adjusted for multiple comparisons P-values are plotted on the top of groups. 
#  P-value adjusting is done using the Holm method (https://en.wikipedia.org/wiki/Holm 
#  You can execute the command ?p.adjust in the R console to see more. 












suppressPackageStartupMessages(
  {library("dndscv")
    library("vcfR")
    #  library("deepSNV")
    library("tidyverse")
    library(stringr)
  }
)

#--------------------------------------------------------------
list.files()

# load all nanoseq muts
nano_muts_all = read.delim("2023-04-19_Sjogrens_CombinedAnnotatedCalls.tsv") # 65634

# fix diagnosis information
for (i in 1:nrow(nano_muts_all)) {
  if (nano_muts_all$donor[i] %in% c("PD45529", "PD45530")) {
    nano_muts_all$tissue[i] <- "Bulk_NonSjogrens"
  }
}

# collapsed by donor
unique_muts_by_donor <- nano_muts_all[which(!(duplicated(paste(nano_muts_all$donor,nano_muts_all$chr,nano_muts_all$pos,nano_muts_all$ref,nano_muts_all$mut, sep = "_")))),]


#--------------------------------------------------------------
# load sample metadata
meta = read.delim("2023-04-19_Sjogrens_SummaryStats_Merged.tsv")

# add columns to metadata
meta$patient = substr(meta$sample, 1, 7)

meta = meta %>% mutate(
  diagnosis = case_when(
    tissue %in% c("Bulk_NonSjogrens") ~ "nonPSS",
    tissue %in% c("Sjogrens_Lymphocyte", "Sjogrens_Epithelium", "Sjogrens_Bulk") ~ "PSS"
  )
)

males = c("PD42761", "PD42764", "PD42766", "PD42770", "PD42773", "PD45531")
females = meta$patient[!meta$patient %in% males]
female_samples = meta$sample[meta$patient %in% females]

meta = meta %>% mutate(
  sex = case_when(
    patient %in% males ~ "male",
    !patient %in% males ~ "female"
  )
)

#--------------------------------------------------------------
# load reference objects
covs = "covariates_20pc_GRCh37-38.epi_strict_outliers.Rdat"
load(covs)

refcds_nano = "RefCDS_GRCh37_vF.Rdat"
load(refcds_nano)

#--------------------------------------------------------------
# dN/dS
setwd("~/Desktop/Code_SjD_paper/dNdS/")


# all muts, collapsed by donor
dnds_all_collapse <- dndscv(unique_muts_by_donor[,1:5], 
                            max_muts_per_gene_per_sample = Inf, 
                            max_coding_muts_per_sample = Inf, 
                            outmats = T, 
                            refdb = refcds_nano, 
                            cv = scores) 

dnds_all_collapse$sel_cv %>% head(10)

write.table(dnds_all_collapse$sel_cv, "Nanoseq2_dnds_all_collapsed_by_donor.tsv", row.names = F, col.names = T, sep = "\t")

#--------------------------------------------------------------
# all Non-PSS samples (doesn't matter if collapsed or not)
dnds_Non_Sjogrens <- dndscv(unique_muts_by_donor[which(unique_muts_by_donor$tissue == "Bulk_NonSjogrens"),1:5], 
                            max_muts_per_gene_per_sample = Inf, 
                            max_coding_muts_per_sample = Inf, 
                            outmats = T, 
                            refdb = refcds_nano, 
                            cv = scores) 

dnds_Non_Sjogrens$sel_cv %>% head(10)

write.table(dnds_Non_Sjogrens$sel_cv, "Nanoseq2_dnds_NonSjogrens_collapsed_by_donor.tsv", row.names = F, col.names = T, sep = "\t")


#--------------------------------------------------------------
# Sjogrens all - collapsed
dnds_Sjogrens_Combined_collapsed <- dndscv(unique_muts_by_donor[which(unique_muts_by_donor$tissue %in% c("Sjogrens_Epithelium","Sjogrens_Bulk","Sjogrens_Lymphocyte")),1:5], 
                                           max_muts_per_gene_per_sample = Inf, 
                                           max_coding_muts_per_sample = Inf, 
                                           outmats = T, 
                                           refdb = refcds_nano,
                                           cv = scores)

dnds_Sjogrens_Combined_collapsed$sel_cv %>% head(20)

write.table(dnds_Sjogrens_Combined_collapsed$sel_cv, "Nanoseq2_dnds_Sjogrens_All_Collapsed_by_donor.tsv", row.names = F, col.names = T, sep = "\t")


#--------------------------------------------------------------

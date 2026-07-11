options(stringsAsFactors = F)
patient=commandArgs(T)[1]
opt=commandArgs(T)[2]
gender=commandArgs(T)[3]

source("germline_exact_binom.R") # See Tim Coorens' github as explained README
source("beta_binom_flt.R") # See Tim Coorens' github as explained README
#source("plot_spectrum.R") # see trinuc_context.R for plotting spectra

library(stringr)
library(dplyr)

#Read in data from cgpVAF:
data = read.table(paste0(patient,"/",patient,".",opt,".tsv"), comment.char="",header=T)
# if sample names contain "-":
colnames(data) = str_replace(colnames(data), "[.]", "-")

Muts = paste(data$Chrom,data$Pos,data$Ref,data$Alt,sep="_")
Genotype = data[,grepl("VAF",colnames(data))&colnames(data)!="PDv37is_VAF"]
NR = data[,grepl("DEP",colnames(data))&colnames(data)!="PDv37is_DEP"]
NV = data[,grepl("MTR",colnames(data))&colnames(data)!="PDv37is_MTR"]

rownames(Genotype)=rownames(NV)=rownames(NR)=Muts
samples = colnames(Genotype)=colnames(NR)=colnames(NV)=gsub("_VAF","",colnames(Genotype))

#lcm_samples=samples[grepl("_lo",samples)]

XY_chromosomal = grepl("X|Y",Muts)
autosomal = !XY_chromosomal
xy_depth=mean(rowMeans(NR[XY_chromosomal,]))
autosomal_depth=mean(rowMeans(NR[autosomal,]))

if(is.na(gender)){
  gender='male'
  if(xy_depth>0.8*autosomal_depth) gender='female'
}

#Filter out germline using exact binomial test
germline=exact.binomial(gender=gender,NV=NV,NR=NR,cutoff = -5) #determine which variants are germline

write.table(Muts[germline],paste0(patient,"/bb_filt/germline_",opt,"_ids.txt"),row.names = F,col.names = F,quote=F)
write.table(Muts[!germline],paste0(patient,"/bb_filt/somatic_",opt,"_ids.txt"),row.names = F,col.names = F,quote=F)

NR_flt=NR[!germline,]
NV_flt=NV[!germline,]

NR_flt_nonzero=NR_flt
NR_flt_nonzero[NR_flt_nonzero==0]=1
shared_muts=rowSums(NV_flt>0)>1

#Use beta-binomial filter on shared muts (unique muts would pass anyway)
rho_est = beta.binom.filter(NR=NR_flt_nonzero[shared_muts,],NV=NV_flt[shared_muts,])
flt_rho=log10(rho_est)<(-1)
rho_filtered_out = rownames(NR_flt_nonzero[shared_muts,])[flt_rho]
write.table(rho_filtered_out,paste0(patient,"/bb_filt/",opt,"_bbinom_filtered_out.txt"))
write.table(rho_est,paste0(patient,"/bb_filt/",opt,"_rho_est.txt"))

NR_flt_2 = NR_flt[!rownames(NR_flt)%in%rho_filtered_out,]
NV_flt_2 = NV_flt[!rownames(NV_flt)%in%rho_filtered_out,]
write.table(NR_flt_2,paste0(patient,"/bb_filt/",opt,"_NR_filtered_all.txt"))
write.table(NV_flt_2,paste0(patient,"/bb_filt/",opt,"_NV_filtered_all.txt"))

NR_flt_nonzero=NR_flt_2
NR_flt_nonzero[NR_flt_nonzero==0]=1

#Convert genotype matrix in binary genotype for mpboot
XY_chromosomal=grepl("X|Y",rownames(NR_flt_2))
autosomal=!XY_chromosomal
genotype_bin=as.matrix(NV_flt_2/NR_flt_nonzero)
if(gender=="male"){
  genotype_bin[autosomal,][genotype_bin[autosomal,]<0.1]=0
  genotype_bin[autosomal,][genotype_bin[autosomal,]>=0.3]=1
  genotype_bin[XY_chromosomal,][genotype_bin[XY_chromosomal,]<0.2]=0
  genotype_bin[XY_chromosomal,][genotype_bin[XY_chromosomal,]>=0.6]=1
  genotype_bin[genotype_bin>0&genotype_bin<1]=0.5
}
if(gender=="female"){
  genotype_bin[genotype_bin<0.1]=0
  genotype_bin[genotype_bin>=0.3]=1
  genotype_bin[genotype_bin>0&genotype_bin<1]=0.5
}

genotype_bin=genotype_bin[,samples]

#Convert genotype matrix to fasta file
dna_strings = list()
muts=as.data.frame(matrix(ncol=4,unlist(strsplit(rownames(genotype_bin),split="_")),byrow = T))
if(opt=="snp"){
  Ref = muts[,3]
  Alt = muts[,4]
  
  dna_strings[1]=paste(Ref,sep="",collapse="")
  
  for (k in 1:ncol(genotype_bin)){
    Mutations = Ref
    Mutations[genotype_bin[,k]==0.5] = '?'
    Mutations[genotype_bin[,k]==1] = Alt[genotype_bin[,k]==1]
    dna_string = paste(Mutations,sep="",collapse="")
    dna_strings[k+1]=dna_string
  }
}
if(opt=="indel"){
  Ref = rep("A",nrow(genotype_bin))
  Alt = rep("T",nrow(genotype_bin))
  
  dna_strings[1]=paste(Ref,sep="",collapse="")
  
  for (n in 1:length(samples)){
    Mutations = Ref
    Mutations[genotype_bin[,n]==0.5] = '?'
    Mutations[genotype_bin[,n]==1] = Alt[genotype_bin[,n]==1]
    dna_string = paste(Mutations,sep="",collapse="")
    dna_strings[n+1]=dna_string
  }
}
bed_muts=cbind(muts[,1],as.numeric(muts[,2])-1,muts[,2])
bed_muts=bed_muts[order(bed_muts[,1],as.numeric(bed_muts[,2])),]
write.table(bed_muts,paste0(patient,"/bb_filt/",opt,"_good.bed"),quote=F,row.names=F,col.names=F,sep="\t")

names(dna_strings)=c("Ancestral",colnames(genotype_bin))
require(seqinr)
write.fasta(dna_strings, names=names(dna_strings),paste0(patient,"/bb_filt/",opt,"_for_MPBoot.fa"))

#Filter original pindel/caveman vcf for good calls, i.e. to look at annotation
if(opt=="snp") algo="caveman_c"
if(opt=="indel") algo="pindel"
for (sample in samples){
  system(paste0("tabix -h -R ",patient,"/bb_filt/",opt,"_good.bed /nfs/cancer_ref01/nst_links/live/2*/",sample,"/",sample,".",algo,".annot.vcf.gz > ",patient,"/bb_filt/",sample,"_good_",opt,".vcf"))
}


for (s in samples){
  if(file.exists(paste0(patient,"/",s,"_complete_final_retained_2.vcf"))){
   Muts_called_by_filters=read.table(paste0(patient, "/",s,"_complete_final_retained_2.vcf"))
    rownames(Muts_called_by_filters)=paste(Muts_called_by_filters$V1,Muts_called_by_filters$V2,Muts_called_by_filters$V4,Muts_called_by_filters$V5,sep="_")
    Muts_called_by_filters=Muts_called_by_filters[rownames(Muts_called_by_filters)%in%rownames(NV_flt_2),]
    plot_spectrum(Muts_called_by_filters[,c(1,2,4,5)],save=paste0(patient,"/bb_filt/",s,"_mut_spectrum.pdf"),add_to_title=paste0(", ",s))
  }
}


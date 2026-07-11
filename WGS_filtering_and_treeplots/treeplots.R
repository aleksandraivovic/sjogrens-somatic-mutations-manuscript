library(ape)
library(ggtree)
library(stringr)
#library(devtools)
#install_git("https://github.com/NickWilliamsSanger/treemut")
library(treemut)

options(stringsAsFactors = F)
#patient=commandArgs(T)[1]
#opt=commandArgs(T)[2]
patient="PD45530"
opt='snp'

# load output from filtering.R step (binomial filters)
NR_flt=read.table(paste0(patient,"/bb_filt/",opt,"_NR_filtered_all.txt"))
NV_flt=read.table(paste0(patient,"/bb_filt/",opt,"_NV_filtered_all.txt"))


# if s~/Desktop/copy_dir/PD42759/# if sample names contain "-":
colnames(NR_flt) = str_replace(colnames(NR_flt), "[.]", "-")
colnames(NV_flt) = str_replace(colnames(NV_flt), "[.]", "-")

tree=read.tree(paste0(patient,"/",opt,"_for_MPBoot.fa.treefile"))

tree$edge.length=rep(1,nrow(tree$edge))
tree=drop.tip(tree,"Ancestral")
NR_flt = as.matrix(NR_flt[,tree$tip.label])
NV_flt = as.matrix(NV_flt[,tree$tip.label])
df=reconstruct_genotype_summary(tree)

#Use treemut package for mutation assignment to trees
res=assign_to_tree(tree=tree, # changed df=df to tree=tree when using treemut R package instead of treemut.R script from farm
                   mtr=NV_flt,
                   dep=NR_flt)

#saveRDS(res,paste0(patient,"/",opt,"_assigned_to_tree.Rdata"))

tree2=tree
#Filter out muts that don't really map to tree:
keep_edge <- res$summary$p_else_where<0.01
# keep_edge <- res$summary$p_else_where<0.1
#keep_edge <- res$summary$p_else_where<=1

edge_length_nonzero = table(res$summary$edge_ml[res$summary$p_else_where<0.01])
#edge_length_nonzero = table(res$summary$edge_ml[keep_edge])
edge_length = rep(0,nrow(tree2$edge))
names(edge_length)=1:nrow(tree2$edge)
edge_length[names(edge_length_nonzero)]=edge_length_nonzero
tree2$edge.length=as.numeric(edge_length)

tree2 = as.polytomy(tree2, feature='branch.length', fun=function(x) as.numeric(x)==0)

# This is fix #1
# Rename the "nodes"
# They're erroneously named according to a percentage derived from mpboot, 
# rather than a unique label. Here we name them sequentially to align with 
# downstream Mapscape expections (if plotting Mapscape)

n_tip <- length(tree2$tip.label)
n_node <- length(tree2$node.label)
tree2$node.label <- paste0(patient,"_",(n_tip+1):(n_tip+n_node))

# New assignment variable based on filtered tree
res2=assign_to_tree(tree=tree2, # changed df=df to tree=tree when using treemut R package instead of treemut.R script from farm
                   mtr=NV_flt,
                   dep=NR_flt)


# write.tree(tree, paste0(patient,"/",opt,"_tree_with_branch_length.tree"))

#Plot trees with and without meaningful branch lengths
p=ggtree(tree2) + geom_tiplab(aes(x=branch),vjust=-0.3)+theme_tree2()+xlim(0,max(fortify(tree2)$x)*1.3)

#pdf(paste0(patient,"/",opt,"_tree_with_branch_length.pdf"))
print(p)
#dev.off()

tree_collapsed=tree2
tree_collapsed$edge.length=rep(1,nrow(tree_collapsed$edge))

#pdf(paste0(patient,"/",opt,"_tree_with_equal_branch_length.pdf"))
ggtree(tree_collapsed) + geom_tiplab(aes(x=branch),vjust=-0.3)+xlim(0,max(fortify(tree_collapsed)$x)*1.3)
#dev.off()

#Write a dataframe with mutations mapped to branches - SampleID will refer to patient_branch
Mutations_per_branch=as.data.frame(matrix(ncol=4,unlist(strsplit(rownames(NR_flt),split="_")),byrow = T))
colnames(Mutations_per_branch)=c("Chr","Pos","Ref","Alt")
# Mutations_per_branch$Branch = tree$edge[res$summary$edge_ml,2]
Mutations_per_branch$Branch = tree2$edge[res2$summary$edge_ml,2]

# Mutations_per_branch=Mutations_per_branch[res$summary$p_else_where<0.01,]
Mutations_per_branch=Mutations_per_branch[keep_edge,]
Mutations_per_branch$Patient = patient
Mutations_per_branch$SampleID = paste(patient,Mutations_per_branch$Branch,sep="_")
# write.table(Mutations_per_branch,paste0(patient,"/",opt,"_assigned_to_branches.txt"),quote=F,row.names=F,sep="\t")

#branch/node numbers can be found using fortify(tree2) 






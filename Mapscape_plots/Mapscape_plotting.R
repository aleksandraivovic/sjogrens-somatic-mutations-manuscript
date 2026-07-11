###### MAPSCAPE plotting ######
######----------------------------------------------------

##
library(ape)
library(ggtree)
library(stringr)
library(plyr)  # ensure imported *before* dplyr
library(dplyr)
library(gplots)
library(funr)
library(treemut)
library(reshape2)
# library(htmltools)
library(wesanderson)
library(htmlwidgets)

## set up directories and scripts
wdir <- "~/Desktop/Somatic_paper/New_mpboot/"
output_dir <- "~/Desktop/Somatic_paper/processed_data/"
dir.create(file.path(output_dir), showWarnings = FALSE)
scripts_dir <-
  funr::get_script_path() # "~/Desktop/Somatic_paper/scripts"
updated_mapscape_path <-
   "~/Desktop/Somatic_paper/Mapscape_Mpboot/mapscape-generator-master/R/mapscape.R"

options(stringsAsFactors = F)
opt = 'snp'

sample_name <- "PD42769"
patient <- "PD42769"
sample_name_ <- "PD42769_"

NV_flt = read.table(paste0(root_dir, "New_mpboot/", sample_name, "/snp_NV_filtered_all.txt"), check.names = F)
NR_flt = read.table(paste0(root_dir, "New_mpboot/", sample_name, "/snp_NR_filtered_all.txt"), check.names = F)
#colnames(NV_flt) = str_replace(colnames(NV_flt), "[.]", "-")
#colnames(NR_flt) = str_replace(colnames(NR_flt), "[.]", "-")

############## 1. Get tree file
tree = read.tree(paste0(root_dir, "New_mpboot/", sample_name, "/", opt, "_for_MPBoot.fa.treefile"))

tree$edge.length = rep(1, nrow(tree$edge))
tree = drop.tip(tree, "Ancestral")
NR_flt = as.matrix(NR_flt[, tree$tip.label])
NV_flt = as.matrix(NV_flt[, tree$tip.label])
df = reconstruct_genotype_summary(tree)

## use treemuts package to assign mutations to branches
res = assign_to_tree(tree = tree, # changed df=df to tree=tree when using treemut R package instead of treemut.R script from farm
                     mtr = NV_flt,
                     dep = NR_flt)

#if (save_rds) {
#  rds_filename <- paste0(patient, "_", opt, "_assigned_to_tree.Rdata")
#  saveRDS(res, paste0(output_dir, rds_filename))
#  rds_filenames[pidx] <- rds_filename
#}

## Filter out muts that don't really map to tree:
keep_edge <- res$summary$p_else_where < 0.01

edge_length_nonzero = table(res$summary$edge_ml[keep_edge])

tree2 <- tree
edge_length = rep(0, nrow(tree2$edge))
names(edge_length) = 1:nrow(tree2$edge)
edge_length[names(edge_length_nonzero)] = edge_length_nonzero
tree2$edge.length = as.numeric(edge_length)

tree2 = as.polytomy(
  tree2,
  feature = 'branch.length',
  fun = function(x)
    as.numeric(x) == 0
)


##  Rename the "nodes"
# They're erroneously named according to a percentage derived from mpboot,
# rather than a unique label. Here we name them sequentially to align with
# mapscape's expectations.

n_tip <- length(tree2$tip.label)
n_node <- length(tree2$node.label)
tree2$node.label <- paste0(patient, "_", (n_tip + 1):(n_tip + n_node))

# New assignment variable based on filtered tree
res = assign_to_tree(tree = tree2, # changed df=df to tree=tree when using treemut R package instead of treemut.R script from farm
                     mtr = NV_flt,
                     dep = NR_flt)

 write.tree(tree, paste0("New_mpboot/", patient,"/",opt,"_tree_with_branch_length.tree"))

 # plot tree with unique branch lengths
 p <- (ggtree(tree2)
       + geom_tiplab(aes(x = branch), vjust = -0.3)
       + theme_tree2()
       + xlim(0, max(fortify(tree2)$x) * 1.3))
 tree_pdf_filename <-
   paste0(patient, "_", opt, "_tree_with_branch_length.pdf")
# pdf(paste0(output_dir, tree_pdf_filename))
# print(p)
# dev.off()
 
 ############## 2. Get mutation table
 #Write a dataframe with mutations mapped to branches - SampleID will refer to
 #patient_branch
 Mutations_per_branch = as.data.frame(matrix(ncol = 4, unlist(strsplit(
   rownames(NR_flt), split = "_"
 )), byrow = T))
 
 colnames(Mutations_per_branch) = c("Chr", "Pos", "Ref", "Alt")
 Mutations_per_branch$Branch = tree2$edge[res$summary$edge_ml, 2]
 
 Mutations_per_branch = Mutations_per_branch[keep_edge,]
 Mutations_per_branch$Patient = patient
 Mutations_per_branch$SampleID = paste(patient, Mutations_per_branch$Branch, sep =
                                         "_")
 
 mutations_per_branch_filename <-
     paste0("New_mpboot/", patient, "_", opt, "_assigned_to_branches.txt")
   
  write.table(
     Mutations_per_branch,
     mutations_per_branch_filename,
     quote = F,
     row.names = F,
     sep = "\t"
   )
   
  #branch/node numbers can be found using fortify(tree2)
  
 
  ######### 3. get clonal prevalence table
  # filter NV_flt to keep only mutations used in filclontered tree2 - to be able to count up muts per sample
  
  NV_flt2 = as.data.frame(NV_flt) # nrow 16,283
  NV_flt2 = NV_flt2 %>% dplyr::filter(keep_edge)
  
  NV_flt2$mutation = rownames(NV_flt2)
  
  Mutations_per_branch_wNV = cbind(Mutations_per_branch, NV_flt2)
  
  clones = Mutations_per_branch_wNV$SampleID %>% unique() # change name to all_clones?
  samples = colnames(Mutations_per_branch_wNV)[grep("PD", colnames(Mutations_per_branch_wNV))]
  
  # number of mutations per clone
  Mutations_per_branch_wNV %>% select(SampleID) %>% table

  ## summary_muts is total mut counts per sample per clone
  summary_muts = data.frame(row.names = samples)
  
  for (clone in clones) {
    a = Mutations_per_branch_wNV[Mutations_per_branch_wNV$SampleID == clone,]
    a = a %>% select(all_of(samples)) %>% colSums() %>% as.data.frame()
    colnames(a) = clone
    summary_muts = cbind(summary_muts, a)
  }

  ## muts_per_sample is counting mutations per sample per clone
  muts_per_sample = data.frame(row.names = samples)
  
  for (clone in clones) {
    a = Mutations_per_branch_wNV[Mutations_per_branch_wNV$SampleID == clone,]
    a = a %>% select(all_of(samples))
    if (dim(a)[1] > 1) {
      a = apply(a > 0, 2, as.integer) %>% colSums() %>% as.data.frame()
    } else {
      a = t(apply(a > 0, 2, as.integer)) %>% colSums() %>% as.data.frame()
    }
    colnames(a) = clone
    muts_per_sample = cbind(muts_per_sample, a)
  }
  
  total_muts_per_sample = rowSums(muts_per_sample)
  total_muts_per_clone = Mutations_per_branch_wNV %>% select(SampleID) %>% table
  
  ## clonal prev table
  clone_muts = t(muts_per_sample) %>% as.data.frame()
  
  rowSums(clone_muts) == colSums(muts_per_sample)
  
  clone_muts$muts_per_clone = rowSums(clone_muts)
  
  ## calculate number of muts per clone per sample divided by total number of muts per clone
  my_clonal_prev = t(muts_per_sample / total_muts_per_sample) %>% as.data.frame()
  
  # melt clonal_prev table
  my_clonal_prev$clone = rownames(my_clonal_prev)
  
  my_clonal_prev_melted = reshape2::melt(
    my_clonal_prev,
    id.vars = "clone",
    variable.name = "sample_id",
    value.name = "clonal_prev"
  )
  
  my_clonal_prev_melted = my_clonal_prev_melted[, c("sample_id", "clone", "clonal_prev")] # reordering
  colnames(my_clonal_prev_melted) = c("sample_id", "clone_id", "clonal_prev")  

  
  ################# 4. get tree edges
  tree2$edge %>% head
  
  tree2$edge %>% dim 
  
  length(clones) 
  length(samples) 
  
  my_tree_edges = tree2$edge
  colnames(my_tree_edges) = c("source", "target")
  
  my_tree_edges = as.data.frame(my_tree_edges)
  my_tree_edges$source = paste0(sample_name_, my_tree_edges$source)
  my_tree_edges$target = paste0(sample_name_, my_tree_edges$target) 

  ################# 5. Get sample locations
  my_sample_loc = read.csv(paste0("New_mpboot/", sample_name,"/",sample_name_,"locations.csv"))
  
  head(my_sample_loc)
  my_sample_loc = my_sample_loc[, c(2, 1, 3:4)]
  colnames(my_sample_loc) = c("sample_id", "location_id", "x", "y")
  
  my_sample_loc <- na.omit(my_sample_loc)
  
#  # If using a half-resolution png due to file size constraints,
#  # adjust the sample x- and y-coordinates accordingly
#  if (grepl('half_res', img_file)){
#    my_sample_loc$x <- my_sample_loc$x %/% 2
#    my_sample_loc$y <- my_sample_loc$y %/% 2
#  }
  
  ## Truncate sample names for legibility
  # truncate_pattern <- ".*_lo00"
   truncate_pattern <- paste0(sample_name, 'a_')
  # truncate_pattern <- ".*_"
  
  my_sample_loc$sample_id <-
    sub(truncate_pattern, "", my_sample_loc$sample_id)
  
  my_clonal_prev_melted$sample_id <-
    sub(truncate_pattern, "", my_clonal_prev_melted$sample_id)

  ################# 6. Adjust clone coloring
  
  # Get major clones
  # top_clones <- my_tree_edges %>%
  #   dplyr::filter(source == "PD")
  # The common ancestor is the only node appearing in "source" but not "target"
  
  common_ancestor <-
    setdiff(my_tree_edges$source, my_tree_edges$target)
  top_clones <- my_tree_edges %>%
    dplyr::filter(source == common_ancestor)
  
  # ID which clones actually present in these samples (will get most distinct colours)
  max_prev <- my_clonal_prev_melted %>%
    dplyr::group_by(clone_id) %>%
    dplyr::summarise(max_prev = max(clonal_prev)) %>%
    dplyr::filter(max_prev >= 0.1)
  clones_present <- max_prev$clone_id
  
  top_clones$present <- top_clones$target %in% clones_present

  # Assign color for top clones, if 8 or fewer, distinct colours chosen,
  # if 8 or fewer are present, distince dark2 pallet used, rest randomly assigned
  # if 12 of fewer are present, disticne Paried pallet used, rest randomly assign
  # otherwise randomly assigned
  if (nrow(top_clones) < 3) {
    top_clone_col <- RColorBrewer::brewer.pal(8, "Dark2")
    top_clone_col <- sample(top_clone_col, nrow(top_clones))
  } else if (nrow(top_clones) >= 3 &&  nrow(top_clones) <= 8) {
    top_clone_col <- RColorBrewer::brewer.pal(nrow(top_clones), "Dark2")
    top_clone_col <- sample(top_clone_col, length(top_clone_col))
  }  else if (sum(top_clones$present) < 3) {
    pres_col <- RColorBrewer::brewer.pal(8, "Dark2")
    pres_col <- sample(pres_col, sum(top_clones$present))
    color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(distinct = T), invert = T)]
    abs_col <- sample(color, nrow(top_clones) - length(pres_col))
  } else if (sum(top_clones$present) >= 3 &&
             sum(top_clones$present) <= 8) {
    pres_col <-
      RColorBrewer::brewer.pal(sum(top_clones$present), "Dark2")
    color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(distinct = T), invert = T)]
    abs_col <- sample(color, nrow(top_clones) - length(pres_col))
  } else if (sum(top_clones$present) <= 12) {
    pres_col <-
      RColorBrewer::brewer.pal(sum(top_clones$present), "Paired")
    color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(distinct = T), invert = T)]
    abs_col <- sample(color, nrow(top_clones) - length(pres_col))
  } else {
    color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(distinct = T), invert = T)]
    top_clone_col <- sample(color, nrow(top_clones))
  }
  
  top_clones$top_col <- ""
  
  if (exists("pres_col")) {
    top_clones[top_clones$present == TRUE, ]$top_col <-
      sample(col2hex(pres_col), length(pres_col), replace = F)
    top_clones[top_clones$present == FALSE, ]$top_col <-
      col2hex(abs_col)
  } else {
    top_clones$top_col <- col2hex(top_clone_col)
  }  
  
  
  library(phylobase)
  tree4 <- as(tree2, "phylo4")
  
  clones <- top_clones$target
  # Get descendants for top clones, ramp colour towards white for descendante
  # For clones with 9 or more descedants second colour is added otherwise too difficult to distinguish
  # For clones with 15 or more descendants, third colour is added
  color_list <- list()
  for (i in 1:length(clones)) {
    this_clone =  clones[i]
    this_clone_shrt = sub(".*_", "", this_clone)
    this_color = top_clones[top_clones$target == this_clone, "top_col"]
    nodes <- phylobase::nodeLabels(tree4)
    node_tbl <- tibble(node = names(nodes), label = nodes)
    new_node <- node_tbl[node_tbl$label == this_clone, ]$node
    new_descendants <-
      descendants(tree4, as.numeric(new_node), type = "ALL")
    if (length(new_descendants) >= 1 &&
        length(new_descendants) <= 9) {
      branch_clones <- names(new_descendants)
      ramp_colors <- colorRampPalette(c(this_color, "white"))
      clone_colors <- ramp_colors(length(branch_clones) + 1)
      clone_colors <- clone_colors[-(length(branch_clones) + 1)]
      color_tbl <-
        tibble(clone_id = branch_clones, colour = clone_colors)
      color_list[[i]] <- color_tbl
    } else if (length(new_descendants) >= 10 &&
               length(new_descendants) < 15) {
      branch_clones <- names(new_descendants)
      # ramp towards white for half of clones, towards black for the other half
      # black ramp is then reveserved and combined with white ramp
      ramp_white <- colorRampPalette(c(this_color, "white"))
      ramp_grey <- colorRampPalette(c(this_color, "black"))
      clone_white <-
        ramp_white(round(length(branch_clones) / 2) + 1)
      clone_white <- clone_white[-(length(clone_white))]
      
      clone_grey <- ramp_grey(round(length(branch_clones) / 2) + 6)
      clone_grey <-
        clone_grey[-c((length(clone_grey) - 3):(length(clone_grey)))]
      
      remain_cols <- length(new_descendants) - length(clone_white)
      
      clone_colors <- c(clone_white, clone_grey[1:remain_cols])
      
      color_tbl <-
        tibble(clone_id = branch_clones, colour = clone_colors)
      color_list[[i]] <- color_tbl
    } else if (length(new_descendants) >= 15) {
      childs <-
        descendants(tree4, as.numeric(new_node), type = "children") # get first level of descendants for new colour scheme
      child_clones <- names(childs)
      if (length(childs) <= 10) {
        pal1 <- wes_palette("Darjeeling1")
        pal2 <- wes_palette("BottleRocket2")
        pals <- c(pal1, pal2)
        child_col <- sample(pals, length(childs))
        child_col <- col2hex(child_col)
      } else {
        color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(distinct = T), invert = T)]
        child_col <- sample(color, length(childs))
        child_col <- col2hex(child_col)
      }
      parent_tbl <-
        tibble(clone_id = this_clone, colour = this_color)
      child_tbl <-
        tibble(clone_id = child_clones, colour = child_col)
      child_tbl <- rbind(parent_tbl, child_tbl)
      child_list <- list()
      for (j in 1:length(child_clones)) {
        this_child = child_clones[j]
        this_color = child_tbl[child_tbl$clone_id == this_child, ]$colour
        this_child_node = node_tbl[node_tbl$label == this_child , ]$node
        if (length(this_child_node) > 0) {
          child_descendants <-
            descendants(tree4, as.numeric(this_child_node), type = "all")
          branch_children <-  names(child_descendants)
          ramp_colors <- colorRampPalette(c(this_color, "white"))
          clone_colors <- ramp_colors(length(branch_children) + 2)
          clone_colors <-
            clone_colors[-c(1, (length(branch_children) + 2))]
          child_tibble <-
            tibble(clone_id = branch_children, colour = clone_colors)
          child_list[[j]] <- child_tibble
        }
        
      }
      child_color_table <- bind_rows(child_list)
      child_color_table <- rbind(child_color_table, child_tbl)
      color_list[[i]] <- child_color_table
    } else {
      color_tbl <-
        tibble(clone_id = as.character(this_clone), colour = this_color)
      color_list[[i]] <- color_tbl
    }
  }
  
  clone_colour_tbl <- bind_rows(color_list)
  PD_tbl <- tibble(clone_id = common_ancestor, colour = "000000")
  clone_colour_tbl <- bind_rows(clone_colour_tbl, PD_tbl)
  clone_colour_tbl$clone_id <-
    as.character(clone_colour_tbl$clone_id)  

  # This is fix #2
  # convert tip names in clone_colour_tbl to follow the convention in my_clonal_prev_melted
  # that is, ensure all nodes/tips are named via the pattern {patient}_{#}
  tips <- phylobase::tipLabels(tree4)
  tip_tbl <- tibble(node = names(tips), label = tips)
  tip_tbl$node_long <- paste0(patient, "_", tip_tbl$node)
  
  # Alternatively, convert tip names in my_clonal_prev_melted and my_tree_edges
  # to follow the convention {patient}_{sample}
  my_clonal_prev_melted$clone_id <-
    mapvalues(my_clonal_prev_melted$clone_id,
              from = tip_tbl$node_long,
              to = tip_tbl$label)
  my_tree_edges$target <-
    mapvalues(my_tree_edges$target,
              from = tip_tbl$node_long,
              to = tip_tbl$label)  

#  # Only necessary for one clone from PD45529. It will throw up
#  # a bunch of warnings about unused values but those can be safely ignored.
#  clone_colour_tbl$clone_id <-
#    mapvalues(clone_colour_tbl$clone_id,
#              from = tip_tbl$node_long,
#              to = tip_tbl$label)
#  
  
  color_tbl = clone_colour_tbl  
  
  ################# 7. Get image file 
  #my_img_ref = "~/Desktop/Somatic_paper/New_mpboot/PD42769/PD42769a-5.png"
  
  ################# run plotting command
  gitdir <- "/Users/aleksandraivovic/Desktop/Somatic_paper/Mapscape_Mpboot/mapscape-generator-master"
  source(paste0(gitdir, "/R/mapscape.R"))
  
  mapscape_path = paste0("~/Desktop/Somatic_paper/New_mpboot/", patient, "/", patient, "_mapscape.html")
  
  mp = mapscape(
    clonal_prev = my_clonal_prev_melted %>% purrr::map_df(rev),
    tree_edges = my_tree_edges %>% purrr::map_df(rev),
    sample_locations = my_sample_loc,
 #   static_tree_image = static_tree_image,
    img_ref = my_img_ref,
    width = 1500,
    height = 1300,
    clone_colours = color_tbl
  )
  
  #################  save HTML file
  htmlwidgets::saveWidget(mp, mapscape_path, selfcontained = T, title = patient)
  
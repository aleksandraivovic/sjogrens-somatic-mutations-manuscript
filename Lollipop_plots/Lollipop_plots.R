## Lollipop plots for all mutations (TGS + WES Nanoseq)

library(ggplot2)
library(ggrepel)
library(maftools)

all_data <- read.delim("Lollipop_plots/annovar_output_select_genes_aa.txt", sep = "\t")

func_map <- c(
  "nonsynonymous SNV"       = "Missense_Mutation",
  "stopgain"                = "Nonsense_Mutation",
  "frameshift deletion"     = "Frame_Shift_Del",
  "frameshift insertion"    = "Frame_Shift_Ins",
  "nonframeshift deletion"  = "In_Frame_Del",
  "nonframeshift insertion" = "In_Frame_Ins",
  "stoploss"                = "Nonstop_Mutation",
  "splicing"                = "Splice_Site"
)

vc_cols <- c(
  "Missense_Mutation" = "#008000",
  "Nonsense_Mutation" = "#000000",
  "Frame_Shift_Del"   = "#FF0000",
  "Frame_Shift_Ins"   = "#FF69B4",
  "In_Frame_Del"      = "#8B4513",
  "In_Frame_Ins"      = "#D2691E",
  "Nonstop_Mutation"  = "#8B0000",
  "Splice_Site"       = "#8B008B"
)

genes <- unique(all_data$Gene.refGene)

for (g in genes) {
  
  gene_data <- all_data[
    all_data$Gene.refGene == g &
      all_data$Func.refGene == "exonic" &
      !is.na(all_data$ExonicFunc.refGene) &
      all_data$ExonicFunc.refGene != "synonymous SNV" &
      !is.na(all_data$AA_Change) &
      all_data$AA_Change != "", ]
  
  if (nrow(gene_data) == 0) next
  
  gene_data$Variant_Classification <- func_map[gene_data$ExonicFunc.refGene]
  
  gene_data$pos <- as.numeric(sapply(
    strsplit(gsub("[^0-9_].*", "", gsub("^[A-Za-z]+", "", gene_data$AA_Change)), "_"),
    function(x) x[1]
  ))
  gene_data <- gene_data[!is.na(gene_data$pos), ]
  if (nrow(gene_data) == 0) next
  
  mut_summary <- as.data.frame(table(
    pos = gene_data$pos,
    conv = gene_data$AA_Change,
    Variant_Classification = gene_data$Variant_Classification
  ))
  mut_summary <- mut_summary[mut_summary$Freq > 0, ]
  colnames(mut_summary)[colnames(mut_summary) == "Freq"] <- "count"
  mut_summary$pos <- as.numeric(as.character(mut_summary$pos))
  mut_summary <- mut_summary[order(mut_summary$pos), ]
  
  # Transcript from data, fall back to maftools default if no domains found
  transcript <- unique(gene_data$Transcript[!is.na(gene_data$Transcript) & gene_data$Transcript != ""])[1]
  
  prot_len <- NA
  domains  <- NULL
  tryCatch({
    prot <- maftools:::.getdomains(geneID = g, refSeqID = transcript)
    if (nrow(prot) == 0) stop("empty")
    prot_len <- max(as.numeric(prot$aa.length), na.rm = TRUE)
    domains  <- as.data.frame(prot[, c("Label", "Start", "End")])
  }, error = function(e) {
    tryCatch({
      prot     <<- maftools:::.getdomains(geneID = g)
      prot_len <<- max(as.numeric(prot$aa.length), na.rm = TRUE)
      domains  <<- as.data.frame(prot[, c("Label", "Start", "End")])
    }, error = function(e2) message("No domain data for ", g))
  })
  
  if (is.na(prot_len)) prot_len <- max(mut_summary$pos) * 1.1
  
  dom_cols <- setNames(
    colorRampPalette(c("#E88C8C","#8CB8E8","#8CE8A0","#F0E68C","#DDA0DD"))(
      max(1, nrow(domains))),
    if (!is.null(domains)) domains$Label else character(0)
  )
  
  # calculate box height of amino acid domain as fration of y-axis (max unique mutations per gene)
  max_count <- max(mut_summary$count)
  box_h     <- max(0.04, 0.025 * max_count)
  
  p <- ggplot() +
    annotate("rect", xmin = 0, xmax = prot_len,
             ymin = 0.58 - box_h, ymax = 0.58,
             fill = "grey85", colour = "black", linewidth = 0.4)
  
  if (!is.null(domains) && nrow(domains) > 0) {
    p <- p +
      geom_rect(data = domains,
                aes(xmin = Start, xmax = End,
                    ymin = 0.58 - box_h * 1.2, ymax = 0.58 + box_h * 0.1,
                    fill = Label),
                colour = "black", linewidth = 0.4) +
      scale_fill_manual(values = dom_cols, guide = "none") +
      geom_text(data = domains,
                aes(x = (Start + End) / 2, y = 0.58 - box_h * 0.55, label = Label),
                size = 3.2, fontface = "italic")
  }
  
  p <- p +
    geom_segment(data = mut_summary,
                 aes(x = pos, xend = pos, y = 0.58, yend = 0.58 + count * 0.35),
                 colour = "grey40", linewidth = 0.4) +
    geom_point(data = mut_summary,
               aes(x = pos, y = 0.58 + count * 0.35,
                   colour = Variant_Classification, size = 3)) +
    scale_colour_manual(name = "Variant class", values = vc_cols, limits = names(vc_cols)) +
 
    geom_text_repel(data = mut_summary,
                    aes(x = pos, y = 0.58 + count * 0.35, label = conv),
                    size = 3.5, angle = 90, nudge_y = 0.15, direction = "x",
                    segment.size = 0.3, segment.color = "grey60",
                    min.segment.length = 0, max.overlaps = Inf,
                    box.padding = 0.25, point.padding = 0) +
    scale_x_continuous(name = "Amino acid position",
                       limits = c(0, prot_len),
                       breaks = pretty(c(0, prot_len), n = 8),
                       expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(name = "Number of mutations",
                       breaks = function(lims) seq(0, floor((lims[2] - 0.58) / 0.35)) * 0.35 + 0.58,
                       labels = function(b) round((b - 0.58) / 0.35)) +
    ggtitle(g) +
    theme_classic(base_size = 12) +
    theme(plot.title         = element_text(face = "bold.italic", hjust = 0.5),
          legend.position    = "right",
          axis.line.x        = element_blank(),
          panel.grid.major.y = element_line(colour = "grey92"))
  
  out_file <- paste0("Lollipop_plots/", g, "_lollipop_all.pdf")
  ggsave(out_file, p, width = 16, height = 8)
  message("Saved: ", out_file)
}
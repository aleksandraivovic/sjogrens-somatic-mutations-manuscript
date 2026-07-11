## Plotting mutations in key genes per patient from Nanoseq Exome only (adapted from Andrew Lawson)
# Load the required libraries
library(ggplot2)
library(ggtext)
# 
# # Recreate the donor_mut_count dataframe
 donor_mut_count_data <- read.csv("~/Desktop/Code_SjD_paper/dNdS/Muts_per_donor_4topgenes.csv")
 
 colnames(donor_mut_count_data) <- c(
   "PD42759", "PD42760", "PD42761", "PD42763", "PD42764", "PD42765", "PD42766", "PD45529", "PD45530",
   "PD42767", "PD42768", "PD42769", "PD42770", "PD42771", "PD42772", "PD42773", "PD42774", "PD45527",
   "PD45528", "PD45531", "PD45532", "PD45533"
 )
 rownames(donor_mut_count_data) <- c("GPR34", "IL1RAPL1", "MTMR1", "RRAGC")
 donor_mut_count <- as.data.frame(donor_mut_count_data)
 
 # Recreate the annotation_df dataframe
 annotation_df <- data.frame(
   Exome_dx = c(272, 216, 230, 215, 221, 189, 222, 26, 99, 82, 204, 200, 112, 146, 64, 39, 110, 207, 66, 36, 49, 37),
   Sex = c("Female", "Female", "Male", "Female", "Male", "Female", "Male", "Female", "Female", "Female", "Female", "Female", "Male", "Female", "Female", "Male", "Female", "Female", "Female", "Male", "Female", "Female"),
   Pathology = c("Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Sjögren", "Sjögren", "Sjögren", "Sjögren", "Sjögren", "Sjögren", "Sjögren", "Sjögren", "Sjögren", "Sjögren", "Sjögren", "Sjögren", "Sjögren")
 )
 rownames(annotation_df) <- colnames(donor_mut_count)

# Calculate total mutation counts
mtmr1_counts <- t(donor_mut_count["MTMR1",])
# mtmr1_counts <- colSums(donor_mut_count)

# Create a dataframe for plotting
plot_df <- as.data.frame(mtmr1_counts)
plot_df$PatientID <- rownames(mtmr1_counts)
plot_df$Count <- plot_df$MTMR1
plot_df$Pathology <- annotation_df[rownames(mtmr1_counts), "Pathology"]
plot_df <- plot_df[,c("PatientID","Count","Pathology")]

# Define colors for pathology
#pathology_colors <- c("Sjögren" = "#E41A1C", "Control" = "#377EB8")
pathology_colors <- c("Sjögren" = "brown2", "Control" = "cornflowerblue")

# To color the axis labels, we create a new column with markdown formatting
# This will be interpreted by element_markdown() in the theme
plot_df$PatientID_colored <- paste0(
  "<span style='color: ",
  ifelse(plot_df$Pathology == "Sjögren", pathology_colors["Sjögren"], pathology_colors["Control"]),
  ";'>",
  plot_df$PatientID,
  "</span>"
)

# Ensure the order of patients is preserved in the plot
plot_df$PatientID_colored <- factor(plot_df$PatientID_colored, levels = plot_df$PatientID_colored)

library(ggpubr)

# Create the plot
pdf(file = "~/Desktop/Code_SjD_paper/dNdS/MTMR1_histogram.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 5.5)


ggplot(plot_df, aes(x = PatientID_colored, y = Count, fill = Pathology)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(
    "Pathology", # Sets the title of the legend
    values = pathology_colors,
    labels = c(
      "Control" = "Biopsy-negative sicca controls",
      "Sjögren" = "Biopsy-positive Sjögren disease cases"
    )) +
  labs(
    x = "",
    y = "Number of Mutations"
  ) +
  theme_pubclean()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text.x = element_markdown(angle = 45, vjust = 1, hjust = 1),
    legend.position = c(0.05, 0.95), # Position legend in the upper-left
    legend.justification = c("left", "top"), # Anchor by top-left corner
    legend.background = element_rect(fill = "white", color = NA) # Give legend a white background
  )

dev.off()

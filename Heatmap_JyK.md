Characterization of extracts from four Amazonian plants by liquid chromatography-mass spectrometry (LC-MS) data analysis:
================
JHESEL ALMEIDA
2025-02-21

## Introduction

This R script aims to record the procedure given for the metabolic profile
of 4 plant species used in Amazonian folk medicine. Each step has a brief explanation, code and graphics.

## HEATMAP_Data

# Install necessary packages (if not already installed)
# The following packages are required for advanced data visualization and manipulation in R.
install.packages("ggpubr")  # For advanced manipulation of ggplot2 graphics
install.packages("cowplot")  # For combining plots and legends
install.packages("grid")     # For low-level graphics in R
install.packages("ggplot2")  # For creating graphics in R

# Install BiocManager if not already installed (required for Bioconductor packages)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")  # For installing Bioconductor packages

# Install ComplexHeatmap from Bioconductor (for advanced heatmaps)
BiocManager::install("ComplexHeatmap")  # For creating advanced heatmaps

# Install additional packages for color management and text handling in graphics
install.packages("circlize")  # For managing colors and scales in heatmaps
install.packages("gridtext")  # For handling text in graphics

# Load necessary libraries
library(dplyr)               # For data manipulation
library(ComplexHeatmap)      # For creating advanced heatmaps
library(circlize)            # For managing colors and scales in heatmaps
library(grid)                # For low-level graphics in R
library(ggplot2)             # For creating graphics in R
library(ggpubr)              # For advanced manipulation of ggplot2 graphics
library(cowplot)             # For combining plots and legends

# Read data from an Excel file
metab_name_hm <- readxl::read_excel("Data/HEATMAP_POS.xlsx", 1)

# Apply logarithmic transformation to the data
hm_scl <- metab_name_hm %>%
  select(Kp_D1:Pa_F3) %>%  # Select data columns
  as.matrix() %>%
  log10()

# Assign row names
rownames(hm_scl) <- metab_name_hm$Metabolite_name

# Handle NA, NaN, and Inf values
hm_scl[is.na(hm_scl) | is.nan(hm_scl) | is.infinite(hm_scl)] <- 0

# Define the color palette for the heatmap (unchanged)
mycol <- colorRamp2(
  breaks = seq(min(hm_scl, na.rm = TRUE), max(hm_scl, na.rm = TRUE), length = 100),
  colors = colorRampPalette(c("#ADD8E6", "white", "red"))(100)
)

# Remove rows with NA in the Superclass column
metab_name_hm <- metab_name_hm[!is.na(metab_name_hm$Superclass), ]

# Define a more professional color palette for Superclass categories
unique_superclass <- unique(metab_name_hm$Superclass)
class_colors <- setNames(
  RColorBrewer::brewer.pal(min(length(unique_superclass), 8), "Set2"),
  unique_superclass
)

# Define a refined color scheme for metabolite identification levels
identification_colors <- c(
  "I" = "#E69F00",  # Orange (Level 1)
  "II" = "#0072B2"  # Blue (Level 2)
)

# Define colors for species with a more professional look
species_colors <- c(
  "Kp" = "#1b9e77",  # Teal
  "Ws" = "#d95f02",  # Orange
  "Gn" = "#7570b3",  # Purple
  "Pa" = "#e7298a"   # Pink
)

# Define colors for conditions (unchanged)
condition_colors <- c("Dry" = "#999999", "Fresh" = "#333333")

# Create row annotations
hm_row_ann <- rowAnnotation(
  Superclass = metab_name_hm$Superclass,
  Identification_level = metab_name_hm$Identification_level,
  col = list(
    Superclass = class_colors,
    Identification_level = identification_colors
  ),
  show_annotation_name = TRUE,
  show_legend = FALSE,
  annotation_name_gp = gpar(fontsize = 8)
)

# Define species information
species_info <- data.frame(
  Species = c(rep("Kp", 6), rep("Ws", 6), rep("Gn", 6), rep("Pa", 6)),
  Condition = rep(c(rep("Dry", 3), rep("Fresh", 3)), 4)
)
rownames(species_info) <- colnames(hm_scl)

# Create top annotation
top_info_ann <- HeatmapAnnotation(
  Species = species_info$Species,
  Condition = species_info$Condition,
  col = list(
    Species = species_colors,
    Condition = condition_colors
  ),
  show_annotation_name = TRUE,
  annotation_name_side = "left",
  gap = unit(2, "mm"),
  border = TRUE
)

# Create the heatmap
column_groups <- factor(
  c(rep("Kp", 6), rep("Ws", 6), rep("Gn", 6), rep("Pa", 6)),
  levels = c("Kp", "Ws", "Gn", "Pa")
)

hm_plot <- Heatmap(
  hm_scl,
  col = mycol,
  border_gp = grid::gpar(col = "black", lty = 0.05),
  rect_gp = grid::gpar(col = "black", lwd = 0.75),
  clustering_distance_columns = "euclidean",
  clustering_method_columns = "complete",
  column_split = column_groups,
  top_annotation = top_info_ann,
  right_annotation = hm_row_ann,
  show_heatmap_legend = TRUE,
  row_names_gp = gpar(fontsize = 8),
  row_names_side = "right",
  row_names_max_width = unit(10, "cm")
)

# Create custom legends
lgd1 <- Legend(
  col_fun = mycol,
  title = "log10 abundance",
  at = seq(floor(min(hm_scl, na.rm = TRUE)), ceiling(max(hm_scl, na.rm = TRUE)), by = 1),
  direction = "horizontal"
)

lgd2 <- Legend(
  labels = names(species_colors),
  legend_gp = gpar(fill = species_colors),
  title = "Species", ncol = 1
)

lgd3 <- Legend(
  labels = names(condition_colors),
  legend_gp = gpar(fill = condition_colors),
  title = "Condition", ncol = 1
)

lgd4 <- Legend(
  labels = unique(metab_name_hm$Superclass),
  legend_gp = gpar(fill = class_colors),
  title = "Metabolite superclass", ncol = 2
)

lgd5 <- Legend(
  labels = names(identification_colors),
  legend_gp = gpar(fill = identification_colors),
  title = "Identification level", ncol = 1
)

# Combine all legends
all_legends <- packLegend(lgd1, lgd2, lgd3, lgd4, lgd5, direction = "horizontal")

# Convert to ggplot objects
gg_heatmap <- grid.grabExpr(draw(hm_plot))
gg_heatmap <- ggpubr::as_ggplot(gg_heatmap)

gg_legend <- grid.grabExpr(draw(all_legends))
gg_legend <- ggpubr::as_ggplot(gg_legend)

# Combine heatmap and legends
gcms_hm <- plot_grid(
  gg_legend,
  gg_heatmap,
  ncol = 1,
  rel_heights = c(0.2, 1)
)

# Display the heatmap
gcms_hm

# Save the heatmap plot
ggsave(filename = "Results/Heatmap_POS.jpeg", plot = gcms_hm,
       width = 9, height = 5, units = "in", dpi = 300, scale = 1.8)


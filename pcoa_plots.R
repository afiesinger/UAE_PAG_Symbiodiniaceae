#!/usr/bin/env Rscript

# Principal Coordinate Analysis (PCoA) of principal coordinates for Cladocopium ITS2 type profiles as output by SymPortal 
# for the overlap dataset of Platygyra daedalea of Howells et al. (2020) (“Year 2012”) and this study (“Year 2022”) separated by sampling location 
# (AA = Al Aqah in the Gulf of Oman, SY = Saadiyat Island in the Persian/Arabian Gulf)
#' ------------------------------------------------------------------------------------------------------------------------------------------------

#' ----------------------- PACKAGES ------------------------
# Load required libraries
library(readr)
library(readxl)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(ggpubr)

#' -------------------------- 2022 DATASET --------------------------
setwd("/path/to/dir")

# A #

# Read the PCoA coordinates CSV file
pcoa_data_A <- read_excel("20230721T154120_braycurtis_profiles_PCoA_coords_A_sqrt_mod.xlsx")

# Read the metadata Excel file
metadata_A <- read_excel("meta_info_A.xlsx")

# Merge PCoA data with metadata
merged_data_A <- left_join(pcoa_data_A, metadata_A, by = c("ITS2_type_profile" = "ITS2_type_profile"))

# Calculate the proportion of variance explained
prop_explained_A <- as.numeric(tail(merged_data_A, 1)[2:6])  # Assuming PCs start from column 2
prop_explained_pc1_A <- prop_explained_A[1]
prop_explained_pc2_A <- prop_explained_A[2]

# Remove the last row (proportion explained)
merged_data_A <- head(merged_data_A, -1)

# C #

# Read the PCoA coordinates CSV file
pcoa_data_C <- read_excel("20230721T154120_braycurtis_profiles_PCoA_coords_C_sqrt_mod.xlsx")

# Read the metadata Excel file
metadata_C <- read_excel("meta_info_C.xlsx")

# Merge PCoA data with metadata
merged_data_C <- left_join(pcoa_data_C, metadata_C, by = c("ITS2_type_profile" = "ITS2_type_profile"))

# Calculate the proportion of variance explained
prop_explained <- as.numeric(tail(merged_data_C, 1)[2:6])  # Assuming PCs start from column 2
prop_explained_pc1 <- prop_explained[1]
prop_explained_pc2 <- prop_explained[2]

# Remove the last row (proportion explained)
merged_data_C <- head(merged_data_C, -1)

# D #

# Read the PCoA coordinates CSV file
pcoa_data_D <- read_excel("20230721T154120_braycurtis_profiles_PCoA_coords_D_sqrt_mod.xlsx")

# Read the metadata Excel file
metadata_D <- read_excel("meta_info_D.xlsx")

# Merge PCoA data with metadata
merged_data_D <- left_join(pcoa_data_D, metadata_D, by = c("ITS2_type_profile" = "ITS2_type_profile"))

# Calculate the proportion of variance explained
prop_explained_D <- as.numeric(tail(merged_data_D, 1)[2:6])  # Assuming PCs start from column 2
prop_explained_pc1_D <- prop_explained_D[1]
prop_explained_pc2_D <- prop_explained_D[2]

# Remove the last row (proportion explained)
merged_data_D <- head(merged_data_D, -1)

### PLOT ###

# Function to create PCoA plot for each genus
create_pcoa_plot <- function(data, prop_explained_pc1, prop_explained_pc2, genus) {
  # Create a color palette for sites
  site_shapes <- c("SA" = 16, "SI" = 15, "SY" = 17, 
                   "SA, SI" = 2, "SA, SY" = 4, "SI, SY" = 5, 
                   "SA, SI, SY" = 8)
  
  # Create a shape palette for species
  species_colors <- c("Phar" = "#825417", "Pdae" = "#FFA000", "Phar, Pdae" = "#737800")
  
  ggplot(data, aes(x = PC1, y = PC2, color = Species, shape = Site)) +
    geom_point(size = 6) +
    geom_text_repel(aes(label = ITS2_type_profile), 
                    size = 6, 
                    max.overlaps = Inf,
                    box.padding = 0.5,
                    point.padding = 0.5,
                    force = 2,
                    hjust = 0.5,
                    vjust = 0.5) +
    scale_color_manual(values = species_colors) +
    scale_shape_manual(values = site_shapes) +
    labs(x = paste0("PC1 (", round(prop_explained_pc1 * 100, 2), "%)"),
         y = paste0("PC2 (", round(prop_explained_pc2 * 100, 2), "%)"),
         title = genus) +
    theme_bw() +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5, face = "italic", size = 25),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_text(size = 18),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 15),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
}

# Create PCoA plots for each genus
plot_A <- create_pcoa_plot(merged_data_A, prop_explained_pc1_A, prop_explained_pc2_A, "Symbiodinium")
plot_C <- create_pcoa_plot(merged_data_C, prop_explained_pc1, prop_explained_pc2, "Cladocopium")
plot_D <- create_pcoa_plot(merged_data_D, prop_explained_pc1_D, prop_explained_pc2_D, "Durusdinium")

# Arrange plots using ggarrange
combined_plot <- ggarrange(plot_A, plot_C, plot_D, 
                           ncol = 3, nrow = 1, 
                           common.legend = TRUE, legend = "right")

combined_plot

# Save the combined plot
ggsave("filename.png", combined_plot, width = 20, height = 10, dpi = 300)

#' ----------------------- 2012 & 2022 COMBINED DATASET ---------------------
setwd("/path/to/dir")

# A #

# Read the PCoA coordinates CSV file
pcoa_data_A <- read_excel("20240207T220214_braycurtis_profiles_PCoA_coords_A_sqrt_modified.xlsx")

# Read the metadata Excel file
metadata_A <- read_excel("meta_info_A.xlsx")

# Merge PCoA data with metadata
merged_data_A <- left_join(pcoa_data_A, metadata_A, by = c("ITS2_type_profile" = "ITS2_type_profile"))

# Calculate the proportion of variance explained
prop_explained_A <- as.numeric(tail(merged_data_A, 1)[2:16])  # Assuming PCs start from column 2
prop_explained_pc1_A <- prop_explained_A[1]
prop_explained_pc2_A <- prop_explained_A[2]

# Remove the last row (proportion explained)
merged_data_A <- head(merged_data_A, -1)

# C #

# Read the PCoA coordinates CSV file
pcoa_data_C <- read_excel("20240207T220214_braycurtis_profiles_PCoA_coords_C_sqrt_modified.xlsx")

# Read the metadata Excel file
metadata_C <- read_excel("meta_info_C.xlsx")

# Merge PCoA data with metadata
merged_data_C <- left_join(pcoa_data_C, metadata_C, by = c("ITS2_type_profile" = "ITS2_type_profile"))

# Calculate the proportion of variance explained
prop_explained <- as.numeric(tail(merged_data_C, 1)[2:54])  # Assuming PCs start from column 2
prop_explained_pc1 <- prop_explained[1]
prop_explained_pc2 <- prop_explained[2]

# Remove the last row (proportion explained)
merged_data_C <- head(merged_data_C, -1)

# D #

# Read the PCoA coordinates CSV file
pcoa_data_D <- read_excel("20240207T220214_braycurtis_profiles_PCoA_coords_D_sqrt_modified.xlsx")

# Read the metadata Excel file
metadata_D <- read_excel("meta_info_D.xlsx")

# Merge PCoA data with metadata
merged_data_D <- left_join(pcoa_data_D, metadata_D, by = c("ITS2_type_profile" = "ITS2_type_profile"))

# Calculate the proportion of variance explained
prop_explained_D <- as.numeric(tail(merged_data_D, 1)[2:33])  # Assuming PCs start from column 2
prop_explained_pc1_D <- prop_explained_D[1]
prop_explained_pc2_D <- prop_explained_D[2]

# Remove the last row (proportion explained)
merged_data_D <- head(merged_data_D, -1)

### PLOT ###

# Function to create PCoA plot for each genus
create_pcoa_plot <- function(data, prop_explained_pc1, prop_explained_pc2, genus) {
  # Create a color palette for sites
  site_shapes <- c("AA" = 15, "SY" = 17, "both" = 8) 
  
  # Create a shape palette for species
  year_colors <- c("2012" = "#ca61b5", "2022" = "#5f278a", "both" = "#D4B9FB")
  
  ggplot(data, aes(x = PC1, y = PC2, color = Year, shape = Site)) +
    geom_point(size = 6) +
    geom_text_repel(aes(label = ITS2_type_profile), 
                    size = 6, 
                    max.overlaps = Inf,
                    box.padding = 0.5,
                    point.padding = 0.5,
                    force = 2,
                    hjust = 0.5,
                    vjust = 0.5) +
    scale_color_manual(values = year_colors) +
    scale_shape_manual(values = site_shapes) +
    labs(x = paste0("PC1 (", round(prop_explained_pc1 * 100, 2), "%)"),
         y = paste0("PC2 (", round(prop_explained_pc2 * 100, 2), "%)"),
         title = genus) +
    theme_bw() +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5, face = "italic", size = 25),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_text(size = 18),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 15),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
}

# Create PCoA plots for each genus
plot_A <- create_pcoa_plot(merged_data_A, prop_explained_pc1_A, prop_explained_pc2_A, "Symbiodinium")
plot_C <- create_pcoa_plot(merged_data_C, prop_explained_pc1, prop_explained_pc2, "Cladocopium")
plot_D <- create_pcoa_plot(merged_data_D, prop_explained_pc1_D, prop_explained_pc2_D, "Durusdinium")

# Arrange plots using ggarrange
combined_plot <- ggarrange(plot_A, plot_C, plot_D, 
                           ncol = 3, nrow = 1, 
                           common.legend = TRUE, legend = "right")

combined_plot

# Save the combined plot
ggsave("filename.png", combined_plot, width = 20, height = 10, dpi = 300)



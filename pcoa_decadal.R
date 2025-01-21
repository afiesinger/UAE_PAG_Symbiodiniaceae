#!/usr/bin/env Rscript

#' ----------------------- PACKAGES ------------------------

# Load required libraries
library(readr)
library(readxl)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(ggpubr)

#' ----------------- PCoA OF CLADOCOPIUM OUTPUT BY SYMPORTAL --------------
setwd("/path/to/dir")

# Read the PCoA coordinates CSV file
pcoa_data <- read_excel("XXX_braycurtis_profiles_PCoA_coords_C_sqrt_modified.xlsx")

# Read the metadata Excel file
metadata <- read_excel("meta_info_C.xlsx")

# Merge PCoA data with metadata
merged_data <- left_join(pcoa_data, metadata, by = c("ITS2_type_profile" = "ITS2_type_profile"))

# Calculate the proportion of variance explained
prop_explained <- as.numeric(tail(merged_data, 1)[2:54])  # Assuming PCs start from column 2
prop_explained_pc1 <- prop_explained[1]
prop_explained_pc2 <- prop_explained[2]

# Remove the last row (proportion explained)
merged_data <- head(merged_data, -1)

# Create the PCoA plot 
ggplot(merged_data, aes(x = PC1, y = PC2, color = factor(Year), shape = Site)) +
  geom_point(size = 5) +
  geom_text_repel(aes(label = ITS2_type_profile), 
                  size = 5, 
                  max.overlaps = Inf,
                  force_pull = 5,
                  color = "black") +
  theme_bw() +
  labs(title = "",
       x = paste0("PC1 (", round(prop_explained_pc1 * 100, 2), "%)"),
       y = paste0("PC2 (", round(prop_explained_pc2 * 100, 2), "%)"),
       color = "Year",
       shape = "Site") +
  scale_color_manual(values = c("2012" = "#ca61b5", "2022" = "#5f278a", "both" = "#D4B9FB")) +
  scale_shape_manual(values = c("AA" = 15, "SY" = 17)) +
  theme(legend.position = "right",
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.title = element_text(size = 15),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank()) 

# Save the plot
ggsave("HF_SYAA_overlap_PCoA_profiles_C.pdf", width = 10, height = 8, dpi = 300)
ggsave("HF_SYAA_overlap_PCoA_profiles_C.png", width = 10, height = 8, dpi = 300)


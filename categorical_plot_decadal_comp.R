#!/usr/bin/env Rscript


#' -------------------------- PACKAGES -----------------------
rm(list = ls())

library(readxl)
library(ggplot2)
library(Cairo)
library(dendextend)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(Redmonder)

#' --------------------------- PREPARE DATA ----------------------------

# read cleaned excel from combined run
setwd("/path/to/dir")

data <- read_excel("filename.xlsx")

# get only profiles without metadata
its2_cols_subset <- names(data)[4:ncol(data)]

# Create a logical matrix indicating presence/absence of ITS2 type profiles
its2_presence_subset <- data[, its2_cols_subset] > 0

# Get the names of columns (ITS2 type profiles) with at least one TRUE entry
its2_names_subset <- its2_cols_subset[colSums(its2_presence_subset) > 0]

# Print the list of ITS2 type profile names
print(its2_names_subset)

# A
# "A1"                                           "A1-A1bw-A1bf-A1bx-A1eb"                       "A1-A1bw-A1bx"                                
# "A1-A1bf-A1bw"                                 "A1/A1bx-A1bw-A1bf"                            "A1-A1nz-A1bw"                                
# "A1-A1eq-A1ep"                                 "A1-A13c-A1lz"                                 "A1-A13a" 

# C
# "C3-C3gulf-C3d-C3i-C115c-C115b-C3ai-C3ak-C18a" "C3/C3c-C3gulf"                                "C3-C3gulf-C3d-C3i-C115c-C3ak-C3al"           
# "C15h-C15o-C15k"                               "C3yd/C3yc-C3-C3gulf"                          "C3-C3gulf-C3ye-C3d-C3i-C115c-C3ak"           
# "C3-C3gulf-C3d-C3i-C115c-C115b-C3ak"           "C39-C39q-C1-C39b-C39p"                        "C3-C3gulf-C3d-C3ai-C3ak-C115b-C3i-C115c"     
# "C3/C3c-C3gulf-C3aq"                           "C3-C3c-C3gulf-C3l-C3m-C3r"                    "C3/C15h-C3gulf"                              
# "C3-C3d-C3i-C3gulf-C115c-C3ai-C3ak"            "C15h-C15k"                                    

# D
# "D5-D5a-D5f-D4"                                 "D5-D5a-D5ai-D4-D5f-D5v"                       "D5/D5a-D4-D4a-D2"                             
# "D5-D5c-D4a-D5b-D4-D5i-D5a"                     "D1/D2"                                        "D5a-D5-D5ah-D4-D4a"                           
# "D5a-D5-D5ah-D4-D5d"                            "D5-D5a-D4-D4a-D4b"                            "D5-D5a-D5f"                                   
# "D1"                                            "D5-D4-D5a"                                    "D1/D4/D2-D4c-D1c-D1h"     

######################################################################################################
###                                          PLOTS                                                 ###
######################################################################################################

#' ---------------------------- DENDROGRAM ----------------------------

# distance matrices
setwd("/path/to/dir")

# Read and process distance matrices and names
Adist <- as.matrix(read.delim("mod_unifrac_profile_distances_A_sqrt_SYAAPdae.dist", sep = "", header = F))
Cdist <- as.matrix(read.delim("mod_unifrac_profile_distances_C_sqrt_SYAAPdae.dist", sep = "", header = F))
Ddist <- as.matrix(read.delim("mod_unifrac_profile_distances_D_sqrt_SYAAPdae.dist", sep = "", header = F))

Anames <- sort(Adist[,1])
Cnames <- sort(Cdist[,1])
Dnames <- sort(Ddist[,1])

# Create distance matrices for hierarchical clustering
dist_mats <- list(dist(Adist[,-1]), dist(Cdist[,-1]), dist(Ddist[,-1]))
name_lists <- list(Anames, Cnames, Dnames)

# Create hierarchical cluster objects and dendrograms
hclusts <- lapply(dist_mats, hclust)
dendrograms <- lapply(hclusts, as.dendrogram)

# Define the dendrogram plotting function
dendrogram_sp <- function(dend, names, orientation = 'left', lwd = 3, 
                          las = 1, color = 'black', size = 5) {
  
  # Set custom labels
  labels(dend) <- names
  
  # Reverse the order of the labels in the dendrogram
  dend <- rev(dend)
  
  # Plot dendrogram
  plot(dend, main = "", ylab = "", xlab = "", 
       sub = "", horiz = (orientation %in% c('left', 'right')), 
       leaflab = "perpendicular", cex.lab = size, 
       las = las, lwd = lwd, 
       edgePar = list(col = color), 
       xaxt = "n", yaxt = "n",
       cex = size)
}

# Set up the plotting area to have 3 rows and save with Cairo
par(mfrow = c(3, 1), mar = c(0, 10, 0, 25) + 0.1, oma = c(2, 0, 2, 0))

# Loop through and plot each dendrogram
for (i in 1:3) {
  dend <- dendrogram_sp(dendrograms[[i]], name_lists[[i]])
}

#' ---------------------- DOTPLOT FOR NUMBER OF PROFILES -----------------

# read data
setwd("/path/to/dir")
data <- read_excel("filename.xlsx")

df = as.data.frame(data)

dfH = df %>% filter(Dataset == "Howells")
dfF = df %>% filter(Dataset == "Fiesinger")

# 2012 #
# transpose data to calculate total no. of profiles in all samples
transposed = t(dfH[, -c(2:4)])
colnames(transposed) = transposed[1,]
transposed = transposed[-1,]

# make logical matrix
transposed <- transposed > 0

# count all TRUE entries in rows 
num_profs_2012 = as.data.frame(rowSums(transposed))
colnames(num_profs_2012) = "2012"

# 2022 #
# transpose data to calculate total no. of profiles in all samples
transposed = t(dfF[, -c(2:4)])
colnames(transposed) = transposed[1,]
transposed = transposed[-1,]

# make logical matrix
transposed <- transposed > 0

# count all TRUE entries in rows 
num_profs_2022 = as.data.frame(rowSums(transposed))
colnames(num_profs_2022) = "2022"

# rehape data
num_profs_all = cbind(rownames(num_profs_2012), num_profs_2012, num_profs_2022)
colnames(num_profs_all) = c("profs", "howells", "fiesinger")
rownames(num_profs_all) = NULL

desired_order <- c("A1", "A1-A13a", "A1-A13c-A1lz", "A1-A1bf-A1bw",  
                   "A1-A1bw-A1bf-A1bx-A1eb", "A1-A1bw-A1bx", "A1-A1eq-A1ep", "A1-A1nz-A1bw", 
                   "A1/A1bx-A1bw-A1bf", "C15h-C15k",  
                   "C15h-C15o-C15k", "C3-C3c-C3gulf-C3l-C3m-C3r", 
                   "C3-C3d-C3i-C3gulf-C115c-C3ai-C3ak",
                   "C3-C3gulf-C3d-C3ai-C3ak-C115b-C3i-C115c", "C3-C3gulf-C3d-C3i-C115c-C115b-C3ai-C3ak-C18a", 
                   "C3-C3gulf-C3d-C3i-C115c-C115b-C3ak", "C3-C3gulf-C3d-C3i-C115c-C3ak-C3al", 
                   "C3-C3gulf-C3ye-C3d-C3i-C115c-C3ak", "C3/C15h-C3gulf", 
                   "C3/C3c-C3gulf", "C3/C3c-C3gulf-C3aq", 
                   "C39-C39q-C1-C39b-C39p", "C3yd/C3yc-C3-C3gulf", "D1", "D1/D2", 
                   "D1/D4/D2-D4c-D1c-D1h", "D5-D4-D5a", "D5-D5a-D4-D4a-D4b", "D5-D5a-D5ai-D4-D5f-D5v", 
                   "D5-D5a-D5f", "D5-D5a-D5f-D4", "D5-D5c-D4a-D5b-D4-D5i-D5a", 
                   "D5a-D5-D5ah-D4-D4a", "D5a-D5-D5ah-D4-D5d")

# Pivot the dataframe to long format
data_long <- num_profs_all %>%
  pivot_longer(cols = -profs, names_to = "variable", values_to = "value")

# Convert value column to numeric
data_long$value <- as.numeric(data_long$value)

# Set the factor levels for variable
data_long$variable <- factor(data_long$variable, levels = c("howells", "fiesinger"))

# Set the factor levels for profs based on desired_order
data_long$profs <- factor(data_long$profs, levels = rev(desired_order))

# Create a new column for non-zero values
data_long$non_zero_value <- ifelse(data_long$value != 0, data_long$value, NA)

# Remove rows where profs is NA
data_long <- data_long %>% filter(!is.na(profs))

# Ensure profs factor levels are updated after filtering
data_long$profs <- factor(data_long$profs, levels = rev(intersect(desired_order, unique(data_long$profs))))

# CREATE THE PLOT #
ggplot(data_long, aes(x = variable, y = profs)) +
  # Add bold "X" for zero or NA values with the same color as the dots
  geom_text(aes(label = ifelse(is.na(value) | value == 0, "X", ""),
                color = variable), 
            size = 3, fontface = "bold") +
  # Add visible points only for non-zero values
  geom_point(aes(size = non_zero_value, color = variable)) +
  scale_size_continuous(range = c(1, 10)) +
  scale_color_manual(values = c("howells" = "#CA61B5", "fiesinger" = "#5F278A")) +
  labs(title = "",
       x = "",
       y = "",
       size = "") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10, hjust = 0),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggsave("dotplot_its2_profs_NEW_nolabs.png", dpi = 300, height = 12, width = 8, unit = "in")
ggsave("dotplot_its2_profs_NEW_nolabs.pdf", dpi = 300, height = 12, width = 8, unit = "in")

#' ------------------------ CATEGORICAL INFORMATION -------------------

# read filtered data from above
setwd("/path/to/dir")

# NAs already removed
filtered_data_categ <- read_excel("filename.xlsx", sheet = 2)
df = as.data.frame(filtered_data_categ)

# Function to get Dataset and Site info for a column
get_info <- function(col_name) {
  info <- df %>%
    filter(!!sym(col_name) > 0) %>%
    select(Dataset, Site, Species) %>%
    distinct()
  
  datasets <- paste(unique(info$Dataset), collapse = ", ")
  sites <- paste(unique(info$Site), collapse = ", ")
  species <- paste(unique(info$Species), collapse = ", ")
  
  c(Datasets = datasets, Sites = sites, Species = species)
}

# Get all column names starting from the 5th column
col_names <- names(df)[5:ncol(df)]

# Create a list to store results
results <- list()

# Loop through each column and get the information
for (col in col_names) {
  info <- get_info(col)
  if (!all(info == "")) {
    results[[col]] <- info
  }
}

# Convert the results list to a data frame
results_df <- data.frame(
  Profile = names(results),
  do.call(rbind, results),
  stringsAsFactors = FALSE
)

# Reorder columns
results_df <- results_df %>%
  select(Profile, Datasets, Sites, Species)
rownames(results_df) = NULL

# View the first few rows of the results
head(results_df)

# Prepare the data for plotting
plot_data <- results_df %>%
  separate_rows(Sites, sep = ", ") %>%
  separate_rows(Species, sep = ", ") %>%
  pivot_longer(cols = c(Sites, Species), names_to = "Category", values_to = "Value") %>%
  mutate(Present = 1) %>%
  # Complete the data correctly
  complete(Profile, Category, Value, fill = list(Present = 0)) %>%
  # Filter out incorrect combinations
  filter((Category == "Sites" & Value %in% c("SY", "SI")) | 
           (Category == "Species" & Value %in% c("Pdae"))) %>%
  # Order factors
  mutate(Value = factor(Value, levels = c("Pdae", "SY", "SI")),
         Category = factor(Category, levels = c("Species", "Sites"))) %>%
  # Assign colors based on Value and Present
  mutate(fill_color = case_when(
    Category == "Species" & Value == "Pdae" & Present == 1 ~ "#65A595",
    Category == "Sites" & Value == "SY" & Present == 1 ~ "#e8744a",
    Category == "Sites" & Value == "SI" & Present == 1 ~ "#fcc762",
    Present == 0 ~ "white"))

plot_data$Profile <- factor(plot_data$Profile, levels = rev(unique(plot_data$Profile)))

#' ------------------------ PLOT CATEGORICAL DIAGRAM --------------------

# Create the plot
ggplot(plot_data, aes(y = Profile, x = Value, fill = fill_color)) +
  geom_tile(color = "white") +
  scale_fill_identity() +
  facet_grid(~ Category, scales = "free_x", space = "free_x") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold")) +
  labs(x = NULL, y = NULL) 

ggsave(filename = "categ_plot_species_site_NEW.png", dpi = 300, width = 8, height = 12, units = "in")
ggsave(filename = "categ_plot_species_site_NEW.pdf", dpi = 300, width = 8, height = 12, units = "in")

#!/usr/bin/env Rscript

# normalize ITS2 type profiles output by SymPortal
# plot ITS2 type profiles as barplots ordered by decreasing order of profile abundance
#' ------------------------------------------------------------------------------------

#' ----------------------------------- PACKAGES AND WORKING ENVIRONMENT ----------------------
rm(list = ls())

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl)
library(gridExtra)

setwd("path/to/dir")

#' ----------------------------- READ & MODIFY ITS2 POST-MED SEQUENCES OUTPUT BY SYMPORTAL -----------------------------

#' ------------------- THRESHOLD <1% SEQUENCE ABUNDANCE ---------------------
# to remove very low abundance sequences

# Read the data
data <- read.delim('XXX.seqs.relative.abund_and_meta.txt')
head(data)

# Select columns to keep
columns_to_keep <- c("sample_name", "noName.Clade.A", "noName.Clade.B", 
                     "noName.Clade.C", "noName.Clade.D", "noName.Clade.E", "noName.Clade.F", 
                     "noName.Clade.G", "noName.Clade.H", "noName.Clade.I")

# Subset the dataset
metadata <- data[, columns_to_keep]
colnames(metadata) <- c("Sample", "A_other", "B_other", "C_other", "D_other", "E_other", "F_other", "G_other", "H_other", "I_other")
seqs <- data[, c(40:length(data))]

data_filt <- cbind(metadata, seqs)

# Function to replace values < 0.01 with 0
replace_small_values <- function(x) {
  ifelse(x < 0.01, 0, x)
}

# Apply the function to sequence columns (not the total "clade" columns)
data_filt[, 2:ncol(data_filt)] <- lapply(data_filt[, 2:ncol(data_filt)], replace_small_values)
final_data <- data_filt

# count total sequences in clades #
# Identify columns for each clade
clade_A_cols <- grep("^A", colnames(final_data), value = TRUE)
clade_B_cols <- grep("^B", colnames(final_data), value = TRUE)
clade_C_cols <- grep("^C", colnames(final_data), value = TRUE)
clade_D_cols <- grep("^D", colnames(final_data), value = TRUE)
clade_E_cols <- grep("^E", colnames(final_data), value = TRUE)
clade_F_cols <- grep("^F", colnames(final_data), value = TRUE)
clade_G_cols <- grep("^G", colnames(final_data), value = TRUE)
clade_H_cols <- grep("^H", colnames(final_data), value = TRUE)
clade_I_cols <- grep("^I", colnames(final_data), value = TRUE)

# Create summary columns
final_data$Symbiodinium_Sum <- rowSums(final_data[, clade_A_cols])
final_data$Brevolium_Sum <- final_data[, clade_B_cols]
final_data$Cladocopium_Sum <- rowSums(final_data[, clade_C_cols])
final_data$Durusdinium_Sum <- rowSums(final_data[, clade_D_cols])
final_data$Effrenium_Sum <- final_data[, clade_E_cols]
final_data$Fugacium_Sum <- rowSums(final_data[, clade_F_cols])
final_data$Gerakladium_Sum <- rowSums(final_data[, clade_G_cols])
final_data$CladeH_Sum <- final_data[, clade_H_cols]
final_data$CladeI_Sum <- final_data[, clade_I_cols]

# View the first few rows of the summarized data
head(final_data)

# Save the summarized data to a new file
write.csv(final_data, "filename.csv")

# COUNT OCCURRENCES -- how many samples have ITS2 sequences from one genus, from two genera, etc?
data <- final_data

# Get the last 9 columns
last_nine_cols <- data[,((length(data)-8):(length(data)))]

# Function to check if a value should be considered empty
is_effectively_empty <- function(x) {
  is.na(x) | x == "" | x == 0 | abs(as.numeric(x)) < 1e-10
}

# Apply the function to each element of the last 9 columns
last_nine_cols_cleaned <- apply(last_nine_cols, 2, function(col) {
  ifelse(is_effectively_empty(col), NA, col)
})

# Count non-empty entries for each row
counts <- rowSums(!is.na(last_nine_cols_cleaned))

# Create a table of counts
result <- table(counts)

# Print the results
print(result)

#' ------------------------------- READ ITS2 TYPE PROFILE DATA OUTPUT FROM SYMPORTAL -----------------------------------

# ITS2 type profile data from SymPortal (https://symportal.org)
# metadata file contains the species name and sampling location of all the samples

# Load data
its2Profs <- read.delim("396_20230721T100004_DBV_20230721T154120.profiles.absolute.abund_and_meta_clean.txt", header = TRUE, check.names = FALSE)
its2MetaData <- read_excel("its2_profiles_metadata.xlsx")

# Combine metadata with abundance data
its2Profs <- cbind(its2MetaData[, c(4,6,8)], its2Profs[, -1]) 

# Rename metadata columns
colnames(its2Profs)[1:3] <- c("Genotype", "Species", "Site")

# Split datasets by Species and Site for plotting

its2ProfsPhar = its2Profs %>% filter(Species == "Phar") %>% group_by(Genotype, Site)
its2ProfsPdae = its2Profs %>% filter(Species == "Pdae") %>% group_by(Genotype, Site)

its2ProfsPharSA = its2ProfsPhar %>% filter(Site == "SA") %>% group_by(Genotype)
its2ProfsPharSY = its2ProfsPhar %>% filter(Site == "SY") %>% group_by(Genotype)
its2ProfsPharSI = its2ProfsPhar %>% filter(Site == "SI") %>% group_by(Genotype)

its2ProfsPharSY = its2ProfsPdae %>% filter(Site == "SY") %>% group_by(Genotype)
its2ProfsPdaeSI = its2ProfsPdae %>% filter(Site == "SI") %>% group_by(Genotype)

# Convert datasets to long format and arrange abundances in decreasing order
process_data <- function(data) {
  data %>%
    pivot_longer(cols = -c(Genotype, Species, Site), names_to = "Profile", values_to = "Abundance") %>%
    group_by(Profile) 
}

its2ProfsPharSA_long <- process_data(its2ProfsPharSA)
its2ProfsPharSY_long <- process_data(its2ProfsPharSY)
its2ProfsPharSI_long <- process_data(its2ProfsPharSI)
its2ProfsPdaeSY_long <- process_data(its2ProfsPdaeSY)
its2ProfsPdaeSI_long <- process_data(its2ProfsPdaeSI)

# modify dataframes for plotting by decreasing abundance of ITS2 type profiles

its2ProfsPharSA_long <- its2ProfsPharSA_long %>%
  arrange(desc(Abundance))

levels <- unique(its2ProfsPharSA_long$Genotype)
its2ProfsPharSA_long$Genotype <- factor(its2ProfsPharSA_long$Genotype, levels = levels)

its2ProfsPharSY <- its2ProfsPharSY %>%
  arrange(desc(count))

levels <- unique(its2ProfsPharSY$Genotype)
its2ProfsPharSY$Genotype <- factor(its2ProfsPharSY$Genotype, levels = levels)

its2ProfsPharSI <- its2ProfsPharSI %>%
  arrange(desc(count))

levels <- unique(its2ProfsPharSI$Genotype)
its2ProfsPharSI$Genotype <- factor(its2ProfsPharSI$Genotype, levels = levels)

its2ProfsPharSY <- its2ProfsPharSY %>%
  arrange(desc(count))

levels <- unique(its2ProfsPharSY$Genotype)
its2ProfsPharSY$Genotype <- factor(its2ProfsPharSY$Genotype, levels = levels)

its2ProfsPharSI <- its2ProfsPharSI %>%
  arrange(desc(count))

levels <- unique(its2ProfsPharSI$Genotype)
its2ProfsPharSI$Genotype <- factor(its2ProfsPharSI$Genotype, levels = levels)

#' ------------------------------------ SET COLOURS FOR PLOTS -----------------------------------

cols_2022 <- c(
  "A1" = "#FF6100", 
  "A1-A13a" = "#FA8F11", 
  "A1-A1bv" = "#F5C372", 
  "A1-A1bw" = "#F6E496", 
  "A1-A1bw-A1bf-A1bx-A1eb" = "#fcaf58", 
  "A1-A1bw-A1bx" = "#E7C703", 
  "A1-A1eq-A1ep" = "#FFFF7A", 
  
  "C15" = "#800080",
  "C15h" = "#8a2be2",
  "C15h-C15o-C15k" = "#b66ee8", 
  "C15/C116" = "#cd34b5",
  "C3-C3cc-C3gulf-C3ye" = "#a01427", 
  "C3-C3d-C3i-C3gulf-C115c-C3ai-C3ak" = "#cc2936", 
  "C3-C3gulf-C3cc" = "#c75c67", 
  "C3-C3gulf-C3d-C3i-C115c-C115b-C3ai-C3ak" = "#e5989b", 
  "C3-C3gulf-C3d-C3i-C115c-C115b-C3ak" = "#ffbfc0", 
  "C3-C3gulf-C3d-C3i-C115c-C3ak-C3al" = "#ffe8f0", 
  "C3-C3gulf-C3ye" = "#e0acd0", 
  "C3/C15h-C3gulf" = "#bf69a2", 
  "C3/C3c-C3gulf" = "#e57f68", 
  "C3/C3gulf" = "#fd7266",
  "C3/C3gulf-C115d" = "#ff0000",
  "C39-C39q-C1-C39b-C39p" = "#a13b42",
  "C3yd" = "#bd9092",
  "C3yd/C3yc-C3-C3gulf" = "#d1c5c5",  
  
  "D1" = "#000080", 
  "D1-D4-D4c-D17ao-D17ap-D17aq" = "#114477",
  "D1-D4-D4c-D17ap-D17aq-D17ao-D17ar" = "#4169e1",
  "D1-D4-D4c-D17d" = "#75a4f8",
  "D1-D4-D4c-D17d-D1r-D17c-D17e" = "#1ab2e5b9",  
  "D1-D4-D4c-D17d-D1r-D17c-D17e-D17j" = "#5e8ebf",
  "D1/D2" = "#5f9ea0", 
  "D1/D2-D4-D4c-D1r" = "#57d4d5", 
  "D1/D4-D4c-D1h" = "#097d6cd0", 
  "D1/D4-D4c-D1r" = "#bcf6f5", 
  "D1/D4/D2-D4c-D1c-D1h" = "#b5d0ff", 
  "D5-D5a-D4-D4a-D4b" = "#bcdbdb", 
  "D5-D5a-D5ai-D4-D5f-D5v" = "#81c784", 
  "D5-D5a-D5f" = "#44aa77", 
  "D5-D5c-D4a-D5b-D4-D5i-D5a" = "#558b2f", 
  "D5/D5a-D4-D4a-D2" = "#1ca62f", 
  "D5a-D5-D5ah-D4" = "#006400"
)

#' ------------------------------- PLOTTING ------------------------

phar_sa = ggplot() +
  geom_bar(aes(y = Abundance, x = Genotype, fill = Profile), data = its2ProfsPharSA_long, stat = "identity", position = "fill")  + 
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) + 
  labs(y = "Percentage", x = "Host genotype", title = expression(paste(italic("N"), " = 40"))) +
  scale_fill_manual(values = cols_2022) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0, size = 6, color = "black"), 
        plot.title = element_text(hjust = 0.5), 
        legend.text = element_text(size = 8, colour = "black"), 
        legend.title = element_text(color = "black", hjust = 0.5, angle = 90), 
        legend.key.size = unit(0.2, "cm"),
        legend.position = "none",
        axis.title.y = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.key = element_blank(),
        legend.background = element_blank(),
        plot.background = element_blank())

phar_sy = ggplot() +
  geom_bar(aes(y = Abundance, x = Genotype, fill = Profile), data = its2ProfsPharSY_long, stat = "identity", position = "fill")  + 
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) + 
  labs(y = "", x = "", title = expression(paste(italic("N"), " = 39"))) +
  scale_fill_manual(values = cols_2022, breaks = p) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0, size = 6, color = "black"), 
        plot.title = element_text(hjust = 0.5), 
        legend.text = element_text(size = 8, colour = "black"), 
        legend.title = element_text(color = "black", hjust = 0.5, angle = 90), 
        legend.key.size = unit(0.2, "cm"),
        legend.position = "none",
        axis.title.y = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.key = element_blank(),
        legend.background = element_blank(),
        plot.background = element_blank())

phar_si = ggplot() +
  geom_bar(aes(y = Abundance, x = Genotype, fill = Profile), data = its2ProfsPharSI_long, stat = "identity", position = "fill")  + 
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) + 
  labs(y = "Percentage", x = "", title = expression(paste(italic("N"), " = 40"))) +
  scale_fill_manual(values = cols_2022, breaks = p) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0, size = 6, color = "black"), 
        plot.title = element_text(hjust = 0.5), 
        legend.text = element_text(size = 8, colour = "black"), 
        legend.title = element_text(color = "black", hjust = 0.5, angle = 90), 
        legend.key.size = unit(0.2, "cm"),
        legend.position = "none",
        axis.title.y = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.key = element_blank(),
        legend.background = element_blank(),
        plot.background = element_blank())

pdae_sy = ggplot() +
  geom_bar(aes(y = Abundance, x = Genotype, fill = Profile), data = its2ProfsPdaeSY_long, stat = "identity", position = "fill")  + 
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) + 
  labs(y = "Percentage", x = "Host genotype", title = expression(paste(italic("N"), " = 39"))) +
  scale_fill_manual(values = cols_2022, breaks = p) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0, size = 6, color = "black"), 
        plot.title = element_text(hjust = 0.5), 
        legend.text = element_text(size = 8, colour = "black"), 
        legend.title = element_text(color = "black", hjust = 0.5, angle = 90), 
        legend.key.size = unit(0.2, "cm"),
        legend.position = "none",
        axis.title.y = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.key = element_blank(),
        legend.background = element_blank(),
        plot.background = element_blank())

pdae_si = ggplot() +
  geom_bar(aes(y = Abundance, x = Genotype, fill = Profile), data = its2ProfsPdaeSI_long, stat = "identity", position = "fill")  + 
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) + 
  labs(y = "", x = "Host genotype", title = expression(paste(italic("N"), " = 40"))) +
  scale_fill_manual(values = cols_2022, breaks = p) +
  guides(fill = guide_legend(ncol = 1, title = expression(paste(italic("ITS2"), " type profile")))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0, size = 6, color = "black"), 
        plot.title = element_text(hjust = 0.5), 
        legend.text = element_text(size = 12, colour = "black"), 
        legend.title = element_text(color = "black", hjust = 0.5), 
        legend.key.size = unit(0.5, "cm"),
        legend.position = "right",
        axis.title.y = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.key = element_blank(),
        legend.background = element_blank(),
        plot.background = element_blank())

g = ggarrange(phar_sa, phar_sy, phar_si, pdae_sy, pdae_si, ncol = 3, nrow = 2)
ggsave(filename = "filename.pdf", plot = g, dpi = 300, width = 18, height = 12, units = "in")

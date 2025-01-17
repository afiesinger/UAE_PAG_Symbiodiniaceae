#!usr/bin/env Rscript

# manipulate tables output by SymPortal
#' ------------------------------------

# ------------------- THRESHOLD <1% SEQUENCE ABUNDANCE ---------------------

rm(list = ls())

# Load required library
library(dplyr)

# setwd
setwd("C:/Users/fiesi/OneDrive/PhD_UKon/Data_analysis/NCCT_MetaB/ITS2/SymPortal_results_fiesinger/post_med_seqs")

# Read the data
data <- read.delim('seqs.relative.abund_and_meta.txt')
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

##### count total sequences in clades

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
write.csv(final_data, "POSTMED_SEQS_Fiesinger_cleaned_1p_THRESHOLD.csv")

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

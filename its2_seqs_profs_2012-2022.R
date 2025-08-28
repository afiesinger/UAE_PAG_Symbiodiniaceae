#!usr/bin/env Rscript

# rarefaction curve of ITS2 post-MED sequences output by SymPortal
# plot ITS2 type profiles as barplots ordered by decreasing order of profile abundance
#' ------------------------------------------------------------------------------------

#' ----------------------------------- PACKAGES AND WORKING ENVIRONMENT ----------------------
rm(list = ls())

library(latex2exp)
library(ggplot2)
library(reshape2)
library(ggforce)
library(vegan)
library(scales)
library(gridExtra)
library(readxl)
library(dplyr)
library(tidyr)
library(ggpubr)

# set working directory
setwd("path/to/dir")

#' ------------- RAREFACTION CURVES OF ITS2 POST-MED SEQUENCES ------------------
# to assess whether sequencing effort was similar and comparable for both sampling campaigns (2012 and 2022) to infer biodiversity patterns

setwd("path/to/dir")
seqs = read.table("seqs_abs_abund_clean.txt", header = TRUE)

seq_mat = as.matrix(seqs[, -1])
seq_mat_howells = as.matrix(seqs[1:79, -1])
seq_mat_fiesinger = as.matrix(seqs[80:158, -1])

# set window
par(mfrow=c(1,3)) 
par(mar = c(5, 5, 5, 5)) 

# plot all rarecurves into one window
rarecurve(x = seq_mat_howells, 
          step = 10, 
          xlab = "No. of sequences", 
          ylab = "Diversity metric", 
          main = "2012 post-MED sequences", 
          label = FALSE,
          cex.lab = 2,   
          cex.main = 2,
          cex.axis = 1.5,
          las = 1)    

rarecurve(x = seq_mat_fiesinger, 
          step = 10, 
          xlab = "No. of sequences", 
          ylab = "", 
          main = "2022 post-MED sequences", 
          label = FALSE,
          cex.lab = 2,
          cex.main = 2,
          cex.axis = 1.5,
          las = 1)

rarecurve(x = seq_mat, 
          step = 10, 
          xlab = "No. of sequences", 
          ylab = "", 
          main = "2012 & 2022 post-MED sequences",
          label = FALSE,
          cex.lab = 2,
          cex.main = 2,
          cex.axis = 1.5,
          las = 1)

#' ------------------------- READ IN ITS2 TYPE PROFILE DATA -----------------------

its2Profs = as.data.frame(read_excel("filtered_data_SY_SI_Pdae_Howells_Fiesinger.xlsx", sheet = 2))
head(its2Profs)

genoshort = read_excel("filtered_data_SY_SI_Pdae_Howells_Fiesinger.xlsx", sheet = 3)

its2Profs = cbind(genoshort[,6], its2Profs[,-1]) # add all together into metadata dataframe
its2Profs = its2Profs[,-3] # remove species column because we only have Pdae
colnames(its2Profs)[1] = "Genotype" # relabel genotype column
head(its2Profs)

# Split dataset for plotting in separate panels
its2ProfsH = its2Profs %>% filter(Dataset == "Howells") %>% group_by(Genotype, Site)
its2ProfsHSY = its2ProfsH %>% filter(Site == "SY") %>% group_by(Genotype)
its2ProfsHSI = its2ProfsH %>% filter(Site == "SI") %>% group_by(Genotype)

its2ProfsF = its2Profs %>% filter(Dataset == "Fiesinger") %>% group_by(Genotype, Site)
its2ProfsFSY = its2ProfsF %>% filter(Site == "SY") %>% group_by(Genotype)
its2ProfsFSI = its2ProfsF %>% filter(Site == "SI") %>% group_by(Genotype)

# convert dataframes to long format for plotting
process_data <- function(data) {
  data_long <- data %>%
    pivot_longer(cols = -c(Genotype, Site, Dataset), names_to = "Profile", values_to = "Abundance")
}

its2Profs_long <- process_data(its2Profs)
p = sort(unique(its2Profs_long$Profile), decreasing = FALSE) # for displaying profiles in legend

its2ProfsHSY_long <- process_data(its2ProfsHSY)
its2ProfsHSI_long <- process_data(its2ProfsHSI)
its2ProfsFSY_long <- process_data(its2ProfsFSY)
its2ProfsFSI_long <- process_data(its2ProfsFSI)

#' ------------------------------------ SET COLOURS FOR PLOTS -----------------------------------
# same colour scheme as 2022 for the profiles that overlap, for better visualization

colors <- c(
  "A1" = "#FF6100",                                                                                                         
  "A1-A13a" = "#FA8F11",
  "A1-A13c-A1lz"=  "#FFA92E",
  "A1-A1bf-A1bw" = "#E1BB86",
  "A1-A1bw-A1bf-A1bx-A1eb"= "#fcaf58",
  "A1-A1bw-A1bx" = "#E7C703",                                                                                        
  "A1-A1eq-A1ep" =  "#FFFF7A",  
  "A1-A1nz-A1bw" = "#FFE92E",                                                                                                  
  "A1/A1bx-A1bw-A1bf"= "#F9F0B3", 
  
  "C15h-C15k" = "#4B0082", 
  "C15h-C15o-C15k" = "#b66ee8",   
  "C3-C3c-C3gulf-C3l-C3m-C3r" = "#771155", 
  "C3-C3d-C3i-C3gulf-C115c-C3ai-C3ak" = "#cc2936", 
  "C3-C3gulf-C3d-C3ai-C3ak-C115b-C3i-C115c" = "#920028", 
  "C3-C3gulf-C3d-C3i-C115c-C115b-C3ai-C3ak-C18a" = "#D18993",
  "C3-C3gulf-C3d-C3i-C115c-C115b-C3ak" = "#ffbfc0",
  "C3-C3gulf-C3d-C3i-C115c-C3ak-C3al" = "#ffe8f0",  
  "C3-C3gulf-C3ye-C3d-C3i-C115c-C3ak" = "#C39B9F", 
  "C3/C15h-C3gulf" = "#bf69a2",   
  "C3/C3c-C3gulf" = "#e57f68",   
  "C3/C3c-C3gulf-C3aq" = "#771122", 
  "C39-C39q-C1-C39b-C39p" = "#a13b42",  
  "C3yd/C3yc-C3-C3gulf" = "#d1c5c5",    
  
  "D1" = "#000080",   
  "D1/D2" = "#5f9ea0", 
  "D1/D4/D2-D4c-D1c-D1h" = "#b5d0ff",
  "D5-D4-D5a" = "#77AADD", 
  "D5-D5a-D4-D4a-D4b" = "#bcdbdb",  
  "D5-D5a-D5ai-D4-D5f-D5v" = "#81c784",
  "D5-D5a-D5f" = "#44aa77",        
  "D5-D5a-D5f-D4" = "#117744", 
  "D5-D5c-D4a-D5b-D4-D5i-D5a" = "#558b2f",
  "D5/D5a-D4-D4a-D2" = "#1ca62f",
  "D5a-D5-D5ah-D4-D4a" = "#3d6f3d", 
  "D5a-D5-D5ah-D4-D5d" = "#105410"
)
                      
#' ------------------- BARPLOT OF ITS2 TYPE PROFILES DECADAL COMPARISON --------------                

hsy = ggplot() +
  geom_bar(aes(y = Abundance, x = Genotype, fill = factor(Profile)), data = its2ProfsHSY_long, stat = "identity", position = "fill")  +  
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) + 
  labs(y = "Percentage", x = "", title = expression(paste(italic("N"), " = 47"))) + 
  scale_fill_manual(values = colors) + 
  #guides(fill = guide_legend(ncol = 1, title = expression(paste(italic("ITS2"), " type profile")))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0, size = 15, color = "black"), 
        axis.title.x = element_text(color = "black", size = 15),
        legend.position = "",
        plot.title = element_text(hjust = 0.5, size = 15), 
        axis.title.y = element_text(color = "black", size = 15),
        axis.text.y = element_text(color = "black", size = 15),
        plot.background = element_blank())

hsi = ggplot() +
  geom_bar(aes(y = Abundance, x = Genotype, fill = factor(Profile)), data = its2ProfsHSI_long, stat = "identity", position = "fill")  +  
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) + 
  labs(y = "", x = "", title = expression(paste(italic("N"), " = 32"))) + 
  scale_fill_manual(values = colors) + 
  #guides(fill = guide_legend(ncol = 1, title = expression(paste(italic("ITS2"), " type profile")))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0, size = 15, color = "black"), 
        axis.title.x = element_text(color = "black", size = 15),
        legend.position = "",
        plot.title = element_text(hjust = 0.5, size = 15), 
        axis.title.y = element_text(color = "black", size = 15),
        axis.text.y = element_text(color = "black", size = 15),
        plot.background = element_blank())

fsy = ggplot() +
  geom_bar(aes(y = Abundance, x = Genotype, fill = factor(Profile)), data = its2ProfsFSY_long, stat = "identity", position = "fill")  +  
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) + 
  labs(y = "Percentage", x = "Host genotype", title = expression(paste(italic("N"), " = 39"))) + 
  scale_fill_manual(values = colors) + 
  #guides(fill = guide_legend(ncol = 1, title = expression(paste(italic("ITS2"), " type profile")))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0, size = 15, color = "black"), 
        axis.title.x = element_text(color = "black", size = 15),
        legend.position = "",
        plot.title = element_text(hjust = 0.5, size = 15), 
        axis.title.y = element_text(color = "black", size = 15),
        axis.text.y = element_text(color = "black", size = 15),
        plot.background = element_blank())

fsi = ggplot() +
  geom_bar(aes(y = Abundance, x = Genotype, fill = factor(Profile)), data = its2ProfsFSI_long, stat = "identity", position = "fill")  +  
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) + 
  labs(y = "", x = "Host genotype", title = expression(paste(italic("N"), " = 40"))) + 
  scale_fill_manual(values = colors) + 
  guides(fill = guide_legend(ncol = 2, title = "ITS2 type profile")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0, size = 15, color = "black"),
        axis.title.x = element_text(color = "black", size = 15),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, size = 15), 
        axis.title.y = element_text(color = "black", size = 15),
        axis.text.y = element_text(color = "black", size = 15),
        plot.background = element_blank())

g = ggarrange(hsy, hsi, fsy, fsi, ncol = 2, nrow = 2)
ggsave(filename = "filename.pdf", plot = g, dpi = 300, width = 18, height = 12, units = "in")




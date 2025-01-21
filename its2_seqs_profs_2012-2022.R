#!usr/bin/env Rscript

#' ----------------------------------- PACKAGES AND WORKING ENVIRONMENT ----------------------
rm(list = ls())

if (!require("pacman")) install.packages("pacman")

pacman::p_load("dplyr", "edgeR", "ggplot2", "MCMC.OTU", "pairwiseAdonis", "rcartocolor", "RColorBrewer", "Redmonder", "reshape2", "vegan", "gridExtra", "scales", "readxl", "ape")

library(latex2exp)
library(ggcorrplot2)
library(ggplot2)
library(reshape2)
library(ggforce)
library(RColorBrewer)
library(vegan)
library(scales)
library(gridExtra)
library(readxl)
library(ape)
library(dplyr)
library(ggpubr)

pacman::p_load_gh("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

if (!require("edgeR")){BiocManager::install("edgeR", update = FALSE)
  library(edgeR)}

options("scipen" = 10)

#' ------------- RAREFACTION CURVES OF ITS2 POST-MED SEQUENCES ------------------
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

#' ------------------------- NORMALIZATION OF ITS2 TYPE PROFILES -----------------------

# set working directory
setwd("path/to/dir")

# data
its2Profs = as.data.frame(read_excel("filename.xlsx))
head(its2Profs)

# normalize the data
its2ProfsTransposed = t(its2Profs[, 5:length(its2Profs[1, ])])
its2ProfsList = DGEList(counts = its2ProfsTransposed)
head(its2ProfsList$samples)

its2ProfsNorm = calcNormFactors(its2ProfsList, method = "TMM")
head(its2ProfsNorm$samples)
its2TMM = t(cpm(its2ProfsNorm, normalized.lib.sizes = TRUE))
its2ProfsNorm = cbind(genoshort, its2Profs[,c(2:4)], its2TMM)
head(its2ProfsNorm)

# plotting the data
colOrder = order(colSums(its2ProfsNorm[5:length(its2ProfsNorm[1,])]), decreasing = TRUE) + 4

its2ProfsPerc = cbind(its2ProfsNorm[,c(1:4)],its2ProfsNorm[,c(colOrder)])

its2ProfsPerc$sum = apply(its2ProfsPerc[, c(5:length(its2ProfsPerc[1,]))], 1, function(x) {
  sum(x, na.rm = T)
})

its2ProfsPerc = cbind(its2ProfsPerc[, c(1:4)], (its2ProfsPerc[, c(5:(ncol(its2ProfsPerc)-1))] / its2ProfsPerc$sum))
head(its2ProfsPerc)

# sanity checking
apply(its2ProfsPerc[, c(5:(ncol(its2ProfsPerc)))], 1, function(x) {
  sum(x, na.rm = TRUE)
})

# prepare data for plotting
gssProf = otuStack(its2ProfsPerc, count.columns = c(5:length(its2ProfsPerc[1, ])), condition.columns = c(1:4))

# inspect the gssProf and remove the summ rows (otu level = summ)
gssProf = otuStack(its2ProfsPerc, count.columns = c(5:length(its2ProfsPerc[1, ])),
                   condition.columns = c(1:4))[1:5530,]

# sort by symbiont type
p = sort(levels(gssProf$otu), decreasing = FALSE)
p # find profiles and match the colours to the colour scheme of the 2022 barplot (see its2_seqs_profs_2022.R)

#' ------------------------------------ SET COLOURS FOR PLOTS -----------------------------------

colors <- c(
  "A1" = "#FF6100", 
  "A1-A13a" = "#FA8F11", 
  "A1-A13c-A1lz" = "#FFA92E", 
  "A1-A1bf-A1bw" = "#FFF27D",
  "A1-A1bw-A1bf-A1bx-A1eb" = "#fcaf58", 
  "A1-A1bw-A1bx" = "#E7C703", 
  "A1-A1eq-A1ep" = "#FFFF7A",
  "A1-A1nz-A1bw" = "#F8F691", 
  "A1/A1bx-A1bw-A1bf" = "#F6E9C4", 
  
  "C15h-C15k" = "#771122",
  "C15h-C15o-C15k" = "#B22222", 
  "C3-C3c-C3gulf-C3l-C3m-C3r" = "#B04456", 
  "C3-C3d-C3i-C3gulf-C115c-C3ai-C3ak" = "#AF6471", 
  "C3-C3gulf-C3d-C3ai-C3ak-C115b-C3i-C115c" = "#771155", 
  "C3-C3gulf-C3d-C3i-C115c-C115b-C3ai-C3ak-C18a" = "#AA4488", 
  "C3-C3gulf-C3d-C3i-C115c-C115b-C3ak" = "#9932CC", 
  "C3-C3gulf-C3d-C3i-C115c-C3ak-C3al" = "#8A2BE2", 
  "C3-C3gulf-C3ye-C3d-C3i-C115c-C3ak" = "#4B0082", 
  "C3/C15h-C3gulf" = "#774411", 
  "C3/C3c-C3gulf" = "#AF513D", 
  "C3/C3c-C3gulf-C3aq" = "#FF7F50", 
  "C39-C39q-C1-C39b-C39p" = "#E70000", 
  "C3yd/C3yc-C3-C3gulf" = "#FFA07A", 
  
  "D1" = "#114477",
  "D1/D2" = "#4169E1", 
  "D1/D4/D2-D4c-D1c-D1h" = "#4477AA", 
  "D5-D4-D5a" = "#77AADD", 
  "D5-D5a-D4-D4a-D4b" = "#117777", 
  "D5-D5a-D5ai-D4-D5f-D5v" = "#44AAAA", 
  "D5-D5a-D5f" = "#57D4D5",
  "D5-D5a-D5f-D4" = "#117744", 
  "D5-D5c-D4a-D5b-D4-D5i-D5a" = "#006400", 
  "D5/D5a-D4-D4a-D2" = "#2E8B57", 
  "D5a-D5-D5ah-D4-D4a" = "#44AA77", 
  "D5a-D5-D5ah-D4-D5d" = "#88CCAA")
)
                      
#' ------------------- BARPLOT OF ITS2 TYPE PROFILES DECADAL COMPARISON --------------                
# split the dataset
gssProfs = gssProf[,-5] # first remove species column
colnames(gssProfs) = c("count", "otu", "sample", "Site", "Dataset")

# make howells (2012) datasets split by site
gssProfH = gssProfs %>% filter(Dataset == "Howells") %>% group_by(sample, otu, count, Site)

## AA
gssProfHAA = gssProfH %>% filter(Site == "AA") %>% group_by(sample, otu, count)
gssProfHAA = gssProfHAA[,-c(4:5)] # remove site and dataset because redundant

# reorder dataset for plot - sorted by its2 type profile
gssProfHAA <- gssProfHAA %>%
  arrange(desc(count))

levels <- unique(gssProfHAA$sample)
gssProfHAA$sample <- factor(gssProfHAA$sample, levels = levels)

## SY
gssProfHSY = gssProfH %>% filter(Site == "SY") %>% group_by(sample, otu, count)
gssProfHSY = gssProfHSY[,-c(4:5)] # remove site and dataset because redundant

gssProfHSY <- gssProfHSY %>%
  arrange(desc(count))

levels <- unique(gssProfHSY$sample)
gssProfHSY$sample <- factor(gssProfHSY$sample, levels = levels)

# make fiesinger (2022) datasets split by site
gssProfF = gssProfs %>% filter(Dataset == "Fiesinger") %>% group_by(sample, otu, count, Site)

## AA
gssProfFAA = gssProfF %>% filter(Site == "AA") %>% group_by(sample, otu, count)
gssProfFAA = gssProfFAA[,-c(4:5)] # remove site and dataset because redundant

# reorder dataset for plot - sorted by its2 type profile
gssProfFAA <- gssProfFSI %>%
  arrange(desc(count))

levels <- unique(gssProfFAA$sample)
gssProfFAA$sample <- factor(gssProfFAA$sample, levels = levels)

## SY
gssProfFSY = gssProfF %>% filter(Site == "SY") %>% group_by(sample, otu, count)
gssProfFSY = gssProfFSY[,-c(4:5)] # remove site and dataset because redundant

gssProfFSY <- gssProfFSY %>%
  arrange(desc(count))

levels <- unique(gssProfFSY$sample)
gssProfFSY$sample <- factor(gssProfFSY$sample, levels = levels)

# plot in 4 panels; plot each ggplot first, then ggarrange

hsy = ggplot() +
  geom_bar(aes(y = count, x = sample, fill = factor(otu)), data = gssProfHSY, stat = "identity", position = "fill")  +  
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) + 
  labs(y = "Percentage", x = "", title = expression(paste(italic("N"), " = 47"))) + 
  scale_fill_manual(values = colors, breaks = p) + 
  #guides(fill = guide_legend(ncol = 1, title = expression(paste(italic("ITS2"), " type profile")))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0, size = 15, color = "black"), 
        axis.title.x = element_text(color = "black", size = 15),
        legend.position = "",
        plot.title = element_text(hjust = 0.5, size = 15), 
        axis.title.y = element_text(color = "black", size = 15),
        axis.text.y = element_text(color = "black", size = 15),
        plot.background = element_blank())

haa = ggplot() +
  geom_bar(aes(y = count, x = sample, fill = factor(otu)), data = gssProfHAA, stat = "identity", position = "fill")  +  
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) + 
  labs(y = "", x = "", title = expression(paste(italic("N"), " = 32"))) + 
  scale_fill_manual(values = colors, breaks = p) + 
  #guides(fill = guide_legend(ncol = 1, title = expression(paste(italic("ITS2"), " type profile")))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0, size = 15, color = "black"), 
        axis.title.x = element_text(color = "black", size = 15),
        legend.position = "",
        plot.title = element_text(hjust = 0.5, size = 15), 
        axis.title.y = element_text(color = "black", size = 15),
        axis.text.y = element_text(color = "black", size = 15),
        plot.background = element_blank())

fsy = ggplot() +
  geom_bar(aes(y = count, x = sample, fill = factor(otu)), data = gssProfFSY, stat = "identity", position = "fill")  +  
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) + 
  labs(y = "Percentage", x = "Host genotype", title = expression(paste(italic("N"), " = 39"))) + 
  scale_fill_manual(values = colors, breaks = p) + 
  #guides(fill = guide_legend(ncol = 1, title = expression(paste(italic("ITS2"), " type profile")))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0, size = 15, color = "black"), 
        axis.title.x = element_text(color = "black", size = 15),
        legend.position = "",
        plot.title = element_text(hjust = 0.5, size = 15), 
        axis.title.y = element_text(color = "black", size = 15),
        axis.text.y = element_text(color = "black", size = 15),
        plot.background = element_blank())

faa = ggplot() +
  geom_bar(aes(y = count, x = sample, fill = factor(otu)), data = gssProfFAA, stat = "identity", position = "fill")  +  
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) + 
  labs(y = "", x = "Host genotype", title = expression(paste(italic("N"), " = 40"))) + 
  scale_fill_manual(values = colors, breaks = p) + 
  #guides(fill = guide_legend(ncol = 1, title = expression(paste(italic("ITS2"), " type profile")))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0, size = 15, color = "black"),
        axis.title.x = element_text(color = "black", size = 15),
        legend.position = "",
        plot.title = element_text(hjust = 0.5, size = 15), 
        axis.title.y = element_text(color = "black", size = 15),
        axis.text.y = element_text(color = "black", size = 15),
        plot.background = element_blank())

ggarrange(hsy, haa, fsy, faa, ncol = 2, nrow = 2)




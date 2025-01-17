#!/usr/bin/env Rscript

#' ----------------------------------- PACKAGES AND WORKING ENVIRONMENT ----------------------
rm(list = ls())

if (!require("pacman")) install.packages("pacman")

pacman::p_load("dplyr", "edgeR", "ggplot2", "MCMC.OTU", "pairwiseAdonis", "rcartocolor", "RColorBrewer", "Redmonder", "reshape2", "vegan", "gridExtra", "scales", "readxl")

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
library(ggpubr)

pacman::p_load_gh("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

if (!require("edgeR")){BiocManager::install("edgeR", update = FALSE)
  library(edgeR)}

options("scipen" = 10)
#' ------------------------------------ READ ITS2 TYPE PROFILE DATA -----------------------------------

# set working directory
setwd("/path/to/dir")

# ITS2 type profile data from SymPortal (https://symportal.org)
# metadata file contains the species name and sampling location of all the samples

its2Profs = read.delim("396_20230721T100004_DBV_20230721T154120.profiles.absolute.abund_and_meta_clean.txt", header = TRUE, check.names = FALSE)
head(its2Profs)

its2MetaData = read_excel("its2_profiles_metadata.xlsx")
head(its2MetaData)

its2Profs = cbind(its2MetaData[,4], its2MetaData[,6], its2MetaData[,8], its2Profs[,c(2:length(its2Profs))])
colnames(its2Profs)[1] = "Genotype"
colnames(its2Profs)[2] = "Species"
colnames(its2Profs)[3] = "Site"
head(its2Profs)

# normalization of the data
its2ProfsTransposed = t(its2Profs[, 4:length(its2Profs[1, ])])
its2ProfsList = DGEList(counts = its2ProfsTransposed)
head(its2ProfsList$samples)

its2ProfsNorm = calcNormFactors(its2ProfsList, method = "TMM")
head(its2ProfsNorm$samples)
its2TMM = t(cpm(its2ProfsNorm, normalized.lib.sizes = TRUE))
its2ProfsNorm = cbind(its2Profs[,c(1:3)], its2TMM)
head(its2ProfsNorm)

# plotting the data
colOrder = order(colSums(its2ProfsNorm[4:length(its2ProfsNorm[1,])]), decreasing = TRUE) + 3

its2ProfsPerc = cbind(its2ProfsNorm[,c(1:3)],its2ProfsNorm[,c(colOrder)])

its2ProfsPerc$sum = apply(its2ProfsPerc[, c(4:length(its2ProfsPerc[1,]))], 1, function(x) {
  sum(x, na.rm = T)
})

its2ProfsPerc = cbind(its2ProfsPerc[, c(1:3)], (its2ProfsPerc[, c(4:(ncol(its2ProfsPerc)-1))] / its2ProfsPerc$sum))
head(its2ProfsPerc)

# sanity checking

apply(its2ProfsPerc[, c(4:(ncol(its2ProfsPerc)))], 1, function(x) {
sum(x, na.rm = TRUE)
})

# prepare data for plotting
gssProf = otuStack(its2ProfsPerc, count.columns = c(4:length(its2ProfsPerc[1, ])),
                   condition.columns = c(1:3))

# inspect the gssProf and remove the summ rows (otu level = summ)
gssProf = otuStack(its2ProfsPerc, count.columns = c(4:length(its2ProfsPerc[1, ])),
                   condition.columns = c(1:3))[1:8316,]

# split dataset by Site to plot different sites
gssProfSY = gssProf %>% filter(Site == "SY") %>% group_by(Genotype, Species, otu, count)
gssProfSI = gssProf %>% filter(Site == "SI") %>% group_by(Genotype, Species, otu, count)
gssProfSA = gssProf %>% filter(Site == "SA") %>% group_by(Genotype, Species, otu, count)

# split dataset by Species to plot the two species separately
gssProfPhar = gssProf %>% filter(Species == "Phar") %>% group_by(Site, Genotype, otu, count)
gssProfPdae = gssProf %>% filter(Species == "Pdae") %>% group_by(Site, Genotype, otu, count)

# split dataset by Species & then Site
gssProfPharSA = gssProfPhar %>% filter(Site == "SA") %>% group_by(Genotype, Species, otu, count)
gssProfPharSY = gssProfPhar %>% filter(Site == "SY") %>% group_by(Genotype, Species, otu, count)
gssProfPharSI = gssProfPhar %>% filter(Site == "SI") %>% group_by(Genotype, Species, otu, count)

gssProfPdaeSY = gssProfPdae %>% filter(Site == "SY") %>% group_by(Genotype, Species, otu, count)
gssProfPdaeSI = gssProfPdae %>% filter(Site == "SI") %>% group_by(Genotype, Species, otu, count)

# sort by symbiont type
p = sort(levels(gssProf$otu), decreasing = F)
p

#' ------------------------------------ SET COLOURS FOR PLOTS -----------------------------------
cols_2022 <- c(
  "A1" = "#FF6100", 
  "A1-A13a" = "#FA8F11", 
  "A1-A1bv" = "#F5C372", 
  "A1-A1bw" = "#F6E496", 
  "A1-A1bw-A1bf-A1bx-A1eb" = "#fcaf58", 
  "A1-A1bw-A1bx" = "#E7C703", 
  "A1-A1eq-A1ep" = "#FFFF7A", 
  
  "C15" = "#8B0000",
  "C15h" = "#DB7093",
  "C15h-C15o-C15k" = "#B22222",
  "C15/C116" = "#FA8072",
  "C3-C3cc-C3gulf-C3ye" = "#A0522D", 
  "C3-C3d-C3i-C3gulf-C115c-C3ai-C3ak" = "#AF6471", 
  "C3-C3gulf-C3cc" = "#EE82EE", 
  "C3-C3gulf-C3d-C3i-C115c-C115b-C3ai-C3ak" = "#771155", 
  "C3-C3gulf-C3d-C3i-C115c-C115b-C3ak" = "#9932CC", 
  "C3-C3gulf-C3d-C3i-C115c-C3ak-C3al" = "#8A2BE2", 
  "C3-C3gulf-C3ye" = "#E79EEB", 
  "C3/C15h-C3gulf" = "#774411", 
  "C3/C3c-C3gulf" = "#AF513D", 
  "C3/C3gulf" = "#E57F68",
  "C3/C3gulf-C115d" = "#E3BAAC",
  "C39-C39q-C1-C39b-C39p" = "#E70000",
  "C3yd" = "#E7CFC0",
  "C3yd/C3yc-C3-C3gulf" = "#FFA07A",  

  "D1" = "#114477", 
  "D1-D4-D4c-D17ao-D17ap-D17aq" = "#BCF6F5",
  "D1-D4-D4c-D17ap-D17aq-D17ao-D17ar" = "#6495ED",
  "D1-D4-D4c-D17d" = "#BCDBDB",
  "D1-D4-D4c-D17d-D1r-D17c-D17e" = "#40E0D0",  
  "D1-D4-D4c-D17d-D1r-D17c-D17e-D17j" = "#5F9EA0",
  "D1/D2" = "#4169E1", 
  "D1/D2-D4-D4c-D1r" = "#00BFFF", 
  "D1/D4-D4c-D1h" = "#000080", 
  "D1/D4-D4c-D1r" = "#008080", 
  "D1/D4/D2-D4c-D1c-D1h" = "#4477AA", 
  "D5-D5a-D4-D4a-D4b" = "#117777", 
  "D5-D5a-D5ai-D4-D5f-D5v" = "#44AAAA", 
  "D5-D5a-D5f" = "#57D4D5", 
  "D5-D5c-D4a-D5b-D4-D5i-D5a" = "#006400", 
  "D5/D5a-D4-D4a-D2" = "#2E8B57", 
  "D5a-D5-D5ah-D4" = "#44AA77"
)

# ------------------------------- PLOTTING ------------------------

# reorder dataset for plot - sorted by its2 type profile

gssProfPharSA <- gssProfPharSA %>%
  arrange(desc(count))

levels <- unique(gssProfPharSA$Genotype)
gssProfPharSA$Genotype <- factor(gssProfPharSA$Genotype, levels = levels)

gssProfPharSY <- gssProfPharSY %>%
  arrange(desc(count))

levels <- unique(gssProfPharSY$Genotype)
gssProfPharSY$Genotype <- factor(gssProfPharSY$Genotype, levels = levels)

gssProfPharSI <- gssProfPharSI %>%
  arrange(desc(count))

levels <- unique(gssProfPharSI$Genotype)
gssProfPharSI$Genotype <- factor(gssProfPharSI$Genotype, levels = levels)

gssProfPdaeSY <- gssProfPdaeSY %>%
  arrange(desc(count))

levels <- unique(gssProfPdaeSY$Genotype)
gssProfPdaeSY$Genotype <- factor(gssProfPdaeSY$Genotype, levels = levels)

gssProfPdaeSI <- gssProfPdaeSI %>%
  arrange(desc(count))

levels <- unique(gssProfPdaeSI$Genotype)
gssProfPdaeSI$Genotype <- factor(gssProfPdaeSI$Genotype, levels = levels)

# plot in 2 panels; plot each ggplot first, then ggarrange

phar_sa = ggplot() +
  geom_bar(aes(y = count, x = Genotype, fill = factor(otu)), data = gssProfPharSA, stat = "identity", position = "fill")  + 
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) + 
  labs(y = "Percentage", x = "Host genotype", title = expression(paste(italic("N"), " = 40"))) +
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

phar_sy = ggplot() +
  geom_bar(aes(y = count, x = Genotype, fill = factor(otu)), data = gssProfPharSY, stat = "identity", position = "fill")  + 
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
  geom_bar(aes(y = count, x = Genotype, fill = factor(otu)), data = gssProfPharSI, stat = "identity", position = "fill")  + 
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
  geom_bar(aes(y = count, x = Genotype, fill = factor(otu)), data = gssProfPdaeSY, stat = "identity", position = "fill")  + 
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) + 
  labs(y = "Percentage", x = "Host genotype", title = expression(paste(italic("N"), " = 39"))) +
  scale_fill_manual(values = cols_2022, breaks = p) +
  #guides(fill = guide_legend(ncol = 4, title = expression(paste(italic("ITS2"), " type profile")))) +
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
  geom_bar(aes(y = count, x = Genotype, fill = factor(otu)), data = gssProfPdaeSI, stat = "identity", position = "fill")  + 
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

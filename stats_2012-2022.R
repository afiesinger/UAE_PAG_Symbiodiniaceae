#!/usr/bin/env Rscript

# statistical testing of decadal community composition changes
# datasets used: May 2022 (Fiesinger et al.) and 2012-2014 (Howells et al. 2020)
#' -----------------------------------------------------------------------------

#' ------------------ PACKAGES -------------------

if (!require("pacman")) install.packages("pacman")

pacman::p_load("vegan", "readxl")
pacman::p_load_gh("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

library(vegan)
library(readxl)
library(pairwiseAdonis)

# We are using an Excel sheet of ITS2 type profiles from the shared dataset of Howells et al. (2020) and our dataset from 2022

#' ----------------------------- ANOSIM --------------------------

# POST-MED SEQUENCES #

# set working directory
setwd("/path/to/dir")

its2seqs = read_excel("filename.xlsx")
head(its2seqs)

data.seqs = its2seqs[, 4:length(its2seqs)]
datas_meta = its2seqs[, 1:3]

datas.dist <- vegdist(data.seqs)
attach(datas_meta)

# 2012 vs. 2022
datas.ano.d <- anosim(datas.dist, Dataset)
summary(datas.ano.d)
plot(datas.ano.d)

# Call:
#   anosim(x = datas.dist, grouping = Dataset) 
# Dissimilarity: bray 
# 
# ANOSIM statistic R: 0.2351 
# Significance: 0.001 
# 
# Permutation: free
# Number of permutations: 999
# 
# Upper quantiles of permutations (null model):
#   90%    95%  97.5%    99% 
#   0.0108 0.0188 0.0260 0.0365 
# 
# Dissimilarity ranks between and within classes:
#   0%  25%  50%   75%    100%    N
# Between   718 4510 6259  9571 12004.5 6241
# Fiesinger   1 1511 6455 11322 12004.5 3081
# Howells     7 1587 3280  7475 12004.5 3081

# The ANOSIM statistic "R" compares the mean of ranked dissimilarities between gr
# oups to the mean of ranked dissimilarities within groups. An R value close to "1.0" su
# ggests dissimilarity between groups while an R value close to "0" suggests an even di
# stribution of high and low ranks within and between groups" 
# In other words, the higher the R value, the more dissimilar your groups are in terms of 
# community composition.

# -> timepoints (dataset; 2012 vs. 2022) are significantly different from each other in terms of ITS2 type profile composition

# Site SY vs Site AA 
datas.ano.s <- anosim(datas.dist, Site)
summary(datas.ano.s)
plot(datas.ano.s)

# Call:
#   anosim(x = datas.dist, grouping = Site) 
# Dissimilarity: bray 
# 
# ANOSIM statistic R: 0.9998 
# Significance: 0.001 
# 
# Permutation: free
# Number of permutations: 999
# 
# Upper quantiles of permutations (null model):
#   90%    95%  97.5%    99% 
#   0.0132 0.0210 0.0283 0.0370 
# 
# Dissimilarity ranks between and within classes:
#   0%     25%    50%      75%    100%    N
# Between 6108 7759.75 9307.5 10855.25 12004.5 6192
# AA         1 2224.75 3084.5  3903.25  6393.0 2556
# SY         2 1336.50 3193.0  5038.00  6075.0 3655

# R = almost 1 -> Sites are very dissimilar!

# PROFILES #

setwd("/path/to/dir")

data.profs = read_excel("filename.xlsx")[, -c(1:4)]
datap_meta = read_excel("filename.xlsx")[, 1:4]

datap.dist <- vegdist(data.profs)
attach(datap_meta)

# 2012 vs. 2022
datap.ano.d <- anosim(datap.dist, Dataset)
summary(datap.ano.d)
plot(datap.ano.d)

# Call:
#   anosim(x = datap.dist, grouping = Dataset) 
# Dissimilarity: bray 
# 
# ANOSIM statistic R: 0.1116 
# Significance: 0.001 
# 
# Permutation: free
# Number of permutations: 999
# 
# Upper quantiles of permutations (null model):
#   90%     95%   97.5%     99% 
#   0.00824 0.01162 0.01714 0.02390 
# 
# Dissimilarity ranks between and within classes:
#   0%  25%  50%  75% 100%    N
# Between   62 7282 7282 7282 7282 6241
# Fiesinger  2 7282 7282 7282 7282 3081
# Howells    1 7282 7282 7282 7282 3081

# Site SY vs. Site AA
datap.ano.s <- anosim(datap.dist, Site)
summary(datap.ano.s)
plot(datap.ano.s)

# Call:
#   anosim(x = datap.dist, grouping = Site) 
# Dissimilarity: bray 
# 
# ANOSIM statistic R: 0.3357 
# Significance: 0.001 
# 
# Permutation: free
# Number of permutations: 999
# 
# Upper quantiles of permutations (null model):
#   90%    95%  97.5%    99% 
#   0.0086 0.0133 0.0174 0.0224 
# 
# Dissimilarity ranks between and within classes:
#   0%    25%  50%  75% 100%    N
# Between 2071 7282.0 7282 7282 7282 6192
# AA         7 7282.0 7282 7282 7282 2556
# SY         1 1407.5 7282 7282 7282 3655

#' -------------- Howells (2012) vs Fiesinger (2022) split by Site -------

# POST-MED SEQUENCES -- SITE AA #
setwd("/path/to/dir")

its2seqs = read_excel("filename.xlsx")
head(its2seqs)

seqs = subset(its2seqs, its2seqs$Site=="AA")
data.seqs = seqs[, 4:length(seqs)]
datas_meta = seqs[, 1:3]

datas.dist <- vegdist(data.seqs)
attach(datas_meta)

datas.ano.d <- anosim(datas.dist, Dataset)
summary(datas.ano.d)
plot(datas.ano.d)

# Call:
#   anosim(x = datas.dist, grouping = Dataset) 
# Dissimilarity: bray 
# 
# ANOSIM statistic R: 0.4532 
# Significance: 0.001 
# 
# Permutation: free
# Number of permutations: 999
# 
# Upper quantiles of permutations (null model):
#   90%    95%  97.5%    99% 
#   0.0290 0.0399 0.0480 0.0680 
# 
# Dissimilarity ranks between and within classes:
#   0%     25%    50%     75% 100%    N
# Between   310 1065.75 1652.5 2065.25 2553 1280
# Fiesinger   1  329.00 1304.5 1923.25 2556  780
# Howells     5  312.50  525.5  823.25 2043  496

# POST-MED SEQUENCES -- SITE SY #

seqs = subset(its2seqs, its2seqs$Site=="SY")
data.seqs = seqs[, 4:length(seqs)]
datas_meta = seqs[, 1:3]

datas.dist <- vegdist(data.seqs)
attach(datas_meta)

datas.ano.d <- anosim(datas.dist, Dataset)
summary(datas.ano.d)
plot(datas.ano.d)

# Call:
#   anosim(x = datas.dist, grouping = Dataset) 
# Dissimilarity: bray 
# 
# ANOSIM statistic R: 0.9974 
# Significance: 0.001 
# 
# Permutation: free
# Number of permutations: 999
# 
# Upper quantiles of permutations (null model):
#   90%    95%  97.5%    99% 
#   0.0246 0.0381 0.0473 0.0696 
# 
# Dissimilarity ranks between and within classes:
#   0%  25%  50%  75% 100%    N
# Between   1669 2281 2739 3197 3655 1833
# Fiesinger    1  293  667 1104 1991  741
# Howells      4  621 1090 1479 2049 1081

# ITS2 TYPE PROFILES -- SITE AA #

# set working directory
setwd("/path/to/dir")

its2norm = read_excel("filename.xlsx")
head(its2norm)

profs = subset(its2norm, its2norm$Site=="AA")
data.profs = profs[, 5:length(profs)]
datap_meta = profs[, 1:4]

datap.dist <- vegdist(data.profs)
attach(datap_meta)

datap.ano.d <- anosim(datap.dist, Dataset)
summary(datap.ano.d)
plot(datap.ano.d)

# Call:
#   anosim(x = datap.dist, grouping = Dataset) 
# Dissimilarity: bray 
# 
# ANOSIM statistic R: 0.3543 
# Significance: 0.001 
# 
# Permutation: free
# Number of permutations: 999
# 
# Upper quantiles of permutations (null model):
#   90%    95%  97.5%    99% 
#   0.0222 0.0326 0.0405 0.0473 
# 
# Dissimilarity ranks between and within classes:
#   0%     25%    50%    75%   100%    N
# Between   13 1536.50 1536.5 1536.5 1536.5 1280
# Fiesinger  1  470.75 1536.5 1536.5 1536.5  780
# Howells    2  304.75  462.5 1536.5 1536.5  496

# ITS2 TYPE PROFILES -- SITE SY #

profs = subset(its2norm, its2norm$Site=="SY")
data.profs = profs[, 5:length(profs)]
datap_meta = profs[, 1:4]

datap.dist <- vegdist(data.profs)
attach(datap_meta)

datap.ano.d <- anosim(datap.dist, Dataset)
summary(datap.ano.d)
plot(datap.ano.d)

# Call:
#   anosim(x = datap.dist, grouping = Dataset) 
# Dissimilarity: bray 
# 
# ANOSIM statistic R: 0.2167 
# Significance: 0.001 
# 
# Permutation: free
# Number of permutations: 999
# 
# Upper quantiles of permutations (null model):
#   90%    95%  97.5%    99% 
#   0.0232 0.0337 0.0428 0.0585 
# 
# Dissimilarity ranks between and within classes:
#   0%  25%    50%    75%   100%    N
# Between   744 1265 2627.5 2627.5 2627.5 1833
# Fiesinger   2  344 1568.0 2627.5 2627.5  741
# Howells     1  563 2627.5 2627.5 2627.5 1081

#' -------------- PERMDISP & PERMUTATION TESTS -----------------

# its2norm is df from above

# Site #
set.seed(694) # setting seed allows repetition of randomized processes
its2dispSite = betadisper(vegdist(its2norm[, c(5:ncol(its2norm))]), its2norm$Site)
anova(its2dispSite)

# Response: Distances
# Df Sum Sq  Mean Sq F value   Pr(>F)   
# Groups      1 0.1840 0.183950  7.1449 0.008317 **
# Residuals 156 4.0163 0.025746

# -- permutation test --
set.seed(694)
its2PermTest = permutest(its2dispSite, permutations = 9999, pairwise = T, model = "full")
its2PermTest

# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 9999
# 
# Response: Distances
# Df Sum Sq  Mean Sq      F N.Perm Pr(>F)   
# Groups      1 0.1840 0.183950 7.1449   9999 0.0077 **
# Residuals 156 4.0163 0.025746                        
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Pairwise comparisons:
#   (Observed p-value below diagonal, permuted p-value above diagonal)
# AA     SY
# AA           0.0074
# SY 0.0083174 

# Dataset (2012 vs. 2022) #
set.seed(694) 
its2dispDat = betadisper(vegdist(its2norm[, c(5:ncol(its2norm))]), its2norm$Dataset)
anova(its2dispDat)

# Response: Distances
# Df  Sum Sq   Mean Sq F value Pr(>F)
# Groups      1 0.00325 0.0032456  0.1775 0.6742
# Residuals 156 2.85327 0.0182902  

# --- permutation test ---

set.seed(694)
its2PermTest2 = permutest(its2dispDat, permutations = 9999, pairwise = T, model = "full")
its2PermTest2

# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 9999
# 
# Response: Distances
# Df  Sum Sq   Mean Sq      F N.Perm Pr(>F)
# Groups      1 0.00325 0.0032456 0.1775   9999 0.6664
# Residuals 156 2.85327 0.0182902                     
# 
# Pairwise comparisons:
#   (Observed p-value below diagonal, permuted p-value above diagonal)
# Fiesinger Howells
# Fiesinger            0.6667
# Howells     0.67415

# --> AA & SY are quite different; timepoints (Howells 2012 vs Fiesinger 2022) not so much


#' ------------------------- PERMANOVA -------------------------

set.seed(694)
its2Adonis = adonis2(its2norm[, c(5:ncol(its2norm))] ~ Site * Dataset, data = its2norm, permutations = 9999, method = "bray")
its2Adonis

# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 9999
# 
# adonis2(formula = its2norm[, c(5:ncol(its2norm))] ~ Site * Dataset, data = its2norm, permutations = 9999, method = "bray")
# Df SumOfSqs      R2      F Pr(>F)    
# Site           1   10.018 0.14652 33.291 0.0001 ***
# Dataset        1    5.967 0.08726 19.827 0.0001 ***
# Site:Dataset   1    6.048 0.08845 20.096 0.0001 ***
# Residual     154   46.344 0.67778                  
# Total        157   68.376 1.00000                  
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# We have both Site and Dataset as significant factors -- continue with PAIRWISE PERMANOVA

#' ----------------- PAIRWISE PERMANOVA -----------------------

# SITE #
set.seed(694)
its2PWAdonisSite = pairwise.adonis(its2norm[, c(5:ncol(its2norm))],
                                   factors = its2norm$Site,
                                   sim.method = "bray", p.adjust.m = "BH", perm = 9999)

its2PWAdonisSite

#   pairs   Df SumsOfSqs  F.Model        R2 p.value p.adjusted sig
# 1 AA vs SY  1  10.01835 26.78056 0.1465176  0.0001     0.0001 ***

# DATASET 2012 vs. 2022 #

set.seed(694)
its2PWAdonisDat = pairwise.adonis(its2norm[, c(5:ncol(its2norm))],
                                  factors = its2norm$Dataset,
                                  sim.method = "bray", p.adjust.m = "BH", perm = 9999)

its2PWAdonisDat

# pairs Df SumsOfSqs  F.Model         R2 p.value p.adjusted sig
# 1 Howells vs Fiesinger  1  5.828557 14.53693 0.08524215  0.0001     0.0001 ***

# Now we will do for each SITE a pairwise PERMANOVA to look for differences in DATASET (i.e. sampling year)

# SITE AA #
profs = subset(its2norm, its2norm$Site=="AA")
data.profs = profs[, 5:length(profs)]

set.seed(694)
pairwise.adonis(data.profs,
                factors = profs$Dataset,
                sim.method = "bray", p.adjust.m = "BH", perm = 9999)

# pairs Df SumsOfSqs  F.Model        R2 p.value p.adjusted sig
# 1 Howells vs Fiesinger  1  6.554576 20.51054 0.2266094  0.0001     0.0001 ***

# SITE SY #
profs = subset(its2norm, its2norm$Site=="SY")
data.profs = profs[, 5:length(profs)]

set.seed(694)
pairwise.adonis(data.profs,
                factors = profs$Dataset,
                sim.method = "bray", p.adjust.m = "BH", perm = 9999)

# pairs Df SumsOfSqs  F.Model        R2 p.value p.adjusted sig
# 1 Howells vs Fiesinger  1  5.459647 19.12955 0.1854905  0.0001     0.0001 ***

#' --------------------- SIMPER TEST ---------------------

# ITS2 TYPE PROFILES -- SITE SY #
profs = subset(its2ProfsNorm, its2ProfsNorm$Site=="SY")

its2Simper1 = simper(sqrt(profs[, c(5:ncol(profs))]), profs$Dataset)
summary(its2Simper1)

# Contrast: Howells_Fiesinger 
# 
# average        sd     ratio       ava       avb cumsum     p    
# C3-C3gulf-C3d-C3i-C115c-C115b-C3ai-C3ak-C18a   0.21845   0.23268   0.93880 558.70000 673.20000  0.359 0.397    
# C3-C3gulf-C3d-C3i-C115c-C3ak-C3al              0.14063   0.21599   0.65110 276.40000  50.30000  0.590 0.115    
# C3-C3gulf-C3ye-C3d-C3i-C115c-C3ak              0.07008   0.16634   0.42130   0.00000 145.80000  0.706 0.006 ** 
# C3-C3gulf-C3d-C3i-C115c-C115b-C3ak             0.05372   0.15099   0.35580  42.10000  77.80000  0.794 0.135    
# C3-C3gulf-C3d-C3ai-C3ak-C115b-C3i-C115c        0.04204   0.13815   0.30430  86.10000   0.00000  0.863 0.947    
# A1                                             0.02486   0.06269   0.39650  41.80000  22.40000  0.904 0.524    
# C3-C3d-C3i-C3gulf-C115c-C3ai-C3ak              0.01150   0.07100   0.16200   0.00000  25.50000  0.923 0.360    
# C3/C3c-C3gulf-C3aq                             0.01051   0.07144   0.14710  21.50000   0.00000  0.940 0.613    
# A1-A1bw-A1bf-A1bx-A1eb                         0.00885   0.03053   0.28980   6.40000  14.20000  0.955 0.275    
# A1-A1bw-A1bx                                   0.00776   0.02269   0.34220   0.00000  17.10000  0.967 0.014 *  
# A1-A1nz-A1bw                                   0.00670   0.02897   0.23140   0.00000  14.30000  0.978 0.148    
# A1-A1bf-A1bw                                   0.00291   0.01979   0.14720   6.70000   0.00000  0.983 0.643    
# A1/A1bx-A1bw-A1bf                              0.00254   0.01725   0.14720   5.80000   0.00000  0.987 0.619    
# D5-D5a-D5f-D4                                  0.00223   0.01512   0.14720   5.00000   0.00000  0.991 0.608    
# D5-D5a-D4-D4a-D4b                              0.00214   0.01452   0.14720   4.80000   0.00000  0.995 0.615    
# A1-A13c-A1lz                                   0.00130   0.00800   0.16190   0.00000   2.80000  0.997 0.323    
# D5-D4-D5a                                      0.00123   0.00838   0.14720   2.70000   0.00000  0.999 0.609    
# A1-A13a                                        0.00081   0.00501   0.16190   0.00000   1.70000  1.000 0.329    
# A1-A1eq-A1ep                                   0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.001 ***
# C3/C3c-C3gulf                                  0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.001 ***
# C15h-C15o-C15k                                 0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.001 ***
# C3yd/C3yc-C3-C3gulf                            0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.001 ***
# C39-C39q-C1-C39b-C39p                          0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.001 ***
# C3-C3c-C3gulf-C3l-C3m-C3r                      0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.001 ***
# C3/C15h-C3gulf                                 0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.001 ***
# C15h-C15k                                      0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.001 ***
# D5-D5a-D5ai-D4-D5f-D5v                         0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.001 ***
# D5/D5a-D4-D4a-D2                               0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.001 ***
# D5-D5c-D4a-D5b-D4-D5i-D5a                      0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.001 ***
# D1/D2                                          0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.001 ***
# D5a-D5-D5ah-D4-D4a                             0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.001 ***
# D5a-D5-D5ah-D4-D5d                             0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.001 ***
# D5-D5a-D5f                                     0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.001 ***
# D1                                             0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.001 ***
# D1/D4/D2-D4c-D1c-D1h                           0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.001 ***
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# Permutation: free
# Number of permutations: 999

# ITS2 TYPE PROFILES -- SITE AA #
profs2 = subset(its2ProfsNorm, its2ProfsNorm$Site=="AA")

its2Simper2 = simper(sqrt(profs2[, c(5:ncol(profs2))]), profs2$Dataset)
summary(its2Simper2)

# Contrast: Howells_Fiesinger 
# 
# average        sd     ratio       ava       avb cumsum     p    
# D5-D5a-D5f-D4                                  0.33290   0.22529   1.47770 694.10000   0.00000  0.341 0.001 ***
# D5-D5a-D5ai-D4-D5f-D5v                         0.23070   0.24307   0.94900   0.00000 478.90000  0.578 0.001 ***
# D5-D5c-D4a-D5b-D4-D5i-D5a                      0.12070   0.20941   0.57640 188.90000 101.10000  0.702 0.213    
# D5a-D5-D5ah-D4-D4a                             0.09720   0.19477   0.49930   0.00000 202.00000  0.801 0.246    
# D5/D5a-D4-D4a-D2                               0.07950   0.18061   0.44000 126.40000  50.50000  0.883 0.149    
# D5a-D5-D5ah-D4-D5d                             0.06030   0.15982   0.37750   0.00000 126.20000  0.945 0.755    
# C39-C39q-C1-C39b-C39p                          0.01670   0.07213   0.23190   0.00000  38.50000  0.962 0.927    
# D1/D2                                          0.01080   0.02632   0.40930  20.80000   3.90000  0.973 0.082 .  
# D1/D4/D2-D4c-D1c-D1h                           0.00980   0.06112   0.16000   0.00000  23.20000  0.983 0.569    
# D1                                             0.00290   0.01791   0.16000   0.00000   6.50000  0.986 0.571    
# A1                                             0.00280   0.01312   0.21270   0.00000   6.10000  0.989 0.779    
# A1-A1eq-A1ep                                   0.00190   0.00873   0.22250   0.00000   4.30000  0.991 0.787    
# D5-D5a-D5f                                     0.00180   0.01127   0.16000   0.00000   4.10000  0.993 0.550    
# C3/C3c-C3gulf                                  0.00160   0.00910   0.17940   3.50000   0.00000  0.994 0.383    
# C3-C3c-C3gulf-C3l-C3m-C3r                      0.00130   0.00741   0.17940   2.90000   0.00000  0.996 0.397    
# C15h-C15o-C15k                                 0.00120   0.00754   0.16000   0.00000   2.60000  0.997 0.554    
# A1-A1bw-A1bx                                   0.00110   0.00684   0.16000   0.00000   2.50000  0.998 0.550    
# C3/C15h-C3gulf                                 0.00080   0.00490   0.16000   0.00000   1.70000  0.999 0.527    
# C3yd/C3yc-C3-C3gulf                            0.00070   0.00465   0.16000   0.00000   1.60000  1.000 0.566    
# C15h-C15k                                      0.00040   0.00224   0.16000   0.00000   0.80000  1.000 0.514    
# A1-A1bw-A1bf-A1bx-A1eb                         0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.001 ***
# A1-A1bf-A1bw                                   0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.001 ***
# A1/A1bx-A1bw-A1bf                              0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.001 ***
# A1-A1nz-A1bw                                   0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.001 ***
# A1-A13c-A1lz                                   0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.001 ***
# A1-A13a                                        0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.001 ***
# C3-C3gulf-C3d-C3i-C115c-C115b-C3ai-C3ak-C18a   0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.001 ***
# C3-C3gulf-C3d-C3i-C115c-C3ak-C3al              0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.001 ***
# C3-C3gulf-C3ye-C3d-C3i-C115c-C3ak              0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.001 ***
# C3-C3gulf-C3d-C3i-C115c-C115b-C3ak             0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.001 ***
# C3-C3gulf-C3d-C3ai-C3ak-C115b-C3i-C115c        0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.001 ***
# C3/C3c-C3gulf-C3aq                             0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.001 ***
# C3-C3d-C3i-C3gulf-C115c-C3ai-C3ak              0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.001 ***
# D5-D5a-D4-D4a-D4b                              0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.001 ***
# D5-D4-D5a                                      0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.001 ***
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# Permutation: free
# Number of permutations: 999



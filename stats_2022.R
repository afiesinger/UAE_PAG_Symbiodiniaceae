#!/usr/bin/env Rscript

# statistical analyses on ITS2 post-med sequences and ITS2 type profiles output by SymPortal
# tests done on the complete dataset from May 2022 collected by Fiesinger et al.
#' -----------------------------------------------------------------------------------------

#' -------------------- PACKAGES & DATA ---------------------
library(vegan)

# We are utilizing the data.frame "its2ProfsNorm" (ITS2 type profiles normalized) from the script "its2_type_profs.R"

#' ------------------------ PERMDISP -----------------------
# Using betadisper() in vegan to look at multivariate homogeneity of dispersion (PERMDISP) between sites and species using Bray-Curtis dissimilarity. 

# PERMDISP is a multivariate extension of Levene's test (Anderson 2006) to examine whether groups differ 
# in plot-to-plot variability. In essence, this technique involves calculating the distance from each 
# data point to its group centroid and then testing whether those distances differ among the groups.

# difference between sites
set.seed(694) # setting seed allows repetition of randomized processes

its2dispSite = betadisper(vegdist(its2ProfsNorm[, c(4:ncol(its2ProfsNorm))]), its2ProfsNorm$Site)

anova(its2dispSite)

# Analysis of Variance Table
# 
# Response: Distances
#             Df  Sum     Sq        Mean Sq F value  Pr(>F)  
# Groups      2   0.5479 0.273941   4.5554          0.01165 *
# Residuals   195 11.7265 0.060136                  
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## significant effect of site (Pr>F = 0.012; p-value = 0.05) on betadiversity

# difference between species
set.seed(694) # setting seed allows repetition of randomized processes

its2dispSpec = betadisper(vegdist(its2ProfsNorm[, c(4:ncol(its2ProfsNorm))]), its2ProfsNorm$Species)

anova(its2dispSpec)

# Analysis of Variance Table
# 
# Response: Distances
#             Df  Sum   Sq        Mean Sq F value  Pr(>F)  
# Groups      1 0.05065 0.050646  3.3306            0.06953 .
# Residuals   196 2.98046 0.015206                  
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# non-significant effect of species on betadiversity

#' ----------------- Permutation test for pairwise comparisons ---------------- 
# Follow up with a permutation test to see where differences occur.

set.seed(694)

its2PermTest = permutest(its2dispSite, permutations = 9999, pairwise = T, model = "full", )
its2PermTest

# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 9999
# 
# Response: Distances
#             Df  Sum   Sq      Mean Sq F N.Perm Pr(>F)  
# Groups      2  0.5479 0.273941 4.5554   9999    0.0124 *
# Residuals   195 11.7265 0.060136                       
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Pairwise comparisons:
# (Observed p-value below diagonal, permuted p-value above diagonal)

#       SA    AA          SY
# SA          0.0134000   0.1284
# AA 0.0141629            0.0046
# SY 0.1205059 0.0044658   

# SA und AA significant (0.014); SY und AA significant (0.0045) -- significant differences between Persian Gulf and Gulf of Oman !

its2PermTest$statistic
# Overall (F)   SA-AA (t)   SA-SY (t)   AA-SY (t) 
# 4.555351   -2.490098   -1.564152    2.885190 

#' ------------------ PERMANOVA --------------- 

# Now let's see how different communities are from each other with PERMANOVA. We will utilize the adonis() function in vegan. 
# We will use Bray-Curtis similarity for our distance matrix and run a total 0f 9,99999 permutations, and test the effects of 
# Site, Species, and the interaction between Site and Species. 

set.seed(694)
its2Adonis = adonis2(its2ProfsNorm[, c(4:ncol(its2ProfsNorm))] ~ Species * Site, data = its2ProfsNorm, permutations = 9999, method = "bray")

its2Adonis
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 9999
# 
# adonis2(formula = its2ProfsNorm[, c(4:ncol(its2ProfsNorm))] ~ Species * Site, data = its2ProfsNorm, permutations = 9999, method = "bray")
#               Df SumOfSqs      R2      F Pr(>F)    
# Species        1    8.922 0.09964 29.646 0.0001 ***
# Site           2   14.570 0.16273 24.207 0.0001 ***
# Species:Site   1    7.962 0.08892 26.456 0.0001 ***
# Residual     193   58.082 0.64871                  
# Total        197   89.535 1.00000                  
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# now Species and Site have a significant effect on ITS2 type profiles; interaction between Species and Site also has a significant effect

#' ------------------- PERMANOVA FOR MULTIPLE COMPARISONS -------------------- 

# pairwise sites
set.seed(694)
its2PWAdonisSite = pairwise.adonis(its2ProfsNorm[, c(4:ncol(its2ProfsNorm))],
                               factors = its2ProfsNorm$Site,
                               sim.method = "bray", p.adjust.m = "BH", perm = 9999)

its2PWAdonisSite
# pairs Df SumsOfSqs  F.Model        R2 p.value p.adjusted sig
# 1 SY vs SA  1  7.080734 20.10952 0.1477452  0.0001     0.0001 ***
# 2 SY vs AA  1  8.172654 20.52398 0.1162674  0.0001     0.0001 ***
# 3 SA vs AA  1  7.341437 18.94598 0.1383464  0.0001     0.0001 ***

# pairwise species
set.seed(694)
its2PWAdonisSpecies = pairwise.adonis(its2ProfsNorm[, c(4:ncol(its2ProfsNorm))],
                                   factors = its2ProfsNorm$Species,
                                   sim.method = "bray", p.adjust.m = "BH", perm = 9999)

its2PWAdonisSpecies
# pairs Df SumsOfSqs  F.Model         R2 p.value p.adjusted sig
# 1 Pdae vs Phar  1  8.921594 21.69164 0.09964388  0.0001     0.0001 ***

# significant effect of all sites - also Persian Gulf samples SY vs. SA
# all sites are significantly different; both species are significantly different; species at sites are significantly different as well

#' ------------------ TEST ROBUSTNESS OF PERMANOVA RESULTS ----------------- 

# Set seed for reproducibility
set.seed(694)

# Calculate Bray-Curtis dissimilarity matrix
its2_dist <- vegdist(its2ProfsNorm[, 4:ncol(its2ProfsNorm)], method = "bray")

# Perform dbRDA using Species and Site as constraints
its2_dbrda <- dbrda(its2_dist ~ Species * Site, data = its2ProfsNorm)

# Summary of the dbRDA
summary(its2_dbrda)

# Check variance explained
R2_adj <- RsquareAdj(its2_dbrda)
R2_adj

# Plot the dbRDA results
plot(its2_dbrda, display = c("sites", "bp"), scaling = 2, main = "dbRDA of ITS2 Profiles")
ordihull(its2_dbrda, groups = its2ProfsNorm$Species, draw = "polygon", col = c("red", "blue"), label = TRUE)
ordihull(its2_dbrda, groups = its2ProfsNorm$Site, draw = "polygon", col = c("green", "purple"), label = TRUE)

# Perform ANOVA on dbRDA to test significance of constraints
anova(its2_dbrda, permutations = 9999)

# Permutation test for dbrda under reduced model
# Permutation: free
# Number of permutations: 9999
# 
# Model: dbrda(formula = its2_dist ~ Species * Site, data = its2ProfsNorm)
# Df SumOfSqs     F Pr(>F)    
# Model      4   31.988 26.82  1e-04 ***
# Residual 193   57.547                 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Variance partitioning: How much each factor contributes individually
anova(its2_dbrda, by = "margin", permutations = 9999)

# Permutation test for dbrda under reduced model
# Marginal effects of terms
# Permutation: free
# Number of permutations: 9999
# 
# Model: dbrda(formula = its2_dist ~ Species * Site, data = its2ProfsNorm)
# Df SumOfSqs      F Pr(>F)    
# Species:Site   1    8.218 27.563  1e-04 ***
# Residual     193   57.547                  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#' -------------- SIMPER TEST BETWEEN PAG AND GO SAMPLES AS WELL AS SPECIES -----------------

# Similarity percentage test (SIMPER) will show us the ITS2 type profiles that contribute the most to the dissimilarity between the different sites and species

# Species
head(its2ProfsNorm)

its2Simper1 = simper(sqrt(its2ProfsNorm[, c(4:ncol(its2ProfsNorm))]), its2ProfsNorm$Species)
summary(its2Simper1)

# Contrast: Pdae_Phar 

#                                           average        sd     ratio       ava       avb cumsum     p    
# C3-C3gulf-C3d-C3i-C115c-C115b-C3ai-C3ak   0.20943   0.26018   0.80490 592.00000   0.00000  0.210 0.001 ***
# D5-D5a-D5ai-D4-D5f-D5v                    0.11909   0.21243   0.56060 233.40000   0.00000  0.330 0.001 ***
# C3-C3gulf-C3cc                            0.10822   0.20067   0.53930   0.00000 229.05000  0.438 1.000    
# C3/C3c-C3gulf                             0.09405   0.17497   0.53750   0.00000 182.54000  0.533 1.000    
# D5a-D5-D5ah-D4                            0.07402   0.17628   0.41990 135.30000  16.36000  0.607 0.001 ***
# D1/D4-D4c-D1h                             0.07265   0.17114   0.42450   0.00000 154.95000  0.680 1.000    
# C3-C3cc-C3gulf-C3ye                       0.05426   0.15227   0.35630   0.00000 114.55000  0.735 1.000    
# A1                                        0.03793   0.12007   0.31590  11.30000  73.99000  0.773 1.000    
# D5-D5c-D4a-D5b-D4-D5i-D5a                 0.02552   0.11077   0.23040  49.30000   0.00000  0.798 0.001 ***
# C3/C3gulf                                 0.02534   0.10313   0.24570   0.00000  56.34000  0.824 1.000    
# C3-C3gulf-C3d-C3i-C115c-C3ak-C3al         0.02258   0.11984   0.18840 148.40000   0.00000  0.846 0.001 ***
# C3-C3gulf-C3ye                            0.01921   0.09384   0.20470   0.00000  40.88000  0.866 1.000    
# C3-C3gulf-C3d-C3i-C115c-C115b-C3ak        0.01914   0.09656   0.19820  37.00000   0.00000  0.885 0.001 ***
# C3yd/C3yc-C3-C3gulf                       0.01901   0.08736   0.21760   0.80000  41.54000  0.904 1.000    
# A1-A1bw-A1bx                              0.01509   0.04747   0.31790  59.60000   9.89000  0.919 0.087 .  
# C3/C3gulf-C115d                           0.01444   0.08124   0.17770  24.60000   4.43000  0.934 0.061 .  
# D5/D5a-D4-D4a-D2                          0.01256   0.07813   0.16080  24.60000   0.00000  0.946 0.001 ***
# C39-C39q-C1-C39b-C39p                     0.01188   0.08023   0.14810  80.30000   0.00000  0.958 0.001 ***
# A1-A1bw-A1bf-A1bx-A1eb                    0.00922   0.03222   0.28620   8.20000  13.17000  0.967 0.797    
# C3-C3d-C3i-C3gulf-C115c-C3ai-C3ak         0.00867   0.07656   0.11320  42.60000   0.00000  0.976 0.001 ***
# D1/D4/D2-D4c-D1c-D1h                      0.00504   0.04455   0.11300  11.30000   0.00000  0.981 0.001 ***
# C15/C116                                  0.00353   0.03916   0.09020   0.00000   8.05000  0.985 1.000    
# C15                                       0.00303   0.03351   0.09030   0.00000   7.69000  0.988 1.000    
# D1                                        0.00148   0.01306   0.11300   3.20000   0.00000  0.989 0.001 ***
# D5-D5a-D5f                                0.00147   0.01300   0.11320  12.30000   0.00000  0.991 0.001 ***
# D1-D4-D4c-D17ap-D17aq-D17ao-D17ar         0.00125   0.01382   0.09020   0.00000   3.01000  0.992 1.000    
# A1-A1eq-A1ep                              0.00101   0.00648   0.15600   2.10000   0.00000  0.993 0.001 ***
# C15h-C15o-C15k                            0.00096   0.00666   0.14400   1.30000   0.74000  0.994 0.160    
# D1-D4-D4c-D17d-D1r-D17c-D17e-D17j         0.00094   0.00606   0.15580   0.00000   2.23000  0.995 1.000    
# D1/D2                                     0.00093   0.00821   0.11300   1.90000   0.00000  0.996 0.001 ***
# D1/D2-D4-D4c-D1r                          0.00083   0.00922   0.09030   0.00000   2.11000  0.997 1.000    
# D1/D4-D4c-D1r                             0.00065   0.00724   0.09020   0.00000   1.49000  0.997 1.000    
# D1-D4-D4c-D17d-D1r-D17c-D17e              0.00065   0.00420   0.15400   0.00000   1.40000  0.998 1.000    
# A1-A13a                                   0.00042   0.00373   0.11300   0.80000   0.00000  0.998 0.001 ***
# C3/C15h-C3gulf                            0.00041   0.00363   0.11300   0.80000   0.00000  0.999 0.001 ***
# C3yd                                      0.00021   0.00230   0.09010   0.00000   0.46000  0.999 1.000    
# C15h                                      0.00018   0.00159   0.11300   0.40000   0.00000  0.999 0.001 ***
# D5-D5a-D4-D4a-D4b                         0.00017   0.00193   0.09010   0.00000   0.38000  0.999 1.000    
# D1-D4-D4c-D17ao-D17ap-D17aq               0.00016   0.00180   0.09010   0.00000   0.35000  1.000 1.000    
# D1-D4-D4c-D17d                            0.00016   0.00175   0.09010   0.00000   0.34000  1.000 1.000    
# A1-A1bw                                   0.00015   0.00169   0.09010   0.00000   0.33000  1.000 1.000    
# A1-A1bv                                   0.00015   0.00164   0.09000   0.00000   0.29000  1.000 1.000    
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# Permutation: free
# Number of permutations: 999

# Sites
head(its2ProfsNorm)

its2Simper2 = simper(sqrt(its2ProfsNorm[, c(4:ncol(its2ProfsNorm))]), its2ProfsNorm$Site)
summary(its2Simper2)

# Contrast: SY_SA 
# 
#                                           average        sd     ratio       ava       avb cumsum     p    
# C3-C3gulf-C3cc                            0.24113   0.25108   0.96040  99.90000 486.70000  0.256 0.001 ***
# C3-C3gulf-C3d-C3i-C115c-C115b-C3ai-C3ak   0.21056   0.25837   0.81500 599.60000   0.00000  0.480 0.001 ***
# C3-C3cc-C3gulf-C3ye                       0.16934   0.23813   0.71110   0.00000 340.80000  0.660 0.001 ***
# C3/C3c-C3gulf                             0.14728   0.20539   0.71710 258.60000  19.40000  0.817 0.001 ***
# C3-C3gulf-C3ye                            0.05994   0.16238   0.36910   0.00000 121.60000  0.881 0.003 ** 
# C3-C3gulf-C3d-C3i-C115c-C3ak-C3al         0.02276   0.12009   0.18960 150.30000   0.00000  0.905 0.290    
# C3-C3gulf-C3d-C3i-C115c-C115b-C3ak        0.01922   0.09613   0.19990  37.40000   0.00000  0.925 0.208    
# C3/C3gulf                                 0.01721   0.08729   0.19720  35.40000   0.00000  0.944 0.917    
# A1-A1bw-A1bx                              0.01715   0.05322   0.32230  67.90000   0.00000  0.962 0.185    
# C3/C3gulf-C115d                           0.01281   0.07901   0.16220  25.00000   0.00000  0.975 0.399    
# A1-A1bw-A1bf-A1bx-A1eb                    0.00895   0.02570   0.34820  15.20000   3.90000  0.985 0.591    
# C3-C3d-C3i-C3gulf-C115c-C3ai-C3ak         0.00876   0.07684   0.11390  43.20000   0.00000  0.994 0.262    
# A1                                        0.00357   0.03133   0.11390   8.40000   0.00000  0.998 1.000    
# D5-D5a-D4-D4a-D4b                         0.00054   0.00345   0.15690   0.00000   1.10000  0.999 0.161    
# A1-A1bw                                   0.00047   0.00302   0.15690   0.00000   1.00000  0.999 0.161    
# A1-A1bv                                   0.00046   0.00296   0.15650   0.00000   0.90000  1.000 0.139    
# A1-A13a                                   0.00042   0.00372   0.11390   0.90000   0.00000  1.000 0.220    
# A1-A1eq-A1ep                              0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.830    
# C3yd/C3yc-C3-C3gulf                       0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.999    
# C39-C39q-C1-C39b-C39p                     0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.966    
# C15h-C15o-C15k                            0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.843    
# C15h                                      0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.596    
# C3/C15h-C3gulf                            0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.635    
# C15                                       0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.602    
# C3yd                                      0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.611    
# C15/C116                                  0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.586    
# D1/D4-D4c-D1h                             0.00000   0.00000       NaN   0.00000   0.00000  1.000 1.000    
# D5-D5a-D5ai-D4-D5f-D5v                    0.00000   0.00000       NaN   0.00000   0.00000  1.000 1.000    
# D5a-D5-D5ah-D4                            0.00000   0.00000       NaN   0.00000   0.00000  1.000 1.000    
# D5-D5c-D4a-D5b-D4-D5i-D5a                 0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.978    
# D1-D4-D4c-D17d-D1r-D17c-D17e-D17j         0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.930    
# D1-D4-D4c-D17d-D1r-D17c-D17e              0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.939    
# D5/D5a-D4-D4a-D2                          0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.833    
# D1-D4-D4c-D17ao-D17ap-D17aq               0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.597    
# D1/D2-D4-D4c-D1r                          0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.602    
# D1/D2                                     0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.589    
# D1-D4-D4c-D17d                            0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.590    
# D5-D5a-D5f                                0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.612    
# D1/D4-D4c-D1r                             0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.586    
# D1/D4/D2-D4c-D1c-D1h                      0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.589    
# D1                                        0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.579    
# D1-D4-D4c-D17ap-D17aq-D17ao-D17ar         0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.558    
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Contrast: SY_AA 
# 
#                                           average        sd     ratio       ava       avb cumsum     p    
# C3-C3gulf-C3d-C3i-C115c-C115b-C3ai-C3ak   0.20447   0.25301   0.80810 599.60000   0.00000  0.206 0.001 ***
# C3/C3c-C3gulf                             0.14106   0.19978   0.70610 258.60000   9.70000  0.348 0.001 ***
# D1/D4-D4c-D1h                             0.11332   0.20876   0.54280   0.00000 230.48000  0.462 0.001 ***
# D5-D5a-D5ai-D4-D5f-D5v                    0.11289   0.20792   0.54300   0.00000 230.48000  0.576 0.001 ***
# D5a-D5-D5ah-D4                            0.07709   0.17940   0.42970   0.00000 157.93000  0.654 0.001 ***
# A1                                        0.05652   0.14854   0.38050   8.40000 113.04000  0.711 0.003 ** 
# C3-C3gulf-C3cc                            0.04960   0.14753   0.33620  99.90000   0.00000  0.760 1.000    
# C3/C3gulf                                 0.03823   0.12691   0.30120  35.40000  49.26000  0.799 0.160    
# C3yd/C3yc-C3-C3gulf                       0.02942   0.11028   0.26680   0.00000  62.59000  0.829 0.053 .  
# D5-D5c-D4a-D5b-D4-D5i-D5a                 0.02419   0.10776   0.22450   0.00000  48.68000  0.853 0.136    
# C3-C3gulf-C3d-C3i-C115c-C3ak-C3al         0.02231   0.11834   0.18850 150.30000   0.00000  0.875 0.045 *  
# C3-C3gulf-C3d-C3i-C115c-C115b-C3ak        0.01860   0.09348   0.19900  37.40000   0.00000  0.894 0.040 *  
# A1-A1bw-A1bx                              0.01748   0.05237   0.33380  67.90000   7.35000  0.912 0.002 ** 
# C3/C3gulf-C115d                           0.01500   0.07997   0.18760  25.00000   6.59000  0.927 0.169    
# D5/D5a-D4-D4a-D2                          0.01191   0.07595   0.15680   0.00000  24.32000  0.939 0.579    
# A1-A1bw-A1bf-A1bx-A1eb                    0.01141   0.03728   0.30610  15.20000  11.01000  0.950 0.192    
# C39-C39q-C1-C39b-C39p                     0.01140   0.07824   0.14580   0.00000  79.33000  0.962 0.454    
# C3-C3d-C3i-C3gulf-C115c-C3ai-C3ak         0.00861   0.07572   0.11370  43.20000   0.00000  0.971 0.190    
# C15/C116                                  0.00549   0.04977   0.11030   0.00000  11.97000  0.976 0.697    
# D1/D4/D2-D4c-D1c-D1h                      0.00478   0.04325   0.11040   0.00000  11.15000  0.981 0.697    
# C15                                       0.00469   0.04240   0.11050   0.00000  11.43000  0.986 0.722    
# D1-D4-D4c-D17ap-D17aq-D17ao-D17ar         0.00194   0.01753   0.11040   0.00000   4.48000  0.988 0.713    
# D1-D4-D4c-D17d-D1r-D17c-D17e-D17j         0.00147   0.00767   0.19140   0.00000   3.32000  0.989 0.400    
# D5-D5a-D5f                                0.00142   0.01268   0.11190   0.00000  12.12000  0.991 0.787    
# D1                                        0.00140   0.01269   0.11040   0.00000   3.14000  0.992 0.704    
# D1/D2-D4-D4c-D1r                          0.00129   0.01166   0.11050   0.00000   3.15000  0.993 0.722    
# C15h-C15o-C15k                            0.00112   0.00718   0.15650   0.00000   2.37000  0.994 0.577    
# D1/D4-D4c-D1r                             0.00101   0.00920   0.11030   0.00000   2.21000  0.995 0.697    
# D1-D4-D4c-D17d-D1r-D17c-D17e              0.00101   0.00533   0.18910   0.00000   2.09000  0.996 0.401    
# A1-A1eq-A1ep                              0.00096   0.00630   0.15220   0.00000   2.07000  0.997 0.601    
# D1/D2                                     0.00088   0.00797   0.11030   0.00000   1.90000  0.998 0.690    
# A1-A13a                                   0.00041   0.00362   0.11340   0.90000   0.00000  0.999 0.302    
# C3/C15h-C3gulf                            0.00039   0.00353   0.11020   0.00000   0.81000  0.999 0.697    
# C3yd                                      0.00032   0.00293   0.11030   0.00000   0.69000  0.999 0.697    
# D1-D4-D4c-D17ao-D17ap-D17aq               0.00025   0.00229   0.11020   0.00000   0.52000  1.000 0.688    
# D1-D4-D4c-D17d                            0.00025   0.00223   0.11020   0.00000   0.51000  1.000 0.675    
# C15h                                      0.00017   0.00154   0.11020   0.00000   0.36000  1.000 0.692    
# A1-A1bw                                   0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.792    
# A1-A1bv                                   0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.790    
# C3-C3cc-C3gulf-C3ye                       0.00000   0.00000       NaN   0.00000   0.00000  1.000 1.000    
# C3-C3gulf-C3ye                            0.00000   0.00000       NaN   0.00000   0.00000  1.000 1.000    
# D5-D5a-D4-D4a-D4b                         0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.792    
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Contrast: SA_AA 
# 
#                                           average        sd     ratio       ava       avb cumsum     p    
# C3-C3gulf-C3cc                            0.24088   0.24327   0.99020 486.70000   0.00000  0.241 0.001 ***
# C3-C3cc-C3gulf-C3ye                       0.16927   0.23242   0.72830 340.80000   0.00000  0.410 0.001 ***
# D1/D4-D4c-D1h                             0.11705   0.21004   0.55730   0.00000 230.48000  0.527 0.016 *  
# D5-D5a-D5ai-D4-D5f-D5v                    0.11661   0.20922   0.55740   0.00000 230.48000  0.644 0.008 ** 
# D5a-D5-D5ah-D4                            0.07963   0.18095   0.44010   0.00000 157.93000  0.724 0.065 .  
# C3-C3gulf-C3ye                            0.05993   0.15949   0.37580 121.60000   0.00000  0.784 0.001 ***
# A1                                        0.05568   0.14888   0.37400   0.00000 113.04000  0.839 0.172    
# C3yd/C3yc-C3-C3gulf                       0.03040   0.11155   0.27250   0.00000  62.59000  0.870 0.217    
# D5-D5c-D4a-D5b-D4-D5i-D5a                 0.02499   0.10895   0.22930   0.00000  48.68000  0.895 0.198    
# C3/C3gulf                                 0.02399   0.10267   0.23360   0.00000  49.26000  0.919 0.744    
# C3/C3c-C3gulf                             0.01573   0.08071   0.19490  19.40000   9.70000  0.935 1.000    
# D5/D5a-D4-D4a-D2                          0.01230   0.07686   0.16010   0.00000  24.32000  0.947 0.228    
# C39-C39q-C1-C39b-C39p                     0.01170   0.07960   0.14700   0.00000  79.33000  0.959 0.305    
# A1-A1bw-A1bf-A1bx-A1eb                    0.00644   0.03103   0.20740   3.90000  11.01000  0.965 0.865    
# C15/C116                                  0.00567   0.05045   0.11250   0.00000  11.97000  0.971 0.215    
# D1/D4/D2-D4c-D1c-D1h                      0.00494   0.04390   0.11250   0.00000  11.15000  0.976 0.224    
# C15                                       0.00484   0.04307   0.11250   0.00000  11.43000  0.981 0.257    
# C3/C3gulf-C115d                           0.00284   0.02524   0.11250   0.00000   6.59000  0.983 0.876    
# D1-D4-D4c-D17ap-D17aq-D17ao-D17ar         0.00200   0.01779   0.11250   0.00000   4.48000  0.985 0.250    
# D1-D4-D4c-D17d-D1r-D17c-D17e-D17j         0.00152   0.00777   0.19510   0.00000   3.32000  0.987 0.301    
# D5-D5a-D5f                                0.00145   0.01290   0.11250   0.00000  12.12000  0.988 0.200    
# D1                                        0.00145   0.01287   0.11250   0.00000   3.14000  0.990 0.234    
# D1/D2-D4-D4c-D1r                          0.00133   0.01185   0.11250   0.00000   3.15000  0.991 0.257    
# C15h-C15o-C15k                            0.00116   0.00727   0.15970   0.00000   2.37000  0.992 0.313    
# D1/D4-D4c-D1r                             0.00105   0.00932   0.11250   0.00000   2.21000  0.993 0.215    
# D1-D4-D4c-D17d-D1r-D17c-D17e              0.00104   0.00540   0.19310   0.00000   2.09000  0.994 0.281    
# A1-A1eq-A1ep                              0.00099   0.00638   0.15520   0.00000   2.07000  0.995 0.305    
# D1/D2                                     0.00091   0.00808   0.11250   0.00000   1.90000  0.996 0.213    
# A1-A1bw-A1bx                              0.00088   0.00783   0.11250   0.00000   7.35000  0.997 0.999    
# D5-D5a-D4-D4a-D4b                         0.00054   0.00340   0.15930   1.10000   0.00000  0.998 0.152    
# A1-A1bw                                   0.00047   0.00297   0.15930   1.00000   0.00000  0.998 0.152    
# A1-A1bv                                   0.00046   0.00290   0.15930   0.90000   0.00000  0.999 0.151    
# C3/C15h-C3gulf                            0.00040   0.00357   0.11250   0.00000   0.81000  0.999 0.212    
# C3yd                                      0.00033   0.00297   0.11250   0.00000   0.69000  0.999 0.232    
# D1-D4-D4c-D17ao-D17ap-D17aq               0.00026   0.00232   0.11250   0.00000   0.52000  1.000 0.236    
# D1-D4-D4c-D17d                            0.00025   0.00226   0.11250   0.00000   0.51000  1.000 0.245    
# C15h                                      0.00018   0.00156   0.11250   0.00000   0.36000  1.000 0.235    
# A1-A13a                                   0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.600    
# C3-C3gulf-C3d-C3i-C115c-C115b-C3ai-C3ak   0.00000   0.00000       NaN   0.00000   0.00000  1.000 1.000    
# C3-C3gulf-C3d-C3i-C115c-C3ak-C3al         0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.951    
# C3-C3gulf-C3d-C3i-C115c-C115b-C3ak        0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.951    
# C3-C3d-C3i-C3gulf-C115c-C3ai-C3ak         0.00000   0.00000       NaN   0.00000   0.00000  1.000 0.627    
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# Permutation: free
# Number of permutations: 999

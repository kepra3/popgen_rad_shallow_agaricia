# Packages
library(vegan)
library(tidyverse)
library(ape)

# Setwd
setwd("~/Dropbox/agaricia_project_2019/shalo_ag/gen_project/3 - Spatial")
all <- read.csv("all_annotations_HORIZ_parallel_XYZadjusted.txt")

# AA1 ####
# Import datasets
#all <- read.csv("oriented_scaled_parallel_xyz_all.csv")
a <- read.delim("3b - Within & between plots/spatial-agaricia/results/AA1/AA1_all.a.txt", header = FALSE)
taxa.pop <- read.csv("3b - Within & between plots/spatial-agaricia/ac_3b_nc_20_AA1_all_.csv", row.names = 1)

## PCOA of genetic distance ####
row1 <- NA
a <- rbind(row1, a)
Vx <- NA
a <- cbind(a, Vx)

a <- as.matrix(a)
diag(a) <- 0
a.dist <- as.dist(a, upper = FALSE)

pcoa.a <- pcoa(a.dist) # cmdscale(a.dist) same thing
plot(pcoa.a$vectors[,1], pcoa.a$vectors[,2])
gen.dist <- as.data.frame(pcoa.a$vectors[,1])
gen.dist[,2] <- pcoa.a$vectors[,2]
colnames(gen.dist) <- c("PCoA1", "PCoA2")

## xyz dataset ####
taxa.geodist <- all[all$Individual %in% taxa.pop$Individual,]
taxa.geodist <- taxa.geodist[order(taxa.geodist$Individual),]

pca <- prcomp(taxa.geodist[,3:4])
PC1xy <- pca$x[,1]
PC2xy <- pca$x[,2]
plot(PC1xy, PC2xy)

plot(gen.dist$PCoA1, PC1xy)
plot(gen.dist$PCoA1, PC2xy)
plot(gen.dist$PCoA1, taxa.geodist$z)

# RDA full ####
rda.f <- rda(gen.dist ~ PC1xy * PC2xy * taxa.geodist$z)
rda.f
RsquareAdj(rda.f)
#$r.squared
#[1] 0.2948992
#$adj.r.squared
#[1] -0.1538013

summary(eigenvals(rda.f,model="constrained"))
#Importance of components:
#RDA1      RDA2
#Eigenvalue            0.00118 0.0001195
#Proportion Explained  0.90807 0.0919309
#Cumulative Proportion 0.90807 1.0000000

summary(eigenvals(rda.f,model="unconstrained"))
#Importance of components:
#PC1       PC2
#Eigenvalue            0.002275 0.0008317
#Proportion Explained  0.732322 0.2676782
#Cumulative Proportion 0.732322 1.0000000

summary(rda.f)
#Partitioning of variance:
#  Inertia Proportion
#Total         0.004406     1.0000
#Constrained   0.001299     0.2949
#Unconstrained 0.003107     0.7051

#                               RDA1     RDA2 PC1 PC2
#PC1xy                       0.94080 -0.12973   0   0
#PC2xy                       0.17076  0.23044   0   0
#taxa.geodist$z              0.07116  0.31512   0   0
#PC1xy:PC2xy                -0.18799 -0.16235   0   0
#PC1xy:taxa.geodist$z       -0.94951  0.11271   0   0
#PC2xy:taxa.geodist$z       -0.30385 -0.06808   0   0
#PC1xy:PC2xy:taxa.geodist$z  0.31779 -0.05682   0   0


#Eigenvalues, and their contribution to the variance 

#Importance of components:
#  RDA1      RDA2      PC1       PC2
#Eigenvalue            0.00108 3.023e-05 0.002368 0.0009276
#Proportion Explained  0.24521 6.861e-03 0.537419 0.2105124
#Cumulative Proportion 0.24521 2.521e-01 0.789488 1.0000000

step(rda.f)
# -AIC scores v. confusing

signif.full.c <- anova.cca(rda.f)
signif.full.c # not significant
#Permutation test for rda under reduced model
#Permutation: free
#Number of permutations: 999
#Model: rda(formula = gen.dist ~ PC1xy + PC2xy + taxa.geodist$z)
#Df  Variance      F Pr(>F)
#Model     3 0.0011107 1.6851  0.172
#Residual 15 0.0032957                

signif.axis.c <- anova.cca(rda.f, by="axis")
signif.axis.c

signif.terms.c <- anova.cca(rda.f, by="terms")
signif.terms.c
# only PC1xy

## RDA 1 ####
rda1 <- rda(gen.dist ~ PC1xy + PC2xy + taxa.geodist$z)
rda1
RsquareAdj(rda1)
#$r.squared
#[1] 0.2520686
#$adj.r.squared
#[1] 0.1024823

summary(eigenvals(rda1,model="constrained"))
#Importance of components:
#RDA1      RDA2
#Eigenvalue            0.00108 3.023e-05
#Proportion Explained  0.97278 2.722e-02
#Cumulative Proportion 0.97278 1.000e+00

summary(eigenvals(rda1,model="unconstrained"))
#Importance of components:
#  PC1       PC2
#Eigenvalue            0.002368 0.0009276
#Proportion Explained  0.718541 0.2814595
#Cumulative Proportion 0.718541 1.0000000

summary(rda1)
#Partitioning of variance:
#  Inertia Proportion
#Total         0.004406     1.0000
#Constrained   0.001111     0.2521
#Unconstrained 0.003296     0.7479

#Eigenvalues, and their contribution to the variance 

#Importance of components:
#  RDA1      RDA2      PC1       PC2
#Eigenvalue            0.00108 3.023e-05 0.002368 0.0009276
#Proportion Explained  0.24521 6.861e-03 0.537419 0.2105124
#Cumulative Proportion 0.24521 2.521e-01 0.789488 1.0000000

step(rda1)
# -AIC scores v. confusing

signif.full.c <- anova.cca(rda1)
signif.full.c # not significant
#Permutation test for rda under reduced model
#Permutation: free
#Number of permutations: 999
#Model: rda(formula = gen.dist ~ PC1xy + PC2xy + taxa.geodist$z)
#Df  Variance      F Pr(>F)
#Model     3 0.0011107 1.6851  0.172
#Residual 15 0.0032957                

signif.axis.c <- anova.cca(rda1, by="axis")
signif.axis.c
# PC2 is doing nothing basically. Makes sense with it explaining 0 variation in the PCA.

signif.terms.c <- anova.cca(rda1, by="terms")
signif.terms.c

## RDA 2 ####
rda2 <- rda(gen.dist ~ PC1xy + taxa.geodist$z)
rda2
RsquareAdj(rda2)
#$r.squared
#[1] 0.2472196
#$adj.r.squared
#[1] 0.1531221

summary(eigenvals(rda2,model="constrained"))
#Importance of components:
#RDA1      RDA2
#Eigenvalue            0.001077 1.235e-05
#Proportion Explained  0.988663 1.134e-02
#Cumulative Proportion 0.988663 1.000e+00
# only RDA1 again. 

summary(eigenvals(rda2,model="unconstrained"))
#Importance of components:
#PC1       PC2
#Eigenvalue            0.002385 0.0009319
#Proportion Explained  0.719068 0.2809324
#Cumulative Proportion 0.719068 1.0000000

summary(rda2)
#Partitioning of variance:
#Inertia Proportion
#Total         0.004406     1.0000
#Constrained   0.001089     0.2472
#Unconstrained 0.003317     0.7528

#Eigenvalues, and their contribution to the variance 

#Importance of components:
#  RDA1      RDA2      PC1       PC2
#Eigenvalue            0.001077 1.235e-05 0.002385 0.0009319
#Proportion Explained  0.244417 2.803e-03 0.541300 0.2114804
#Cumulative Proportion 0.244417 2.472e-01 0.788520 1.0000000
#RDA1    RDA2 PC1 PC2
#PC1xy          -0.98554  0.1695   0   0
#taxa.geodist$z -0.07179 -0.9974   0   0
# biplot interesting PC1 and PC2 are 0?

step(rda2)
# -AIC scores v. confusing again

signif.full.c <- anova.cca(rda2)
signif.full.c # almost significant
#Model: rda(formula = gen.dist ~ PC1xy + taxa.geodist$z)
#Df  Variance      F Pr(>F)  
#Model     2 0.0010893 2.6273  0.076 .
#Residual 16 0.0033170               

signif.axis.c <- anova.cca(rda2, by="axis")
signif.axis.c
# RDA2 not explaining much, take out depth then?
signif.terms.c <- anova.cca(rda2, by="terms")
signif.terms.c

## RDA 3 ####
rda3 <- rda(gen.dist ~ PC1xy)
rda3
RsquareAdj(rda3)
#$r.squared
#[1] 0.2374778
#$adj.r.squared
#[1] 0.1926236

summary(eigenvals(rda3,model="constrained"))
summary(eigenvals(rda3,model="unconstrained"))
#Importance of components:
#PC1       PC2
#Eigenvalue            0.002399 0.0009606
#Proportion Explained  0.714101 0.2858993
#Cumulative Proportion 0.714101 1.0000000

summary(rda3)
#Partitioning of variance:
#Inertia Proportion
#Total         0.004406     1.0000
#Constrained   0.001046     0.2375
#Unconstrained 0.003360     0.7625

#Eigenvalues, and their contribution to the variance 

#Importance of components:
#  RDA1      PC1       PC2
#Eigenvalue            0.001046 0.002399 0.0009606
#Proportion Explained  0.237478 0.544518 0.2180046
#Cumulative Proportion 0.237478 0.781995 1.0000000
#Cumulative Proportion 0.244417 2.472e-01 0.788520 1.0000000
#Biplot scores for constraining variables
#RDA1 PC1 PC2
#PC1xy    1   0   0
# biplot interesting PC1 and PC2 are 0?
# PC1 > RDA1


signif.full.c <- anova.cca(rda3)
signif.full.c # almost significant
#Model: rda(formula = gen.dist ~ PC1xy)
#Df  Variance      F Pr(>F)   
#Model     1 0.0010464 5.2944  0.005 **
#  Residual 17 0.0033599                 

signif.axis.c <- anova.cca(rda3, by="axis")
signif.axis.c
#Model: rda(formula = gen.dist ~ PC1xy)
#Df  Variance      F Pr(>F)   
#RDA1      1 0.0010464 5.2944   0.01 **
#  Residual 17 0.0033599  

signif.terms.c <- anova.cca(rda3, by="terms")
signif.terms.c

# AA2 ####
## Import datasets
a <- read.delim("3b - Within & between plots/spatial-agaricia/results/AA2/AA2_all.a.txt", header = FALSE)
taxa.pop <- read.csv("3b - Within & between plots/spatial-agaricia/ac_3b_nc_20_AA2_all_.csv", row.names = 1)

## PCOA of genetic distance ####
row1 <- NA
a <- rbind(row1, a)
Vx <- NA
a <- cbind(a, Vx)

a <- as.matrix(a)
diag(a) <- 0
a.dist <- as.dist(a, upper = FALSE)

pcoa.a <- pcoa(a.dist) # cmdscale(a.dist) same thing
plot(pcoa.a$vectors[,1], pcoa.a$vectors[,2])
gen.dist <- as.data.frame(pcoa.a$vectors[,1])
gen.dist[,2] <- pcoa.a$vectors[,2]
colnames(gen.dist) <- c("PCoA1", "PCoA2")

## xyz dataset ####
taxa.geodist <- all[all$Individual %in% taxa.pop$Individual,]
taxa.geodist <- taxa.geodist[order(taxa.geodist$Individual),]

pca <- prcomp(taxa.geodist[,5:6])
PC1xy <- pca$x[,1]
PC2xy <- pca$x[,2]

plot(PC1xy, PC2xy)

plot(gen.dist$PCoA1, PC1xy)
plot(gen.dist$PCoA2, PC1xy)
plot(gen.dist$PCoA1, PC2xy)
plot(gen.dist$PCoA2, PC2xy)
plot(gen.dist$PCoA1, taxa.geodist$z)
plot(gen.dist$PCoA2, taxa.geodist$z)

# RDA full ####
rda.f <- rda(gen.dist ~ PC1xy * PC2xy * taxa.geodist$z)
rda.f
RsquareAdj(rda.f)
#$r.squared
#[1] 0.5023841
#$adj.r.squared
#[1] 0.4786881


summary(eigenvals(rda.f,model="constrained"))
#Importance of components:
#RDA1      RDA2
#Eigenvalue            0.001159 1.447e-05
#Proportion Explained  0.987665 1.234e-02
#Cumulative Proportion 0.987665 1.000e+00

summary(eigenvals(rda.f,model="unconstrained"))
#Importance of components:
#PC1       PC2
#Eigenvalue            0.0007027 0.0004595
#Proportion Explained  0.6046061 0.3953939
#Cumulative Proportion 0.6046061 1.0000000

summary(rda.f)
#Partitioning of variance:
#Inertia Proportion
#Total         0.002336     1.0000
#Constrained   0.001173     0.5024
#Unconstrained 0.001162     0.4976

#                               RDA1    RDA2 PC1 PC2
#PC1xy                       0.98771  0.1085   0   0
#PC2xy                       0.11607 -0.4933   0   0
#taxa.geodist$z             -0.06755  0.5512   0   0
#PC1xy:PC2xy                -0.01830  0.2184   0   0
#PC1xy:taxa.geodist$z       -0.92676 -0.1769   0   0
#PC2xy:taxa.geodist$z       -0.20890  0.5781   0   0
#PC1xy:PC2xy:taxa.geodist$z -0.32241 -0.2555   0   0


#Eigenvalues, and their contribution to the variance 

#Importance of components:
#  RDA1      RDA2      PC1       PC2
#Eigenvalue            0.00108 3.023e-05 0.002368 0.0009276
#Proportion Explained  0.24521 6.861e-03 0.537419 0.2105124
#Cumulative Proportion 0.24521 2.521e-01 0.789488 1.0000000

step(rda.f)
# rda(formula = gen.dist ~ PC1xy * PC2xy * taxa.geodist$z)

signif.full.c <- anova.cca(rda.f)
signif.full.c # not significant
#Model: rda(formula = gen.dist ~ PC1xy * PC2xy * taxa.geodist$z)
#Df  Variance      F Pr(>F)    
#Model      7 0.0011734 21.201  0.001 ***
#  Residual 147 0.0011623  
signif.axis.c <- anova.cca(rda.f, by="axis")
signif.axis.c
# RDA1       1 0.00115891 151.5636  0.001 ***

signif.terms.c <- anova.cca(rda.f, by="terms")
signif.terms.c
#PC1xy                        1 0.00069941 91.0238  0.001 ***
#PC2xy                        1 0.00003390  4.4118  0.025 *  
#  taxa.geodist$z               1 0.00000105  0.1364  0.883    
#PC1xy:PC2xy                  1 0.00000828  1.0777  0.312    
#PC1xy:taxa.geodist$z         1 0.00000248  0.3229  0.729    
#PC2xy:taxa.geodist$z         1 0.00005595  7.2811  0.003 ** 
#  PC1xy:PC2xy:taxa.geodist$z   1 0.00005674  7.3841  0.002 ** 

## RDA 1 ####
rda1 <- rda(gen.dist ~ PC1xy + PC2xy + PC2xy:taxa.geodist$z + PC1xy:PC2xy:taxa.geodist$z)
rda1
#                Inertia Proportion Rank
#Total         0.0025098  1.0000000     
#Constrained   0.0007717  0.3074792    2
#Unconstrained 0.0017381  0.6925208    2

RsquareAdj(rda1)
#$r.squared
#[1] 0.3074792
#$adj.r.squared
#[1] 0.2947724

summary(eigenvals(rda1,model="constrained"))
#Importance of components:
#RDA1      RDA2
#Eigenvalue            0.001146 4.523e-06
#Proportion Explained  0.996069 3.931e-03
#Cumulative Proportion 0.996069 1.000e+00

summary(eigenvals(rda1,model="unconstrained"))
#Importance of components:
#PC1       PC2
#Eigenvalue            0.0007116 0.0004732
#Proportion Explained  0.6006163 0.3993837
#Cumulative Proportion 0.6006163 1.0000000

summary(rda1)
#                   RDA1    RDA2 PC1 PC2
#PC1xy           0.99319 -0.1022   0   0
#PC2xy           0.11638  0.8932   0   0
#taxa.geodist$z -0.06756 -0.9923   0   0

step(rda1)
# -AIC scores
# step reckons: gen.dist ~ PC1xy + PC2xy + PC2xy:taxa.geodist$z

signif.full.c <- anova.cca(rda1)
signif.full.c # significant
#          Df   Variance      F Pr(>F)    
#Model      4 0.00077171 24.198  0.001 ***
#  Residual 218 0.00173809                 

signif.axis.c <- anova.cca(rda1, by="axis")
signif.axis.c
#          Df   Variance      F Pr(>F)    
#RDA1       1 0.00074681 94.528  0.001 ***
#  RDA2       1 0.00002490  3.152  0.379    
#Residual 220 0.00173809  
signif.terms.c <- anova.cca(rda1, by="terms")
signif.terms.c
#Df   Variance       F Pr(>F)    
#PC1xy                        1 0.00069941 86.9479  0.001 ***
#  PC2xy                        1 0.00003390  4.2142  0.021 *  
#  PC2xy:taxa.geodist$z         1 0.00001765  2.1948  0.101    
#PC1xy:PC2xy:taxa.geodist$z   1 0.00000526  0.6544  0.547    

## RDA 2 ####
rda2 <- rda(gen.dist ~ PC1xy + PC2xy)
rda2
#                Inertia Proportion Rank
#Total         0.0025098  1.0000000     
#Constrained   0.0007704  0.3069532    2
#Unconstrained 0.0017394  0.6930468    2
RsquareAdj(rda2)
#$r.squared
#[1] 0.3069532
#$adj.r.squared
#[1] 0.2974594

summary(eigenvals(rda2,model="constrained"))
#Importance of components:
#RDA1      RDA2
#Eigenvalue            0.001146 3.656e-06
#Proportion Explained  0.996820 3.180e-03
#Cumulative Proportion 0.996820 1.000e+00
# only RDA1 again. 

summary(eigenvals(rda2,model="unconstrained"))
#Importance of components:
#PC1       PC2
#Eigenvalue            0.002385 0.0009319
#Proportion Explained  0.719068 0.2809324
#Cumulative Proportion 0.719068 1.0000000

summary(rda2)
#       RDA1    RDA2 PC1 PC2
#PC1xy 0.9778 -0.2096   0   0
#PC2xy 0.2096  0.9778   0   0

#Eigenvalues, and their contribution to the variance 

#Importance of components:
#                        RDA1    RDA2 PC1 PC2
#PC1xy                 0.9679  0.1942   0   0
#PC2xy                 0.2321 -0.4491   0   0
#PC2xy:taxa.geodist$z -0.2846  0.5393   0   0

step(rda2)

signif.full.c <- anova.cca(rda2)
signif.full.c
#Model: rda(formula = gen.dist ~ PC1xy + PC2xy)
#Df  Variance      F Pr(>F)    
#Model      2 0.0011499 73.704  0.001 ***
#  Residual 152 0.0011857               

signif.axis.c <- anova.cca(rda2, by="axis")
signif.axis.c
# RDA2 not explaining much
signif.axis.c <- anova.cca(rda2, by="axis")
signif.axis.c

signif.terms.c <- anova.cca(rda2, by="terms")
signif.terms.c
#Df   Variance       F Pr(>F)    
#PC1xy                  1 0.00069941 88.0583  0.001 ***
#PC2xy                  1 0.00004519  5.6902  0.006 ** 
#PC2xy:taxa.geodist$z   1 0.00002579  3.2475  0.058 .  

## RDA 3 ####
rda3 <- rda(gen.dist ~ PC1xy)
rda3
RsquareAdj(rda3)
#$r.squared
#[1] 0.4841378
#$adj.r.squared
#[1] 0.4807661

summary(eigenvals(rda3,model="constrained"))
summary(eigenvals(rda3,model="unconstrained"))
#Importance of components:
#PC1       PC2
#Eigenvalue            0.0007153 0.0004895
#Proportion Explained  0.5937117 0.4062883
#Cumulative Proportion 0.5937117 1.0000000

summary(rda3)
#Importance of components:
#  RDA1       PC1       PC2
#Eigenvalue            0.001131 0.0007153 0.0004895
#Proportion Explained  0.484138 0.3062734 0.2095888
#Cumulative Proportion 0.484138 0.7904112 1.0000000
# RDA1 > PC1

signif.full.c <- anova.cca(rda3)
signif.full.c
#Model: rda(formula = gen.dist ~ PC1xy)
#Df  Variance      F Pr(>F)    
#Model      1 0.0011308 143.59  0.001 ***
#  Residual 153 0.0012049 
signif.axis.c <- anova.cca(rda3, by="axis")
signif.axis.c
#Model: rda(formula = gen.dist ~ PC1xy)
#Df  Variance      F Pr(>F)    
#RDA1       1 0.0011308 143.59  0.001 ***
#  Residual 153 0.0012049  

step(rda3) #AIC improves by heaps if you remove it?

## RDA 4 ####
rda4 <- rda(gen.dist ~ taxa.geodist$z | PC1xy + PC2xy)
summary(rda4)

# AH1 ####
# Import datasets
all <- read.csv("oriented_scaled_parallel_xyz_all.csv")
a <- read.delim("3b - Within & between plots/spatial-agaricia/results/AH1/AH1_all.a.txt", header = FALSE)
taxa.pop <- read.csv("3b - Within & between plots/spatial-agaricia/hu_3b_nc_20_AH1_all_.csv", row.names = 1)

## PCOA of genetic distance ####
row1 <- NA
a <- rbind(row1, a)
Vx <- NA
a <- cbind(a, Vx)

a <- as.matrix(a)
diag(a) <- 0
a.dist <- as.dist(a, upper = FALSE)

pcoa.a <- pcoa(a.dist) # cmdscale(a.dist) same thing
plot(pcoa.a$vectors[,1], pcoa.a$vectors[,2])
gen.dist <- as.data.frame(pcoa.a$vectors[,1])
gen.dist[,2] <- pcoa.a$vectors[,2]
colnames(gen.dist) <- c("PCoA1", "PCoA2")

## xyz dataset ####
taxa.geodist <- all[all$Individual %in% taxa.pop$Individual,]
taxa.geodist <- taxa.geodist[order(taxa.geodist$Individual),]

pca <- prcomp(taxa.geodist[,5:6])
PC1xy <- pca$x[,1]
PC2xy <- pca$x[,2]

plot(PC1xy, PC2xy)

plot(gen.dist$PCoA1, PC1xy)
plot(gen.dist$PCoA2, PC1xy)
plot(gen.dist$PCoA1, PC2xy)
plot(gen.dist$PCoA2, PC2xy)
plot(gen.dist$PCoA1, taxa.geodist$z)
plot(gen.dist$PCoA2, taxa.geodist$z)

# RDA full ####
rda.f <- rda(gen.dist ~ PC1xy * PC2xy * taxa.geodist$z)
rda.f
#Inertia Proportion Rank
#Total          0.8585     1.0000     
#Constrained    0.5507     0.6415    2
#Unconstrained  0.3078     0.3585    2
RsquareAdj(rda.f)
#$r.squared
#[1] 0.6414528
#$adj.r.squared
#[1] 0.3625828


summary(eigenvals(rda.f,model="constrained"))
#Importance of components:
#RDA1   RDA2
#Eigenvalue            0.3740 0.1766
#Proportion Explained  0.6793 0.3207
#Cumulative Proportion 0.6793 1.0000

summary(eigenvals(rda.f,model="unconstrained"))
#Importance of components:
#PC1     PC2
#Eigenvalue            0.2082 0.09956
#Proportion Explained  0.6765 0.32347
#Cumulative Proportion 0.6765 1.00000

summary(rda.f)
#                              RDA1     RDA2 PC1 PC2
#PC1xy                       0.6008  0.15526   0   0
#PC2xy                      -0.1915  0.09222   0   0
#taxa.geodist$z             -0.4012  0.19183   0   0
#PC1xy:PC2xy                -0.2200  0.64747   0   0
#PC1xy:taxa.geodist$z       -0.6510  0.01659   0   0
#PC2xy:taxa.geodist$z        0.2936 -0.16831   0   0
#PC1xy:PC2xy:taxa.geodist$z  0.4260 -0.58325   0   0

step(rda.f)
# wants to keep 3way interaction!
#gen.dist ~ PC1xy * PC2xy *taxa.geodist$z

signif.full.c <- anova.cca(rda.f)
signif.full.c # almost significant
#Model: rda(formula = gen.dist ~ PC1xy * PC2xy * taxa.geodist$z)
#Df Variance      F Pr(>F)  
#Model     7  0.55066 2.3002  0.062 .
#Residual  9  0.30780  
signif.axis.c <- anova.cca(rda.f, by="axis")
signif.axis.c
#         Df Variance       F Pr(>F)  
#RDA1      1  0.37404 17.0129  0.040 *
#  RDA2      1  0.17662  8.0336  0.581  
#Residual 14  0.30780                 

signif.terms.c <- anova.cca(rda.f, by="terms")
signif.terms.c
#PC1xy                       1 0.139274 4.0724  0.035 *
#PC1xy:taxa.geodist$z        1 0.111051 3.2471  0.076 .
#PC1xy:PC2xy:taxa.geodist$z  1 0.152346 4.4546  0.030 *

## RDA 1 ####
rda1 <- rda(gen.dist ~ PC1xy + PC1xy:taxa.geodist$z + PC1xy:PC2xy:taxa.geodist$z)
rda1
#              Inertia Proportion Rank
#Total          0.8585     1.0000     
#Constrained    0.2531     0.2948    2
#Unconstrained  0.6054     0.7052    2
RsquareAdj(rda1)
#$r.squared
#[1] 0.2948077
#$adj.r.squared
#[1] 0.132071


summary(eigenvals(rda1,model="constrained"))
#Importance of components:
#RDA1    RDA2
#Eigenvalue            0.1674 0.08572
#Proportion Explained  0.6613 0.33871
#Cumulative Proportion 0.6613 1.00000


summary(eigenvals(rda1,model="unconstrained"))
#Importance of components:
#PC1    PC2
#Eigenvalue            0.3731 0.2323
#Proportion Explained  0.6162 0.3838
#Cumulative Proportion 0.6162 1.0000

summary(rda1)
#                              RDA1    RDA2 PC1 PC2
#PC1xy                       0.8270 -0.5381   0   0
#PC1xy:taxa.geodist$z       -0.9449  0.3267   0   0
#PC1xy:taxa.geodist$z:PC2xy  0.7695  0.5802   0   0

step(rda1)
# very small changes in AIC

signif.full.c <- anova.cca(rda1)
signif.full.c # not significant
#         Df Variance      F Pr(>F)
#Model     3  0.25308 1.8116   0.14
#Residual 13  0.60538  
signif.axis.c <- anova.cca(rda1, by="axis")
signif.axis.c
#        Df Variance      F Pr(>F)
#RDA1      1  0.16736 3.8704  0.227
#RDA2      1  0.08572 1.9824  0.482           

signif.terms.c <- anova.cca(rda1, by="terms")
signif.terms.c
#PC1xy                       1  0.13927 2.9908  0.068 .
#PC1xy:taxa.geodist$z        1  0.07883 1.6929  0.197  
#PC1xy:taxa.geodist$z:PC2xy  1  0.03497 0.7510  0.490  
# Not sure why this one isn't now.

# AH2 ####

# Import datasets
all <- read.csv("oriented_scaled_parallel_xyz_all.csv")
a <- read.delim("3b - Within & between plots/spatial-agaricia/results/AH2/AH2_all.a.txt", header = FALSE)
taxa.pop <- read.csv("3b - Within & between plots/spatial-agaricia/hu_3b_nc_20_AH2_all_.csv", row.names = 1)

## PCOA of genetic distance ####
row1 <- NA
a <- rbind(row1, a)
Vx <- NA
a <- cbind(a, Vx)

a <- as.matrix(a)
diag(a) <- 0
a.dist <- as.dist(a, upper = FALSE)

pcoa.a <- pcoa(a.dist) # cmdscale(a.dist) same thing
plot(pcoa.a$vectors[,1], pcoa.a$vectors[,2])
gen.dist <- as.data.frame(pcoa.a$vectors[,1])
gen.dist[,2] <- pcoa.a$vectors[,2]
colnames(gen.dist) <- c("PCoA1", "PCoA2")

## xyz dataset ####
taxa.geodist <- all[all$Individual %in% taxa.pop$Individual,]
taxa.geodist <- taxa.geodist[order(taxa.geodist$Individual),]

pca <- prcomp(taxa.geodist[,5:6])
PC1xy <- pca$x[,1]
PC2xy <- pca$x[,2]

plot(PC1xy, PC2xy)

plot(gen.dist$PCoA1, PC1xy)
plot(gen.dist$PCoA2, PC1xy)
plot(gen.dist$PCoA1, PC2xy)
plot(gen.dist$PCoA2, PC2xy)
plot(gen.dist$PCoA1, taxa.geodist$z)
plot(gen.dist$PCoA2, taxa.geodist$z)

# RDA full ####
rda.f <- rda(gen.dist ~ PC1xy * PC2xy * taxa.geodist$z)
rda.f
#              Inertia Proportion Rank
#Total         0.52045    1.00000     
#Constrained   0.44211    0.84948    2
#Unconstrained 0.07834    0.15052    2
RsquareAdj(rda.f)
#$r.squared
#[1] 0.8494792
#$adj.r.squared
#[1] 0.6387502


summary(eigenvals(rda.f,model="constrained"))
#Importance of components:
#RDA1   RDA2
#Eigenvalue            0.2853 0.1569
#Proportion Explained  0.6452 0.3548
#Cumulative Proportion 0.6452 1.0000

summary(eigenvals(rda.f,model="unconstrained"))
#Importance of components:
#PC1     PC2
#Eigenvalue            0.05809 0.02025
#Proportion Explained  0.74157 0.25843
#Cumulative Proportion 0.74157 1.00000

summary(rda.f)
#                              RDA1     RDA2 PC1 PC2
#PC1xy                      -0.1353  0.61741   0   0
#PC2xy                      -0.3854 -0.03816   0   0
#taxa.geodist$z              0.3514  0.27693   0   0
#PC1xy:PC2xy                 0.6566  0.11437   0   0
#PC1xy:taxa.geodist$z       -0.0219 -0.66757   0   0
#PC2xy:taxa.geodist$z        0.3910  0.08329   0   0
#PC1xy:PC2xy:taxa.geodist$z -0.6172 -0.27202   0   0

step(rda.f)
# wants to keep 3way interaction!
#gen.dist ~ PC1xy * PC2xy * taxa.geodist$z)

signif.full.c <- anova.cca(rda.f)
signif.full.c # almost significant
#Model: rda(formula = gen.dist ~ PC1xy * PC2xy * taxa.geodist$z)
#Df Variance      F Pr(>F)  
#Model     7  0.44211 4.0311  0.027 *
#  Residual  5  0.07834  
signif.axis.c <- anova.cca(rda.f, by="axis")
signif.axis.c
#        Df Variance      F Pr(>F)  
#RDA1      1 0.285259 36.413  0.030 *
#  RDA2      1 0.156854 20.023  0.268  
#Residual 10 0.078339               

signif.terms.c <- anova.cca(rda.f, by="terms")
signif.terms.c
#Model: rda(formula = gen.dist ~ PC1xy * PC2xy * taxa.geodist$z)
#Df Variance      F Pr(>F)  
#PC1xy                       1 0.065010 4.1493  0.052 .
#PC2xy                       1 0.042605 2.7193  0.124  
#taxa.geodist$z              1 0.063198 4.0336  0.066 .
#PC1xy:PC2xy                 1 0.116102 7.4103  0.016 *
#  PC1xy:taxa.geodist$z        1 0.100715 6.4282  0.012 *
#  PC2xy:taxa.geodist$z        1 0.021317 1.3606  0.317  
#PC1xy:PC2xy:taxa.geodist$z  1 0.033166 2.1168  0.181  

## RDA 1 ####
rda1 <- rda(gen.dist ~ PC1xy + PC2xy + taxa.geodist$z + PC1xy:PC2xy + PC1xy:taxa.geodist$z)
rda1
#              Inertia Proportion Rank
#Total          0.5205     1.0000     
#Constrained    0.3876     0.7448    2
#Unconstrained  0.1328     0.2552    2
RsquareAdj(rda1)
#$r.squared
#[1] 0.7447959
#$adj.r.squared
#[1] 0.5625073


summary(eigenvals(rda1,model="constrained"))
#Importance of components:
#RDA1   RDA2
#Eigenvalue            0.2619 0.1258
#Proportion Explained  0.6756 0.3244
#Cumulative Proportion 0.6756 1.0000


summary(eigenvals(rda1,model="unconstrained"))
#Importance of components:
#PC1     PC2
#Eigenvalue            0.0909 0.04193
#Proportion Explained  0.6843 0.31566
#Cumulative Proportion 0.6843 1.00000

summary(rda1)
#                          RDA1     RDA2 PC1 PC2
#PC1xy                -0.165383  0.67824   0   0
#PC2xy                -0.400232 -0.07221   0   0
#taxa.geodist$z        0.355363  0.33592   0   0
#PC1xy:PC2xy           0.679909  0.17809   0   0
#PC1xy:taxa.geodist$z  0.003566 -0.74628   0   0

step(rda1)
# PC1xy + PC2xy + taxa.geodist$z + PC1xy:PC2xy + PC1xy:taxa.geodist$z
# not very small but better to keep them all!

signif.full.c <- anova.cca(rda1)
signif.full.c # significant
#Model: rda(formula = gen.dist ~ PC1xy + PC2xy + taxa.geodist$z + PC1xy:PC2xy + PC1xy:taxa.geodist$z)
#Df Variance      F Pr(>F)  
#Model     5  0.38763 4.0858   0.02 *
#  Residual  7  0.13282 
signif.axis.c <- anova.cca(rda1, by="axis")
signif.axis.c
#         Df Variance       F Pr(>F)  
#RDA1      1  0.26188 19.7166  0.012 *
#  RDA2      1  0.12575  9.4677  0.220  
#Residual 10  0.13282     

signif.terms.c <- anova.cca(rda1, by="terms")
signif.terms.c
#                     Df Variance      F Pr(>F)  
#PC1xy                 1 0.065010 3.4262  0.051 .
#PC2xy                 1 0.042605 2.2454  0.152  
#taxa.geodist$z        1 0.063198 3.3307  0.069 .
#PC1xy:PC2xy           1 0.116102 6.1188  0.012 *
#  PC1xy:taxa.geodist$z  1 0.100715 5.3079  0.011 *
#  Residual              7 0.132822   
# very interesting but only three individuals found at 10 m.... over parameterised?
# AH3 ####

# Import datasets
all <- read.csv("oriented_scaled_parallel_xyz_all.csv")
a <- read.delim("3b - Within & between plots/spatial-agaricia/results/AH3/AH3_all.a.txt", header = FALSE)
taxa.pop <- read.csv("3b - Within & between plots/spatial-agaricia/hu_3b_nc_20_AH3_all_.csv", row.names = 1)

## PCOA of genetic distance ####
row1 <- NA
a <- rbind(row1, a)
Vx <- NA
a <- cbind(a, Vx)

a <- as.matrix(a)
diag(a) <- 0
a.dist <- as.dist(a, upper = FALSE)

pcoa.a <- pcoa(a.dist) # cmdscale(a.dist) same thing
plot(pcoa.a$vectors[,1], pcoa.a$vectors[,2])
gen.dist <- as.data.frame(pcoa.a$vectors[,1])
gen.dist[,2] <- pcoa.a$vectors[,2]
colnames(gen.dist) <- c("PCoA1", "PCoA2")

## xyz dataset ####
taxa.geodist <- all[all$Individual %in% taxa.pop$Individual,]
taxa.geodist <- taxa.geodist[order(taxa.geodist$Individual),]

pca <- prcomp(taxa.geodist[,5:6])
PC1xy <- pca$x[,1]
PC2xy <- pca$x[,2]

plot(PC1xy, PC2xy)

plot(gen.dist$PCoA1, PC1xy)
plot(gen.dist$PCoA2, PC1xy)
plot(gen.dist$PCoA1, PC2xy)
plot(gen.dist$PCoA2, PC2xy)
plot(gen.dist$PCoA1, taxa.geodist$z)
plot(gen.dist$PCoA2, taxa.geodist$z)

# RDA full ####
rda.f <- rda(gen.dist ~ PC1xy * PC2xy * taxa.geodist$z) # too many parameters
rda.f
#              Inertia Proportion Rank
#Total          9.3524     1.0000     
#Constrained    7.2053     0.7704    2
#Unconstrained  2.1470     0.2296    2
RsquareAdj(rda.f)
#$r.squared
#[1] 0.7704273
#$adj.r.squared
#[1] -0.03307735


summary(eigenvals(rda.f,model="constrained"))
#Importance of components:
#RDA1   RDA2
#Eigenvalue            4.3914 2.8140
#Proportion Explained  0.6095 0.3905
#Cumulative Proportion 0.6095 1.0000

summary(eigenvals(rda.f,model="unconstrained"))
#                         PC1    PC2
#Eigenvalue            1.7320 0.4151
#Proportion Explained  0.8067 0.1933
#Cumulative Proportion 0.8067 1.0000

summary(rda.f)
#                               RDA1    RDA2 PC1 PC2
#PC1xy                      -0.07563 -0.4198   0   0
#PC2xy                      -0.22732  0.2363   0   0
#taxa.geodist$z             -0.17984  0.3569   0   0
#PC1xy:PC2xy                -0.22851  0.2366   0   0
#PC1xy:taxa.geodist$z        0.04281  0.4538   0   0
#PC2xy:taxa.geodist$z        0.20147 -0.2907   0   0
#PC1xy:PC2xy:taxa.geodist$z  0.20215 -0.2907   0   0


step(rda.f)
# wants to keep 3way interaction! maybe an artefact of too many parameters?
#gen.dist ~ PC1xy * PC2xy * taxa.geodist$z)

signif.full.c <- anova.cca(rda.f)
signif.full.c # not significant
#Model: rda(formula = gen.dist ~ PC1xy * PC2xy * taxa.geodist$z)
#Df Variance      F Pr(>F)
#Model     7   7.2053 0.9588    0.6
#Residual  2   2.1470   
signif.axis.c <- anova.cca(rda.f, by="axis")
signif.axis.c
#         Df Variance       F Pr(>F)
#RDA1      1   4.3914 14.3171  0.623
#RDA2      1   2.8140  9.1743  0.814
#Residual  7   2.1470               

signif.terms.c <- anova.cca(rda.f, by="terms")
signif.terms.c
#None

## RDA 1 ####
rda1 <- rda(gen.dist ~ PC1xy + PC2xy + taxa.geodist$z)
rda1
#              Inertia Proportion Rank
#Total          9.3524     1.0000     
#Constrained    1.0223     0.1093    2
#Unconstrained  8.3300     0.8907    2
RsquareAdj(rda1)
#$r.squared
#[1] 0.1093108
#$adj.r.squared
#[1] -0.3360338 # still negative


summary(eigenvals(rda1,model="constrained"))
#                        RDA1   RDA2
#Eigenvalue            0.7073 0.3150
#Proportion Explained  0.6919 0.3081
#Cumulative Proportion 0.6919 1.0000


summary(eigenvals(rda1,model="unconstrained"))
#                         PC1    PC2
#Eigenvalue            5.5329 2.7971
#Proportion Explained  0.6642 0.3358
#Cumulative Proportion 0.6642 1.0000

summary(rda1)
#                   RDA1   RDA2 PC1 PC2
#PC1xy           0.8274 0.3423   0   0
#PC2xy          -0.4979 0.8139   0   0
#taxa.geodist$z -0.7326 0.6195   0   0

step(rda1)
# remove all

signif.full.c <- anova.cca(rda1)
signif.full.c # not significant

signif.axis.c <- anova.cca(rda1, by="axis")
signif.axis.c # not


signif.terms.c <- anova.cca(rda1, by="terms")
signif.terms.c # not

# AL1 ####

# Import datasets
all <- read.csv("oriented_scaled_parallel_xyz_all.csv")
a <- read.delim("3b - Within & between plots/spatial-agaricia/results/AL1/AL1_all.a.txt", header = FALSE)
taxa.pop <- read.csv("3b - Within & between plots/spatial-agaricia/lm_3b_nc_20_AL1_all_.csv", row.names = 1)

## PCOA of genetic distance ####
row1 <- NA
a <- rbind(row1, a)
Vx <- NA
a <- cbind(a, Vx)

a <- as.matrix(a)
diag(a) <- 0
a.dist <- as.dist(a, upper = FALSE)

pcoa.a <- pcoa(a.dist) # cmdscale(a.dist) same thing
plot(pcoa.a$vectors[,1], pcoa.a$vectors[,2])
gen.dist <- as.data.frame(pcoa.a$vectors[,1])
gen.dist[,2] <- pcoa.a$vectors[,2]
colnames(gen.dist) <- c("PCoA1", "PCoA2")

## xyz dataset ####
taxa.geodist <- all[all$Individual %in% taxa.pop$Individual,]
taxa.geodist <- taxa.geodist[order(taxa.geodist$Individual),]

pca <- prcomp(taxa.geodist[,5:6])
PC1xy <- pca$x[,1]
PC2xy <- pca$x[,2]

plot(PC1xy, PC2xy)

plot(gen.dist$PCoA1, PC1xy)
plot(gen.dist$PCoA2, PC1xy)
plot(gen.dist$PCoA1, PC2xy)
plot(gen.dist$PCoA2, PC2xy)
plot(gen.dist$PCoA1, taxa.geodist$z)
plot(gen.dist$PCoA2, taxa.geodist$z)

# RDA full ####
rda.f <- rda(gen.dist ~ PC1xy * PC2xy * taxa.geodist$z) # too many parameters
rda.f
#               Inertia Proportion Rank
#Total         0.005011   1.000000     
#Constrained   0.002813   0.561405    2
#Unconstrained 0.002198   0.438595    2
RsquareAdj(rda.f)
#$r.squared
#[1] 0.5614052
#$adj.r.squared
#[1] 0.04971136


summary(eigenvals(rda.f,model="constrained"))
#Importance of components:
#RDA1     RDA2
#Eigenvalue            0.001794 0.001019
#Proportion Explained  0.637739 0.362261
#Cumulative Proportion 0.637739 1.000000

summary(eigenvals(rda.f,model="unconstrained"))
#Importance of components:
#PC1       PC2
#Eigenvalue            0.001502 0.0006957
#Proportion Explained  0.683469 0.3165314
#Cumulative Proportion 0.683469 1.0000000

summary(rda.f)
#                               RDA1     RDA2 PC1 PC2
#PC1xy                       0.69387 -0.06184   0   0
#PC2xy                       0.03872 -0.24394   0   0
#taxa.geodist$z              0.15033 -0.07058   0   0
#PC1xy:PC2xy                -0.38191 -0.34351   0   0
#PC1xy:taxa.geodist$z       -0.68565  0.07499   0   0
#PC2xy:taxa.geodist$z       -0.04553  0.25092   0   0
#PC1xy:PC2xy:taxa.geodist$z  0.37323  0.34711   0   0

step(rda.f)
# wants to keep 3way interaction! maybe an artefact of too many parameters?
#gen.dist ~ PC1xy * PC2xy * taxa.geodist$z)

signif.full.c <- anova.cca(rda.f)
signif.full.c # not significant
#         Df  Variance      F Pr(>F)
#Model     7 0.0028133 1.0972  0.471
#Residual  6 0.0021978   
signif.axis.c <- anova.cca(rda.f, by="axis")
signif.axis.c

signif.terms.c <- anova.cca(rda.f, by="terms")
signif.terms.c
#None

## RDA 1 ####
rda1 <- rda(gen.dist ~ PC1xy + PC2xy + taxa.geodist$z)
rda1
#               Inertia Proportion Rank
#Total         0.005011   1.000000     
#Constrained   0.001554   0.310027    2
#Unconstrained 0.003458   0.689973    2
RsquareAdj(rda1)
#$r.squared
#[1] 0.3100274
#$adj.r.squared
#[1] 0.1030357

summary(eigenvals(rda1,model="constrained"))
#                          RDA1      RDA2
#Eigenvalue            0.001393 0.0001604
#Proportion Explained  0.896771 0.1032292
#Cumulative Proportion 0.896771 1.0000000


summary(eigenvals(rda1,model="unconstrained"))
#                           PC1      PC2
#Eigenvalue            0.002083 0.001375
#Proportion Explained  0.602459 0.397541
#Cumulative Proportion 0.602459 1.000000

summary(rda1)
#                 RDA1     RDA2 PC1 PC2
#PC1xy          0.7642 -0.58084   0   0
#PC2xy          0.1072  0.54315   0   0
#taxa.geodist$z 0.1809  0.01101   0   0

step(rda1)
# only 2AIC in it all none of the models are that diff from each other

signif.full.c <- anova.cca(rda1)
signif.full.c # not significant

signif.axis.c <- anova.cca(rda1, by="axis")
signif.axis.c # not


signif.terms.c <- anova.cca(rda1, by="terms")
signif.terms.c # xy a little
#PC1xy           1 0.0008677 2.5096  0.097 .

## RDA2 ####
rda2 <- rda(gen.dist ~ PC1xy + taxa.geodist$z)
rda2
#               Inertia Proportion Rank
#Total         0.005011   1.000000     
#Constrained   0.000989   0.197355    2
#Unconstrained 0.004022   0.802645    2
RsquareAdj(rda2)
#$r.squared
#[1] 0.1973553
#$adj.r.squared
#[1] 0.05141995

summary(eigenvals(rda2,model="constrained"))
#                          RDA1      RDA2
#Eigenvalue            0.000986 2.981e-06
#Proportion Explained  0.996985 3.015e-03
#Cumulative Proportion 0.996985 1.000e+00


summary(eigenvals(rda2,model="unconstrained"))
#                           PC1      PC2
#Eigenvalue            0.002425 0.001597
#Proportion Explained  0.602930 0.397070
#Cumulative Proportion 0.602930 1.000000

summary(rda2)
#                 RDA1    RDA2 PC1 PC2
#PC1xy          0.9379  0.3469   0   0
#taxa.geodist$z 0.2083 -0.9781   0   0

step(rda2)
# Call: rda(formula = gen.dist ~ PC1xy) but only small diff in AIC

signif.full.c <- anova.cca(rda2)
signif.full.c # not significant

signif.axis.c <- anova.cca(rda2, by="axis")
signif.axis.c # not


signif.terms.c <- anova.cca(rda2, by="terms")
signif.terms.c # not
# AL2 ####

# Import datasets
all <- read.csv("oriented_scaled_parallel_xyz_all.csv")
a <- read.delim("3b - Within & between plots/spatial-agaricia/results/AL2/AL2_all.a.txt", header = FALSE)
taxa.pop <- read.csv("3b - Within & between plots/spatial-agaricia/lm_3b_nc_20_AL2_all_.csv", row.names = 1)

## PCOA of genetic distance ####
row1 <- NA
a <- rbind(row1, a)
Vx <- NA
a <- cbind(a, Vx)

a <- as.matrix(a)
diag(a) <- 0
a.dist <- as.dist(a, upper = FALSE)

pcoa.a <- pcoa(a.dist) # cmdscale(a.dist) same thing
plot(pcoa.a$vectors[,1], pcoa.a$vectors[,2])
gen.dist <- as.data.frame(pcoa.a$vectors[,1])
gen.dist[,2] <- pcoa.a$vectors[,2]
colnames(gen.dist) <- c("PCoA1", "PCoA2")

## xyz dataset ####
taxa.geodist <- all[all$Individual %in% taxa.pop$Individual,]
taxa.geodist <- taxa.geodist[order(taxa.geodist$Individual),]

pca <- prcomp(taxa.geodist[,5:6])
PC1xy <- pca$x[,1]
PC2xy <- pca$x[,2]

plot(PC1xy, PC2xy)

plot(gen.dist$PCoA1, PC1xy)
plot(gen.dist$PCoA2, PC1xy)
plot(gen.dist$PCoA1, PC2xy)
plot(gen.dist$PCoA2, PC2xy)
plot(gen.dist$PCoA1, taxa.geodist$z)
plot(gen.dist$PCoA2, taxa.geodist$z)

# RDA full ####
rda.f <- rda(gen.dist ~ PC1xy * PC2xy * taxa.geodist$z) # too many parameters
rda.f
#               Inertia Proportion Rank
#Total         0.005372   1.000000     
#Constrained   0.000692   0.128798    2
#Unconstrained 0.004680   0.871202    2
RsquareAdj(rda.f)
#$r.squared
#[1] 0.128798
#$adj.r.squared
#[1] -0.1151386


summary(eigenvals(rda.f,model="constrained"))
#                           RDA1      RDA2
#Eigenvalue            0.0005606 0.0001313
#Proportion Explained  0.8102120 0.1897880
#Cumulative Proportion 0.8102120 1.0000000

summary(eigenvals(rda.f,model="unconstrained"))
#                           PC1      PC2
#Eigenvalue            0.003435 0.001246
#Proportion Explained  0.733879 0.266121
#Cumulative Proportion 0.733879 1.000000

summary(rda.f)
#                               RDA1     RDA2 PC1 PC2
#PC1xy                      -0.33488 -0.22478   0   0
#PC2xy                       0.11890 -0.03158   0   0
#taxa.geodist$z             -0.12107  0.07276   0   0
#PC1xy:PC2xy                -0.03137  0.75111   0   0
#PC1xy:taxa.geodist$z        0.33569  0.11645   0   0
#PC2xy:taxa.geodist$z        0.01453 -0.03927   0   0
#PC1xy:PC2xy:taxa.geodist$z -0.01390 -0.65854   0   0

step(rda.f)
# Call: rda(formula = gen.dist ~ PC2xy + taxa.geodist$z + PC2xy:taxa.geodist$z)
# wants to keep interaction! with PC2 and z

signif.full.c <- anova.cca(rda.f)
signif.full.c # not significant

signif.axis.c <- anova.cca(rda.f, by="axis")
signif.axis.c

signif.terms.c <- anova.cca(rda.f, by="terms")
signif.terms.c
#None

## RDA 1 ####
rda1 <- rda(gen.dist ~ PC2xy * taxa.geodist$z)
rda1
#                Inertia Proportion Rank
#Total         0.0053724  1.0000000     
#Constrained   0.0003811  0.0709332    2
#Unconstrained 0.0049914  0.9290668    2
RsquareAdj(rda1)
#$r.squared
#[1] 0.07093321
#$adj.r.squared
#[1] -0.02517715


summary(eigenvals(rda1,model="constrained"))
#                           RDA1      RDA2
#Eigenvalue            0.0003732 7.888e-06
#Proportion Explained  0.9793000 2.070e-02
#Cumulative Proportion 0.9793000 1.000e+00

summary(eigenvals(rda1,model="unconstrained"))
#                           PC1      PC2
#Eigenvalue            0.003604 0.001388
#Proportion Explained  0.721979 0.278021
#Cumulative Proportion 0.721979 1.000000

summary(rda1)
#                         RDA1     RDA2 PC1 PC2
#PC2xy                 0.14631 -0.09267   0   0
#taxa.geodist$z       -0.15424 -0.06729   0   0
#PC2xy:taxa.geodist$z  0.02246  0.12969   0   0

step(rda1)
# Call: rda(formula = gen.dist ~ PC2xy + taxa.geodist$z + PC2xy:taxa.geodist$z)
# not really different... 

signif.full.c <- anova.cca(rda1)
signif.full.c # not significant

signif.axis.c <- anova.cca(rda1, by="axis")
signif.axis.c

signif.terms.c <- anova.cca(rda1, by="terms")
signif.terms.c
#None

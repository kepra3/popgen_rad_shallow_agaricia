# Packages
library(vegan)
library(tidyverse)
library(ape)
library(tidyr)
library(dplyr)
library(corrplot)


# Arguments
args <- commandArgs(TRUE)
taxa <- args[1]

print("*****************************************")
print(paste0("Performing Redundancy Analysis for ", taxa))
print("*****************************************")

all <- read.csv("../data/all_annotations_X_DEPTH_parallel_XYZ_adjusted.txt")
all <- all %>% separate(Individual, into = c("Sample", "Species", "Site"), remove = FALSE) %>% 
  unite(Individual, Sample, Species, Site, remove = FALSE) %>% 
  dplyr::select(Individual, Sample, Site, Loc, Depth, x, y, z)

# TAXA ####
## Import datasets
a <- read.delim(paste0("../results/ibd/", taxa, "/", taxa, "_all.a.txt"), header = FALSE)
if (taxa == "AA1" | taxa == "AA2") {
  taxa.pop <- read.csv(paste0("../data/", "ac_3b_nc_20_", taxa, "_all_.csv"), row.names = 1)
} else if (taxa == "AH1" | taxa == "AH2" | taxa == "AH3") {
  taxa.pop <- read.csv(paste0("../data/", "hu_3b_nc_20_", taxa, "_all_.csv"), row.names = 1)
} else if (taxa == "AL1" | taxa == "AL2") {
  taxa.pop <- read.csv(paste0("../data/", "lm_3b_nc_20_", taxa, "_all_.csv"), row.names = 1)
}


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

# PCA ####
pca <- prcomp(taxa.geodist[,6:7])
PC1xy <- pca$x[,1]
PC2xy <- pca$x[,2]
dat <- cbind(taxa.geodist$Individual, PC1xy, PC2xy, taxa.geodist[,6:8])
write.csv(dat, paste0("../results/rda/rda_data_", taxa, ".csv"), quote = FALSE)

# RDAs #####
# Full model
rda1 <- rda(gen.dist ~ PC1xy * PC2xy * taxa.geodist$z)
rda1
RsquareAdj(rda1)
signif.full.c <- anova.cca(rda1)
signif.full.c
signif.terms.c <- anova.cca(rda1, by="terms")
signif.terms.c
step(rda1)

write.csv(signif.full.c, paste0("../results/rda/fullmodeloverall_", taxa, ".csv"), quote = FALSE)
write.csv(signif.terms.c, paste0("../results/rda/fullmodelterms_", taxa, ".csv"), quote = FALSE)

if (taxa == "AA1" | taxa == "AA2" | taxa == "AH2" | taxa == "AH3" | taxa == "AL1" | taxa == "AL2") {
  
  if (taxa == "AA2") {
    print("Three way interaction significant within some sites & plots, Z or PC2xy significant,
        but only small effect")
  }
  print("Remove all non-significant terms")
  rda2 <- rda(gen.dist ~ PC1xy)
  rda2
  RsquareAdj(rda2)
  signif.full.c <- anova.cca(rda2)
  signif.full.c
  signif.terms.c <- anova.cca(rda2, by="terms")
  signif.terms.c
  write.csv(signif.full.c, paste0("../results/rda/finalmodeloverall_", taxa, ".csv"), quote = FALSE)
  write.csv(signif.terms.c, paste0("../results/rda/finalmodelterms_", taxa, ".csv"), quote = FALSE)
  
  rda3 <- rda(gen.dist ~ PC1xy + PC2xy + taxa.geodist$z)
  rda3
  RsquareAdj(rda3)
  signif.full.c <- anova.cca(rda3)
  signif.full.c
  signif.terms.c <- anova.cca(rda3, by="terms")
  signif.terms.c
  
  
  rda4 <- rda(gen.dist ~ PC1xy + PC2xy)
  rda4
  RsquareAdj(rda4)
  signif.full.c <- anova.cca(rda4)
  signif.full.c
  signif.terms.c <- anova.cca(rda4, by="terms")
  signif.terms.c
  
  
  rda5 <- rda(gen.dist ~ PC1xy + taxa.geodist$z)
  rda5
  RsquareAdj(rda5)
  signif.full.c <- anova.cca(rda5)
  signif.full.c
  signif.terms.c <- anova.cca(rda5, by="terms")
  signif.terms.c
  
} else if (taxa == "AH1") {
  print("Remove all non-significant terms")
  rda2 <- rda(gen.dist ~ PC1xy + PC2xy + taxa.geodist$z + PC2xy:taxa.geodist$z)
  rda2
  RsquareAdj(rda2)
  signif.full.c <- anova.cca(rda2)
  signif.full.c
  signif.terms.c <- anova.cca(rda2, by="terms")
  signif.terms.c
  step(rda2)
  print("Interaction significant within some plots Z or PC2xy significant,
        but only small effect") 
  print("PC2xy significant, but PC2xy and Z very similar (the order of the model may be causing this effect)")

  rda3 <- rda(gen.dist ~ PC1xy + PC2xy)
  rda3
  RsquareAdj(rda3)
  signif.full.c <- anova.cca(rda3)
  signif.full.c
  signif.terms.c <- anova.cca(rda3, by="terms")
  signif.terms.c
  step(rda3)
  # PC2xy has a slight effect
  write.csv(signif.full.c, paste0("../results/rda/finalmodeloverallPC2xy_", taxa, ".csv"), quote = FALSE)
  write.csv(signif.terms.c, paste0("../results/rda/finalmodeltermsPC2xy_", taxa, ".csv"), quote = FALSE)
  
  rda4 <- rda(gen.dist ~ PC1xy + taxa.geodist$z)
  rda4
  RsquareAdj(rda4)
  signif.full.c <- anova.cca(rda4)
  signif.full.c
  signif.terms.c <- anova.cca(rda4, by="terms")
  signif.terms.c
  step(rda4)
  # Z explains the same amount of variation
  write.csv(signif.full.c, paste0("../results/rda/finalmodeloverallZ_", taxa, ".csv"), quote = FALSE)
  write.csv(signif.terms.c, paste0("../results/rda/finalmodeltermsZ_", taxa, ".csv"), quote = FALSE)
  
  # Conditioning out variables to see individual variation explained
  rda5 <- rda(gen.dist ~ PC1xy  + Condition(PC2xy + taxa.geodist$z))
  rda5
  RsquareAdj(rda5)
  signif.full.c <- anova.cca(rda5)
  signif.full.c
  signif.terms.c <- anova.cca(rda5, by="terms")
  signif.terms.c
  
  rda6 <- rda(gen.dist ~ PC2xy  + Condition(PC1xy + taxa.geodist$z))
  rda6
  RsquareAdj(rda6)
  signif.full.c <- anova.cca(rda6)
  signif.full.c
  signif.terms.c <- anova.cca(rda6, by="terms")
  signif.terms.c
  
  rda7 <- rda(gen.dist ~ taxa.geodist$z + Condition(PC1xy + PC2xy))
  rda7
  RsquareAdj(rda7)
  signif.full.c <- anova.cca(rda7)
  signif.full.c
  signif.terms.c <- anova.cca(rda7, by="terms")
  signif.terms.c
  
  rda8 <- rda(gen.dist ~ PC1xy + PC2xy + Condition(taxa.geodist$z))
  rda8
  RsquareAdj(rda8)
  signif.full.c <- anova.cca(rda8)
  signif.full.c
  signif.terms.c <- anova.cca(rda8, by="terms")
  signif.terms.c
  
  rda9 <- rda(gen.dist ~ PC1xy + taxa.geodist$z + Condition(PC2xy))
  rda9
  RsquareAdj(rda9)
  signif.full.c <- anova.cca(rda9)
  signif.full.c
  signif.terms.c <- anova.cca(rda9, by="terms")
  signif.terms.c
  
  rda10 <- rda(gen.dist ~ PC2xy + taxa.geodist$z + Condition(PC1xy))
  rda10
  RsquareAdj(rda10)
  signif.full.c <- anova.cca(rda10)
  signif.full.c
  signif.terms.c <- anova.cca(rda10, by="terms")
  signif.terms.c
  
}

print("*****************************************")
print(paste0("Finished RDA analysis for ", taxa))
print("*****************************************")
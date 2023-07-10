# Quick PCA separate species
# Author: Katharine Prata
# Description: multiple types of PCAs
# Last edit: 28/06/23

## Dependencies ================================================
library(ggplot2)
library(adegenet)
library(dplyr)
library(tidyr)
library(stringr)
library(vcfR)

## Functions ===================================================
sort_pca <- function(pca, pop, nf) {
  pca_info <- cbind(pop, pca$scores[, 1:nf])
  return(pca_info)
}
pca_plot_species <- function(info, eig, title, int1, int2, colours, n.loc){
  ggplot(info, aes_(x = as.name(paste0("PC",int1)), y = as.name(paste0("PC",int2)),
                    fill= info$Species)) + 
    geom_point(size = 2, shape = 21, colour = "grey") +
    scale_fill_manual(name = "Species", values = colours, breaks = c("AA", "AH", "AL", "AG")) +
    xlab(paste0("PC", int1, " ", eig[int1], "%")) +
    ylab(paste0("PC", int2, " ", eig[int2], "%")) +
    #ggtitle(paste0("PCA of ", DATA_NAME, ", n = ", count(info), ", SNPs = ", n.loc), 
    #        subtitle = title) +
    theme(plot.background = element_blank(),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          legend.background = element_blank(),
          legend.key = element_blank(),
          panel.border = element_rect(colour = "grey", fill = NA, size = 1),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 8)) + facet_wrap(~Depth, nrow = 6, ncol = 1)
}
pca_plot_depth <- function(info, eig, title, int1, int2, colours, shapes, n.loc){
  ggplot(info, aes_(x = as.name(paste0("PC", int1)), y = as.name(paste0("PC", int2)),
                    colour = info$Loc, shape = info$Depth)) + 
    geom_point(alpha = 1/2, size = 5) +
    scale_color_manual(name = "Location", values = colours, breaks = c("WP", "CA", "SB", "SQ")) +
    scale_shape_manual(name = "Depth", values = shapes, breaks = c("5", "10", "12", "20")) +
    xlab(paste0("PC", int1, " ", eig[int1], "%")) +
    ylab(paste0("PC", int2, " ", eig[int2], "%")) +
    #ggtitle(paste0("PCA of ", DATA_NAME, ", n = ", count(info), ", SNPs = ", n.loc), 
    #        subtitle = title) +
    theme(panel.background = element_blank(),
          plot.background = element_blank(),
          panel.grid = element_blank(),
          legend.background = element_blank(),
          legend.key = element_blank(),
          panel.border = element_rect(colour = "grey", fill = NA, size = 1),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 8))
}
pca_plot_depth_all <- function(info, eig, title, int1, int2, colours, shapes, n.loc){
  ggplot(info, aes_(x = as.name(paste0("PC", int1)), y = as.name(paste0("PC", int2)),
                    colour = info$Loc, shape = info$Depth)) + 
    geom_point(alpha = 1/2, size = 5) +
    scale_color_manual(name = "Location", values = colours, breaks = c("WP", "CA", "SB", "CR", "SQ", "BR")) +
    scale_shape_manual(name = "Depth", values = shapes) +
    xlab(paste0("PC", int1, " ", eig[int1], "%")) +
    ylab(paste0("PC", int2, " ", eig[int2], "%")) +
    #ggtitle(paste0("PCA of ", DATA_NAME, ", n = ", count(info), ", SNPs = ", n.loc), 
    #        subtitle = title) +
    theme(plot.background = element_blank(),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          legend.background = element_blank(),
          legend.key = element_blank(),
          axis.line.x = element_line(colour = 'black', size = 0.5, linetype = 'solid'),
          axis.line.y = element_line(colour = 'black', size = 0.5, linetype = 'solid'))
}
pca_plot_depth_facet <- function(info, eig, title, int1, int2, colours){
  ggplot(info, aes_(x = as.name(paste0("PC", int1)), y = as.name(paste0("PC", int2)),
                    fill = info$Site)) + 
    geom_point(size = 2, colour = "black", shape = 21) +
    scale_fill_manual(name = "Site", values = colours, breaks = c("WP05", "WP10", "WP20",
                                                                  "CA05", "CA10", "CA20", 
                                                                  "SB05", "SB10", "SB20",
                                                                  "SQ12", "SQ20")) +
    xlab(paste0("PC", int1, " ", eig[int1], "%")) +
    ylab(paste0("PC", int2, " ", eig[int2], "%")) +
    theme(panel.background = element_blank(),
          plot.background = element_blank(),
          panel.grid = element_blank(),
          legend.background = element_blank(),
          legend.key = element_blank(),
          panel.border = element_rect(colour = "grey", fill = NA, size = 1),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 8)) + facet_wrap(~ Depth, nrow = 3, ncol = 1)
}
# Separate lm_pure fucntions
pca_plot_depth_lm_pure <- function(info, eig, title, int1, int2, colours, shapes, n.loc){
  ggplot(info, aes_(x = as.name(paste0("PC", int1)), y = as.name(paste0("PC", int2)),
                    colour = info$Loc, shape = info$Depth)) + 
    geom_point(alpha = 1/2, size = 5) +
    scale_color_manual(name = "Location", values = colours, breaks = c("WP", "CA", "SB", "SQ", "AL1", "AL2")) +
    scale_shape_manual(name = "Depth", values = shapes) +
    xlab(paste0("PC", int1, " ", eig[int1], "%")) +
    ylab(paste0("PC", int2, " ", eig[int2], "%")) +
    #ggtitle(paste0("PCA of ", DATA_NAME, ", n = ", count(info), ", SNPs = ", n.loc), 
    #        subtitle = title) +
    theme(plot.background = element_blank(),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          legend.background = element_blank(),
          legend.key = element_blank(),
          panel.border = element_rect(colour = "grey", fill = NA, size = 1),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 8))
}
pca_plot_depth_lm_pure_facet <- function(info, eig, title, int1, int2, colours, shapes, n.loc){
  ggplot(info, aes_(x = as.name(paste0("PC", int1)), y = as.name(paste0("PC", int2)),
                    fill = info$Site, shape = info$Depth)) + 
    geom_point(size = 2, shape = 21) +
    scale_fill_manual(name = "Location", values = colours, breaks = c("WP10", "WP20",
                                                                      "CA10", "CA20", 
                                                                      "SB10", "SB20",
                                                                      "SQ12", "SQ20",
                                                                      "AL1-15", "AL2-15", "AL1-50")) +
    scale_shape_manual(name = "Depth", values = shapes) +
    xlab(paste0("PC", int1, " ", eig[int1], "%")) +
    ylab(paste0("PC", int2, " ", eig[int2], "%")) +
    #ggtitle(paste0("PCA of ", DATA_NAME, ", n = ", count(info), ", SNPs = ", n.loc), 
    #        subtitle = title) +
    theme(plot.background = element_blank(),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          legend.background = element_blank(),
          legend.key = element_blank(),
          panel.border = element_rect(colour = "grey", fill = NA, size = 1),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 8)) + facet_wrap(~Depth, nrow = 4, ncol = 1)
}
# Admixture
pca_plot_admixture <- function(info, eig, title, int1, int2, n.loc) {
  ggplot(info, aes_(x = as.name(paste0("PC", int1)), y = as.name(paste0("PC", int2)),
                    fill = info$Clusters, shape = info$Loc)) + 
    geom_point(size = 2, colour = 'black') + 
    scale_fill_manual(name = 'Assignment',
                      values =  colours) +
    scale_shape_manual(name = "Loc", values = 21:(20 + length(levels(pca_info$Loc)))) +
    xlab(paste0("PC", int1, " ", eig[int1], "%")) +
    ylab(paste0("PC", int2, " ", eig[int2], "%")) +
    theme(plot.background = element_blank(),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          legend.background = element_blank(),
          legend.key = element_blank(),
          panel.border = element_rect(colour = "grey", fill = NA, size = 1),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 8)) + 
    facet_wrap(~Depth, nrow = length(levels(pca_info$Depth)), ncol = 1)
} 
pca_plot_admixture2 <- function(info, eig, title, int1, int2, n.loc) {
  ggplot(info, aes_(x = as.name(paste0("PC", int1)), y = as.name(paste0("PC", int2)),
                    fill = info$Clusters, shape = info$Loc)) + 
    geom_point(size = 2, colour = 'black') + 
    scale_fill_manual(name = 'Assignment',
                      values =  colours) +
    scale_shape_manual(name = "Loc", values = 21:(20 + length(levels(pca_info$Loc)))) +
    xlab(paste0("PC", int1, " ", eig[int1], "%")) +
    ylab(paste0("PC", int2, " ", eig[int2], "%")) +
    theme(plot.background = element_blank(),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          legend.background = element_blank(),
          legend.key = element_blank(),
          panel.border = element_rect(colour = "grey", fill = NA, size = 1),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 8)) + 
    facet_wrap(~Depth+Loc, nrow = length(levels(pca_info$Depth)), ncol = length(levels(pca_info$Loc)))
} 

## Arguments ==================================================
args = commandArgs(TRUE)
# defining the vcf file, the variable that is specified on the command line will be defined here.
# args[1] is the name of vcf file, with 4 parameters, first is the dataset "KP01-8", OR "ac"
# second is the step of the pipeline "1d", third is the max missing sites per indv used when 
# removing individuals "7" which means 70% or another dataset attritube, i.e., "nc" meaning no clones,
# and forth is max-missing data across sites "10", 10% max-missing.
VCF_NAME <-  as.character(args[1])

# number of PC axis to use, it can be any *even* integer, preferably 2-10
nf <- as.numeric(args[2]) 

# plot type, one of three: species, depth and admixture
plot_type <- args[3]

# Data variables to be used when naming files =====
DATA_PARAMS <- str_split_fixed(VCF_NAME, "_", n = 4) # split all of the vcf name
DATA_NAME <- str_c(DATA_PARAMS[1:4], collapse = "_")
# path for outputs
PATH = paste0("../results/pca/",DATA_PARAMS[1], "_stats/", VCF_NAME)

# Import data and formatting  =================================
# GENLIGHT object
GENLIGHT <- vcfR2genlight(read.vcfR(paste0("../data/", VCF_NAME, ".vcf")))
#GENLIGHT <- as.matrix(GENLIGHT[,1:100]) # sampling for quick troubleshooting
#GENLIGHT <- as.genlight(GENLIGHT) # sampling for quick troubleshooting

# popfile
pop <- read.table(paste0('../data/pop_', str_c(DATA_PARAMS[1:3], collapse = "_"), '.txt'), header = FALSE)
# popfile formatting
if (DATA_PARAMS[1] == "all-aga") {
  colnames(pop) <- c("Individual","Population")
pop  <- pop %>% 
  separate(Population, sep = "_", into = c("Species","Site"), remove = FALSE) %>% 
  separate(Site, sep = 2, into = c("Loc","Depth"), remove = FALSE)
pop <- pop %>% 
  unite(Site, Loc, Depth, sep = "", remove = FALSE) %>% 
  unite(Population, Species, Site, remove = FALSE)
} else {
  colnames(pop) <- c("Individual","Site")
  pop <- pop %>% 
    separate(Site, sep = 2, into = c("Loc", "Depth"), remove = FALSE)
}
pop$Depth[pop$Depth == "05"] = "5"

# PCA =============================
pca <- glPca(GENLIGHT, nf = nf)
n.loc <- GENLIGHT$n.loc
eig <- pca$eig
write.csv(eig, file = paste0(PATH, "_eig.csv"), quote = FALSE, row.names = FALSE)
pca_info <- sort_pca(pca, pop, nf)
title <- paste0(DATA_PARAMS[4],"% missing data")
clusters <- c("Clust1")

# Using PCA scores and fixing mislabelled samples ####
if (DATA_PARAMS[1] == "all-aga") {
  pca_info$Species[pca_info$Species == "AC"] = "AA"
  pca_info$Species[pca_info$Species == "LM"] = "AL"
  pca_info$Species[pca_info$Species == "HU"] = "AH"
  pca_info$Species[pca_info$Species == "L1"] = "AL"
  pca_info$Species[pca_info$Species == "L2"] = "AL"
  pca_info$Species[pca_info$Species == "G1"] = "AG"
  pca_info$Species[pca_info$Species == "G2"] = "AG"
  pca_info$Depth[pca_info$Individual == "KP0406_NA_NANA"] = "20"
  pca_info$Loc[pca_info$Individual == "KP0406_NA_NANA"] = "WP"
  pca_info$Depth[pca_info$Individual == "KP0673_NA_NANA"] = "20"
  pca_info$Loc[pca_info$Individual == "KP0673_NA_NANA"] = "SB"
  pca_info$Loc[pca_info$Individual == "KP0673_NA_NANA"] = "WP"
  pca_info$Depth[pca_info$Individual == "KP0844_UN_NANA"] = "10"
  pca_info$Loc[pca_info$Individual == "KP0844_UN_NANA"] = "WP"
  pca_info$Species[pca_info$Individual == "KP0844_UN_NANA"] = "AA"
  pca_info$Depth[pca_info$Individual == "KP0723_NA_NANA"] = "20"
  pca_info$Loc[pca_info$Individual == "KP0723_NA_NANA"] = "SB"
  pca_info$Species[pca_info$Individual == "KP0723_NA_NANA"] = "AL"
  pca_info$Species[pca_info$Individual == "KP1036_UN_CA10"] = "AH"
  pca_info$Species[pca_info$Individual == "KP1153_UN_CA20"] = "AH"
  pca_info$Species[pca_info$Individual == "KP1159_UN_CA20"] = "AH"
  pca_info$Species[pca_info$Individual == "KP1156_UN_CA20"] = "AA"
  pca_info$Species[pca_info$Individual == "KP1160_UN_CA20"] = "AL"
  pca_info$Loc[pca_info$Loc == "CK"] = "WP"
  pca_info$Loc[pca_info$Loc == "CS"] = "SQ"
  pca_info$Species[pca_info$Species == "CA"] = "AA"
  # mislabels agaricites, humilis, lamarcki
  mislabels_AA = c("KP1139_LM_CA20", "KP0361_LM_WP20", "KP0349_LM_WP20","KP0611_LM_SB20",
                   "KP0399_LM_WP20", "KP1146_LM_CA20", "KP0109_LM_SB10", "KP0583_LM_WP20",
                   "KP0244_LM_SQ20", "KP0565_LM_WP20", "KP0406_NA_NANA", "KP0673_NA_NANA")
  write.csv(mislabels_AA, file = "metadata/mislabels_AA.csv")
  for (a in 1:length(mislabels_AA)) {
    pca_info$Species[pca_info$Individual == mislabels_AA[a]] = "AA"}
  mislabels_AH = c("KP0639_LM_SB20", "KP0899_LM_WP10")
  write.csv(mislabels_AH, file = "metadata/mislabels_AH.csv")
  for (h in 1:length(mislabels_AH)) {
    pca_info$Species[pca_info$Individual == mislabels_AH[h]] = "AH"}
  mislabels_AL = c("KP1124_AC_CA20", "KP0557_AC_WP20", "KP0203_AC_SQ12", "KP1095_HU_CA20",
                   "KP0186_AC_SQ12", "KP0851_AC_WP10", "KP0605_AC_SB20", "KP0609_AC_SB20",
                   "KP0096_AC_SB10", "KP0765_AC_WP10", "KP0192_AC_SQ12", "KY0578_AC_WP20",
                   "KP0671_AC_SB20", "KP0108_AC_SB10")
  write.csv(mislabels_AL, file = "metadata/mislabels_AL.csv")
  for (l in 1:length(mislabels_AL)) {
    pca_info$Species[pca_info$Individual == mislabels_AL[l]] = "AL"}
  pca_info$Species <- as.factor(pca_info$Species)
  pca_info$Depth[pca_info$Depth == "12"] = "10-12"
  pca_info$Depth[pca_info$Depth == "10"] = "10-12"
} else if (DATA_PARAMS[1] == "ac") {
  pca_info$Depth[pca_info$Individual == "KP0406_NA_NANA"] = "20"
  pca_info$Loc[pca_info$Individual == "KP0406_NA_NANA"] = "WP"
  pca_info$Site[pca_info$Individual == "KP0406_NA_NANA"] = "WP20"
# pca_info$Depth[pca_info$Individual == "KP0673_NA_NANA"] = "20"
# pca_info$Loc[pca_info$Individual == "KP0673_NA_NANA"] = "SB" Not assigned to AC?
# pca_info$Site[pca_info$Individual == "KP0673_NA_NANA"] = "SB20"
  pca_info$Depth[pca_info$Individual == "KP0844_UN_NANA"] = "10"
  pca_info$Loc[pca_info$Individual == "KP0844_UN_NANA"] = "WP"
  pca_info$Site[pca_info$Individual == "KP0844_UN_NANA"] = "WP10"
  
  pca_info$Depth[pca_info$Depth == "12"] = "10-12"
  pca_info$Depth[pca_info$Depth == "10"] = "10-12"
  shapes = c(2, 6, 17, 4)
} else if (DATA_PARAMS[1] == "lm" & DATA_PARAMS[3] == "nc") {
  pca_info$Depth[pca_info$Individual == "KP0723_NA_NANA"] = "20"
  pca_info$Loc[pca_info$Individual == "KP0723_NA_NANA"] = "SB"
  pca_info$Site[pca_info$Individual == "KP0723_NA_NANA"] <- "SB20"
  pca_info$Depth[pca_info$Depth == "12"] = "10-12"
  pca_info$Depth[pca_info$Depth == "10"] = "10-12"
  shapes = c(6, 6, 17)
} else if (DATA_PARAMS[1] == "hu") {
  pca_info$Depth[pca_info$Depth == "12"] = "10-12"
  pca_info$Depth[pca_info$Depth == "10"] = "10-12"
} else if (DATA_PARAMS[1] == "lm" & (DATA_PARAMS[2] == "1d-pure" | DATA_PARAMS[3] == "nc-wnr")) {
  pca_info$Depth[pca_info$Individual == "KP0723_NA_NANA"] = "20"
  pca_info$Loc[pca_info$Individual == "KP0723_NA_NANA"] = "SB"
  pca_info$Site[pca_info$Individual == "KP0723_NA_NANA"] <- "SB20"
  pca_info$Loc[pca_info$Loc == "CK"] = "WP"
  pca_info$Loc[pca_info$Loc == "CS"] = "SQ"
  shapes = c(6, 6, 5, 17, 15)
  pca_info$Depth[pca_info$Depth == "12"] = "10-12"
  pca_info$Depth[pca_info$Depth == "10"] = "10-12"
  pca_info$Site[pca_info$Individual == "NR4671_L1_CS50"] = "AL1-50"
  pca_info$Site[pca_info$Individual == "NR5050_L2_CS15"] = "AL2-15"
  pca_info$Site[pca_info$Individual == "NR5053_L2_CS15"] = "AL2-15"
  pca_info$Site[pca_info$Individual == "NR5242_L1_CK50"] = "AL1-50"
  pca_info$Site[pca_info$Individual == "NR5533_L2_CR15"] = "AL2-15"
  pca_info$Site[pca_info$Individual == "NR5534_L1_CR15"] = "AL1-15"
}

write.csv(pca_info, file = paste0(PATH, "_pcscores", nf, ".csv"), quote = FALSE, row.names = TRUE)

# Save eigenvectors
eig <- pca$eig
eig <- round(eig[1:length(eig)]/sum(eig)*100, 2)
dat <- as.data.frame(eig[1:20])
colnames(dat) <- "eig"

ggplot(dat, aes(1:20, y = eig)) + geom_point() + geom_line() +
  scale_y_continuous(name = "% variation explained", breaks = c(seq(1,  (max(dat$eig) + 2), 2))) +
  scale_x_continuous(name = "PC axes", breaks = c(seq(1,  20, 2))) + theme_bw()
ggsave(paste0(PATH, "_eig_plot.pdf"), height = 5, width = 5, units = "cm", dpi = 400)

# Making plots ####
for (i in 1:(nf/2)) {
  j = seq(1, nf, 2)
  k = seq(2, nf, 2)
  if (plot_type == "depth") {
    if (DATA_PARAMS[1] == "all-aga") {
      colours <- c("#EC9C9D", "#C4AABB", "#9BB7DB", "#FFE1A6", "#8BD3D7", "#88C762")
      shapes = c(2, 6, 6, 5, 17, 15)
      pca_info$Depth = factor(pca_info$Depth, levels = c("5", "10-12", "20", "15", "50"))
      pca_info$Loc <- as.factor(pca_info$Loc)
      pca_info$Site <- as.factor(pca_info$Site)
      pca.plot <- pca_plot_depth_all(pca_info, eig, title, j[i], k[i], colours, shapes, n.loc)
      ggsave(pca.plot, file = paste0(PATH,"_depth_facet_pca-x", j[i], k[i], ".pdf"), width = 8, height = 12, units = "cm", dpi = 400) # TODO: might need to alter this!!!
    } else if (DATA_PARAMS[1] == "lm" & (DATA_PARAMS[2] == "1d-pure" | DATA_PARAMS[3] == "nc-wnr")) {
      colours <- c("#EC9C9D", "#AA5459",
                   "#C4AABB", "#8C6C82",
                   "#9BB7DB", "#5179A6",
                   "#8BD3D7", "#74A8B6",
                   "#6F1242", "#6652A2",
                   "#4E1136")
      pca_info$Loc[pca_info$Individual == "NR4671_L1_CS50"] = "AL1"
      pca_info$Loc[pca_info$Individual == "NR5050_L2_CS15"] = "AL2"
      pca_info$Loc[pca_info$Individual == "NR5053_L2_CS15"] = "AL2"
      pca_info$Loc[pca_info$Individual == "NR5242_L1_CK50"] = "AL1"
      pca_info$Loc[pca_info$Individual == "NR5533_L2_CR15"] = "AL2"
      pca_info$Loc[pca_info$Individual == "NR5534_L1_CR15"] = "AL1"
      pca_info$Depth = factor(pca_info$Depth, levels = c("10-12", "20", "15", "50"))
      pca_info$Loc <- as.factor(pca_info$Loc)
      pca_info$Site <- as.factor(pca_info$Site)
      #pca.plot <- pca_plot_depth_lm_pure(pca_info, eig, title, j[i], k[i], colours, shapes, n.loc)
      pca.plot <- pca_plot_depth_lm_pure_facet(pca_info, eig, title, j[i], k[i], colours)
      ggsave(pca.plot, file = paste0(PATH,"_depth_facet_pca-x", j[i], k[i], ".pdf"), width = 8, height = 12, units = "cm", dpi = 400)
      } else {
      colours <- c("#FBD1D7", "#EC9C9D", "#AA5459",
                   "#E2DCE0", "#C4AABB", "#8C6C82",
                   "#C0DFF5", "#9BB7DB", "#5179A6",
                   "#8BD3D7", "#74A8B6")
      pca_info$Depth <- factor(pca_info$Depth, levels = c("5", "10-12", "20"))
      pca_info$Loc <- as.factor(pca_info$Loc)
      pca_info$Site <- as.factor(pca_info$Site)
      #pca.plot <- pca_plot_depth(pca_info, eig, title, j[i], k[i], colours, shapes, n.loc)
      pca.plot <- pca_plot_depth_facet(pca_info, eig, title, j[i], k[i], colours)
      ggsave(pca.plot, file = paste0(PATH,"_depth_facet_pca-x", j[i], k[i], ".pdf"), width = 8, height = 12, units = "cm", dpi = 400)
    }
    plotPATH = paste0(PATH,"_depth_pca", j[i], k[i], ".pdf")
    
    # Admixture ============================
  } else if (plot_type == "admixture") {
    if (DATA_PARAMS[1] == "all-aga") {
      admix_k <- 4
      taxa_k <- 3
      colours <- c("#2D4E24", #agaricites1
                   "#6F1242", #lamarcki1
                   "#B06327", #humilis1
                   "#6652A2") #lamarcki2
      threshold <- 0.8
    } else if (DATA_PARAMS[1] == "ac") {
      admix_k <- nf
      taxa_k <- 2
      colours <- c("#2D4E24", "#8FC73E", "#43BB93", "#347F66")
      loc_levels <- c("WP", "CA", "SB", "SQ")
      depth_levels <- c("5", "10-12", "20")
      threshold <- 0.90
    } else if (DATA_PARAMS[1] == "hu") {
      admix_k <- 4
      taxa_k <- 3
      colours <- c("#B06327", "#FAA41A", "#F05123", "#821F27")
      loc_levels <- c("WP", "CA", "SB", "SQ")
      depth_levels <- c("5", "10-12", "20")
      threshold <- 0.8
    } else if (DATA_PARAMS[1] == "lm" & DATA_PARAMS[3] == "nc") {
      admix_k <- 2
      taxa_k <- 2
      colours <- c("#6F1242", "#6652A2")
      loc_levels <- c("WP", "CA", "SB", "SQ")
      depth_levels <- c("10-12", "20")
      threshold <- 0.8
    } else if (DATA_PARAMS[1] == "lm" & (DATA_PARAMS[2] == "1d-pure" | DATA_PARAMS[3] == "nc-wnr")) {
      admix_k <- 2
      taxa_k <- 2
      colours <- c("#6F1242", "#6652A2")
      loc_levels <- c("WP", "CA", "SB", "SQ", "CR")
      depth_levels <- c("10-12", "20", "15", "50")
      threshold <- 0.8
    } else {
      admix_k <- 2
      taxa_k <- 2
      colours <- c("red", "blue")
    }
    pca_info$Loc <- factor(pca_info$Loc, levels = loc_levels)
    pca_info$Depth <- factor(pca_info$Depth, levels = depth_levels)
    clusters <- append(clusters, paste0("Clust", 2:admix_k))
    print(paste0(DATA_PARAMS[1], " dataset set to k=", admix_k, " for admixture assignment and taking from sortedQ"))
    Q <- read.table(paste0("../2b - Admixture/", DATA_PARAMS[4], "percent_", DATA_PARAMS[1], "_",
                           DATA_PARAMS[2], "_", DATA_PARAMS[3], "/sortedQ/", DATA_PARAMS[1], "_2bi_", DATA_PARAMS[2],
                           "_", DATA_PARAMS[3], "_",DATA_PARAMS[4], ".", admix_k, ".Q"), header = TRUE)
    taxaQ <- read.table(paste0("../2b - Admixture/", DATA_PARAMS[4], "percent_", DATA_PARAMS[1], "_",
                               DATA_PARAMS[2], "_", DATA_PARAMS[3], "/sortedQ/", DATA_PARAMS[1], "_2bi_", DATA_PARAMS[2],
                               "_", DATA_PARAMS[3], "_",DATA_PARAMS[4], ".", taxa_k, ".Q"), header = TRUE)
    source('scripts/func-k.R')
    admix_data <- cbind(pop, taxaQ)
    sorted_clusters <- sort_clusters(Q, admix_k)
    admix_data <- choose_Qprop(admix_data, taxa_k, threshold)
    pca_info$Clusters <- sorted_clusters
    pca_info$Taxa <- admix_data$Taxa
    for (column in c("Site", "Clusters", "Taxa")) {
      pca_info[, column] <- as.factor(pca_info[, column])
    }
    pca_info$Q <- Q
    print("Writing sorted admix assignments to data ;)")
    write.csv(pca_info, paste0("../../data/", VCF_NAME, "_", admix_k, ".csv"), quote = FALSE,
              row.names = FALSE)
    pca.plot <- pca_plot_admixture(pca_info, eig, title, j[i], k[i], n.loc)
    plotPATH = paste0(PATH,"_admixture_pca", j[i], k[i], ".pdf")
    ggsave(pca.plot, file = plotPATH, width = 8, height = 12, units = "cm", dpi = 400)
  } else {
    # species
    colours = c("#2D4E24", "#B06327", "#6F1242", "#066991", "grey")
    pca_info$Depth = factor(pca_info$Depth, levels = c("5", "10-12", "20", "15", "50"))
    pca_info$Loc <- as.factor(pca_info$Loc)
    pca_info$Site <- as.factor(pca_info$Site)
    pca.plot <- pca_plot_species(pca_info, eig, title, j[i], k[i], colours, n.loc)
    plotPATH = paste0(PATH,"_species_pca-x", j[i], k[i], ".pdf")
    ggsave(pca.plot, file = plotPATH, width = 7, height = 20, units = "cm", dpi = 400)
  }
}


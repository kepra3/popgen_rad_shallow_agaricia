# Title: Convert vcf to genepop and use gene pop
# Author: Katharine Prata
# Date created: 19/4/22

# Required packages
library(vcfR)
library(adegenet)
library(genepop)
library(tidyr)
library(dplyr)

# Functions  ==================================================
standardise_plot_dist <- function(metadata, plot) {
  metadata$x[metadata$Site == plot] <- metadata$x[metadata$Site == plot] + abs(min(metadata$x[metadata$Site == plot]))
  metadata$y[metadata$Site == plot] <- metadata$y[metadata$Site == plot] + abs(min(metadata$y[metadata$Site == plot]))
  return(metadata)
}
genind2genepop <- function(taxa.genind, taxa.pop, category) {
  if (nchar(category) < 3) {
    taxa.genind@pop <- as.factor(taxa.pop$Depth)
    taxa.pop <- taxa.pop[taxa.pop$Depth == category, ]
    taxa.genind.pop <- taxa.genind[taxa.genind@pop == category]
  } else if (nchar(category) == 4) {
    taxa.genind@pop <- as.factor(taxa.pop$Site)
    taxa.pop <- taxa.pop[taxa.pop$Site == category, ]
    taxa.genind.pop <- taxa.genind[taxa.genind@pop == category]
  } else {
    taxa.genind.pop <- taxa.genind
  }

  write.csv(taxa.pop, file = paste(vcf_name, taxa, category, ".csv", sep = "_"), quote = FALSE)
  
  # Turn taxa.genind into dataframe
  df <- genind2df(taxa.genind.pop, usepop = FALSE)
  
  # Convert 0 to 01 and 1 to 02 and NA's to 0000
  mat <- as.matrix(df)
  mat <- apply(mat, FUN = gsub, MARGIN = 2, pattern = "0", replacement = "A")
  mat <- apply(mat, FUN = gsub, MARGIN = 2, pattern = "1", replacement = "T")
  mat <- apply(mat, FUN = gsub, MARGIN = 2, pattern = "A", replacement = "01")
  mat <- apply(mat, FUN = gsub, MARGIN = 2, pattern = "T", replacement = "02")
  mat[is.na(mat)] = "0000"
  
  # Add the distances and individual name, and make distances positive
  xy <- paste(taxa.pop$x, taxa.pop$y, paste0(indNames(taxa.genind.pop), ","))
  mat <- cbind(xy, mat)
  
  # Insert a Pop row between each population
  # Double all rows
  mat <- mat[rep(1:nrow(mat), 1, each = 2), ]
  
  # Replace all duplicates with blank cells
  mat[c(seq(2, dim(mat)[1], by = 2)), ] <- ""
  mat[c(seq(2, dim(mat)[1] - 1, by = 2)), 1] <- "Pop"
  mat <- mat[-nrow(mat),]
  
  # Genepop header
  file_date <- format(Sys.time(), "%Y%m%d@%H%M") # date and time
  header <- c(paste("Genepop file format", file_date), rep("", ncol(mat) - 1))
  loc_names <- c(paste(locNames(taxa.genind), collapse = ","), rep("", ncol(mat) - 1))
  popline <- c("Pop", rep("", ncol(mat) - 1))
  first_lines <- rbind(header, loc_names, popline)
  mat <- rbind(first_lines, mat)
  return(mat)
}

# Arguments ==================================================
args = commandArgs(TRUE)
vcf_name <- args[1]
taxa <- args[2]
category <- args[3]
scale <- args[4]

## Example use
#vcf_name <- "ac_3b_nc_20"
#taxa <- "AA1"
#category <- "all"
#scale <- "within"

# Import data  ==================================================
setwd("~/Dropbox/agaricia_project_2019/shalo_ag/gen_project/3 - Spatial/3b - Within & between plots/spatial-agaricia")
genind <- vcfR2genind(read.vcfR(paste0("data/", vcf_name, ".vcf")))
metadata <- read.csv("../../all_annotations_X_HORIZ_parallel.txt", sep = "\t")

# Renaming clusters #### (need to check whether these files are the best check quick.pca)

if (vcf_name == "ac_3b_nc_20") {
  clusters <- read.csv("data/ac_1div_nc_20_4.csv")
  clusters$Taxa[clusters$Taxa == "Clust2"] = "AA2"
  clusters$Taxa[clusters$Taxa == "Clust1"]  = "AA1"
} else if (vcf_name == "hu_3b_nc_20") {
  clusters <- read.csv("data/hu_1div_nc_20_4.csv")
  clusters$Taxa[clusters$Taxa == "Clust1"] = "AH1"
  clusters$Taxa[clusters$Taxa == "Clust2"]  = "AH2"
  clusters$Taxa[clusters$Taxa == "Clust3"] = "AH3"
} else if (vcf_name == "lm_3b_nc_20") {
  clusters <- read.csv("data/lm_1div_nc-wnr_20_2.csv")
  clusters$Taxa[clusters$Taxa == "Clust1"] = "AL1"
  clusters$Taxa[clusters$Taxa == "Clust2"]  = "AL2"
}

# sort metadata
colnames(metadata)[colnames(metadata) == "X"] <- "Individual"
metadata <- metadata %>% separate(Individual, into = c("Sample", "Species", "Site"), sep = "_", remove = FALSE) %>% 
    separate(Site, into = c("Loc", "Depth"), sep = 2, remove = FALSE) %>% 
  unite(Individual, Sample, Species, Site, remove = FALSE) %>% 
  select(Individual, Sample, Species, Site, Depth, x, y, z)

metadata$Depth[metadata$Depth == "12"] = "10"
metadata$Site[metadata$Site == "SQ12"] = "SQ10"

# Standardise distances  ==================================================
if (nchar(category) <= 3) {
  # distances from https://latlongdata.com/distance-calculator/ the Haversine formula
  # distance between WP and CA
  dist_CA <- 17828
  # distance between WP and SB
  dist_SB <- 31509
  # distance between WP and SQ
  dist_SQ <- 42880 }
if (nchar(category) == 2) {
  metadata <- standardise_plot_dist(metadata, paste0("WP", category))
  metadata <- standardise_plot_dist(metadata, paste0("SB", category))
  metadata <- standardise_plot_dist(metadata, paste0("CA", category))
  metadata$y[metadata$Site == paste0("SB", category)] <- metadata$y[metadata$Site == paste0("SB", category)] + dist_SB
  metadata$y[metadata$Site == paste0("CA", category)] <- metadata$y[metadata$Site == paste0("CA", category)] + dist_CA
  if (category != "05") {
  metadata <- standardise_plot_dist(metadata, paste0("SQ", category))
  metadata$y[metadata$Site == paste0("SQ", category)] <- metadata$y[metadata$Site == paste0("SQ", category)] + dist_SQ
  }
}

#  For all dataset ========================================================
if (category == "all") {
  ## Add distances b/w depths & sites #####
  dist_5_10 <- 30
  # 1 - 6 cattle tags
  SBdist_10_20 <- sum(c(13.678, 14.800, 14.777, 15.191, 17.111, 17.083))/6
  WPdist_10_20 <- sum(c(21.12, 19.27, 18.45,18.65, 18.85, 16.88))/6
  dist_10_20 <- sum(SBdist_10_20, WPdist_10_20)/2
  for (depth in c("05", "10", "20")) {
    # can do this in a simpler way but won't worry about it for now.
  metadata <- standardise_plot_dist(metadata, paste0("WP", depth))
  metadata <- standardise_plot_dist(metadata, paste0("SB", depth))
  metadata <- standardise_plot_dist(metadata, paste0("CA", depth))
  metadata$x[metadata$Site == paste0("SB", depth)] <- metadata$x[metadata$Site == paste0("SB", depth)] + dist_SB
  metadata$x[metadata$Site == paste0("CA", depth)] <- metadata$x[metadata$Site == paste0("CA", depth)] + dist_CA
  if (depth != "05") {
    metadata <- standardise_plot_dist(metadata, paste0("SQ", depth))
    metadata$x[metadata$Site == paste0("SQ", depth)] <- metadata$x[metadata$Site == paste0("SQ", depth)] + dist_SQ
  }}
  # Depth distances
  metadata$y[metadata$Depth == "10"] = metadata$y[metadata$Depth == "10"] - dist_5_10
  metadata$y[metadata$Site == "WP20"] = metadata$y[metadata$Site == "WP20"] - dist_5_10 - WPdist_10_20
  metadata$y[metadata$Site == "SB20"] = metadata$y[metadata$Site == "SB20"] - dist_5_10 - SBdist_10_20
  metadata$y[metadata$Site == "CA20"] = metadata$y[metadata$Site == "CA20"] - dist_5_10 - dist_10_20
  metadata$y[metadata$Site == "SQ20"] = metadata$y[metadata$Site == "SQ20"] - dist_5_10 - dist_10_20
  }

# Match and subset datasets  ================================================== 
pop <- metadata[metadata$Individual %in% indNames(genind),]
pop <- pop[order(pop$Individual),]
subset.clusters <- clusters[clusters$Individual %in% pop$Individual,]
identical(subset.clusters$Individual, pop$Individual)
pop$Clusters <- subset.clusters$Taxa
pop <- pop[!is.na(pop$Clusters), ]
subset.genind <- genind[indNames(genind) %in% pop$Individual,]
identical(indNames(subset.genind), pop$Individual)
subset.genind@pop <- as.factor(pop$Clusters)
taxa.genind <- subset.genind[subset.genind@pop == taxa]
taxa.pop <- pop[pop$Clusters == taxa,]
identical(indNames(taxa.genind), taxa.pop$Individual)

taxa.pop$Site <- as.factor(taxa.pop$Site)
summary(taxa.pop$Site)

# Convert 2 genepop and export file  ==================================================
genepop.mat <- genind2genepop(taxa.genind, taxa.pop, category)
write.table(genepop.mat, file = paste0("data/", taxa,"_", category, ".genepop.txt"), quote = FALSE, col.names = FALSE, row.names = FALSE)

# Gene pop  ==================================================
if (scale == "all") {
  minD = 1e-03
  maxD = 1e+05
} else if (scale == "within") {
  minD = 1e-03
  maxD = 1e+02
} else if (scale == "between") {
  minD = 1e+02
  maxD = 1e+05
}
ibd(paste0("data/", taxa, "_", category, ".genepop.txt"),
  outputFile = paste0("results/", taxa, "/", taxa, "_", category, ".", scale, ".results.txt"),
  settingsFile = "",
  dataType = "Diploid",
  statistic = "a",
  geographicScale = "2D",
  CIcoverage = 0.95,
  testPoint = 0,
  minimalDistance = minD,
  maximalDistance = maxD,
  mantelPermutations = 1000,
  mantelRankTest = FALSE,
  verbose = interactive())


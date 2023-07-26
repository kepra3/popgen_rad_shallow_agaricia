# Title: Convert vcf to genepop and use gene pop
# Author: Katharine Prata
# Date created: 19/4/22

# Required packages
# R v4.2.0
library(vcfR) # v1.12.0
library(adegenet) # v2.1.7
library(genepop) # v1.1.7
library(tidyr) # v1.2.0
library(dplyr) # v1.0.9

# Functions  ==================================================
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

  write.csv(taxa.pop, file = paste0("../data/", paste(vcf_name, taxa, category, ".csv", sep = "_")), quote = FALSE)
  
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
taxa <- args[1]
category <- args[2]
scale <- args[3]

## Example use
#taxa <- "AA1"
#category <- "all" #(depth category)
#scale <- "within"

# Import data  ==================================================
if (taxa == "AA1" | taxa == "AA2") {
  vcf_name = "ac_3b_nc_20"
} else if (taxa == "AH1" | taxa == "AH2" | taxa == "AH3") {
  vcf_name = "hu_3b_nc_20"
} else if (taxa == "AL1" | taxa == "AL2") {
  vcf_name = "lm_3b_nc_20"
}

genind <- vcfR2genind(read.vcfR(paste0("../data/", vcf_name, ".vcf")))
metadata <- read.csv("../data/all_annotations_X_HORIZ_parallel_XYZ_adjusted.txt", sep = ",")

# Renaming clusters #### (need to check whether these files are the best check quick.pca)

if (vcf_name == "ac_3b_nc_20") {
  clusters <- read.csv("../data/ac_1div_nc_20_4.csv")
  clusters$Taxa[clusters$Taxa == "Clust2"] = "AA2"
  clusters$Taxa[clusters$Taxa == "Clust1"]  = "AA1"
} else if (vcf_name == "hu_3b_nc_20") {
  clusters <- read.csv("../data/hu_1div_nc_20_4.csv")
  clusters$Taxa[clusters$Taxa == "Clust1"] = "AH1"
  clusters$Taxa[clusters$Taxa == "Clust2"]  = "AH2"
  clusters$Taxa[clusters$Taxa == "Clust3"] = "AH3"
} else if (vcf_name == "lm_3b_nc_20") {
  clusters <- read.csv("../data/lm_1div_nc-wnr_20_2.csv")
  clusters$Taxa[clusters$Taxa == "Clust1"] = "AL1"
  clusters$Taxa[clusters$Taxa == "Clust2"]  = "AL2"
}

# sort metadata
metadata <- metadata %>% separate(Individual, into = c("Sample", "Species", "Site"), sep = "_", remove = FALSE) %>% 
    separate(Site, into = c("Loc", "Depth"), sep = 2, remove = FALSE) %>% 
  unite(Individual, Sample, Species, Site, remove = FALSE) %>% 
  select(Individual, Sample, Species, Site, Depth, x, y, z)

metadata$Depth[metadata$Depth == "12"] = "10"
metadata$Site[metadata$Site == "SQ12"] = "SQ10"


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
write.table(genepop.mat, file = paste0("../data/", taxa,"_", category, ".genepop_copy.txt"), quote = FALSE, col.names = FALSE, row.names = FALSE)

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
ibd(paste0("../data/", taxa, "_", category, ".genepop_copy.txt"),
  outputFile = paste0("../results/ibd/", taxa, "/", taxa, "_", category, ".", scale, ".results_copy.txt"),
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


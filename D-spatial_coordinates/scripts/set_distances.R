# Title: Set distances between plots
# Created by: Katharine Prata

# Packages
# R v4.2.0
library(tidyverse) # v1.3.1

# Functions
clean_df <- function(data) {
  data <- data %>% separate(X, into = c("Sample", "Sp", "Site"), sep = "_", remove = FALSE) %>% 
    dplyr::select(X, Site, x, y, z)
  columns <- c("Individual", "Site", "x", "y", "z")
  colnames(data) <- columns
  return(data)
}
standardise_plot_dist <- function(metadata) {
  for (site in c("WP", "CA", "SB")) {
    for (depth in c("05", "10", "20")) {
    metadata$x[metadata$Site == paste0(site, depth)] <- metadata$x[metadata$Site == paste0(site, depth)] + abs(min(metadata$x[metadata$Site == paste0(site, depth)]))
    metadata$y[metadata$Site == paste0(site, depth)] <- metadata$y[metadata$Site == paste0(site, depth)] + abs(min(metadata$y[metadata$Site == paste0(site, depth)]))
    metadata$z[metadata$Site == paste0(site, depth)] <- metadata$z[metadata$Site == paste0(site, depth)] + (-as.numeric(depth) - mean(metadata$z[metadata$Site == paste0(site, depth)]))
    }}
  for (depth in c("12", "20")) {
    site = "SQ"
    metadata$x[metadata$Site == paste0(site, depth)] <- metadata$x[metadata$Site == paste0(site, depth)] + abs(min(metadata$x[metadata$Site == paste0(site, depth)]))
    metadata$y[metadata$Site == paste0(site, depth)] <- metadata$y[metadata$Site == paste0(site, depth)] + abs(min(metadata$y[metadata$Site == paste0(site, depth)]))
    metadata$z[metadata$Site == paste0(site, depth)] <- metadata$z[metadata$Site == paste0(site, depth)] + (-as.numeric(depth) - mean(metadata$z[metadata$Site == paste0(site, depth)]))
  }
  return(metadata)
}

# Arguments
args <- commandArgs(TRUE)
annotation_file <- args[1]


# Read data
metadata <- read.csv(paste0("../results/", annotation_file, ".txt"), sep = ",")
metadata <- clean_df(metadata)
metadata <- na.omit(metadata)
metadata <- standardise_plot_dist(metadata)

metadata <- metadata %>% separate(Site, into = c("Loc", "Depth"), sep = 2, remove = FALSE)
metadata$Depth[metadata$Depth == "12"] = "10"
metadata$Site[metadata$Site == "SQ12"] = "SQ10"
# Standardise distances  ==================================================
# distances from https://latlongdata.com/distance-calculator/ the Haversine formula
# distance between WP and CA
dist_CA <- 17828
# distance between WP and SB
dist_SB <- 31509
# distance between WP and SQ
dist_SQ <- 42880

#  For all dataset ========================================================
  ## Add distances b/w depths & sites #####
dist_5_10 <- 30
# 1 - 6 cattle tags
SBdist_10_20 <- sum(c(13.678, 14.800, 14.777, 15.191, 17.111, 17.083))/6
WPdist_10_20 <- sum(c(21.12, 19.27, 18.45, 18.65, 18.85, 16.88))/6
dist_10_20 <- sum(SBdist_10_20, WPdist_10_20)/2
 for (depth in c("05", "10", "20")) {
    metadata$x[metadata$Site == paste0("SB", depth)] <- metadata$x[metadata$Site == paste0("SB", depth)] + dist_SB
    metadata$x[metadata$Site == paste0("CA", depth)] <- metadata$x[metadata$Site == paste0("CA", depth)] + dist_CA
    if (depth != "05") {
      metadata$x[metadata$Site == paste0("SQ", depth)] <- metadata$x[metadata$Site == paste0("SQ", depth)] + dist_SQ
    }}
# Depth distances
metadata$y[metadata$Depth == "10"] = metadata$y[metadata$Depth == "10"] - dist_5_10
metadata$y[metadata$Site == "WP20"] = metadata$y[metadata$Site == "WP20"] - dist_5_10 - WPdist_10_20
metadata$y[metadata$Site == "SB20"] = metadata$y[metadata$Site == "SB20"] - dist_5_10 - SBdist_10_20
metadata$y[metadata$Site == "CA20"] = metadata$y[metadata$Site == "CA20"] - dist_5_10 - dist_10_20
metadata$y[metadata$Site == "SQ20"] = metadata$y[metadata$Site == "SQ20"] - dist_5_10 - dist_10_20

write.csv(metadata, file = paste0("../results/", annotation_file, "_XYZ_adjusted_copy.txt"), row.names = FALSE, quote = FALSE)

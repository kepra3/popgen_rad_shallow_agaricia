# Title: Genotypes vs. clone distances distributions
# Author: Katharine Prata

# Packages
# R v4.2.0
library(ggplot2) # v3.4.0
library(spaa) # v0.2.2
library(tidyr) # v1.2.0

# Functions
list_df <- function(metadata, cluster) {
  clust.metadata <- metadata[metadata$Clusters == cluster,]
  list <- dist2list(dist(clust.metadata[,c(7:9)]))
  list <- list[!duplicated(list$value), ]
  list$col <- as.character(list$col)
  list$row <- as.character(list$row)
  list$col <- clust.metadata$Individual[match(list$col, rownames(clust.metadata))]
  list$row <- clust.metadata$Individual[match(list$row, rownames(clust.metadata))]
  list$col_group <- clust.metadata$clone_group[match(list$col, clust.metadata$Individual)]
  list$row_group <- clust.metadata$clone_group[match(list$row, clust.metadata$Individual)]
  list <- list[list$value != 0,]
  list <- list %>% unite(group_pair, col_group, row_group, sep = "-", remove = FALSE)
  list$group_pair <- as.factor(list$group_pair)
  list$group_stat <- "NA"
  list[list$col_group == list$row_group, 7] <- "within"
  list[list$col_group != list$row_group, 7] <- "between"
  # Sort list from 0 - high
  list <- list[order(list$value), ]
  return(list)
}

# Variables
args = commandArgs(TRUE)
taxa <- args[1]

# Colours
colours <- c("#274e13", "#8fce00", "#b06100", "#ffa500", "#ff4e00", "#6f0641", "#7223f0")

if (taxa == "AA1") {
  colour <- colours[1]
} else if (taxa == "AA2") {
  colour <- colours[2]
} else if (taxa == "AH1") {
  colour <- colours[3]
} else if (taxa == "AH2") {
  colour <- colours[4]
} else if (taxa == "AH3") {
  colour <- colours[5] 
} else if (taxa == "AL1") {
  colour <- colours[6]
} else if (taxa == "AL2") {
  colour <- colours[7]
} else {
  print("taxa not chosen")
}

# Files
metadata <- read.table("../results/clone_dist/all-aga_spatgeo.txt", header = TRUE, stringsAsFactors = TRUE)

list <- list_df(metadata, taxa)
list[list$group_pair == "non-clone-non-clone",]$group_stat <- "between"

list_wloc <- list[list$value < 50, ]

list_wloc <- list_wloc %>% 
  separate(col = col, into = c("x", "col_loc"), remove = FALSE, sep = -4) %>%
  separate(col = row, into = c("y", "row_loc"), remove = FALSE, sep = -4)

list_wloc <- list_wloc[,c(-2,-5)]

list_wloc$match <- FALSE # differences between different plots
list_wloc[list_wloc$col_loc == list_wloc$row_loc,]$match <- TRUE
list_wloc <- list_wloc[list_wloc$match == TRUE,]

# add n = for each site
ggplot(list_wloc, aes(x = group_stat, y = value, colour = taxa)) +
  geom_jitter(size = 2, shape = 21) + theme_classic() +
  scale_color_manual(values = colour) +
  coord_flip()

ggplot(list_wloc, aes(x = group_stat, y = value)) +
  geom_boxplot() + theme_classic() + geom_jitter(fill = colour, shape = 21,
                                                 colour = "black") +
  scale_fill_manual(values = colour)  +
  ylab("distances (m)") +
  xlab("Genotypes") +
  facet_wrap(~col_loc)

### same format as the clone distances but not including clones

if (taxa == "AA2") {
  ggplot(list_wloc[list_wloc$group_stat == "between",], aes(value)) +
    geom_histogram(breaks = seq(0, 30, 1), fill = colour, colour = "black") +
    theme_classic() + xlab("distance between genotypes (m)") + ylab("Frequnecy") +
    xlim(c(0,30)) +
    scale_y_continuous(limits = c(0, 1250), expand = c(0, 0), breaks = seq(0, 1500, 250))
  ggsave(paste0("../results/clone_dist/", taxa, "_genotype_distrib.pdf"), height = 10, width = 5, units = "cm")
} else {
  ggplot(list_wloc[list_wloc$group_stat == "between",], aes(value)) +
    geom_histogram(breaks = seq(0, 30, 1), fill = colour, colour = "black") +
    theme_classic() + xlab("distance between genotypes (m)") + ylab("Frequnecy") +
    xlim(c(0,30)) +
    scale_y_continuous(limits = c(0, 50), expand = c(0, 0), breaks = seq(0, 50, 10))
  ggsave(paste0("../results/clone_dist/", taxa, "_genotype_distrib.pdf"), height = 5, width = 5, units = "cm")
}



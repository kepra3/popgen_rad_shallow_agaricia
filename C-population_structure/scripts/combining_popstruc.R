# Combining Plots
# Author: Katharine Prata
# Date created: 21/7/21
# Last edit: 20/2/23
# Description: Combination population structure plots

## Dependencies  ===================================================
# R v4.2.0
library(tidyverse) # v1.3.1
library(adegenet) # v2.1.7
library(vcfR) # v1.12.0
library(ggtree) # v3.4.0
library(treeio) # v1.20.0
library(ggstance) # v0.3.5
library(RColorBrewer) # v1.1-3
library(phytools) # v1.0-3
library(ggplotify) # v0.1.0

# Functions ==========================================
organise_data <- function(qfile, popfile, clusters){
  # where q file is your STRUCURE or Admixture output "Q" file
  # where popfile is first columns of your sample names and next column your population sites
  # clusters is a list of your clusters, i.e., clusters = c("Clust1", "Clust2, ... etc)
  data <- as.data.frame(cbind(popfile,qfile))
  clusters <- expand.grid(popfile[,1], clusters)
  n <- ncol(data)
  Assignment_t <- as.matrix(data[,3:n])
  Assignment <- matrix(Assignment_t, ncol = 1)
  plotdata <- cbind(clusters,Assignment)
  plotdata <- cbind(plotdata,popfile[,2])
  colnames(plotdata) <- c("Individuals","Clusters","Assignment","pop")
  return(plotdata)
}
basic_tree <- function(tree, pop_meta, k_selection) {
  rooted_tree <- midpoint.root(tree)
  # Set groups by using tibble
  tibble <- as_tibble(rooted_tree)
  col = which(colnames(pop_meta) == k_selection)
  group_tibble <- tibble(label = pop_meta$Individual,
                         group = pop_meta[,col])
  # Set tree info
  tree_tibble_group <- full_join(tibble, group_tibble, by = 'label')
  tree_tibble_group
  tree_x <- as.phylo(tree_tibble_group)
  # Basic tree plot
  tree_plot <- ggtree(tree_x) + 
    geom_tippoint(aes(colour = tree_tibble_group$group), size = 1) #+ 
  #geom_tiplab(aes(colour = tree_tibble_group$group), size = 1, align = TRUE)
  return(tree_plot)
}
advanced_tree <- function(tree_plot, pop_meta, colours, clusters, k) {
  tree_plot + theme_tree2() +
    geom_facet(panel = "Admixture", data = admix, geom = geom_barh,
               mapping = aes(x = Assignment, fill = Clusters), stat = 'identity') + 
    scale_fill_manual(name = "Cluster",
                      values = colours) +
    # Note: Breaks here controls the order of the colour values
    # before it was in alphanumeric order...
    scale_color_manual(name = "Assignment", values = colours,
                       label = clusters, breaks = clusters) +
    theme(legend.position = "none",
          strip.background.x = element_blank(),
          strip.text.x = element_text(size = 12),
          panel.border = element_rect(fill = "NA", linewidth = 0.2),
          panel.background = element_blank())
}
advanced_tree_loc_depth <- function(tree_plot, pop_meta, colours, clusters, k) {
  tree_plot + theme_tree2() +
    geom_facet(panel = "Admixture", data = admix, geom = geom_barh,
               mapping = aes(x = Assignment, fill = Clusters), stat = 'identity') + 
    geom_facet(panel = "Location", data = pop_meta, geom = geom_point,
               mapping = aes(x = Loc_value, colour = Loc), size = 1) +
    geom_facet(panel = "Depth", data = pop_meta, geom = geom_point,
               mapping = aes(x = Depth_values, colour = Site), size = 1) +
    scale_fill_manual(name = "Cluster",
                      values = colours[(length(levels(pop_meta$Site)) + 1):(length(levels(pop_meta$Site)) + k)]) +
    # Note: Breaks here controls the order of the colour values
    # before it was in alphanumeric order...
    scale_color_manual(name = "Assignment", values = colours,
                       label = c(levels(pop_meta$Site),
                                 clusters, levels(pop_meta$Loc)), breaks = c(levels(pop_meta$Site), clusters,
                                                                   levels(pop_meta$Loc))) +
    theme(legend.position = "none",
          
          strip.background.x = element_blank(),
          strip.text.x = element_text(size = 12),
          panel.border = element_rect(fill = "NA", linewidth = 0.2),
          panel.background = element_blank())
}

circle_tree <- function()

# Arguments ====
args = commandArgs(TRUE)
VCF_NAME = "ac" # "ac" "hu" ""lm"
STAGE = "1div" # "1div"
PARAM = "nc" # "nc" or "wc"
k_selection = "Taxa" # "Taxa" or "Cluster"
include.het = "no" # "yes" or "no"

# plot_widths = list(as.numeric(strsplit(args[6], ",")[[1]])) # how spaces the plot widths are (given how many panels)
# TODO: plot widths


# Set specific values here ====
if (VCF_NAME == "ac") {
  cluster_k <- 4
  taxa_k <- 2
  colours_clust <- c("#2D4E24", "#8FC73E", "#43BB93", "#347F66")
  levels <- c("WP", "CA", "SB", "SQ")
  site_colours <- c("#FBD1D7", "#EC9C9D", "#AA5459",
                    "#E2DCE0", "#C4AABB", "#8C6C82",
                    "#C0DFF5", "#9BB7DB", "#5179A6",
                    "#8BD3D7", "#74A8B6")
  site_labels = c("WP05", "WP10", "WP20",
                  "CA05", "CA10", "CA20", 
                  "SB05", "SB10", "SB20",
                  "SQ12", "SQ20")
} else if (VCF_NAME == "hu") {
  cluster_k <- 4
  taxa_k <- 3
  colours_clust <- c("#B06327", "#FAA41A", "#F05123", "#821F27")
  levels <- c("WP", "CA", "SB", "SQ")
  site_colours <- c("#FBD1D7", "#EC9C9D", "#AA5459",
                    "#E2DCE0", "#C4AABB", "#8C6C82",
                    "#C0DFF5", "#9BB7DB", "#5179A6",
                    "#8BD3D7", "#74A8B6")
  site_labels = c("WP05", "WP10", "WP20",
                  "CA05", "CA10", "CA20", 
                  "SB05", "SB10", "SB20",
                  "SQ12", "SQ20")
} else if (VCF_NAME == "lm") {
  cluster_k <- 2
  taxa_k <- 2
  colours_clust <- c("#6F1242", "#6652A2")
  site_colours <- c("#EC9C9D", "#AA5459",
                    "#C4AABB", "#8C6C82",
                    "#9BB7DB", "#5179A6",
                    "#8BD3D7", "#74A8B6")
  site_labels <- c("WP10", "WP20",
                   "CA10", "CA20", 
                   "SB10", "SB20",
                   "SQ12", "SQ20")
} else {
  cluster_k <- 5
  taxa_k <- 3
  colours_clust <- rainbow(cluster_k, alpha = 1)
  # need to add arbitrary site_colours and site labels - but don't know how many sites,
  # so will have to come later
}

# clusters names
if (k_selection == "Taxa") {
  k <- taxa_k
} else if (k_selection == "Clusters") {
  k <- cluster_k
} else {
  print("Incorrect entry for k_selection")
}

clusters <- rep(k, k)
for (i in 1:k) {
  clusters[i] <- paste0("Clust", i)}

# Data needed ====
pop_meta <- read.csv(paste0("../results/", VCF_NAME, "_", STAGE, "_", PARAM, "_20_", cluster_k, ".csv"), header = TRUE, stringsAsFactors = FALSE)
# miss_dat <- read.table(paste0("../data/miss-indiv/", VCF_NAME, "_", STAGE, "_", PARAM, "_20_miss-indiv.txt"), header = TRUE)
tree <- read.nexus(paste0("../results/", VCF_NAME,  "_", STAGE, "_", PARAM, "_20_tree.nex"))

# Adjustments to data files =====
# Clusters & colours
# For using less k than k all, we first set all entries to be the first colour
pop_meta$colours <- colours_clust[1]
for (i in 1:k) {
  pop_meta$colours[pop_meta$Clusters == paste0("Clust", i)] <- colours_clust[i]
}

# Levels, factors, values & colours =====
pop_meta$Site <- factor(pop_meta$Site, levels = site_labels)
pop_meta$Loc <- factor(pop_meta$Loc, levels = c("WP", "CA", "SB", "SQ"))
for (i in 1:length(levels(pop_meta$Loc))) {
  pop_meta$Loc_value[pop_meta$Loc == levels(pop_meta$Loc)[i]] = i + k
  pop_meta$Loc_value[is.na(pop_meta$Loc)] = length(levels(pop_meta$Loc)) + k + 1
}
pop_meta$Depth[pop_meta$Depth == "10-12"] <- "10"
pop_meta$Depth <- factor(pop_meta$Depth, levels = c("5", "10", "20"))
pop_meta$Depth_values <- as.numeric(pop_meta$Depth)
loc_colours <- c("#EC9C9D", "#C4AABB", "#9BB7DB", "#8BD3D7")
colours <- c(site_colours, colours_clust[1:k], loc_colours)

# Admixture ====
if (k_selection == "Taxa") {
  if (VCF_NAME == "lm") {
    admix_column_start = which(colnames(pop_meta) == "Q.Clust1")
    admix_column_end = admix_column_start + (k - 1) 
  } else {
  admix_column_start = which(colnames(pop_meta) == "taxaQ.Clust1")
  admix_column_end = admix_column_start + (k - 1)
  }
  admix <- organise_data(pop_meta[,admix_column_start:admix_column_end], pop_meta[,1:2], clusters)
} else if (k_selection == "Clusters") {
  admix_column_start = which(colnames(pop_meta) == "Q.Clust1")
  admix_column_end = admix_column_start + (k - 1)
  admix <- organise_data(pop_meta[,(11 + (k - 1)):(11 + (k - 1))], pop_meta[,1:2], clusters)
}

# Heterozygosity ====
if (include.het == "yes") {
  het <- read.csv(paste0("../../E-intraspecific/results/fstat/", VCF_NAME, "_", STAGE, "_", PARAM, "_indiv.het.csv"), row.names = 1) # need to redo heterozygosity for without clones
  het <- het[order(het$X),]
  het$Loc <- pop_meta$Loc[match(het$X, pop_meta$Individual)]
  het$Clust <- pop_meta$Clusters[match(het$X, pop_meta$Individual)]
} else if (include.het == "no") {
  print("not including heterozygosity")
} else {
  print("do you want to include heterozygosity?")
}

# Basic tree ====
tree_plot <- basic_tree(tree, pop_meta, k_selection)
tree_plot

# Advanced trees ====
if (PARAM == "nc") {
  complot_admixture <- advanced_tree(tree_plot, pop_meta, colours[(length(levels(pop_meta$Site)) + 1):(length(levels(pop_meta$Site)) + k)], clusters, k)
  ggsave(paste0("../results/", VCF_NAME, "_", STAGE, "_nc_20_complot_admixture.pdf"), complot_admixture,
         dpi = 300, width = 3, height = 10, units = "cm")
  complot_loc_depth <- advanced_tree_loc_depth(tree_plot, pop_meta, colours, clusters, k)
  ggsave(paste0("../results/", VCF_NAME, "_", STAGE, "_nc_20_complot_admixture_final.pdf"), complot_loc_depth,
         dpi = 300, width = 3, height = 10, units = "cm")
}

# TODO: Circular plot
# Use basic plot details here...
p <- ggtree(tree_x, layout = "fan", open.angle = 20) + 
  geom_tippoint(aes(colour = tree_tibble_group$group), size = 2) +
  scale_colour_manual(values = colours_clust[1:2])

heatmap_data <- as.data.frame(pop_meta[,1])
for (sites in c("WP", "CA", "SB", "SQ")) {
  heatmap_data[[sites]] <- as.character(pop_meta$Loc)
  heatmap_data[[sites]][heatmap_data[[sites]] != sites] <- 0
}
heatmap_data$WP[heatmap_data$WP == "WP"] <- 3
heatmap_data$CA[heatmap_data$CA == "CA"] <- 4
heatmap_data$SB[heatmap_data$SB == "SB"] <- 5
heatmap_data$SQ[heatmap_data$SQ == "SQ"] <- 6

rownames(heatmap_data) <- heatmap_data[,1]
heatmap_data <- heatmap_data[,-1]
rn <- rownames(heatmap_data)
heatmap_data <- as.data.frame(sapply(heatmap_data, as.character))
rownames(heatmap_data) <- rn

heatmap.colours <- c("transparent", loc_colours)
names(heatmap.colours) <- c(0, 3:6)

gheatmap(p, heatmap_data, color = NULL, offset = 0.01, width = 0.2,
         colnames_angle = 90, colnames_offset_y = -1) +
  scale_fill_manual(values = heatmap.colours,
                    breaks = names(heatmap.colours)) +
  theme(legend.position = "none",
        plot.background = element_blank(),
        panel.background = element_blank())
ggsave("../results/ac_1div_nc_circ_plot.pdf", units = "cm", height = 15, width = 15)

### Depth
heatmap_data <- as.data.frame(pop_meta[,1])
for (i in 1:3) {
  heatmap_data[, 1 + i] <- as.character(pop_meta$Depth_values)
  heatmap_data[, 1 + i][heatmap_data[, 1 + i] != i] <- 0
}
rownames(heatmap_data) <- heatmap_data[,1]
heatmap_data <- heatmap_data[,-1]
rn <- rownames(heatmap_data)
heatmap_data <- as.data.frame(sapply(heatmap_data, as.character))
rownames(heatmap_data) <- rn

heatmap.colours <- c("transparent", "lightblue", "blue", "darkblue")
names(heatmap.colours) <- c(0, 1:3)

gheatmap(p, heatmap_data, color = NULL, offset = 0.01, width = 0.2,
         colnames_angle = 90, colnames_offset_y = -1) +
  scale_fill_manual(values = heatmap.colours,
                    breaks = names(heatmap.colours)) +
  theme(legend.position = "none",
        plot.background = element_blank(),
        panel.background = element_blank())
ggsave("../results/ac_1div_nc_circ_plot_depth.pdf", units = "cm", height = 15, width = 15)

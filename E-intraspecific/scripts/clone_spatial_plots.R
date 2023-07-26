# Investigation of spatial genetic structure and clone distributions
# Author: Katharine Prata
# Last edit: 15/4/22

# Packages
# R v4.2.0
library(tidyverse) # v1.3.1
library(adegenet) # v2.1.7
library(vcfR) # v1.12.0
library(spaa) # v0.2.2
library(vegan) # v2.6-2
library(ape) # v5.6-3
library(reshape2) # v1.4.4

# Functions
make_clone_distribution <- function(clonegroups, clusters, site) {
  if (clusters == "AA1") {
    color = colours[1]
    site_colours <- site_colours[c(2, 3, 5, 6, 9, 11)]
  } else if (clusters == "AA2") {
    color = colours[2]
  } else if (clusters == "AH1") {
    color = colours[3]
    site_colours <- site_colours[-11]
  } else if (clusters == "AH2") {
    color = colours[4]
    site_colours <- site_colours[c(1,2,7,8)]
  } else if (clusters == "AH3") {
    color = colours[5]
    site_colours <- site_colours[c(1,2,4,7)]
  } else if (clusters == "AL1") {
    color = colours[6]
    site_colours <- site_colours[c(2,3,5,6,8,9,11)]
  } else if (clusters == "AL2") {
    color = colours[7]
    site_colours <- site_colours[c(2,3,6,8,9,10,11)]
    }
  
  if (site != "all sites") {
  clonegroups.sub <- all.clonegroups[all.clonegroups$Clusters == clusters | all.clonegroups$Site == site,] 
  clone_distrib <- as.data.frame(table(table(clonegroups.sub$Groups)))
  ggplot(clone_distrib, aes(x = Var1, y = Freq)) + geom_bar(stat = "identity", fill = color) +
    theme_classic() + ylab("Frequency of clone group") + xlab("Number of colonies within a group") +
    scale_x_discrete(labels = paste0(clone_distrib$Var1, "\n n = ", clone_distrib$Freq)) + ggtitle(paste(clusters, site))
  } else {
    clonegroups.sub <- all.clonegroups[all.clonegroups$Clusters == clusters,]
    clone_distrib <- as.data.frame(table(table(clonegroups.sub$Groups)))
    print("Using all sites")
    ggplot(clonegroups.sub, aes(x = Freq, fill = Site)) +
      geom_bar(stat = "count") + theme_classic() + scale_fill_manual(values = site_colours) +
      scale_x_discrete(labels = paste0(clone_distrib$Var1, "\n n = ", clone_distrib$Freq)) +
      ggtitle(paste(clusters))
    ggsave(paste0("3a - Clone distances/", clusters,"_clonegroups.pdf"))
  }
}
make_clone_freq_table <- function(clonegroups) {
  freq_table <- data.frame(Group = NA, Freq_group = NA, Site = NA, Taxa = NA)
  for (taxa in c("AA1", "AA2", "AH1", "AH2", "AH3", "AL1", "AL2")) {
    subset.clonegroups <- all.clonegroups[all.clonegroups$Clusters == taxa,]
    subset.clonegroups$Site <- droplevels(subset.clonegroups$Site)
    for (site in levels(subset.clonegroups$Site)) {
      print(site)
      dat <- as.data.frame(table(table(subset.clonegroups$Groups[subset.clonegroups$Site == site])))
      if (length(dat) > 1) {
        dat$Site <- site
        dat$Taxa <- taxa
        colnames(dat) <- c("Group", "Freq_group", "Site", "Taxa")
        freq_table <- rbind(freq_table, dat)
      }}
    }
    freq_table <- na.omit(freq_table)
    freq_table$Site <- as.factor(freq_table$Site)
    freq_table$Group <- as.factor(freq_table$Group)
  return(freq_table)}
make_clone_freq_distribution <- function(freq_table, clusters) {
  if (clusters == "AA1") {
    color = colours[1]
    site_colours <- site_colours[c(2, 3, 5, 6, 9, 11)]
  } else if (clusters == "AA2") {
    color = colours[2]
  } else if (clusters == "AH1") {
    color = colours[3]
    site_colours <- site_colours[-11]
  } else if (clusters == "AH2") {
    color = colours[4]
    site_colours <- site_colours[c(1,2,7,8)]
  } else if (clusters == "AH3") {
    color = colours[5]
    site_colours <- site_colours[c(1,2,4,7)]
  } else if (clusters == "AL1") {
    color = colours[6]
    site_colours <- site_colours[c(2,3,5,6,8,9,11)]
  } else if (clusters == "AL2") {
    color = colours[7]
    site_colours <- site_colours[c(2,3,6,8,9,10,11)]
  }
  
  subset.freq_table <- freq_table[freq_table$Taxa == clusters, ]

  ggplot(subset.freq_table, aes(x = Group, y = Freq_group, fill = Site)) +
    geom_bar(position = "stack", stat = "identity") + theme_classic() + scale_fill_manual(values = site_colours) +
    xlab("Number of colonies in a group") +
    ylab("Frequency of groups")
  ggsave(paste0("../results/clone_dist/", clusters,"_freqclonegroups.pdf"), height = 10, width = 7, units = "cm")
}
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
sep_plot <- function(metadata, plot) {
  plot.data <- metadata[metadata$Site == plot,]
  plot.data$clone_status = "clone"
  plot.data[plot.data$Individual %in% ac.mlgs$V1, length(plot.data)] <- "ac.mlg"
  plot.data[plot.data$Individual %in% lm.mlgs$V1, length(plot.data)] <- "lm.mlg"
  plot.data[plot.data$Individual %in% hu.mlgs$V1, length(plot.data)] <- "hu.mlg"
  return(plot.data)
}
match_clones <- function(plot, clonegroups) {
  plot$clone_group <- clonegroups[match(plot$Individual, clonegroups$Sample), 3]
  clone_list <- unique(plot$clone_group[plot$clone_status == "clone"])
  for (i in clone_list) {
    plot$clone_status[plot$clone_group == i] <- paste0("clone", i)
  }
  print(clone_list)
  return(plot)
}

# Import  ####
metadata <- read.csv('../data/all_annotations_X_DEPTH_parallel_XYZ_adjusted.txt')
ac.mlgs <- read.csv("../data/ac_1d_wc_20_colours.txt", header = FALSE)
hu.mlgs <- read.csv("../data/hu_1d_wc_20_colours.txt", header = FALSE)
lm.mlgs <- read.csv("../data/lm_1d_wc_20_colours.txt", header = FALSE)
ac.clonegroups <- read.csv("../data/clone_groups_ac.csv")
hu.clonegroups <- read.csv("../data/clone_groups_hu.csv")
lm.clonegroups <- read.csv("../data/clone_groups_lm.csv")
pop.wc.ac <- read.csv("../data/ac_1d_wc_20_4.csv")
pop.wc.hu <- read.csv("../data/hu_1d_wc_20_4.csv")
pop.wc.lm <- read.csv("../data/lm_1d_wc_20_2.csv")

# Colours
              #AA1        #AA2        #AH1       #AH2      #AH3      #AL1      #AL2
colours <- c("#274e13", "#8fce00", "#b06100", "#ffa500", "#ff4e00", "#6f0641", "#7223f0")
site_colours <- c("#FBD1D7", "#EC9C9D", "#AA5459",
                  "#E2DCE0", "#C4AABB", "#8C6C82",
                  "#C0DFF5", "#9BB7DB", "#5179A6",
                  "#8BD3D7", "#74A8B6")

# Organising clone groups ####
# Fixing NAs
ac.clonegroups[ac.clonegroups$Pop == "NA_NANA",]
ac.clonegroups$Pop[ac.clonegroups$Sample == "KP0406_NA_NANA"] <- "AC_WP20"
ac.clonegroups$Pop[ac.clonegroups$Sample == "KP0673_NA_NANA"] <- "AC_SB20"
ac.clonegroups[ac.clonegroups$Pop == "UN_NANA",]
ac.clonegroups$Pop[ac.clonegroups$Sample == "KP0844_UN_NANA"] <- "AC_WP10"
lm.clonegroups[lm.clonegroups$Pop == "NA_NANA",]
lm.clonegroups$Pop[lm.clonegroups$Pop == "NA_NANA"] <- "LM_SB20"


# Remove technical replicates
# ac
ac.clonegroups[ac.clonegroups$Sample == "KX0394_AC_WP20",]
ac.clonegroups[ac.clonegroups$Sample == "KX0396_AC_WP20",]
ac.clonegroups[ac.clonegroups$Sample == "KX0604_AC_SB20",]
ac.clonegroups[ac.clonegroups$Sample == "KX0743_AC_WP10",]

ac.clonegroups <- ac.clonegroups[-c(508:511),]

# lm
lm.clonegroups[lm.clonegroups$Sample == "KX0330_LM_WP20",]
lm.clonegroups[lm.clonegroups$Sample == "KX0657_LM_SB20",]
lm.clonegroups <- lm.clonegroups[-c(103:104),]

hu.clonegroups[hu.clonegroups$Sample == "KX0075_HU_SB10",]
hu.clonegroups <- hu.clonegroups[-c(142),]

# Combining datasets and counting clone groups ####
# When unequal sample sizes have to have smaller list first and use "match" not %in%
ac.clonegroups$Clusters <- pop.wc.ac$Clusters[match(ac.clonegroups$Sample, pop.wc.ac$Individual)]
ac.clonegroups$Clusters[ac.clonegroups$Clusters == "Clust1"] = "AA1"
ac.clonegroups$Clusters[ac.clonegroups$Clusters == "Clust2" |
                          ac.clonegroups$Clusters == "Clust3"  | ac.clonegroups$Clusters == "Clust4"] = "AA2"
ac.clonegroups <- ac.clonegroups %>% unite(all_groups, Groups, Clusters, remove = FALSE)

hu.clonegroups$Clusters <- pop.wc.hu$Clusters[match(hu.clonegroups$Sample, pop.wc.hu$Individual)]
hu.clonegroups$Clusters[hu.clonegroups$Clusters == "Clust1" | hu.clonegroups$Clusters == "Clust4"] = "AH1"
hu.clonegroups$Clusters[hu.clonegroups$Clusters == "Clust2"] = "AH2"
hu.clonegroups$Clusters[hu.clonegroups$Clusters == "Clust3"] = "AH3"
hu.clonegroups <- hu.clonegroups %>% unite(all_groups, Groups, Clusters, remove = FALSE)

lm.clonegroups$Clusters <- pop.wc.lm$Clusters[match(lm.clonegroups$Sample, pop.wc.lm$Individual)]
lm.clonegroups$Clusters[lm.clonegroups$Clusters == "Clust1"] = "AL1"
lm.clonegroups$Clusters[lm.clonegroups$Clusters == "Clust2"] = "AL2"
lm.clonegroups <- lm.clonegroups %>% unite(all_groups, Groups, Clusters, remove = FALSE)

all.clonegroups <- rbind(ac.clonegroups, hu.clonegroups, lm.clonegroups)
rm(lm.clonegroups, ac.clonegroups, hu.clonegroups)

str(all.clonegroups)
all.clonegroups$all_groups <- as.factor(all.clonegroups$all_groups)
all.clonegroups$Clusters <- as.factor(all.clonegroups$Clusters)
all.clonegroups$Pop <- as.factor(all.clonegroups$Pop)
all.clonegroups <- all.clonegroups %>% separate(Pop, into = c("Species", "Site")) %>%
  separate(Site, into = c("Loc", "Depth"), sep = 2, remove = FALSE)
x <- as.data.frame(table(all.clonegroups$all_groups))
all.clonegroups$Freq <- x$Freq[match(all.clonegroups$all_groups, x$Var1)]
all.clonegroups$Freq <- as.factor(all.clonegroups$Freq)

# Clonal richness ####
table(all.clonegroups$Clusters)
# AA1
length(unique(all.clonegroups$Groups[all.clonegroups$Clusters == "AA1"]))
table(all.clonegroups$Groups[all.clonegroups$Clusters == "AA1"])
table(table(all.clonegroups$Groups[all.clonegroups$Clusters == "AA1"]))
aa1.R <- (length(unique(all.clonegroups$Groups[all.clonegroups$Clusters == "AA1"])) -1) / (46 - 1)
print(paste0("Clonal richness for AA1 is ", aa1.R))
#[1] 0.7333333

# AA2
length(unique(all.clonegroups$Groups[all.clonegroups$Clusters == "AA2"]))
table(all.clonegroups$Groups[all.clonegroups$Clusters == "AA2"])
table(table(all.clonegroups$Groups[all.clonegroups$Clusters == "AA2"]))
aa2.R <- (length(unique(all.clonegroups$Groups[all.clonegroups$Clusters == "AA2"])) -1) / (462 - 1)
print(paste0("Clonal richness for AA2 is ", aa2.R))
#[1] 0.6464208

# AH1
length(unique(all.clonegroups$Groups[all.clonegroups$Clusters == "AH1"]))
table(all.clonegroups$Groups[all.clonegroups$Clusters == "AH1"])
table(table(all.clonegroups$Groups[all.clonegroups$Clusters == "AH1"]))
ah1.R <- (length(unique(all.clonegroups$Groups[all.clonegroups$Clusters == "AH1"])) -1) / (88 - 1)
print(paste0("Clonal richness for AH1 is ", ah1.R))
#[1] 0.816092

# AH2
length(unique(all.clonegroups$Groups[all.clonegroups$Clusters == "AH2"]))
table(all.clonegroups$Groups[all.clonegroups$Clusters == "AH2"])
table(table(all.clonegroups$Groups[all.clonegroups$Clusters == "AH2"]))
ah2.R <- (length(unique(all.clonegroups$Groups[all.clonegroups$Clusters == "AH2"])) -1) / (31 - 1)
print(paste0("Clonal richness for AH2 is ", ah2.R))
#[1] 0.9333333

# AH3
length(unique(all.clonegroups$Groups[all.clonegroups$Clusters == "AH3"]))
table(all.clonegroups$Groups[all.clonegroups$Clusters == "AH3"])
table(table(all.clonegroups$Groups[all.clonegroups$Clusters == "AH3"]))
ah3.R <- (length(unique(all.clonegroups$Groups[all.clonegroups$Clusters == "AH3"])) -1) / (22 - 1)
print(paste0("Clonal richness for AH3 is ", ah3.R))
#[1] 0.9047619

# AL1
length(unique(all.clonegroups$Groups[all.clonegroups$Clusters == "AL1"]))
table(all.clonegroups$Groups[all.clonegroups$Clusters == "AL1"])
table(table(all.clonegroups$Groups[all.clonegroups$Clusters == "AL1"]))
al1.R <- (length(unique(all.clonegroups$Groups[all.clonegroups$Clusters == "AL1"])) -1) / (34 - 1)
print(paste0("Clonal richness for AL1 is ", al1.R))
#[1] 0.9090909

# AL2
length(unique(all.clonegroups$Groups[all.clonegroups$Clusters == "AL2"]))
table(all.clonegroups$Groups[all.clonegroups$Clusters == "AL2"])
table(table(all.clonegroups$Groups[all.clonegroups$Clusters == "AL2"]))
al2.R <- (length(unique(all.clonegroups$Groups[all.clonegroups$Clusters == "AL2"])) -1) / (69 - 1)
print(paste0("Clonal richness for AL2 is ", al2.R))
#[1] 0.8970588

# Clone distrib plot ####
all.clonegroups$Site <- factor(all.clonegroups$Site, levels = c("WP05", "WP10", "WP20",
                                                                "CA05", "CA10", "CA20",
                                                                "SB05", "SB10", "SB20",
                                                                "SQ12", "SQ20"))

#make_clone_distribution(all.clonegroups, "AA1", site = "all sites")
#make_clone_distribution(all.clonegroups, "AA2", site = "all sites")
#make_clone_distribution(all.clonegroups, "AH1", site = "all sites")
#make_clone_distribution(all.clonegroups, "AH2", site = "all sites")
#make_clone_distribution(all.clonegroups, "AH3", site = "all sites")
#make_clone_distribution(all.clonegroups, "AL1", site = "all sites")
#make_clone_distribution(all.clonegroups, "AL2", site = "all sites")
#make_clone_distribution(all.clonegroups, "AA1", site = "WP20")
#make_clone_distribution(all.clonegroups, "AA1", site = "SB20")

ggplot(all.clonegroups, aes(x = Freq, fill = Site)) +
  geom_bar(stat = "count") + theme_classic() + scale_fill_manual(values = site_colours) +
  facet_wrap(~Clusters)
ggsave(paste0("../results/clone_dist/all_clonegroups.pdf"))


freq_table <- make_clone_freq_table(all.clonegroups)
freq_table$Site <- factor(freq_table$Site, levels = c("WP05", "WP10", "WP20",
                                                      "CA05", "CA10", "CA20",
                                                      "SB05", "SB10", "SB20",
                                                      "SQ12", "SQ20"))



make_clone_freq_distribution(freq_table, "AA2")

freq_table_noAA2 <- freq_table[freq_table$Taxa != "AA2",]

ggplot(freq_table_noAA2, aes(x = Group, y = Freq_group, fill = Site)) +
  geom_bar(position = "stack", stat = "identity") + theme_classic() + scale_fill_manual(values = site_colours) +
  xlab("Number of colonies in a group") +
  ylab("Frequency of groups") + facet_wrap(~Taxa)
ggsave(paste0("../results/clone_dist/_freqclonegroups_noAA2.pdf"), height = 10, width = 10, units = "cm")

# Within singles ####
freq_table2 <- freq_table[freq_table$Group != "1",]

ggplot(freq_table2[freq_table2$Taxa == "AA2",], aes(x = Group, y = Freq_group, fill = Site)) +
  geom_bar(position = "stack", stat = "identity") + theme_classic() + scale_fill_manual(values = site_colours[c(2,3,5,6,8:11)]) +
  facet_wrap(~Taxa) +
  theme(legend.position = "none",
        axis.title = element_blank())
ggsave(paste0("../results/clone_dist/freqclonegroups_AA2.pdf"), height = 10, width = 3, units = "cm")

freq_table2_noAA2 <- freq_table2[freq_table2$Taxa != "AA2",]

row1 <- c(7, 0, "SB05", "AA1")
freq_table2_noAA2$Group
freq_table2_noAA2 <- rbind(freq_table2_noAA2, row1)
freq_table2_noAA2$Freq_group <- as.numeric(freq_table2_noAA2$Freq_group)
ggplot(freq_table2_noAA2, aes(x = Group, y = Freq_group, fill = Site)) +
  geom_bar(position = "stack", stat = "identity") + theme_classic() + scale_fill_manual(values = site_colours[2:11]) +
  xlab("Number of colonies in a group") +
  ylab("Frequency of groups") + facet_wrap(~Taxa)
ggsave(paste0("../results/clone_dist/freqclonegroups_noAA2.pdf"), height = 10, width = 10, units = "cm")

# Organising spatial metadata ####
metadata <- metadata %>% separate(Individual, into=c("Sample", "Species", "Site")) %>% 
  unite(Individual, Sample, Species, Site, remove = FALSE) %>% 
  dplyr::select(Individual, Sample, Species, Site, Loc, Depth, x, y, z)
metadata$clone = "clone"
metadata[metadata$Individual %in% ac.mlgs$V1, length(metadata)] <- "ac.mlg"
metadata[metadata$Individual %in% lm.mlgs$V1, length(metadata)] <- "lm.mlg"
metadata[metadata$Individual %in% hu.mlgs$V1, length(metadata)] <- "hu.mlg"

all.clonegroups$clones <- all.clonegroups$all_groups
levels(all.clonegroups$clones) <- c(levels(all.clonegroups$clones), "non-clone")
all.clonegroups$clones[all.clonegroups$Freq == "1"] = "non-clone"
all.clonegroups$clones <- droplevels(all.clonegroups$clones)
metadata$clone_group <- all.clonegroups$clones[match(metadata$Individual, all.clonegroups$Sample)]
# remove instances no genotypes colonies
metadata <- na.omit(metadata) # 636
metadata$Clusters <- all.clonegroups$Clusters[match(metadata$Individual, all.clonegroups$Sample)]

# Writing metadata ####
write.table(metadata, file = "../results/clone_dist/all-aga_spatgeo.txt", quote = FALSE, row.names = FALSE)
write.table(all.clonegroups, file = "../results/clone_dist/all-aga_clones.txt", quote = FALSE, row.names = FALSE)

# Euclidean spatial distances between clones ####
# aa1
list_aa1 <- list_df(metadata, "AA1")
ggplot(list_aa1, aes(x = group_pair, y = value, color = group_stat)) + geom_jitter() + theme_classic() + coord_flip()
list.small <- list_aa1[list_aa1$col_group == list_aa1$row_group & list_aa1$group_pair != "non-clone-non-clone", ]
ggplot(list.small, aes(x = group_pair, y = value, fill = group_pair)) + geom_jitter(size = 2, shape = 21) + theme_classic() + coord_flip()
ggsave("../results/clone_dist/AA1_clone_geo_dist.pdf")

ggplot(list.small, aes(value)) + 
  geom_histogram(breaks = seq(0, 0.6, 0.1), fill = colours[1], colour = "black") +
  theme_classic() + xlab("distance between clones (m)") + ylab("Frequnecy") +
  xlim(c(0,3)) +
  scale_y_continuous(limits = c(0, 10), expand = c(0, 0), breaks = seq(0, 10, 2))
ggsave("../results/clone_dist/AA1_clone_geo_distancedistrib.pdf", height = 5, width = 5, units = "cm")

# list small aa2
list_aa2 <- list_df(metadata, "AA2")
list.small <- list_aa2[list_aa2$col_group == list_aa2$row_group & list_aa2$group_pair != "non-clone-non-clone", ]
list.small <- list.small[list.small$value < 100,]
ggplot(list.small, aes(x = group_pair, y = value, fill = group_pair)) + geom_jitter(size = 2, shape = 21) + theme_classic() + coord_flip()
ggsave("../results/clone_dist/AA2_clone_geo_dist.pdf", height = 30, width = 50, units = "cm")

ggplot(list.small, aes(value)) + 
  geom_histogram(breaks = seq(0, 3, 0.1), fill = colours[2], colour = "black") +
  theme_classic() + xlab("distance between clones (m)") + ylab("Frequnecy") +
  scale_y_continuous(limits = c(0, 80), expand = c(0, 0), breaks = seq(0, 70, 10))
ggsave("../results/clone_dist/AA2_clone_geo_distancedistrib.pdf", height = 10, width = 5, units = "cm")

# al1
list_al1 <- list_df(metadata, "AL1")
list.small <- list_al1[list_al1$col_group == list_al1$row_group & list_al1$group_pair != "non-clone-non-clone", ]
ggplot(list.small, aes(x = group_pair, y = value, fill = group_pair)) + geom_jitter(size = 2, shape = 21) + theme_classic() + coord_flip()

ggplot(list.small, aes(value)) + 
  geom_histogram(breaks = seq(0, 3, 0.1), fill = colours[6], colour = "black") +
  theme_classic() + xlab("distance between clones (m)") + ylab("Frequnecy") +
  xlim(c(0,3)) +
  scale_y_continuous(limits = c(0, 10), expand = c(0, 0), breaks = seq(0, 10, 2))
ggsave("../results/clone_dist/AL1_clone_geo_distancedistrib.pdf", height = 5, width = 5, units = "cm")

# al2
list_al2 <- list_df(metadata, "AL2")
list.small <- list_al2[list_al2$col_group == list_al2$row_group & list_al2$group_pair != "non-clone-non-clone", ]

# mislabelled samples
list.small <- list.small[list.small$value < 100, ]

ggplot(list.small, aes(x = group_pair, y = value, fill = group_pair)) + geom_jitter(size = 2, shape = 21) +
  theme_classic() + coord_flip()
ggsave("../results/clone_dist/AL2_clone_geo_dist.pdf")

ggplot(list.small, aes(value)) + 
  geom_histogram(breaks = seq(0, 1, 0.1), fill = colours[7], colour = "black") +
  theme_classic() + xlab("distance between clones (m)") + ylab("Frequnecy")  +
  xlim(c(0,3)) +
  scale_y_continuous(limits = c(0, 10), expand = c(0, 0), breaks = seq(0, 10, 2))
ggsave("../results/clone_dist/AL2_clone_geo_distancedistrib.pdf", height = 5, width = 5, units = "cm")

# ah1
list_ah1 <- list_df(metadata, "AH1")
list.small <- list_ah1[list_ah1$col_group == list_ah1$row_group & list_ah1$group_pair != "non-clone-non-clone", ]
ggplot(list.small, aes(x = group_pair, y = value, fill = group_pair)) + geom_jitter(size = 2, shape = 21) + theme_classic() + coord_flip()
ggsave("../results/clone_dist/AH1_clone_geo_dist.pdf")

# there are clones between WP 10 - 20m
ggplot(list.small, aes(value)) + 
  geom_histogram(breaks = seq(0, 30, 1), fill = colours[3], colour = "black") +
  theme_classic() + xlab("distance between clones (m)") + ylab("Frequnecy") +
  xlim(c(0,30)) +
  scale_y_continuous(limits = c(0, 10), expand = c(0, 0), breaks = seq(0, 10, 2))
ggsave("../results/clone_dist/AH1_clone_geo_distancedistrib.pdf", height = 5, width = 10, units = "cm")

# ah2
list_ah2 <- list_df(metadata, "AH2")
list.small <- list_ah2[list_ah2$col_group == list_ah2$row_group & list_ah2$group_pair != "non-clone-non-clone", ]
ggplot(list.small, aes(x = group_pair, y = value, fill = group_pair)) + geom_jitter(size = 2, shape = 21) + theme_classic() + coord_flip()
ggsave("../results/clone_dist/AH2_clone_geo_dist.pdf")

ggplot(list.small, aes(value)) + 
  geom_histogram(breaks = seq(0, 1, 0.1), fill = colours[4], colour = "black") +
  theme_classic() + xlab("distance between clones (m)") + ylab("Frequnecy") +
  xlim(c(0,3)) +
  scale_y_continuous(limits = c(0, 10), expand = c(0, 0), breaks = seq(0, 10, 2))
ggsave("../results/clone_dist/AH2_clone_geo_distancedistrib.pdf", height = 5, width = 5, units = "cm")

# ah3
list_ah3 <- list_df(metadata, "AH3")
list.small <- list_ah3[list_ah3$col_group == list_ah3$row_group & list_ah3$group_pair != "non-clone-non-clone", ]
ggplot(list.small, aes(x = group_pair, y = value, fill = group_pair)) + geom_jitter(size = 2, shape = 21) + theme_classic() + coord_flip()
ggsave("../results/clone_dist/AH3_clone_geo_dist.pdf")

# there are clones between WP05 - WP10
ggplot(list.small, aes(value)) + 
  geom_histogram(breaks = seq(0, 30, 1), fill = colours[5], colour = "black") +
  theme_classic() + xlab("distance between clones (m)") + ylab("Frequnecy")  +
  xlim(c(0,30)) +
  scale_y_continuous(limits = c(0, 10), expand = c(0, 0), breaks = seq(0, 10, 2))
ggsave("../results/clone_dist/AH3_clone_geo_distancedistrib.pdf", height = 5, width = 10, units = "cm")

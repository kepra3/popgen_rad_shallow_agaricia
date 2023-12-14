# Title: Plotting kinship distances
# Author: Katharine Prata

# Packages  ####
library(dplyr)
library(tidyr)
library(ggplot2)

# Functions  ####
kinship_distance <- function(sample1, sample2) {
  # First sample
  one <- metadata[metadata$Individual == sample1,]
  # Second sample
  two <- metadata[metadata$Individual == sample2,]
  pair <- rbind(one, two)
  # distance
  pair.dist <- dist(pair[,7:9])
  return(pair.dist)
}
kinship_df <- function(kinship_txt, relation, taxa) {
  df.length <- length(kinship_txt[,1])
  dists <- data.frame(1:df.length)
  for (i in 1:df.length) {
    dists[i] <- as.numeric(kinship_distance(kinship_txt[i, 1],	kinship_txt[i, 2]))
    }
  kinship.df <- data.frame(Sample1 = kinship_txt[1:df.length, 1],
                           Sample2 = kinship_txt[1:df.length,2],
                           Status = rep(relation, df.length), 
                           Taxa = c(rep(taxa, df.length)),
                           distances = t(dists)[1:df.length])
  return(kinship.df)
                                     
}
kinship_plot <- function(kinship.df, max.dist, max.count, colours, width) {
  ggplot(kinship.df, aes(distances, fill = Taxa)) + 
  geom_histogram(breaks = seq(0, max.dist, 1), colour = "black") +
  theme_classic() + xlab(paste("distance between", unique(kinship.df$Status), "(m)")) + ylab("Frequency")  +
  scale_fill_manual(values = colours) +
  xlim(c(0, max.dist)) +
  scale_y_continuous(limits = c(0, max.count), expand = c(0, 0), breaks = seq(0, max.count, 2)) +
  theme(legend.position = "none")
  status <- paste(unique(kinship.df$Status))
  taxa <- paste0(unique(kinship.df$Taxa), collapse = "_")
  
ggsave(paste0("3c - Kinship/plots/", taxa, "_", status, "_distances.pdf"), height = 5, width = width, units = "cm")
}

# Set working directory
setwd("~/Dropbox/agaricia_project_2019/shalo_ag/gen_project/3 - Spatial/")

# Import data  ####
metadata <- read.csv('all_annotations_X_HORIZ_parallel_XYZadjusted.txt')

# aa1
parents.aa1 <-  read.table("3c - Kinship/Colony2_Mac_01_02_2022/aa1/Parents_aa1_1div_nc_50_0.2.txt", sep = ",",
                       header = TRUE)
FullSib.aa1 <-  read.table("3c - Kinship/Colony2_Mac_01_02_2022/aa1/FullSib_aa1_1div_nc_50_0.2.txt", sep = ",",
                       header = TRUE)
HalfSib.aa1 <-  read.table("3c - Kinship/Colony2_Mac_01_02_2022/aa1/HalfSib_aa1_1div_nc_50_0.2.txt", sep = ",",
                       header = TRUE)

# aa2 ####
# ca - only half and fullsibs
FullSib.aa2.ca <-  read.table("3c - Kinship/Colony2_Mac_01_02_2022/aa2-ca/FullSib_aa2-ca_3c_nc_50_0.2.txt", sep = ",",
                           header = TRUE)
HalfSib.aa2.ca <-  read.table("3c - Kinship/Colony2_Mac_01_02_2022/aa2-ca/HalfSib_aa2-ca_3c_nc_50_0.2.txt", sep = ",",
                           header = TRUE)
# sb - only half and fullsibs
FullSib.aa2.sb <-  read.table("3c - Kinship/Colony2_Mac_01_02_2022/aa2-sb/FullSib_aa2-sb_3c_nc_50_0.2.txt", sep = ",",
                              header = TRUE)
HalfSib.aa2.sb <-  read.table("3c - Kinship/Colony2_Mac_01_02_2022/aa2-sb/HalfSib_aa2-sb_3c_nc_50_0.2.txt", sep = ",",
                              header = TRUE)
# sq - only half and fullsibs
FullSib.aa2.sq <-  read.table("3c - Kinship/Colony2_Mac_01_02_2022/aa2-sq/FullSib_aa2-sq_3c_nc_50_0.2.txt", sep = ",",
                              header = TRUE)
HalfSib.aa2.sq <-  read.table("3c - Kinship/Colony2_Mac_01_02_2022/aa2-sq/HalfSib_aa2-sq_3c_nc_50_0.2.txt", sep = ",",
                              header = TRUE)
# wp - all
parents.aa2.wp <-  read.table("3c - Kinship/Colony2_Mac_01_02_2022/aa2-wp/Parents_aa2-wp_3c_nc_50_0.2.txt", sep = ",",
                           header = TRUE)
FullSib.aa2.wp <-  read.table("3c - Kinship/Colony2_Mac_01_02_2022/aa2-wp/FullSib_aa2-wp_3c_nc_50_0.2.txt", sep = ",",
                              header = TRUE)
HalfSib.aa2.wp <-  read.table("3c - Kinship/Colony2_Mac_01_02_2022/aa2-wp/HalfSib_aa2-wp_3c_nc_50_0.2.txt", sep = ",",
                              header = TRUE)

# ah ####
# ah1 - half/full
FullSib.ah1 <-  read.table("3c - Kinship/Colony2_Mac_01_02_2022/ah1/FullSib_ah1_1div_nc_50_0.1.txt", sep = ",",
                              header = TRUE)
HalfSib.ah1 <-  read.table("3c - Kinship/Colony2_Mac_01_02_2022/ah1/HalfSib_ah1_1div_nc_50_0.1.txt", sep = ",",
                              header = TRUE)

# ah2 - half/full
FullSib.ah2 <-  read.table("3c - Kinship/Colony2_Mac_01_02_2022/ah2/FullSib_ah2_1div_nc_50_0.1.txt", sep = ",",
                              header = TRUE)
HalfSib.ah2 <-  read.table("3c - Kinship/Colony2_Mac_01_02_2022/ah2/HalfSib_ah2_1div_nc_50_0.1.txt", sep = ",",
                              header = TRUE)

# ah3 - half
HalfSib.ah3 <-  read.table("3c - Kinship/Colony2_Mac_01_02_2022/ah3/HalfSib_ah3_1div_nc_50_0.1.txt", sep = ",",
                              header = TRUE)

# Colours
#AA1        #AA2        #AH1       #AH2      #AH3      #AL1      #AL2
colours <- c("#274e13", "#8fce00", "#b06100", "#ffa500", "#ff4e00", "#6f0641", "#7223f0")

# Sort data  ####
metadata <- metadata %>% separate(Individual, into=c("Sample", "Species", "Site")) %>% 
  unite(Individual, Sample, Species, Site, remove = FALSE) %>% 
  dplyr::select(Individual, Sample, Species, Site, Loc, Depth, x, y, z)


# Calculate distances ####
## Parent Off plots ####
aa1.parent.df <- kinship_df(parents.aa1, "Parent-Offspring", "AA1")
aa2.parent.df <- kinship_df(parents.aa2.wp, "Parent-Offspring", "AA2")

parents <- rbind(aa1.parent.df, aa2.parent.df)

kinship_plot(parents, max.dist = 5, max.count = 10, colours = colours[1:2], width = 5)


## Full sib plots ####
# aa1.fullsib.df <- kinship_df(FullSib.aa1, "Full-siblings", "AA1")
# KP0769_AC_WP10 KP1132_AC_CA20
# KP0769_AC_WP10 KP1132_AC_CA20
# Fullsiblings for AA1 messed up...
aa2.ca.fullsib.df <- kinship_df(FullSib.aa2.ca, "Full-siblings", "AA2")
aa2.sb.fullsib.df <- kinship_df(FullSib.aa2.sb, "Full-siblings", "AA2")
aa2.sq.fullsib.df <- kinship_df(FullSib.aa2.sq, "Full-siblings", "AA2")
aa2.wp.fullsib.df <- kinship_df(FullSib.aa2.wp, "Full-siblings", "AA2")
ah1.fullsib.df <- kinship_df(FullSib.ah1, "Full-siblings", "AH1")
# two missing from AH1
ah2.fullsib.df <- kinship_df(FullSib.ah2, "Full-siblings", "AH2")

fullsibs <- rbind(aa2.ca.fullsib.df, aa2.sb.fullsib.df, aa2.sq.fullsib.df, aa2.wp.fullsib.df, ah1.fullsib.df, ah2.fullsib.df)
fullsibs <- na.omit(fullsibs)
kinship_plot(fullsibs, max.dist = 30, max.count = 10, colours = colours[c(2,4)], width = 10)

## Half-sib plots ####
aa1.halfsib.df <- kinship_df(HalfSib.aa1, "Half-siblings", "AA1")
aa2.ca.halfsib.df <- kinship_df(HalfSib.aa2.ca, "Half-siblings", "AA2")
aa2.sb.halfsib.df <- kinship_df(HalfSib.aa2.sb, "Half-siblings", "AA2")
aa2.sq.halfsib.df <- kinship_df(HalfSib.aa2.sq, "Half-siblings", "AA2")
aa2.wp.halfsib.df <- kinship_df(HalfSib.aa2.wp, "Half-siblings", "AA2")
ah1.halfsib.df <- kinship_df(HalfSib.ah1, "Half-siblings", "AH1")
ah2.halfsib.df <- kinship_df(HalfSib.ah2, "Half-siblings", "AH2")
ah3.halfsib.df <- kinship_df(HalfSib.ah3, "Half-siblings", "AH3")



halfsibs <- rbind(aa1.halfsib.df, aa2.ca.halfsib.df, aa2.sb.halfsib.df, aa2.sq.halfsib.df, aa2.wp.halfsib.df, ah1.halfsib.df,
                  ah2.halfsib.df, ah3.halfsib.df)
halfsibs <- na.omit(halfsibs)
kinship_plot(halfsibs, max.dist = 30, max.count = 10, colours = colours[1:4], width = 10)


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
  
ggsave(paste0("../results/kinship/", taxa, "_", status, "_distances.pdf"), height = 5, width = width, units = "cm")
}

# Import data  ####
metadata <- read.csv('../data/all_annotations_X_DEPTH_parallel_XYZ_adjusted.txt')

# aa1
kin.aa1.wp <-  read.table("../results/kinship/aa1-wp/aa1-wp_1div_nc_20_0.2.kin.txt", sep = ",",
                       header = TRUE)
kin.aa1.wp$Relationship <- as.factor(kin.aa1.wp$Relationship)
summary(kin.aa1.wp$Relationship)
# aa2 ####
kin.aa2.wp <- read.table("../results/kinship/aa2-wp/aa2-wp_1div_nc_20_0.2.kin.txt", sep = ",",
                        header = TRUE)
kin.aa2.wp$Relationship <- as.factor(kin.aa2.wp$Relationship)
summary(kin.aa2.wp$Relationship)
kin.aa2.ca <- read.table("../results/kinship/aa2-ca/aa2-ca_1div_nc_20_0.2.kin.txt", sep = ",",
                         header = TRUE)
kin.aa2.ca$Relationship <- as.factor(kin.aa2.ca$Relationship)
summary(kin.aa2.ca$Relationship)
kin.aa2.sb <- read.table("../results/kinship/aa2-sb/aa2-sb_1div_nc_20_0.2.kin.txt", sep = ",",
                         header = TRUE)
kin.aa2.sb$Relationship <- as.factor(kin.aa2.sb$Relationship)
summary(kin.aa2.sb$Relationship)
kin.aa2.sq <- read.table("../results/kinship/aa2-sq/aa2-sq_1div_nc_20_0.2.kin.txt", sep = ",",
                         header = TRUE)
kin.aa2.sq$Relationship <- as.factor(kin.aa2.sq$Relationship)
summary(kin.aa2.sq$Relationship)

# ah ####
# ah1 - half
kin.ah1.ca <- read.table("../results/kinship/ah1-ca/ah1-ca_1div_nc_20_0.2.kin.txt", sep = ",",
                         header = TRUE)
# ah2 - half
kin.ah2.wp <- read.table("../results/kinship/ah2-wp/ah2-wp_1div_nc_20_0.2.kin.txt", sep = ",",
                         header = TRUE)
kin.ah2.sb <- read.table("../results/kinship/ah2-sb/ah2-sb_1div_nc_20_0.2.kin.txt", sep = ",",
                         header = TRUE)

# ah3 - half
kin.ah3.wp <- read.table("../results/kinship/ah3-wp/ah3-wp_1div_nc_20_0.2.kin.txt", sep = ",",
                         header = TRUE)

# Colours
#AA1        #AA2        #AH1       #AH2      #AH3      #AL1      #AL2
colours <- c("#274e13", "#8fce00", "#b06100", "#ffa500", "#ff4e00", "#6f0641", "#7223f0")

# Sort data  ####
metadata <- metadata %>% separate(Individual, into = c("Sample", "Species", "Site")) %>% 
  unite(Individual, Sample, Species, Site, remove = FALSE) %>% 
  dplyr::select(Individual, Sample, Species, Site, Loc, Depth, x, y, z)


# Calculate distances ####
## Parent Off plots ####
aa1.parent.df <- kinship_df(kin.aa1.wp[1:2,1:2], "Parent-Offspring", "AA1")
aa2.parent.df <- kinship_df(kin.aa2.wp[1:4,1:2], "Parent-Offspring", "AA2")

parents <- rbind(aa1.parent.df, aa2.parent.df)

kinship_plot(parents, max.dist = 5, max.count = 10, colours = colours[1:2], width = 5)


## Full sib plots ####
aa1.fullsib.df <- kinship_df(kin.aa1.wp[kin.aa1.wp$Relationship == "FS",1:2], "Full-siblings", "AA1")
aa2.wp.fullsib.df <- kinship_df(kin.aa2.wp[kin.aa2.wp$Relationship == "FS",1:2], "Full-siblings", "AA2")
aa2.ca.fullsib.df <- kinship_df(kin.aa2.ca[kin.aa2.ca$Relationship == "FS",1:2], "Full-siblings", "AA2")
aa2.sb.fullsib.df <- kinship_df(kin.aa2.sb[kin.aa2.sb$Relationship == "FS",1:2], "Full-siblings", "AA2")
aa2.sq.fullsib.df <- kinship_df(kin.aa2.sq[kin.aa2.sq$Relationship == "FS",1:2], "Full-siblings", "AA2")

fullsibs <- rbind(aa1.fullsib.df, aa2.ca.fullsib.df, aa2.sb.fullsib.df, aa2.sq.fullsib.df, aa2.wp.fullsib.df)
fullsibs.o <- na.omit(fullsibs)
kinship_plot(fullsibs.o, max.dist = 30, max.count = 10, colours = colours[c(2,4)], width = 10)

## Half-sib plots ####
aa1.halfsib.df <- kinship_df(kin.aa1.wp[kin.aa1.wp$Relationship == "HS",1:2], "Half-siblings", "AA1")
aa2.wp.halfsib.df <- kinship_df(kin.aa2.wp[kin.aa2.wp$Relationship == "HS",1:2], "Half-siblings", "AA2")
aa2.ca.halfsib.df <- kinship_df(kin.aa2.ca[kin.aa2.ca$Relationship == "HS",1:2], "Half-siblings", "AA2")
aa2.sb.halfsib.df <- kinship_df(kin.aa2.sb[kin.aa2.sb$Relationship == "HS",1:2], "Half-siblings", "AA2")
aa2.sq.halfsib.df <- kinship_df(kin.aa2.sq[kin.aa2.sq$Relationship == "HS",1:2], "Half-siblings", "AA2")
ah1.halfsib.df <- kinship_df(kin.ah1.ca[kin.ah1.ca$Relationship == "HS",1:2], "Half-siblings", "AH1")
ah2.wp.halfsib.df <- kinship_df(kin.ah2.wp[kin.ah2.wp$Relationship == "HS",1:2], "Half-siblings", "AH2")
ah2.sb.halfsib.df <- kinship_df(kin.ah2.sb[kin.ah2.sb$Relationship == "HS",1:2], "Half-siblings", "AH2")
ah3.halfsib.df <- kinship_df(kin.ah3.wp[kin.ah3.wp$Relationship == "HS",1:2], "Half-siblings", "AH3")

halfsibs <- rbind(aa1.halfsib.df, aa2.ca.halfsib.df, aa2.sb.halfsib.df, aa2.sq.halfsib.df, aa2.wp.halfsib.df, ah1.halfsib.df,
                  ah2.wp.halfsib.df, ah2.sb.halfsib.df, ah3.halfsib.df)
halfsibs.o <- na.omit(halfsibs)
kinship_plot(halfsibs.o, max.dist = 30, max.count = 10, colours = colours[1:4], width = 10)


# @Title: Admixture plots
# @Author: Katharine Prata
# @Date created: 07/04/21
# @Last edit: 06/08/21
# @Name of file: admix_plots2.R
# @Description: Takes results from admixture and plot
# @TODO: Make work with custom colour file!

## Dependencies ===================================================
library(ggplot2) # version 3.4.0
library(dplyr) # version 1.0.9
library(tidyr) # version 1.2.0

## Functions ===================================================
organise_data <- function(qfile,popfile,clusters){
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
admixture_plot <- function(plotdata,title){
  plot <- ggplot(plotdata, aes(fill = Clusters, y = Assignment, x = Individuals)) + 
    geom_bar(position = "fill", stat = "identity", width = 1) +
    facet_grid(~ pop, space = "free_x", 
               scales = "free_x", switch = "x") +
    guides(fill = "none") +
    ggtitle(paste(title)) +
    scale_fill_manual(name = "Cluster", values = colours) +
    theme(plot.title = element_text(size = 6),
          plot.background = element_blank(),
          legend.background = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 6),
          axis.ticks = element_blank(),
          panel.spacing.x = unit(0.1, "lines"),
          panel.spacing.y = unit(-2, "lines"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(fill = NA, 
                                      colour = "black", 
                                      size = 0.2, 
                                      linetype = "solid"),
          strip.text.x = element_text(size = 4, angle = 90))
  return(plot)
}
admixture_plot_species <- function(plotdata, title) {
  plot <- ggplot(plotdata, aes(fill = Clusters, y = Assignment, x = Individuals)) + 
    geom_bar(position = "fill", stat = "identity", width = 1) +
    facet_grid(~ pop, space = "free_x", 
               scales = "free_x", switch = "x") +
    guides(fill = "none") +
    ggtitle(paste(title)) +
    scale_fill_manual(name = "Cluster", values = colours) +
    theme(plot.title = element_text(size = 6),
          plot.background = element_blank(),
          legend.background = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 6),
          axis.ticks = element_blank(),
          panel.spacing.x = unit(-0.05, "lines"),
          panel.spacing.y = unit(-2, "lines"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          strip.background = element_blank(),
          panel.border = element_blank(),
          strip.text.x = element_text(size = 4, angle = 90))
  return(plot)
}

## Arguments ===================================================
args = commandArgs(TRUE)
VCF_NAME <- as.character(args[1]) 
STAGE <- args[2]
PARAMS <- args[3]
MISS_DAT <- args[4]
NUM_K <- as.numeric(args[5])
sorted <- args[6]
  
# Start cluster list
Clusters <- "Clust1"

# Import data ===================================================
popfile <- read.table(paste0("../data/pop_", VCF_NAME, "_", STAGE, "_", PARAMS, ".txt"), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(popfile) <- c("Individual", "pop")
# Note: altered individual species popfiles to be only at sites and changed NAs to appropriate location
colours <- rainbow(NUM_K, alpha = 0.9)

# Specify colours ands levels ===================================================
if (VCF_NAME == "all-aga") {
  colours <- c("#2D4E24", #agaricites
               "#B06327", #humilis
               "#6F1242", #lamarcki
               "#6652A2", #lamarcki2
               "#8FC73E", #agarcites2 
               "#821F27") #humilis2)
  if (STAGE == '1div' & PARAMS == "nc-wnr") {
    colours <- c("#2D4E24", #agaricites
                 "#6F1242", #lamarcki
                 "#B06327", #humilis
                 "#6652A2", #lamarcki2
                 "#8FC73E", #agarcites2 
                 "#821F27") #humilis2)
  }
  popfile <- popfile %>% separate(pop, into = c("Species", "Site"), remove = TRUE)
  popfile$Species[popfile$Species == "AC"] = "AA"
  popfile$Species[popfile$Species == "HU"] = "AH"
  popfile$Species[popfile$Species == "LM"] = "AL"
  # Unknowns
  popfile$Species[popfile$Individual == "KP0844_UN_NANA"] = "AA"
  popfile$Species[popfile$Individual == "KP0723_NA_NANA"] = "AL"
  popfile$Species[popfile$Individual == "KP1036_UN_CA10"] = "AH"
  popfile$Species[popfile$Individual == "KP1153_UN_CA20"] = "AH"
  popfile$Species[popfile$Individual == "KP1159_UN_CA20"] = "AH"
  popfile$Species[popfile$Individual == "KP1156_UN_CA20"] = "AA"
  popfile$Species[popfile$Individual == "KP1160_UN_CA20"] = "AL"
  popfile$Site[popfile$Individual == "KP0406_NA_NANA"] = "WP20"
  popfile$Site[popfile$Individual == "KP0673_NA_NANA"] = "SB20"
  popfile$Site[popfile$Individual == "KP0844_UN_NANA"] = "WP10"
  popfile$Site[popfile$Individual == "KP0723_NA_NANA"] = "SB20"
  # mislabels agaricites, humilis, lamarcki
  mislabels_AA = c("KP1139_LM_CA20", "KP0361_LM_WP20", "KP0349_LM_WP20","KP0611_LM_SB20",
                   "KP0399_LM_WP20", "KP1146_LM_CA20", "KP0109_LM_SB10", "KP0583_LM_WP20",
                   "KP0244_LM_SQ20", "KP0565_LM_WP20", "KP0406_NA_NANA", "KP0673_NA_NANA")
  for (a in 1:length(mislabels_AA)) {
    popfile$Species[popfile$Individual == mislabels_AA[a]] = "AA"}
  mislabels_AH = c("KP0639_LM_SB20", "KP0899_LM_WP10")
  for (h in 1:length(mislabels_AH)) {
    popfile$Species[popfile$Individual == mislabels_AH[h]] = "AH"}
  mislabels_AL = c("KP1124_AC_CA20", "KP0557_AC_WP20", "KP0203_AC_SQ12", "KP1095_HU_CA20",
                   "KP0186_AC_SQ12", "KP0851_AC_WP10", "KP0605_AC_SB20", "KP0609_AC_SB20",
                   "KP0096_AC_SB10", "KP0765_AC_WP10", "KP0192_AC_SQ12", "KY0578_AC_WP20",
                   "KP0671_AC_SB20", "KP0108_AC_SB10")
  for (l in 1:length(mislabels_AL)) {
    popfile$Species[popfile$Individual == mislabels_AL[l]] = "AL"}
  
  
  popfile <- popfile %>% unite(pop, Species, Site, sep = "_", remove = TRUE)
  
  if (PARAMS == "nc-wnr") {
   popfile$pop[popfile$pop == "G1_CS50"] = "AG1_50"
   popfile$pop[popfile$pop == "G2_BR50"] = "AG2_50"
   
   popfile$pop[popfile$pop == "L1_CS50"] = "AL1_50"
   popfile$pop[popfile$pop == "L2_CS15"] = "AL2_15"
   popfile$pop[popfile$pop == "L1_CK50"] = "AL1_50"
   popfile$pop[popfile$pop == "L2_CR15"] = "AL2_15"
   popfile$pop[popfile$pop == "L1_CR15"] = "AL1_15"
  
   levels <- c("AA_WP05", "AA_WP10", "AA_WP20", "AA_CA05", "AA_CA10", "AA_CA20", "AA_SB05", "AA_SB10", "AA_SB20", "AA_SQ12", "AA_SQ20",
               "AH_WP05", "AH_WP10", "AH_WP20", "AH_CA05", "AH_CA10", "AH_CA20", "AH_SB05", "AH_SB10", "AH_SB20", "AH_SQ12",
               "AL_WP10", "AL_WP20", "AL_CA10", "AL_CA20", "AL_SB10", "AL_SB20", "AL_SQ12", "AL_SQ20", "AL1_15", "AL1_50",
               "AL2_15", "AG1_50", "AG2_50")
  } else {
  levels <- c("AA_WP05", "AA_WP10", "AA_WP20", "AA_CA05", "AA_CA10", "AA_CA20", "AA_SB05", "AA_SB10", "AA_SB20", "AA_SQ12", "AA_SQ20",
              "AH_WP05", "AH_WP10", "AH_WP20", "AH_CA05", "AH_CA10", "AH_CA20", "AH_SB05", "AH_SB10", "AH_SB20", "AH_SQ12",
              "AL_WP10", "AL_WP20", "AL_CA10", "AL_CA20", "AL_SB10", "AL_SB20", "AL_SQ12", "AL_SQ20")
  }
  
  

  
  
} else if (VCF_NAME == "ac") {
  # green-aqua colour pallete
  colours <- c("#2D4E24", "#8FC73E", "#43BB93", "#347F66")
  levels <- c("WP05", "WP10", "WP20", "CA05", "CA10", "CA20", "SB05", "SB10", "SB20", "SQ12", "SQ20")
  popfile$pop[popfile$Individual == "KP0844_UN_NANA"] = "WP10"
  popfile$pop[popfile$Individual == "KP0406_NA_NANA"] = "WP20"
  } else if (VCF_NAME == "hu") {
  # brown-red colour pallete
  colours <- c("#B06327", "#FAA41A", "#F05123", "#821F27")
  levels <- c("WP05", "WP10", "WP20", "CA05", "CA10", "CA20", "SB05", "SB10", "SB20", "SQ12")
} else if (VCF_NAME == "lm-gr") {
  # lamarcki colour palette OG purple pink
  colours <-  c("#6F1242", "#6652A2")
} else if (VCF_NAME == "lm" & (STAGE == "1d-pure" | PARAMS == "nc-wnr")) {
  colours <- c("#6F1242", "#6652A2", "red", "purple")
  popfile$pop[popfile$Individual == "NR5242_L1_CK50" | 	popfile$Individual == "NR4671_L1_CS50"] = "AL1-50"
  popfile$pop[popfile$Individual == "NR5534_L1_CR15"] = "AL1-15"
  popfile$pop[popfile$Individual == "NR5533_L2_CR15" | 	popfile$Individual == "NR5050_L2_CS15" | popfile$Individual == "NR5053_L2_CS15"] = "AL2-15"
  popfile$pop[popfile$Individual == "KP0723_NA_NANA"] = "SB20"
  levels <- c("WP10", "WP20", "CA10", "CA20", "SB10", "SB20", "SQ12", "SQ20", "AL1-15", "AL1-50", "AL2-15")
} else if (VCF_NAME == "lm" & STAGE != "1d-pure") {
  # lamarcki colour OG purple pink
  colours <- c("#6F1242", "#6652A2", "red", "purple")
  levels <- c("WP10", "WP20", "CA10", "CA20", "SB10", "SB20", "SQ12", "SQ20")
  popfile$pop[popfile$Individual == "KP0723_NA_NANA"] = "SB20"
} else {
  print("colour not chosen")
  if (VCF_NAME == "ah3") {
    levels <- c("WP05", "WP10", "CA05", "SB05")
  } else if (VCF_NAME == "ah1") {
    levels <- c("WP05", "WP10", "WP20", "CA05", "CA10", "CA20", "SB05", "SB10", "SB20", "SQ12")
  } else if (VCF_NAME == "ah2") {
    levels <- c("WP05", "WP10", "SB05", "SB10")
    
  }
  
}
popfile[,2] <- factor(popfile[,2], levels = levels)


for (i in 2:NUM_K) {
  if (sorted == "yes") {
    Q <- read.table(paste0("../results/admix_runs/", MISS_DAT, "percent_", VCF_NAME, "_", STAGE, "_", PARAMS, '/sortedQ/', VCF_NAME, "_2bi_", STAGE, "_", PARAMS, "_",
                           MISS_DAT, ".", i, ".Q"), header = TRUE)
  } else {
    Q <- read.table(paste0("../results/admix_runs/", MISS_DAT, "percent_", VCF_NAME, "_", STAGE, "_", PARAMS, '/', VCF_NAME, "_2bi_", STAGE, "_", PARAMS, "_",
                           MISS_DAT, ".", i, ".Q"), header = FALSE)
  }
# Specify clusters ===================================================
Clusters <- append(Clusters, paste0("Clust", i))



# Data organisation ===================================================
data <- organise_data(Q, popfile, Clusters)
write.csv(data, paste0("../data/", VCF_NAME, "_", STAGE, "_", PARAMS, "_", NUM_K, "_admix_data.csv"), quote = FALSE, row.names = FALSE)
title <- paste0(MISS_DAT, "% missing data, K = ", i)

# Plot data ===================================================
plot <- admixture_plot(data, title)
ggsave(paste0(VCF_NAME, "_", STAGE, "_", PARAMS, "_", i, "_admixture.pdf"), plot = plot, width = 22, height = 4,
              units = "cm",path = "../results/admix_plots", limitsize = FALSE)
if (VCF_NAME == "all-aga") {
  plot2 <- admixture_plot_species(data, title)
  ggsave(paste0(VCF_NAME, "_", STAGE, "_", PARAMS, "_", i, "_admixture_species.pdf"), plot = plot2, width = 22, height = 4,
       units = "cm", path = "../results/admix_plots", limitsize = FALSE)
}
}


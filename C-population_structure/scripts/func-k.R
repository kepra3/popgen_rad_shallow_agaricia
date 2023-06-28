# Title: Functions to make stuff easier
# Author: Katharine Prata
# Date created: 23/02/2021
# Last edit: 23/03/20201
# Description: Script pronounced "funk" "kay" ;)

# Functions:

# As factor a column  or multiple
k_factor <- function(df, col) {
  df$col <- as.factor(df$col)
  return(df)
}

# renaming sample names!
rename_samples <- function(pop_names) {
  pop_names <- pop_names %>% separate(V1, into = c("Sample_no", "Species", "Site"), sep = "_",
                            remove = TRUE, extra = "merge") %>%
  separate(Site, into = c("X", "Site"), fill = "left", remove = TRUE) %>% 
  separate(Site, into = c("Loc", "Depth"), sep = 2, remove = TRUE)  
  pop_names$Depth[pop_names$Depth == "5"] = "05"
  pop_names$Species[pop_names$Species == "AL1"] = "L1"
  pop_names$Species[pop_names$Species == "AL2"] = "L2"
  pop_names$Species[pop_names$Species == "CA"] = "AC"
  pop_names <- pop_names %>% unite(Sample_name, Sample_no, Species, Loc, Depth, sep = "_",
                                 remove = TRUE) %>% 
  select(Sample_name)
  return(pop_names)
} 

# Organise admixture Q file data into long format for plotting with barplot
organise_data <- function(qfile,popfile,clusters){
  data <- as.data.frame(cbind(popfile,qfile))
  clusters <- expand.grid(popfile[,1], clusters)
  n <- ncol(data)
  Assignment_t <- as.matrix(data[,3:n])
  Assignment <- matrix(Assignment_t, ncol = 1)
  plot.data <- cbind(clusters,Assignment)
  plot.data <- cbind(plot.data,popfile[,2])
  colnames(plot.data) <- c("Individuals","Clusters","Assignment","pop")
  return(plot.data)
}

# Using an organised admixture file choosing highest admixture assignment as the cluster
# Fave function <3
sort_clusters <- function(Q, k) {
  sorted = rep(1, length(Q[,1]))
  names = colnames(Q)
  for (i in 1:length(Q[,1])) {
      number = max(Q[i,1:k])
      sorted[i] = names[Q[i,] == number]}
  return(sorted)}


choose_Qprop <- function(popQ, k, threshold) {
  popQ$Taxa <- NA
  for (i in 1:k) {
  num <- length(popQ) - k - 1
  popQ$Taxa[popQ[, num + i] > threshold] <- paste0("Clust", i)
  } 
  return(popQ)
}

# Title: Functions to make stuff easier
# Author: Katharine Prata
# Date created: 23/02/2021
# Last edit: 23/03/2021
# Description: Script pronounced "funk" "kay" ;)

# Functions:

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

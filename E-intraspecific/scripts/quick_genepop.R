# Title: quick gene pop
# Author: Katharine Prata
# Date created: 15/12/23

# Required packages
# R v4.2.0
library(genepop) # v1.1.7

# Arguments ==================================================
args = commandArgs(TRUE)
taxa <- args[1]
category <- args[2]
scale <- args[3]

# Gene pop  ==================================================
if (scale == "all") {
  minD = 1e-03
  maxD = 1e+05
  dimension = "1D"
} else if (scale == "within") {
  minD = 1e-03
  maxD = 1e+02
  dimension = "2D"
} else if (scale == "between") {
  minD = 1e+02
  maxD = 1e+05
  dimension = "1D"
}
ibd(paste0("../data/", taxa, "_", category, ".genepop.txt"),
    outputFile = paste0("../results/ibd/", taxa, "/", taxa, "_", category, ".",
                        scale, ".", dimension, ".results.txt"),
    settingsFile = "",
    dataType = "Diploid",
    statistic = "a",
    geographicScale = dimension,
    CIcoverage = 0.95,
    testPoint = 0,
    minimalDistance = minD,
    maximalDistance = maxD,
    mantelPermutations = 1000,
    mantelRankTest = FALSE,
    verbose = interactive())

# Ancestry proportions plots
# Created by: Katharine Prata

# Packages
# R v4.2.0
library(stringr) # v1.4.0

# Arguments
args = commandArgs(TRUE)
VCF_NAME = args[1]
K = args[2]

DATA_PARAMS <- str_split_fixed(VCF_NAME, "_", n = 5)

if (DATA_PARAMS[1] == "hu" & DATA_PARAMS[3] == "1div") {
  colour_palette <- c("#B06327","#FAA41A", "#F05123", "#821F27")
  species = c("AH1", "AH2", "AH3")
  main = "A. humilis"
  threshold = 0.8
  x_pos = 40
} else if (DATA_PARAMS[1] == "ac" & DATA_PARAMS[3] == "1div") {
  colour_palette <- c("#2D4E24", "#8FC73E", "#43BB93", "#347F66")
  species = c("AA1", "AA2", "AA3", "AA4")
  main = "A. agaricites"
  threshold = 0.8
  x_pos = 170
} else if (DATA_PARAMS[1] == "lm" & DATA_PARAMS[3] == "1div") {
  colour_palette <- c("#6F1242", "#6652A2", "red")
  species = c("AL1", "AL2")
  main = "A. lamarcki"
  threshold = 0.8
  x_pos = 40
} else if (DATA_PARAMS[1] == "all-aga") {
  colour_palette <- c("#2D4E24", #agaricites
                      "#6F1242", #lamarcki
                      "#B06327") #humilis
  species = c("AA", "AL", "AH")
  main = "All individuals"
  threshold = 0.8
  x_pos = 250
}

y_pos = threshold + 0.05

for (i in 2:K) {
  Qfile <- read.delim(paste0("../results/admix_runs/", DATA_PARAMS[5], "percent", "_", DATA_PARAMS[1], "_", DATA_PARAMS[3], "_", 
                         DATA_PARAMS[4], "/sortedQ/", DATA_PARAMS[1], "_", DATA_PARAMS[2], "_",
                         DATA_PARAMS[3], "_", DATA_PARAMS[4], "_", DATA_PARAMS[5], ".", i, ".Q"))
  for (j in 1:i) {
    pdf(file = paste0("../results/admix_plots/", VCF_NAME, "_K", i, "_Clust", j, "Q.pdf"))
    plot(sort(Qfile[, j]), xlab = "Sample",
         ylab = paste("Admixture proportion for", species[j]),
         main = paste0(main, ", K = ", i),
         col = colour_palette[j])
    dev.off()
}}

for (j in 1:K) {
  pdf(file = paste0("../results/admix_plots/", VCF_NAME, "_K", K, "_", species[j], "Q.pdf"), height = 4, width = 4)
  par(mgp = c(1.8, 0.5, 0),mar = c(5, 4, 4, 2) + 0.1)
  plot(sort(Qfile[, j]), xlab = "Sample",
       ylab = paste("Admixture proportion for", species[j]),
       main = paste0(main, ", K = ", K),
       col = colour_palette[j])
  abline(h = threshold, col = 'red', lwd = 3, lty = 2)
  text(x_pos, y_pos, paste("Taxa threshold =", threshold))
  dev.off()}

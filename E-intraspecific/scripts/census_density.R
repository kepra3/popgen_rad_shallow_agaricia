# Title: Census per plot
# Author: Katharine Prata
# Date created: 29/04/24

# In loop of taxa
for (taxa in c("AA1", "AA2", "AH1", "AH2", "AH3", "AL1", "AL2")) {
  print(taxa)
  category = "all"
  scale = "within"
  
  # Choose species
  if (taxa == "AA1" | taxa == "AA2") {
    vcf = "ac_3b_nc_20"
  } else if (taxa == "AL1" | taxa == "AL2") {
    vcf = "lm_3b_nc_20"
  } else if (taxa == "AH1" | taxa == "AH2" | taxa == "AH3") {
    vcf = "hu_3b_nc_20"
  }
  
  # Import data
  taxa.depth.pop <- read.csv(paste0("../data/", vcf, "_", taxa, "_",
                                    category, "_.csv"),
                             stringsAsFactors = TRUE)
  
  # Number of samples per site that have been genotyped and annotated
  taxa.depth.pop$Site <- as.factor(taxa.depth.pop$Site)
  summary(taxa.depth.pop$Site)
  #total_samples = sum(summary(taxa.depth.pop$Site))
  summary_samples <- table(taxa.depth.pop$Site)
  # genotyping rate = 2/3
  N <- summary_samples / (2/3)
  # Annotation rate & De effective density - see density_calcs.R
  if (taxa == "AA1") {
    annotation_rate = 0.8571429
  } else if (taxa == "AA2") {
    annotation_rate = 0.8233333
  } else if (taxa == "AH1") {
    annotation_rate = 0.5211268
  } else if (taxa == "AH2") {
    annotation_rate = 0.7241379
  } else if (taxa == "AH3") {
    annotation_rate = 0.7647059
  } else if (taxa == "AL1") {
    annotation_rate = 0.6896552
  } else if (taxa == "AL2") {
    annotation_rate = 0.7580645
  }
  # Adjust by annotation rate
  N <- N / annotation_rate
  exhaustive_sample <- c(N[names(N) == "SQ20"], N[names(N) == "SQ10"],
                           N[names(N) == "WP05"], N[names(N) == "CA05"],
                           N[names(N) == "SB05"])
  # Adjust by discovery rate 70%
  N[names(N) != names(exhaustive_sample)] <- N[names(N) != names(exhaustive_sample)]/0.7
  
  # D individuals per plot within 25 x 2 = 50m2
  D <- N / 50 # individuals per m2
  D
  print(mean(D))
  write.csv(D, paste0("../data/", taxa,
                      ".census.txt"), quote = FALSE, row.names = FALSE)
  print(round(N, 1))
  print(round(D, 2))
}

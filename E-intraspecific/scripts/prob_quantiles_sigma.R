# Title: Probability quantiles of sigma
# Author: Katharine Prata
# Date created: 8/05/24


for (taxa in c("AA1", "AA2", "AH1", "AH3")) {
  # loiselle
  l.sigma.dne <- read.csv(paste0("../results/distributions/",
                                 taxa, ".sigma_DNE-loiselle.txt"))
  sigma_quantiles <- as.data.frame(round(quantile(l.sigma.dne$x,
                                          probs = c(0.10, 0.25, 0.5, 0.75, 0.9)),
                                 2))
  colnames(sigma_quantiles) <- taxa
  write.csv(sigma_quantiles, paste0("../results/distributions/",
            taxa, ".sigma_DNE-loiselle_quantiles.txt"), quote = FALSE)
  
  l.sigma.dnc <- read.csv(paste0("../results/distributions/",
                                 taxa, ".sigma_DNC-loiselle.txt"))
  sigma_quantiles <- as.data.frame(round(quantile(l.sigma.dnc$V1,
                                                  probs = c(0.10, 0.25, 0.5, 0.75, 0.9)),
                                         2))
  colnames(sigma_quantiles) <- taxa
  write.csv(sigma_quantiles, paste0("../results/distributions/",
                                    taxa, ".sigma_DNC-loiselle_quantiles.txt"), quote = FALSE)
  
  if (taxa == "AA1" | taxa == "AA2") {
    # rousset
    r.sigma.dne <- read.csv(paste0("../results/distributions/",
                                   taxa, ".sigma_DNE-rousset.txt"))
    sigma_quantiles <- as.data.frame(round(quantile(l.sigma.dne$x,
                                                    probs = c(0.10, 0.25, 0.5, 0.75, 0.9)),
                                           2))
    colnames(sigma_quantiles) <- taxa
    write.csv(sigma_quantiles, paste0("../results/distributions/",
                                      taxa, ".sigma_DNE-rousset_quantiles.txt"), quote = FALSE)
    r.sigma.dnc <- read.csv(paste0("../results/distributions/",
                                   taxa, ".sigma_DNC-rousset.txt"))
    sigma_quantiles <- as.data.frame(round(quantile(l.sigma.dnc$V1,
                                                    probs = c(0.10, 0.25, 0.5, 0.75, 0.9)),
                                           2))
    colnames(sigma_quantiles) <- taxa
    write.csv(sigma_quantiles, paste0("../results/distributions/",
                                      taxa, ".sigma_DNC-rousset_quantiles.txt"), quote = FALSE)
  }}
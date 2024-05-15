# Title: effective densities
# Author: Katharine Prata
# Date created: 29/04/24

# Kinship results summarised in Supplementry Table 10
aa1 <- c("AA1", "pop1", 32, 56, 159, 528)
aa2.1 <- c("AA2", "pop1", 547, 707, 955, 1328)
aa2.2 <- c("AA2", "pop2", 222, 341, 667, 1328)
aa2.3 <- c("AA2", "pop3", 328, 483, 8320, 1328)
aa2.4 <- c("AA2", "pop4", 106, 193, 603, 528)
ah1 <- c("AH1", "pop1", 154, 316, 5415, 1328) #census equiv 485
# For AH2 and AH3 upper limit = Inf,
# Let's just go with AH1 upper limit incase taxa are similar
ah2.1 <- c("AH2", "pop1", 121, 291, 5415, 850) #census equiv 234
ah2.2 <- c("AH2", "pop2", 110, 338, 5415, 850) #census equiv 234
ah3 <- c("AH3", "pop1", 103, 334, 5415, 850) #census equiv 155

effective.densities <- rbind(aa1, aa2.1, aa2.2, aa2.3,
                             aa2.4, ah1, ah2.1, ah2.2, ah3)

effective.densities <- as.data.frame(effective.densities)

colnames(effective.densities) <- c("taxa",
                                  "pop",
                                  "Ne_low",
                                  "Ne",
                                  "Ne_high",
                                  "Area")
write.csv(effective.densities,
          "../data/agaricia_effective_density.txt",
          quote = FALSE, row.names = FALSE)

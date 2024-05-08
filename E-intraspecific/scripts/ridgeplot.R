
# Packages
library(ggplot2)
library(ggridges)
library(tidyverse)

combined_data <- data.frame(V1 = numeric(0),
                            Taxa = numeric(0),
                            Density = numeric(0))

for (taxa in c("AA1", "AA2", "AH1", "AH3")) {
  sigma_DNC <- read.csv(paste0("~/git/CalcluatingSigmaWithError/agaricia/",
                               taxa, ".sigma_DNC.txt"))
  sigma_DNC$Taxa <- taxa
  sigma_DNC$Density <- "census"
  sigma_DNE <- read.csv(paste0("~/git/CalcluatingSigmaWithError/agaricia/",
                              taxa, ".sigma_DNE.txt"))
  sigma_DNE$Taxa <- taxa
  sigma_DNE$Density <- "effective"
  colnames(sigma_DNE) <- c("V1", "Taxa", "Density") 
  combined_data <- rbind(combined_data, sigma_DNC, sigma_DNE)
  rm(sigma_DNC, sigma_DNE)
}
str(combined_data)
combined_data$Taxa <- factor(combined_data$Taxa, levels = c("AH3", "AH1", "AA2", "AA1"))
combined_data$Density <- as.factor(combined_data$Density)

# Calculate quantiles
quantiles <- combined_data %>%
  group_by(Taxa, Density) %>%
  summarize(q10 = quantile(log10(V1), probs = 0.10),
            q25 = quantile(log10(V1), probs = 0.25),
            q50 = quantile(log10(V1), probs = 0.50),
            q75 = quantile(log10(V1), probs = 0.75),
            q90 = quantile(log10(V1), probs = 0.90))

# Create ridge plot
p <- ggplot(combined_data, aes(x = log10(V1), y = Taxa, fill = Density)) +
  geom_density_ridges(alpha = 0.6) +
  theme_ridges() +
  scale_fill_manual(values = c("census" = "blue", "effective" = "red")) +
  xlab(expression(log[10](sigma) ~ "m"))

quantiles <- quantiles[quantiles$Density == "census",]

# Add vertical lines for quantiles
p + geom_point(data = quantiles, aes(x = q10, y = Taxa), color = "black") +
  geom_point(data = quantiles, aes(x = q50, y = Taxa), color = "black") +
  geom_point(data = quantiles, aes(x = q90, y = Taxa), color = "black") +
  geom_text(data = quantiles, aes(x = q10, label = sprintf("%.0f", 10^q10),
                                  y = Taxa), color = "black",
            vjust = 1.5, hjust = 1.2, size = 2, fontface = "bold") +
  geom_text(data = quantiles, aes(x = q50, label = sprintf("%.0f", 10^q50),
                                  y = Taxa), color = "black",
            vjust = 1.5, hjust = -0.1, size = 2, fontface = "bold") +
  geom_text(data = quantiles, aes(x = q90, label = sprintf("%.0f", 10^q90),
                                  y = Taxa), color = "black",
            vjust = 1.5, hjust = -0.1, size = 2, fontface = "bold") +
  theme(text = element_text(size = 8), # Set default text size to 8 pt
        axis.title = element_text(size = 10),  # Set axis title size to 10 pt
        axis.text.y = element_text(size = 10),  # Set y-axis label size to 10 pt
        axis.text.x = element_text(size = 10),  # Set x-axis label size to 10 pt
        strip.text = element_text(size = 10),  # Set facet labels (taxa labels) size to 10 pt
        legend.title = element_text(size = 10),  # Set legend title size to 10 pt
        )
ggsave("~/git/CalcluatingSigmaWithError/agaricia_ridge.pdf",
       units = "cm", height = 10, width = 16, dpi = 400)

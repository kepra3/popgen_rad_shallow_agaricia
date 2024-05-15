
# Packages
library(ggplot2)
library(ggridges)
library(tidyverse)

combined_data <- data.frame(V1 = numeric(0),
                            Taxa = numeric(0),
                            Density = numeric(0),
                            TDensity = numeric(0))

for (taxa in c("AA1", "AA2", "AH1", "AH3")) {
  sigma_DNC <- read.csv(paste0("../results/distributions/",
                               taxa, ".sigma_DNC-loiselle.txt"))
  sigma_DNC$Taxa <- taxa
  sigma_DNC$Density <- "census"
  sigma_DNC$TDensity <- paste0(taxa, "census")
  sigma_DNE <- read.csv(paste0("../results/distributions/",
                              taxa, ".sigma_DNE-loiselle.txt"))
  sigma_DNE$Taxa <- taxa
  sigma_DNE$Density <- "effective"
  sigma_DNE$TDensity <- paste0(taxa, "effective")
  colnames(sigma_DNE) <- c("V1", "Taxa", "Density", "TDensity") 
  combined_data <- rbind(combined_data, sigma_DNC, sigma_DNE)
  rm(sigma_DNC, sigma_DNE)
}
str(combined_data)
combined_data$Taxa <- factor(combined_data$Taxa, levels = c("AH3", "AH1", "AA2", "AA1"))
combined_data$Density <- as.factor(combined_data$Density)
combined_data$TDensity <- as.factor(combined_data$TDensity)

# Calculate quantiles
quantiles <- combined_data %>%
  group_by(Taxa, Density, TDensity) %>%
  summarize(q10 = quantile(log10(V1), probs = 0.10),
            q25 = quantile(log10(V1), probs = 0.25),
            q50 = quantile(log10(V1), probs = 0.50),
            q75 = quantile(log10(V1), probs = 0.75),
            q90 = quantile(log10(V1), probs = 0.90))

# Colours
              #AA1        #AA2        #AH1       #AH3
colours <- c("#274e13", "#8fce00", "#b06100", "#ff4e00")
colours <- adjustcolor(colours, alpha.f = 0.9)
adjusted_colours <- adjustcolor(colours, alpha.f = 0.5)

# Create a vector to store both original and adjusted colours
all_colours <- rep(NA, length(colours) * 2)
all_colours[seq(1, length(all_colours), by = 2)] <- colours
all_colours[seq(2, length(all_colours), by = 2)] <- adjusted_colours

# Create ridge plot
p <- ggplot(combined_data, aes(x = log10(V1), y = Taxa, fill = TDensity)) +
  geom_density_ridges() +
  theme_ridges() +
  #scale_fill_manual(values = c("census" = "black", "effective" = "lightgray")) +
  scale_fill_manual(values = all_colours) +
  xlab(expression(log[10](sigma) ~ "m"))

quantiles <- quantiles[quantiles$Density == "census",]

# Add vertical lines for quantiles
p + geom_point(data = quantiles, aes(x = q10, y = Taxa), color = "black") +
  geom_point(data = quantiles, aes(x = q50, y = Taxa), color = "black") +
  geom_point(data = quantiles, aes(x = q90, y = Taxa), color = "black") +
  geom_text(data = quantiles, aes(x = q10, label = sprintf("%.0f", 10^q10),
                                  y = Taxa), color = "black",
            vjust = -0.8, hjust = 0.8, size = 2, fontface = "bold") +
  geom_text(data = quantiles, aes(x = q50, label = sprintf("%.0f", 10^q50),
                                  y = Taxa), color = "black",
            vjust = -0.8, hjust = -0.1, size = 2, fontface = "bold") +
  geom_text(data = quantiles, aes(x = q90, label = sprintf("%.0f", 10^q90),
                                  y = Taxa), color = "black",
            vjust = -0.8, hjust = -0.1, size = 2, fontface = "bold") +
  theme(text = element_text(size = 9), # Set default text size to 8 pt
        axis.title = element_text(size = 10),  # Set axis title size to 10 pt
        axis.text.y = element_text(size = 10),  # Set y-axis label size to 10 pt
        axis.text.x = element_text(size = 10),  # Set x-axis label size to 10 pt
        strip.text = element_text(size = 10),  # Set facet labels (taxa labels) size to 10 pt
        legend.title = element_text(size = 10),  # Set legend title size to 10 pt
        )
ggsave("../results/distributions/agaricia_ridge_colour_cp.pdf",
       units = "cm", height = 10, width = 16, dpi = 400)

# Re-calculate quantiles
quantiles <- combined_data %>%
  group_by(Taxa, Density, TDensity) %>%
  summarize(q10 = quantile(V1, probs = 0.10),
            q25 = quantile(V1, probs = 0.25),
            q50 = quantile(V1, probs = 0.50),
            q75 = quantile(V1, probs = 0.75),
            q90 = quantile(V1, probs = 0.90))

# Metre scale
combined_data <- combined_data[combined_data$V1 < 10^2,]

# Create ridge plot
p <- ggplot(combined_data, aes(x = V1, y = Taxa, fill = TDensity)) +
  geom_density_ridges() +
  theme_ridges() +
  #scale_fill_manual(values = c("census" = "black", "effective" = "lightgray")) +
  scale_fill_manual(values = all_colours) +
  xlab(expression(sigma ~ "(m)"))

quantiles <- quantiles[quantiles$Density == "census",]
quantiles$q90[quantiles$Taxa == "AA2"] <- NA

# Add vertical lines for quantiles
p + geom_point(data = quantiles, aes(x = q10, y = Taxa), color = "black") +
  geom_point(data = quantiles, aes(x = q50, y = Taxa), color = "black") +
  geom_point(data = quantiles, aes(x = q90, y = Taxa), color = "black") +
  geom_text(data = quantiles, aes(x = q10, label = round(q10, 0),
                                  y = Taxa), color = "black",
            vjust = -0.8, hjust = 0.8, size = 2, fontface = "bold") +
  geom_text(data = quantiles, aes(x = q50, label = round(q50, 0),
                                  y = Taxa), color = "black",
            vjust = -0.8, hjust = -0.1, size = 2, fontface = "bold") +
  geom_text(data = quantiles, aes(x = q90, label = round(q90, 0),
                                  y = Taxa), color = "black",
            vjust = -0.8, hjust = -0.1, size = 2, fontface = "bold") +
  theme(text = element_text(size = 9), # Set default text size to 8 pt
        axis.title = element_text(size = 10),  # Set axis title size to 10 pt
        axis.text.y = element_text(size = 10),  # Set y-axis label size to 10 pt
        axis.text.x = element_text(size = 10),  # Set x-axis label size to 10 pt
        strip.text = element_text(size = 10),  # Set facet labels (taxa labels) size to 10 pt
        legend.title = element_text(size = 10),  # Set legend title size to 10 pt
  )
ggsave("../results/distributions/agaricia_ridge_colour_metres.pdf",
       units = "cm", height = 10, width = 16, dpi = 400)

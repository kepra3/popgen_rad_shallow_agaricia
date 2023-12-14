# Title: Isolation-by-distance plots
# Author: Katharine Prata

# Required packages
# R v4.2.0
library(ggplot2) # v3.4.0
library(spaa) # v0.2.2

# Arguments  ==================================================
args = commandArgs(TRUE)
#taxa = args[1] #"AA2"
#category = args[2] #"WP05"
#scale = args[3] # "within"
taxa = "AH3"
category = "all"
scale = "within"

# Import data ==================================================
a <- read.delim(paste0("../results/ibd/", taxa, "/", taxa, "_", category, ".a.txt"), header = FALSE)
lndist <- read.delim(paste0("../results/ibd/", taxa, "/", taxa, "_", category, ".lndist.txt"), header = FALSE)
results <- read.csv(paste0("../results/ibd/", taxa, "/",taxa, "_", category, ".", scale, ".results_short.txt"))
if (taxa == "AA1" | taxa == "AA2") {
  vcf = "ac_3b_nc_20"
} else if (taxa == "AL1" | taxa == "AL2") {
  vcf = "lm_3b_nc_20"
} else if (taxa == "AH1" | taxa == "AH2" | taxa == "AH3") {
  vcf = "hu_3b_nc_20"
}
taxa.depth.pop <- read.csv(paste0("../data/",vcf, "_", taxa, "_", category, "_.csv"), stringsAsFactors = TRUE)
table(taxa.depth.pop$Site)

results.round <- results
for (i in 1:7) {
  results.round[,i] <- round(results[,i], 3)
}

row1 <- rep(NA, nrow(a))
a <- rbind(row1, a)
lndist <- rbind(row1, lndist)

b = results$slope
b.low = results$s.lowCI
b.high = results$s.highCI
regression <- c(results$intercept, b)
confidence_interval_lower <- c(results$i.lowCI, b.low)
confidence_interval_upper <- c(results$i.highCI, b.high)
pvalue = results$p.slope
sign = "="

# Re-organise data ==================================================
x <- lndist[lower.tri(lndist)]
y <- a[lower.tri(a)]
dat <- as.data.frame(cbind(x, y))
colnames(dat) <- c("distance", "a")

## setting axes & scales ####
if (scale == "within") {
  dat <- dat[dat$distance < 8,]
  xlim <- c(-4, 4)
  if (taxa == "AA1" | taxa == "AA2") {
    xlim <- c(-4, 4.1)
  }
} else if (scale == "all") {
  xlim <- c(-4, 11)
}

if (taxa == "AH2" | taxa == "AH1") {
  ylim <- c(0, 4)
} else if (taxa == "AH3") {
  ylim <- c(0, 11)
} else {
    ylim <- c(-0.45, 0.45)
  }

regression <- as.data.frame(cbind(regression[1], regression[2]))
confidence_interval_lower <- as.data.frame(cbind(confidence_interval_lower[1], confidence_interval_lower[2]))
confidence_interval_upper <- as.data.frame(cbind(confidence_interval_upper[1], confidence_interval_upper[2]))

y_lower <- dat$distance * confidence_interval_lower$V2 + confidence_interval_lower$V1
y_mean <- dat$distance * regression$V2 + regression$V1
y_upper <- dat$distance * confidence_interval_upper$V2 + confidence_interval_upper$V1

# Make plots  ========================
slope <- round(regression[2], 3)
res <- data.frame(cond1 = "regression",
                  x = dat$distance,
                  y = y_mean,
                  ymin = y_lower,
                  ymax = y_upper)

rib <- geom_ribbon(data = res, aes(x = x, y = y, ymin = ymin, ymax = ymax,
                                   fill = cond1), fill = 'blue', alpha = 0.2)

p <- ggplot(dat, aes(distance, a)) +
  geom_point(shape = 21) +
  ylab("â") +
  xlab("log distance") +
  geom_abline(data = regression, aes(intercept = V1, slope = V2), colour = 'red') +
  geom_line(data = res, aes(x, ymin), colour = 'blue', linetype = "dashed") + 
  geom_line(data = res, aes(x, ymax), colour = 'blue', linetype = "dashed") + rib + theme_bw() + 
  ggtitle(paste0("Slope = ", slope, " p ", sign, " ", round(pvalue, 4))) +
  ylim(ylim) +
  xlim(xlim) +
  theme(plot.title = element_text(size = 10),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8)) 
p

ggsave(paste0("../results/ibd/plots/", taxa, "_", category, "_", scale, "_genvdist_copy.pdf"), height = 4, width = 4, units = "cm", dpi = 400)


# Dispersal kernal ==================================================
# Density estimate
if (pvalue <= 0.5) {
  taxa.depth.pop$Site <- as.factor(taxa.depth.pop$Site)
  summary(taxa.depth.pop$Site)
  total_samples = sum(summary(taxa.depth.pop$Site))
  # genotyping rate = 2/3
  N <- total_samples / (2/3)
  # Annotation rate
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
  N <- N / annotation_rate
  if (any(taxa.depth.pop$Site == "SQ20") | any(taxa.depth.pop$Site == "SQ10") | any(taxa.depth.pop$Depth == "5")) {
    exhuastive_sample <- sum(taxa.depth.pop$Site == "SQ20", taxa.depth.pop$Site == "SQ10", taxa.depth.pop$Depth == "5")
    N <- (N - exhuastive_sample) / 0.7 # discovery rate.. 70%
    N <- N + exhuastive_sample
  }
  # per site
  plots <- length(levels(taxa.depth.pop$Site))
  N_plot <- N/plots
  # D individuals within 25 x 2 = 50m2
  D <- N_plot / 50 # individuals per m2
  Neighbourhood <- 1 / b
  Neighbourhood.low <- 1 / b.low
  Neighbourhood.high <- 1 / b.high
  
  sigma2 <- (1 / (4*D*pi*b))
  sigma <- sqrt((1 / (4*D*pi*b)))
  sigma.low <-  sqrt((1 / (4*D*pi*b.low)))
  sigma.high <-  sqrt((1 / (4*D*pi*b.high)))
  
  
  distance <- seq(0, 50, 5)
  p.d <- (1/sigma) * exp(-distance/sigma)
  p.d.low <- (1/sigma.low) * exp(-distance/sigma.low)
  p.d.high <- (1/sigma.high) * exp(-distance/sigma.high)
  plot(p.d ~ distance)
  kernel <- as.data.frame(cbind(distance,p.d,p.d.low,p.d.high))
  
  Lab <- expression(P(d) == paste(frac(1, sigma),
                                  " ", e^{frac(-d, sigma)}))
  
  
  res <- data.frame(cond1 = "regression",
                    x = distance,
                    y = p.d,
                    ymin = p.d.low,
                    ymax = p.d.high)
  
  rib <- geom_ribbon(data = res, aes(x = x, y = y, ymin = ymin, ymax = ymax,
                                     fill = cond1), fill = 'darkgrey', alpha = 0.8)
  
  
  kernelplot <- ggplot(kernel, aes(distance, p.d)) + geom_line() + geom_line(aes(distance, p.d.low), linetype = "dashed", colour = "blue") + 
    geom_line(aes(distance, p.d.high), linetype = "dashed", colour = "red") + 
    theme_bw() +
    annotate("text", x = 25, y = max(kernel$p.d.high) - 0.08,
             label = Lab, parse = T, size = 2) +
    annotate("text", x = 18, y = max(kernel$p.d.high) + 0.01,
             label = "Upper 95% CI", colour = "red", size = 2) +
    annotate("text", x = 30, y = mean(kernel$p.d.high),
             label = "Lower 95% CI", colour = "blue", size = 2) +
    ggtitle(paste0("Neighbourhood = ", round(Neighbourhood, 0))) +
    annotate("text", x = 30, y = min(kernel$p.d) + sd(kernel$p.d*1.5),
             label = expression(paste(sigma, " = ")), size = 2, parse = T) +
    annotate("text", x = 36, y = min(kernel$p.d) + sd(kernel$p.d*1.5),
             label = paste0(round(sigma, 1)), size = 2) +
    ylab("P(d)") +
    xlab("d (m)") + rib + theme(plot.title = element_text(size = 6),
                                axis.title = element_text(size = 8),
                                axis.text = element_text(size = 6)) 
  kernelplot
  ggsave(paste0("../results/ibd/plots/dispersal_kernel_", taxa, "_", category, "_", scale, ".pdf"), units = "cm", height = 4, width = 4, dpi = 400)
  
  dispersal_data <- cbind(taxa, Neighbourhood, Neighbourhood.low, Neighbourhood.high,
                          D, sigma2, sigma)
  write.csv(dispersal_data, paste0("../results/ibd/", taxa, "_", category, "_", scale, "_dispersal_results_new.csv"), quote = FALSE,
            row.names = FALSE)
}


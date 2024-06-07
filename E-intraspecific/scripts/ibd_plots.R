# Title: Isolation-by-distance plots
# Author: Katharine Prata

# Required packages
# R v4.2.0
library(ggplot2) # v3.4.0
library(spaa) # v0.2.2

# Arguments  ==================================================
args = commandArgs(TRUE)
taxa = args[1] #"AA2"
category = args[2] #"WP05"
scale = args[3] # "within"
dimension = args[4] # "2D

# Import data ==================================================
a <- read.delim(paste0("../results/ibd/", taxa, "/", taxa, "_", "a.txt"), header = FALSE, sep = "\t")
dist <- read.delim(paste0("../results/ibd/", taxa, "/", taxa, "_", dimension,".txt"), header = FALSE, sep = "\t")
results <- read.csv(paste0("../results/ibd/", taxa, "/",taxa, "_", category,
                           ".", scale, ".", dimension, ".results_short.txt"))
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
  results.round[,i] <- signif(results[,i], 3)
}

row1 <- rep(NA, nrow(a))
a <- rbind(row1, a)
dist <- rbind(row1, dist)

b = results$slope
b.low = results$s.lowCI
b.high = results$s.highCI
regression <- c(results$intercept, b)
confidence_interval_lower <- c(results$i.lowCI, b.low)
confidence_interval_upper <- c(results$i.highCI, b.high)
pvalue = results$p.slope
sign = "="

# Re-organise data ==================================================
x <- dist[lower.tri(dist)]
y <- a[lower.tri(a)]
dat <- as.data.frame(cbind(x, y))
colnames(dat) <- c("distance", "a")

## setting axes & scales ####
if (scale == "within") {
  dat <- dat[dat$distance < 8,]
  xlim <- c(-4, 4.01)
} else if (scale == "all") {
  xlim <- c(-4, 11)
} else if (scale == "between") {
  xlim <- c(10000, 43000)
}

if (taxa == "AH2" | taxa == "AH1") {
  ylim <- c(0, 4)
} else if (taxa == "AH3") {
  ylim <- c(0, 10)
} else {
    ylim <- c(-0.25, 0.50)
  }

regression <- as.data.frame(cbind(regression[1], regression[2]))
confidence_interval_lower <- as.data.frame(cbind(confidence_interval_lower[1], confidence_interval_lower[2]))
confidence_interval_upper <- as.data.frame(cbind(confidence_interval_upper[1], confidence_interval_upper[2]))

y_lower <- dat$distance * confidence_interval_lower$V2 + confidence_interval_lower$V1
y_mean <- dat$distance * regression$V2 + regression$V1
y_upper <- dat$distance * confidence_interval_upper$V2 + confidence_interval_upper$V1

# Make plots  ========================
slope <- signif(regression[2], 3)
res <- data.frame(cond1 = "regression",
                  x = dat$distance,
                  y = y_mean,
                  ymin = y_lower,
                  ymax = y_upper)

rib <- geom_ribbon(data = res, aes(x = x, y = y, ymin = ymin, ymax = ymax,
                                   fill = cond1), fill = 'blue', alpha = 0.2)

scaleFUN <- function(x) {sprintf("%.2f", x)}
if (dimension == "1D") {
  lab <- "distance (km)"
} else if (dimension == "2D") {
  lab <- "log distance (m)"
}

p <- ggplot(dat, aes(distance, a)) +
  geom_point(shape = 21) +
  ylab("â") +
  xlab(lab) +
  geom_abline(data = regression, aes(intercept = V1, slope = V2), colour = 'red') +
  geom_line(data = res, aes(x, ymin), colour = 'blue', linetype = "dashed") + 
  geom_line(data = res, aes(x, ymax), colour = 'blue', linetype = "dashed") + rib + theme_bw() + 
  ggtitle(paste0("b = ", slope, " p ", sign, " ", signif(pvalue, 3))) +
  scale_y_continuous(labels = scaleFUN, limits = ylim) +
  xlim(xlim) +
  theme(plot.title = element_text(size = 10),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        plot.background = element_blank()) 
p

ggsave(paste0("../results/ibd/plots/", taxa, "_", category, "_", scale, "_", dimension, "_genvdist.pdf"),
       height = 4, width = 4, units = "cm", dpi = 400)

Neighbourhood <- 1 / b
Neighbourhood.low <- 1 / b.low
Neighbourhood.high <- 1 / b.high

dispersal_data <- data.frame(taxa = taxa, 
                             min.dist = min(dat$distance[dat$distance > xlim[1]]),
                             max.dist = max(dat$distance[dat$distance < xlim[2]]),
                             Nb = signif(Neighbourhood, 2), 
                             Nb.low = signif(Neighbourhood.low, 2), 
                             Nb.high = signif(Neighbourhood.high, 2))
write.csv(dispersal_data,
          paste0("../results/ibd/", taxa, "_", category, "_", scale, "_dispersal_results.csv"),
          quote = FALSE, row.names = FALSE)
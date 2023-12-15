# Title: Get results from spagedi output
# Author: Katharine Prata
# Description: Using spagedi import and data handling functions based on https://github.com/lukembrowne/rSpagedi

# Packages
library(ggplot2)
library(tidyr)
library(stringr)

# Functions
readSpagediTable <- function(path_to_out, type){
  out_lines <- readLines(path_to_out, warn = FALSE)
  
  if (type == "perm") {
    start <- grep("LOCATIONS, INDIVIDUALS and/or GENES PERMUTATION TESTS",
                  out_lines) + 3}
  
  if (type == "kin") {
    start <- grep("Genetic analyses at INDIVIDUAL level",
                  out_lines) + 10}
  
  if (type == "dist") {
    start <- grep("Genetic analyses at INDIVIDUAL level",
                  out_lines) + 1}
  
  if (type == "diversity") {
    start <- grep("GENE DIVERSITY and ALLELE FREQUENCIES",
                  out_lines) + 1}
  
  if (type == "spatialdist") {
    start <- grep("PAIRWISE SPATIAL AND GENETIC DISTANCES written as matrices",
                  out_lines) + 3}
  
  if (type == "kindist") {
    all <- grep("Pairwise KINSHIP coefficients ",
                  out_lines)
    start <- all[length(all) - 1] + 2
  }
  
  if (type == "sigma_iter") {
    start <- grep("Estimated gene dispersal parameters for an assumed effective pop density",
                  out_lines) + 1}

  if (length(start) == 0) {
    cat("Date of type", type, "not found in Spagedi output...\n")
    return(NULL)
  }
  
  # Find the last line of the table by looking for the next blank line
  if (type == "kin") {
  end <- min(grep("Estimated gene dispersal", out_lines)) - 1
  } else {
  end = min(which(out_lines[start:length(out_lines)] == "")) + start - 2
  }
  
  raw <- out_lines[start:end] 
  
  split <- strsplit(raw, split = "\t")
  
  if (type == "sigma_iter") {
    split <- lapply(split, function(x) x[!x %in% ""])
  }
  
  maxcol <- max(sapply(split, length))
  
  # Make empty dataframe to store info and set row names
  if (type == "diversity") { # Need to cut out duplicate row names
    split <- split[-c(seq(2, length(split), by = 2))]
    tab <- data.frame(matrix(ncol = maxcol - 1, nrow = length(split)))
    row.names(tab) <- sapply(split, "[[", 1)
  } else if (type == "sigma_iter") {
    tab <- data.frame(matrix(ncol = maxcol, nrow = length(raw) - 1))
  } else { 
    tab <- data.frame(matrix(ncol = maxcol - 1, nrow = length(raw)))
    row.names(tab) <- sapply(split, "[[", 1)
  }
  
  # Set column names
  if (type == "perm") {
    labels  <-  2
    names(tab) <- split[[labels]][-1]
  } else if (type == "sigma_iter") {
    names(tab) <- c(split[[1]][c(1,3)], split[[2]][c(1,3)])
  } else {
    labels  <- 1
    names(tab) <- split[[labels]][-1]
  }
  
  # For diversity - Cut out allele frequency data 
  if (type == "diversity") {
    # Find first column that is NA (blank column) and cut everything after
    tab <- tab[, -c(min(which(is.na(names(tab)))):ncol(tab))]
  }
  
  if (type != "sigma_iter") {
  # Pick out list elements to fill in with empty data to reach 12 columns
  fill_in <- which(sapply(split, length) < maxcol)
  
  for (x in fill_in) {
    split[[x]][!(1:length(split[[2]]) %in% 1:length(split[[x]]))] <- NA
  }
  
  # Fill in columns with actual data
  for (j in 1:ncol(tab)) {
    tab[, j] <- sapply(split, "[[", (j + 1))
  }} else {
   tab[1] <- sapply(split, "[[", 2)[1]
   tab[2] <- sapply(split, "[[", 4)[1]
   tab[3] <- sapply(split, "[[", 2)[2]
   tab[4] <- sapply(split, "[[", 4)[2]
  }

  if (type == "perm") {tab  <- tab[-c(1,2,3), ]}
  if (type == "kin") {tab <- tab[-1, ]}
  if (type == "dist") {tab <- tab[, -1]}
  if (type == "diversity") {tab <- tab[-1, -1]}
  if (type == "spatialdist") {tab <- tab[-1, ]}
  if (type == "kindist") {tab <- tab[-1, ]}
  
  # Convert columns to numeric
  if (type == "diversity") {tab[, -c(2,3)] <- apply(tab[, -c(2,3)], 2, as.numeric)}
  if (type != "diversity") {tab[,] <- apply(tab, 2, as.numeric)}
  return(tab)
}
plotAutoCor <- function(spagediList, overlay = FALSE, color = "black", max_dist){
  
  symbols <- 19
  
  perm <- spagediList$perm
  dist <- spagediList$dist
  dist <- as.numeric(spagediList$dist["Mean ln(distance)", ])
  
  if (!is.null(perm)) {
    slope <- perm[, (ncol(perm) - 1):ncol(perm)]
    perm <- perm[, c(2:(ncol(perm) - 4))] # chop off last 3 columns and 1st col
    
    obs_slope <- as.numeric(slope["Obs val", ])
    conf_hi_slope <- as.numeric(slope["95%CI-sup", ])
    conf_low_slope <- as.numeric(slope["95%CI-inf", ])
    
    obs <- as.numeric(perm["Obs val", ])
    conf_hi <- as.numeric(perm["95%CI-sup", ])
    conf_low <- as.numeric(perm["95%CI-inf", ])
    
    # Closed symbol if permutation was significant
    sig_slope <- apply(slope[c(8,9,10), ], 2, FUN = function(x) {any(x < 0.05)})
    sig_slope
    sig <- apply(perm[c(8,9,10), ], 2, FUN = function(x) {any(x < 0.05)})
    sig
    symbols <- rep(1, ncol(perm))
    symbols[sig] <- 19
    
  } else {
    kin <- spagediList$kin
    kin <- kin[, c(2:(ncol(kin) - 4))]
    obs <- kin["ALL LOCI", ]
    conf_hi = 0
    conf_low = 0
  }
  
  if (overlay) {
    par(new = TRUE)
    points(dist, obs, type = "b", pch = symbols, xlab = "", ylab = "",
           yaxt = "n", xaxt = "n", col = color, lwd = 2, lty = 1)
  if (!is.null(perm)) {
      lines(dist, conf_hi, lty = 4, col = "red")
      lines(dist, conf_low, lty = 4, col = "red")
    }
  } else {
    
    dat <- data.frame(distance = dist, kinship = obs, upper = conf_hi, lower = conf_low,
                      significant = symbols)
    dat$significant <- as.factor(dat$significant)
    
    p <- ggplot(dat, aes(distance, kinship)) +
    geom_hline(yintercept = 0, colour = "black", size = 0.3) +
    geom_line(size = 1, colour = "black", linetype = "longdash") + 
    geom_point(aes(fill = significant), size = 4, shape = 21) + 
    scale_fill_manual(values = c("white", "black")) +
    ylim(c(min(obs, conf_low, na.rm = TRUE) * 1.1, 
           max(obs, conf_hi, na.rm = TRUE) * 1.1)) +
    xlim(c(min(dat$distance), log(max_dist))) + 
    ylab("Kinship") +
    xlab("Log distance (m)") +
    geom_line(aes(dist, upper), colour = "red", linetype = "dotdash", size = 0.5) +
    geom_line(aes(dist, lower), colour = "red", linetype = "dotdash", size = 0.5) + 
    theme_bw() +
    theme(legend.position = "none")
  }
  return(p)
}
calcSp <- function(kin, row_name){
  return( ( -kin[row_name, (length(kin) - 7)]) / (1 - kin[row_name, 2]) )
}
SpSummary <- function(spagediList){
  kin <- spagediList$kin
  mean <- calcSp(kin, "ALL LOCI")
  mean_jck <- calcSp(kin, "Mean")
  se_jck <- calcSp(kin, "SE")
  
  cat("\n\n")
  cat("Mean Sp across all loci       == ", mean, "\n")
  cat("-------\n")
  
  if (!is.na(mean_jck)) {
    cat("Mean Sp across jacknifed loci == ", mean_jck, "\n")
    cat("-------\n")
    cat("SE Sp across jacknifed loci == ", -se_jck, "\n")
    cat("-------\n")
  }
  # Select rows that have individual loci
  if (!is.na(mean_jck)) {
    loci_rows <- (grep("ALL LOCI", row.names(kin)) + 1):(grep("Jack", row.names(kin)) - 1)
  } else{
    loci_rows <- 2:nrow(kin)
  }
  sp_loci <- calcSp(kin, loci_rows )
  names(sp_loci) <- row.names(kin)[loci_rows]
  # ignore all loci estimates...
  
  sp_stats <- data.frame(mean_sp = mean,
                         mean_jack_sp = mean_jck,
                         se_jack_sp = se_jck)
  return(sp_stats)
}
calc_sigma <- function(nb, density, dimension) {
  if (dimension == "1d") {
    sigma <- sqrt((nb/(4*density)))
  } else if (dimension == "2d") {
    sigma <- sqrt((nb/(4*pi*density)))
  }
  return(sigma)
}

# spagedi data
# Make outfile an argument
args <- commandArgs(TRUE)
taxa <- args[1]
rep <- as.numeric(args[2])

if (taxa == "aa1") {
  TAXA <- "AA1"
} else if (taxa == "aa2") {
  TAXA <- "AA2"
} else if (taxa == "ah1") {
  TAXA <- "AH1"
} else if (taxa == "ah2") {
  TAXA <- "AH2"
} else if (taxa == "ah3") {
  TAXA <- "AH3"
} else if (taxa == "al1") {
  TAXA <- "AL1"
} else if (taxa == "al2") {
  TAXA <- "AL2"
}

# Change which outfile you use for the different taxa!
outfile = paste0("../results/ibd/", TAXA, "/", taxa, ".spagedi-results", rep, ".txt")
name <- str_split(str_split(outfile, "/")[[1]][length(str_split(outfile, "/")[[1]])], "\\.")[[1]][1]

spagediList <- list()
spagediList$perm <- readSpagediTable(outfile, "perm")
spagediList$dist <- readSpagediTable(outfile, "dist")
spagediList$kin <- readSpagediTable(outfile, "kin")
spagediList$spatialdist <- readSpagediTable(outfile, "spatialdist")
spagediList$kindist <- readSpagediTable(outfile, "kindist")
spagediList$sigma_iter <- readSpagediTable(outfile, "sigma_iter")

sum <- SpSummary(spagediList)

nb <- 1/sum$mean_jack_sp
nb_upp <- 1/(sum$mean_jack_sp + sum$se_jack_sp)
nb_low <- 1/(sum$mean_jack_sp - sum$se_jack_sp)

print(paste("Neighbourhood size:", round(nb,0) , "upp:", round(nb_upp, 0), "low:", round(nb_low, 0)))

if (name == "aa1") {
  dc <- 0.24
  de <- 0.27
} else if (name == "aa2") {
  dc <- 1.14
  de <- 3.28
} else if (name == "ah1") {
  dc <- 0.36
  de <- 0.45
} else if (name == "ah2") {
  dc <- 0.27
  de <- 0.79
} else if (name == "ah3") {
  dc <- 0.16
  de <- NA
} else if (name == "al1") {
  dc <- 0.25
  de <- NA
} else if (name == "al2") {
  dc <- 0.37
  de <- NA
}

sigma_dc <- calc_sigma(nb, dc, "2d")
print(paste("Census simga is:", round(sigma_dc, 2), "using census of", dc))

sigma_de <- calc_sigma(nb, de, "2d")
print(paste("Effective simga is:", round(sigma_de, 2), "using effective of", de))

print("Sigma iterations using spagedi")
spagediList$sigma_iter
# doesn't work if it cycles between - TODO

# reorg
spatialdist <- spagediList$spatialdist
loiselle <- spagediList$kindist
x <- spatialdist[lower.tri(spatialdist)]
y <- loiselle[lower.tri(loiselle)]
dat <- as.data.frame(cbind(x, y))
colnames(dat) <- c("distance", "F")

if (rep == 2 & !is.na(sigma_de)) {
  dat <- dat[dat$distance < sigma_de*20,] # restrict by 20sigma
  dat <- dat[dat$distance > sigma_de,] # restrict by sigma
} else if (rep == 2 & is.na(sigma_de)) {
  dat <- dat[dat$distance < sigma_dc*20,] # restrict by 20sigma
  dat <- dat[dat$distance > sigma_dc,] # restrict by sigma
}


results <- data.frame(Taxon = name,
                      Mindist = min(dat$distance),
                      Maxdist = max(dat$distance),
                      Nb = nb,
                      Nb_upp = nb_upp,
                      Nb_low = nb_low,
                      Sigma_dc = sigma_dc,
                      Sigma_de = sigma_de)

write.table(results,
            paste0("../results/ibd/", taxa, ".spagedi-summary", rep, ".txt"),
            quote = FALSE, row.names = FALSE)


## setting axes & scales ####
if (rep == 2) { 
  #ggplot(dat, aes(log_distance, r)) + geom_point() + theme_bw()
  #ggplot(dat, aes(distance, r)) + geom_point() + theme_bw()
  
  dat$log_distance <- log(dat[,1])
  
  # Setting up regression details
  b <- spagediList$perm[3,length(spagediList$perm)]
  b.low <- spagediList$perm[3,length(spagediList$perm)] + spagediList$perm[6,length(spagediList$perm)]
  b.high <- spagediList$perm[3,length(spagediList$perm)] + spagediList$perm[7,length(spagediList$perm)]
  regression <- data.frame(V1 = spagediList$kin[1,(length(spagediList$kin) - 6)], V2 = b)
  confidence_interval_lower <- c(spagediList$kin[1,(length(spagediList$kin) - 6)], b.low)
  confidence_interval_upper <- c(spagediList$kin[1,(length(spagediList$kin) - 6)], b.high)
  pvalue = spagediList$perm[10,length(spagediList$perm)]
  sign = "="
  y_lower <- dat$log_distance * confidence_interval_lower[2] + confidence_interval_lower[1]
  y_mean <- dat$log_distance * regression$V2 + regression$V1
  y_upper <- dat$log_distance * confidence_interval_upper[2] + confidence_interval_upper[1]
  
  
  # Make plot
  slope <- round(regression[2], 3)
  res <- data.frame(cond1 = "regression",
                    x = dat$log_distance,
                    y = y_mean,
                    ymin = y_lower,
                    ymax = y_upper)
  
  rib <- geom_ribbon(data = res, aes(x = x, y = y, ymin = ymin, ymax = ymax,
                                     fill = cond1), fill = 'blue', alpha = 0.2)
  xlim <- c(0, 4.01)
  ylim <- c(-0.25, 0.4)
  
  p <- ggplot(dat, aes(log_distance, F)) +
    geom_point(shape = 21) +
    ylab("kinship (F)") +
    xlab("log distance") +
    geom_abline(data = regression, aes(intercept = V1, slope = V2), colour = 'red') +
    geom_line(data = res, aes(x, ymin), colour = 'blue', linetype = "dashed") + 
    geom_line(data = res, aes(x, ymax), colour = 'blue', linetype = "dashed") + rib + theme_bw() + 
    ggtitle(paste0("Slope = ", slope, " p ", sign, " ", signif(pvalue, 3))) +
    ylim(ylim) +
    xlim(xlim) +
    theme(plot.title = element_text(size = 10),
          axis.title = element_text(size = 8),
          axis.text = element_text(size = 8)) 
  p
  
  ggsave(paste0("../results/ibd/plots/", TAXA, "_loisellevdist.pdf"),
         height = 4, width = 4, units = "cm", dpi = 400) }

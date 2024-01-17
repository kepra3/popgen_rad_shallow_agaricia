# Title: Density calculations
# Author: Katharine Prata
# Description: Calculating effective density for each taxa based on area covered by individuals used for analysis
# i.e., per location, thus including area of plots at different depths within locations and area between plots

# Functions
calc_area <- function(plots) {
  area <- 50 * length(plots)
  # Average across sites, known cattetag distances for WP and SB averaged these for CA and SQ,
  # then averaged across all three estimates
  dist_10_20 <- sum(17.155, 15.44, 18.87)/3
  # Rough average distance measure as seen from google maps
  dist_5_10 <- 30
  
  if (("5" %in% plots) & ("10" %in% plots)) {
    area <- area + (dist_5_10 * 25)
  } else {
    print("No 5 and 10")
  }
  if (("10" %in% plots) & ("20" %in% plots)) {
    area <- area + (dist_10_20 * 25)
  } else {
    print("No 10 and 20")
  }
  return(area)
}

# aa1
aa1.ne <- 56
aa1.plots <- c("10", "20")
aa1.area <- calc_area(aa1.plots)
aa1.density <- aa1.ne/aa1.area
print(paste("aa1 density =", round(aa1.density, 2)))

# aa2
aa2.ne1 <- sum(707, 341, 483)/3
aa2.ne2 <- 193
aa2.plots1 <- c("5", "10", "20")
aa2.plots2 <- c("10", "20")
aa2.area1 <- calc_area(aa2.plots1)
aa2.area2 <- calc_area(aa2.plots2)
aa2.density1 <- aa2.ne1/aa2.area1
aa2.density2 <- aa2.ne2/aa2.area2
aa2.density <- sum(aa2.density1, aa2.density2)/2
print(paste("aa2 density =", round(aa2.density, 2)))

# 1d harmonic mean
aa2.pop1.Ne = 349.9 #(prcit = 0.01)
aa2.pop2.Ne = 243.9 #(pcrit = 0.02)
aa2.pop3.Ne = 287.8 #(pcrit = 0.02)
aa2.pop4.Ne = 25.6 #(pcrit= 0.03)
n = 4
harm_mean <- n / (sum(1/aa2.pop1.Ne, 1/aa2.pop2.Ne,
                  1/aa2.pop3.Ne, 1/aa2.pop4.Ne))
harm_mean

# ah1
ah1.ne <- 316
ah1.plots <- c("5", "10", "20")
ah1.area <- calc_area(ah1.plots)
ah1.density <- ah1.ne/ah1.area
print(paste("ah1 density =", round(ah1.density, 2)))

# ah2
ah2.ne <- sum(291, 338)/2
ah2.plots <- c("5", "10")
ah2.area <- calc_area(ah2.plots)
ah2.density <- ah2.ne/ah2.area
print(paste("ah2 density =", round(ah2.density, 2)))

# ah3
ah3.ne <- 344
ah3.plots <- c("5", "10")
ah3.area <- calc_area(ah3.plots)
ah3.density <- ah3.ne/ah3.area
print(paste("ah3 density =", round(ah3.density, 2)))

# Title: Convert vcf to spagedi format
# Author: Katharine Prata
# Last edit: 14/12/23

# Packages ####
library(dplyr)
library(tidyr)
library(vcfR) # 1.12.0
library(adegenet) # 2.1.7

# Functions
convert2snpdf <- function(genind) {
  gendata <- genind2df(genind)
  samples <- rownames(gendata)
  gendata <- cbind(samples, gendata)
  return(gendata)
}
addxy <- function(gendata, geodata) {
  geodata <- geodata %>% separate(Individual, into = c("Sample", "Spp", "Pop", "X"), sep = "_") %>% 
    unite(Individual, Sample, Spp, Pop, sep = "_") %>% 
    select(Individual, x, y, z)
  gendata$x <- geodata$x[match(gendata[,1], geodata[,1])]
  gendata$y <- geodata$y[match(gendata[,1], geodata[,1])]
  gendata <- gendata[!is.na(gendata$x),]
  return(gendata)
}
gen2spagedi <- function(gendata, distance_classes) {
  colnames(gendata)[colnames(gendata) == colnames(gendata)[1]] <- "sample"
  # Check format of data
  #gendata[,2:(length(gendata) - 2)] <- lapply(gendata[,2:(length(gendata) - 2)], factor)
  #list <- c(0)
  #for (column in 2:(length(gendata) - 2)) {
  #  list <- append(list, length(levels(gendata[,column])))
  #}
  mat <- as.matrix(gendata[,2:(length(gendata) - 2)])
  mat[is.na(mat)] = FALSE
  if (any(mat == 44) | any(mat == 33)) {
    print("Alleles in nucleotide bases, 1,2,3,4 with 0 as missing data")
    # convert to nucleotides
    mat <- apply(mat, FUN = gsub, MARGIN = 2, pattern = "1", replacement = "A")
    mat <- apply(mat, FUN = gsub, MARGIN = 2, pattern = "2", replacement = "T")
    mat <- apply(mat, FUN = gsub, MARGIN = 2, pattern = "3", replacement = "G")
    mat <- apply(mat, FUN = gsub, MARGIN = 2, pattern = "4", replacement = "C")
    # convert to genepop SNPs (e.g., allele1 or allele2)
    mat <- apply(mat, FUN = gsub, MARGIN = 2, pattern = "0", replacement = "00")
    mat <- apply(mat, FUN = gsub, MARGIN = 2, pattern = "A", replacement = "01")
    mat <- apply(mat, FUN = gsub, MARGIN = 2, pattern = "T", replacement = "02")
    mat <- apply(mat, FUN = gsub, MARGIN = 2, pattern = "G", replacement = "01")
    mat <- apply(mat, FUN = gsub, MARGIN = 2, pattern = "C", replacement = "02")
  } else if (any(mat == 01) | any(mat == 10)) {
    print("Alleles in 0 and 1, with NA is missing data")
    mat <- apply(mat, FUN = gsub, MARGIN = 2, pattern = "0", replacement = "G")
    mat <- apply(mat, FUN = gsub, MARGIN = 2, pattern = "1", replacement = "C")
    mat <- apply(mat, FUN = gsub, MARGIN = 2, pattern = "G", replacement = "01")
    mat <- apply(mat, FUN = gsub, MARGIN = 2, pattern = "C", replacement = "02")
    mat[mat == FALSE] = "0000"
  } else if (any(nchar(mat) == 6)) {
    print("Microsatt alleles in number of repeats, 00 is missing data")
    mat <- apply(mat, FUN = gsub, MARGIN = 2, pattern = "-9", replacement = "0000")
  } else if (any(mat != 44) & any(mat != 33) & any(mat != 01) & any(mat != 10) & !any(nchar(mat) == 6)) {
    print("Alleles in 1 and 2, with 0 and -9 missing data")
    mat <- apply(mat, FUN = gsub, MARGIN = 2, pattern = "1", replacement = "G")
    mat <- apply(mat, FUN = gsub, MARGIN = 2, pattern = "2", replacement = "C")
    mat <- apply(mat, FUN = gsub, MARGIN = 2, pattern = "0", replacement = "00")
    mat <- apply(mat, FUN = gsub, MARGIN = 2, pattern = "G", replacement = "01")
    mat <- apply(mat, FUN = gsub, MARGIN = 2, pattern = "C", replacement = "02")
    mat <- apply(mat, FUN = gsub, MARGIN = 2, pattern = "-9", replacement = "00")
  } else {
    print("Different format")
  }
  # genotype matrix (very similar to genepop)
  gen_lines <- cbind(gendata$sample, gendata$x, gendata$y, mat)
  
  # header
  header <- c(paste("// Spagedi file format", format(Sys.time(), "%Y%m%d@%H%M")), rep("", ncol(gen_lines) - 1))
  # first line, 6 format numbers separated by a tab
  #1. # of indiviudals, 2. # of categories(?), 3. # of spatial coordinates in a cartesian coordinate system (0 to 3),
  # or put -2 for latitude & longitude, 4. # of loci, 5. # of digits used to code one allele (1 to 3), 6. ploidy
  first_line <- c(nrow(gendata), 0, 2, ncol(gendata) - 3, 2, 2, rep("", ncol(gen_lines) - 6))
  # second line, number of distance intervals, preceded by a negative sign to choose ideal number based on samples
  second_line <- c(length(distance_classes), distance_classes, rep("", ncol(gen_lines) - (length(distance_classes) + 1)))
  # the names used as column labels
  # a generic name for individuals(?), spatial coordinates, name of each locus
  third_line <- c("Ind", "X", "Y", colnames(gendata)[2:(length(gendata) - 2)])
  end_line <- c("END", rep("", ncol(gen_lines) - 1))
  spagedi <- rbind(header, first_line, second_line, third_line, gen_lines, end_line)
  return(spagedi)
}

# Arguments
args = commandArgs(TRUE)
taxa <- args[1]

# Load data
gendata <- vcfR2genind(read.vcfR(paste0("../data/", taxa, "_1div_nc_20.vcf")))
geodata <- read.csv("../data/all_annotations_X_HORIZ_parallel_XYZ_adjusted.txt")

# Convert and combine
gendata <- convert2snpdf(gendata)
gendata <- addxy(gendata, geodata)

# Distances
all.distances <- dist(c(gendata$x, gendata$y))
within.distances <- all.distances[all.distances < 200]
max(within.distances) 
hist(within.distances)

if (taxa == "aa1") {
  distances <- seq(5, 75, 5)
} else if (taxa == "aa2") {
  distances <- seq(5, 75, 5)
} else if (taxa == "ah1" | any(grepl("^ah1-c", taxa))) {
  distances <- seq(10, 70, 10)
  distances <- c(5, distances)
} else if (taxa == "ah2") {
  distances <- seq(10, 60, 10)
} else if (taxa == "ah3" | any(grepl("^ah3-c", taxa))) {
  distances <- seq(10, 50, 10)
} else if (taxa == "al1") {
  distances <- seq(5, 75, 5)
} else if (taxa == "al2") {
  distances <- seq(5, 75, 5)
} else {
  distances <- seq(5, 100, 5)
}

spag <- gen2spagedi(gendata, distances)
write.table(spag, file = paste0("../data/", taxa, ".spagedi.txt"),
              quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

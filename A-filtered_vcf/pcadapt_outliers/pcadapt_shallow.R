# Title: PCAdapt
# Author: Katharine Prata
# Date created: 7/3/22
# Last edit: 19/4/22

# Packages
# R version 4.2.0
library("pcadapt") # version 4.3.3, https://doi.org/10.1093/molbev/msaa053
library("tidyr") # version 1.2.0, https://CRAN.R-project.org/package=tidyr
library("qvalue") # version 2.28.0, http://github.com/jdstorey/qvalue
library("vcfR") # version 1.12.0, http://dx.doi.org/10.1111/1755-0998.12549 
library("adegenet") # version 2.1.7, 10.1093/bioinformatics/btn129, 10.1093/bioinformatics/btr521
library("dplyr") # version 1.0.9, https://CRAN.R-project.org/package=dplyr

# functions ####
find.snps <- function(loc.names,snp_pc){
  snps <- as.data.frame(snp_pc[,1])
  snps <- as.numeric(levels(snps[,1]))[snps[,1]]
  pcs <- as.data.frame(snp_pc[,2])
  outlier_names <- as.data.frame(matrix(nrow = length(snps), ncol = 1))
  for (i in seq(1:length(snps))) {
    outlier_names[i,1] <- paste(loc.names[snps[i],1])
    outlier_names[i,2] <- pcs[i,1]
  }
  return(outlier_names)
}

# all-aga ####
VCF <- "all-aga_1d_nc_20.vcf"
VCF.adapt <- read.pcadapt(VCF, type = "vcf")
genlight <- vcfR2genlight(read.vcfR(VCF))
genlight$ind.names

x <- pcadapt(VCF.adapt, K = 100)
plot(x, option = "screeplot")
plot(x, option = "screeplot", K = 10)
pop <- read.csv("all-aga_1d_nc_20_6.csv")
pop.ordered <- pop[order(pop[,1], decreasing = FALSE),]
poplist <- pop.ordered$Clusters
poplist[poplist == 1] <- "AA1"
poplist[poplist == 2] <- "AL2"
poplist[poplist == 3] <- "AH1"
poplist[poplist == 4] <- "AA2"
poplist[poplist == 5] <- "AL1"
poplist[poplist == 6] <- "AH2"
plot(x, option = "scores", pop = poplist)
plot(x, option = "scores", i = 1, j = 3, pop = poplist)
plot(x, option = "scores", i = 2 , j = 3, pop = poplist)
plot(x, option = "scores", i = 3, j = 4, pop = poplist)
plot(x, option = "scores", i = 4, j = 5, pop = poplist)
plot(x, option = "scores", i = 5, j = 6, pop = poplist)
x <- pcadapt(VCF.adapt, K = 6)
summary(x) # 8864 passed MAF out of 15536

plot(x, option = "manhattan")
plot(x, option = "qqplot")

# histograms of the test statistic and of the p-values
hist(x$pvalues, xlab = "p-values", main = NULL, breaks = 120, col = "pink")

# presence of outlier visible when plotting a histogram of the test statistic
plot(x, option = "stat.distribution")

# choosing a cutoff for outlier detection

# q-values
qval <- qvalue(x$pvalues)$qvalues
alpha <- 0.1
outliers1 <- which(qval < alpha) # expected false discovery rate lower than 10%
length(outliers1) # 336 outliers!
2931/15536 # 18.9%

# Benjamin-Hochberg Procedure
padj <- p.adjust(x$pvalues, method = "BH")
alpha <- 0.1
outliers2 <- which(padj < alpha)
length(outliers2) # same again

# Bonferroni correction
padj <- p.adjust(x$pvalues, method = "bonferroni")
outliers3 <- which(padj < alpha) 
length(outliers3) # 2139
2139/15536 # 13.8%

# Association between PCs and outliers
snp_pc <- get.pc(x, outliers1)
snp_pc <- as.data.frame(snp_pc)
snp_pc$PC <- as.factor(snp_pc$PC)
snp_pc$SNP <- as.factor(snp_pc$SNP)
str(snp_pc)
summary(snp_pc)
#1. 373 association with separation b/w AA/AH and AL
#2. 2306 association with separation b/w/ AA and AH
#3. 26 separate AL1 and AL2
#4. 73 separate AH1 and AH2
#5. 122 separate AA1 and AA2
#6. 31 separate AA1 geographically

loc.names <- as.data.frame(genlight@loc.names)

outlier_snp_list <- find.snps(loc.names,snp_pc)

first_spp_outliers <- subset(outlier_snp_list, V2 == 1)
length(first_spp_outliers$V1) # 373

second_spp_outliers <- subset(outlier_snp_list, V2 == 2)
length(second_spp_outliers$V1) # 2306

lam_outliers <- subset(outlier_snp_list, V2 == 3)
length(lam_outliers$V1) # 26

hum_outliers <- subset(outlier_snp_list, V2 == 4)
length(hum_outliers$V1) # 73

aga_outliers <- subset(outlier_snp_list, V2 == 5)
length(aga_outliers$V1) # 122

geo_outliers <- subset(outlier_snp_list, V2 == 6)
length(geo_outliers$V1) # 31

# exporting snp list to filtering folder
all_outlier_list <- as.data.frame(outlier_snp_list$V1)
colnames(all_outlier_list) <- "locus_pos"
all_outlier_list <- all_outlier_list %>% 
  separate(locus_pos, sep = "_", into = c("loc", "position"))
write.table(all_outlier_list,
            file = "all-aga_pcadapt_outliers.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# agaricites ####
VCF.ac <- "ac_1d-m_nc_20.vcf"
ac.adapt <- read.pcadapt(VCF.ac, type = "vcf")
ac.genlight <- vcfR2genlight(read.vcfR(VCF.ac))

x <- pcadapt(input = ac.adapt, K = 10)
# the scree plot displays in decreasing order the percentage of variance by each PC.
# It corresponds to the eigenvalues in decreading order. The ideal pattern
# is a steep curve followed by a bend and a straight line.
# The eigenvalues that correspone to random variation lie on a straight line
# whereas the ones that correspond to population structure lie on a steep curve
# Keep the PCs that correspond to eigenvalues to the left of the straight line (Cattell's rule)
plot(x, option = "screeplot")
# 4 or 5?
# add a popfile
ac.pop <- read.csv("ac_1d_nc_20_4.csv", header = TRUE)
ac.pop.ordered <- ac.pop[order(ac.pop[,1], decreasing = FALSE),]
ac.poplist <- ac.pop.ordered$Clusters

plot(x, option = "scores", pop = ac.poplist) # i = 1, j = 2
plot(x, option = "scores", i = 3, j = 4, pop = ac.poplist) 
plot(x, option = "scores", i = 5, j = 6, pop = ac.poplist) # nothing

# computing the test statistic based on PCA
x <- pcadapt(ac.adapt, K = 4)
# the user can also set the paramater min.maf that corresponds to a threshold
# of minor allele frequency. By default it is set to 5%.
summary(x)
# 3708 loci/snps were used because others did not pass the MAF=0.05

# manhattan plot of outliers
plot(x, option = "manhattan")

# expected uniform distribution of the p-values
plot(x, option = "qqplot")
# the smallest p-values are smaller than expected confirming the presence
# of outliers

# histograms of the test statistic and of the p-values
hist(x$pvalues, xlab = "p-values", main = NULL, breaks = 120, col = "pink")

# presence of outlier visible when plotting a histogram of the test statistic
plot(x, option = "stat.distribution")

# choosing a cutoff for outlier detection

# q-values
qval <- qvalue(x$pvalues)$qvalues
alpha <- 0.05
outliers1 <- which(qval < alpha) # expected false discovery rate lower than 10%
length(outliers1) # 347 outliers!
347/10137
# ~ 3%

# Benjamin-Hochberg Procedure
padj <- p.adjust(x$pvalues, method = "BH")
alpha <- 0.05
outliers2 <- which(padj < alpha)
length(outliers2) # 347

# Bonferroni correction
padj <- p.adjust(x$pvalues, method = "bonferroni")
outliers3 <- which(padj < alpha)
length(outliers3) # 170

# Association between PCs and outliers
snp_pc <- get.pc(x, outliers1)
snp_pc <- as.data.frame(snp_pc)
snp_pc$PC <- as.factor(snp_pc$PC)
snp_pc$SNP <- as.factor(snp_pc$SNP)
str(snp_pc)
summary(snp_pc)

ac.loc.names <- as.data.frame(ac.genlight@loc.names)

ac_outliers <- find.snps(ac.loc.names, snp_pc)
length(ac_outliers$V1) # 347

ac_outliers_list <- as.data.frame(ac_outliers$V1)
colnames(ac_outliers_list) <- "locus_pos"
ac_outliers_list <- ac_outliers_list %>% 
  separate(locus_pos, sep = "_", into = c("loc","position")) %>% 
  separate(loc, sep = 3, into = c("X", "locusnum")) %>% 
  separate(position, sep = 3, into = c("Y", "pos"))
ac_outliers_list$Z <- "RAD"
ac_outliers_list <- ac_outliers_list %>% 
  unite("RADNAME", "Z", "locusnum")

ac_outliers_list <- ac_outliers_list[,c(2, 4)]
  

write.table(ac_outliers_list,
            file = "ac_pcadapt_outliers.txt",
            row.names = FALSE, col.names = FALSE , quote = FALSE)


# humilis ####
VCF.hu <- "hu_1d-m_nc_20.vcf"
hu.adapt <- read.pcadapt(VCF.hu, type = "vcf")
# 3651
hu.genlight <- vcfR2genlight(read.vcfR(VCF.hu))

x <- pcadapt(input = hu.adapt, K = 100)
# the scree plot displays in decreasing order the percentage of variance by each PC.
# It corresponds to the eigenvalues in decreading order. The ideal pattern
# is a steep curve followed by a bend and a straight line.
# The eigenvalues that correspone to random variation lie on a straight line
# whereas the ones that correspond to population structure lie on a steep curve
# Keep the PCs that correspond to eigenvalues to the left of the straight line (Cattell's rule)
plot(x, option = "screeplot")
plot(x, option = "screeplot", K = 10)
# 3 or 4
# add a popfile
hu.pop <- read.csv("hu_1d_nc_20_4.csv", header = TRUE)
hu.pop.ordered <- hu.pop[order(hu.pop[,1], decreasing = FALSE),]
hu.poplist <- hu.pop.ordered$Clusters

plot(x, option = "scores", pop = hu.poplist) # 1. splits clust 3 and 1/4&2 and 2 splits Clust 2
plot(x, option = "scores", i = 3, j = 4, pop = hu.poplist) # splits 4 and 1, 4 doesn't do much a few outlier indivs?
plot(x, option = "scores", i = 5, j = 6, pop = hu.poplist) # nothing

# computing the test statistic based on PCA
x <- pcadapt(hu.adapt, K = 4)
# the user can also set the paramater min.maf that corresponds to a threshold
# of minor allele frequncy. By default it is set to 5%.
summary(x)
# 2370 loci/snps were used because others did not pass the MAF=0.05

# manhattan plot of outliers
plot(x, option = "manhattan")

# expected uniform distribution of the p-values
plot(x, option = "qqplot")
# the smallest p-values are smaller than expected confirming the presence
# of outliers

# histograms of the test statistic and of the p-values
hist(x$pvalues, xlab = "p-values", main = NULL, breaks = 120, col = "pink")

# presence of outlier visible when plotting a histogram of the test statistic
plot(x, option = "stat.distribution")

# choosing a cutoff for outlier detection

# q-values
qval <- qvalue(x$pvalues)$qvalues
alpha <- 0.05
outliers1 <- which(qval < alpha) # expected false discovery rate lower than 10%
length(outliers1) # 362 outliers!
221/3651
# ~ 6%

# Benjamin-Hochberg Procedure
padj <- p.adjust(x$pvalues, method = "BH")
alpha <- 0.05
outliers2 <- which(padj < alpha)
length(outliers2) # 221

# Bonferroni correction
padj <- p.adjust(x$pvalues, method = "bonferroni")
outliers3 <- which(padj < alpha)
length(outliers3) # 98


# Association between PCs and outliers
snp_pc <- get.pc(x, outliers1)
snp_pc <- as.data.frame(snp_pc)
snp_pc$PC <- as.factor(snp_pc$PC)
snp_pc$SNP <- as.factor(snp_pc$SNP)
str(snp_pc)
summary(snp_pc)

hu.loc.names <- as.data.frame(hu.genlight@loc.names)

hu_outliers <- find.snps(hu.loc.names, snp_pc)
length(hu_outliers$V1) # 221

hu_outliers_list <- as.data.frame(hu_outliers$V1)
colnames(hu_outliers_list) <- "locus_pos"
hu_outliers_list <- hu_outliers_list %>% 
  separate(locus_pos, sep = "_", into = c("loc","position")) %>% 
  separate(loc, sep = 3, into = c("X", "locusnum")) %>% 
  separate(position, sep = 3, into = c("Y", "pos"))
hu_outliers_list$Z <- "RAD"
hu_outliers_list <- hu_outliers_list %>% 
  unite("RADNAME", "Z", "locusnum")

hu_outliers_list <- hu_outliers_list[,c(2, 4)]

write.table(hu_outliers_list,
            file = "hu_pcadapt_outliers.txt",
            row.names = FALSE, col.names = FALSE , quote = FALSE)



# lamarcki ####
VCF.lm <- "lm_1d-m_nc_20.vcf"
lm.adapt <- read.pcadapt(VCF.lm, type = "vcf")
lm.genlight <- vcfR2genlight(read.vcfR(VCF.lm))

x <- pcadapt(input = lm.adapt, K = 10)
# the scree plot displays in decreasing order the percentage of variance by elmh PC.
# It corresponds to the eigenvalues in decreading order. The ideal pattern
# is a steep curve followed by a bend and a straight line.
# The eigenvalues that correspone to random variation lie on a straight line
# whereas the ones that correspond to population structure lie on a steep curve
# Keep the PCs that correspond to eigenvalues to the left of the straight line (Cattell's rule)
plot(x, option = "screeplot")
# 4 or 5?
# add a popfile
lm.pop <- read.csv("lm_1d_nc_20_2.csv", header = TRUE)
lm.pop.ordered <- lm.pop[order(lm.pop[,1], decreasing = FALSE),]
lm.poplist <- lm.pop.ordered$Clusters

plot(x, option = "scores", pop = lm.poplist) # i = 1, j = 2
plot(x, option = "scores", i = 3, j = 4, pop = lm.poplist) 
plot(x, option = "scores", i = 5, j = 6, pop = lm.poplist) # nothing

# computing the test statistic based on PCA
x <- pcadapt(lm.adapt, K = 2)
# the user can also set the paramater min.maf that corresponds to a threshold
# of minor allele frequncy. By default it is set to 5%.
summary(x)
# 3380 loci/snps were used because others did not pass the MAF=0.05

# manhattan plot of outliers
plot(x, option = "manhattan")

# expected uniform distribution of the p-values
plot(x, option = "qqplot")
# the smallest p-values are smaller than expected confirming the presence
# of outliers

# histograms of the test statistic and of the p-values
hist(x$pvalues, xlab = "p-values", main = NULL, breaks = 120, col = "pink")

# presence of outlier visible when plotting a histogram of the test statistic
plot(x, option = "stat.distribution")

# choosing a cutoff for outlier detection

# q-values
qval <- qvalue(x$pvalues)$qvalues
alpha <- 0.05
outliers1 <- which(qval < alpha) # expected false discovery rate lower than 10%
length(outliers1) # 189 outliers!
132/10554
# ~ 1%

# Benjamin-Hochberg Procedure
padj <- p.adjust(x$pvalues, method = "BH")
alpha <- 0.05
outliers2 <- which(padj < alpha)
length(outliers2) # 132

# Bonferroni correction
padj <- p.adjust(x$pvalues, method = "bonferroni")
outliers3 <- which(padj < alpha)
length(outliers3) # 43

# Association between PCs and outliers
snp_pc <- get.pc(x, outliers1)
snp_pc <- as.data.frame(snp_pc)
snp_pc$PC <- as.factor(snp_pc$PC)
snp_pc$SNP <- as.factor(snp_pc$SNP)
str(snp_pc)
summary(snp_pc)

lm.loc.names <- as.data.frame(lm.genlight@loc.names)

lm_outliers <- find.snps(lm.loc.names, snp_pc)
length(lm_outliers$V1) # 132

lm_outliers_list <- as.data.frame(lm_outliers$V1)
colnames(lm_outliers_list) <- "locus_pos"
lm_outliers_list <- lm_outliers_list %>% 
  separate(locus_pos, sep = "_", into = c("loc","position")) %>% 
  separate(loc, sep = 3, into = c("X", "locusnum")) %>% 
  separate(position, sep = 3, into = c("Y", "pos"))
lm_outliers_list$Z <- "RAD"
lm_outliers_list <- lm_outliers_list %>% 
  unite("RADNAME", "Z", "locusnum")

lm_outliers_list <- lm_outliers_list[,c(2, 4)]

write.table(lm_outliers_list,
            file = "lm_pcadapt_outliers.txt",
            row.names = FALSE, col.names = FALSE , quote = FALSE)


# Packages ####
library("vcfR")
library("adegenet")
library("dplyr")
library("tidyr")
library("hierfstat")
library("sjmisc")
library("pegas")

# Arguments
args <- commandArgs(TRUE)
taxa <- args[1]
missing.dat <- args[2]

# Importing VCF and popfile =====
pop <- read.csv(paste0("../data/pop_", taxa, "_1div_nc.txt"), header = FALSE, sep = "\t")
colnames(pop) <- c("Individual", "Site")
# mislabels
pop$Site[pop$Individual == "KP0406_NA_NANA"] = "WP20"
pop$Site[pop$Individual == "KP0844_UN_NANA"] = "WP10"
pop$Site[pop$Individual == "KP0723_NA_NANA"] = "SB20"
#pop$Site[pop$Individual == "KP0673_UN_NANA"] = "SB20"
pop <- pop %>% separate(Site, into = c("Loc", "Depth"), sep = 2, remove = FALSE)
pop$Site <- as.factor(pop$Site)
pop$Loc <- as.factor(pop$Loc)
pop$Depth <- as.factor(pop$Depth)
genind <- vcfR2genind(read.vcfR(paste0("../data/", taxa, "_1div_nc_", missing.dat, ".vcf")))
genind@pop <- pop$Loc
genlight <- vcfR2genlight(read.vcfR(paste0("../data/", taxa, "_1div_nc_", missing.dat, ".vcf")))

# Build a result table for fstatistics
results <- data.frame(Taxa = taxa,
                      HoL = NA, Ho = NA, HoU = NA,
                      HeL = NA, He = NA, HeU = NA,
                      fstL = NA, fst = NA, fstU = NA,
                      fitL = NA, fit = NA, fitU = NA,
                      fis = NA, fis = NA, fisU = NA)

# Heirfstat ####
g_h <- genind2hierfstat(genind)

statistics <- as.data.frame(boot.vc(g_h[1], g_h[-1])$ci)
results[,2:4] <- round(statistics$`Hobs`, 3)
results[, 5:7] <- round(statistics$`H-Total`, 3)
results[, 8:10] <- round(statistics$`F-pop/Total`, 3)
results[, 11:13] <- round(statistics$`F-Ind/Total`, 3)
results[, 14:16] <- round(statistics$`F-Ind/pop`, 3)

write.csv(results, paste0("../results/fstat/", taxa, "_", missing.dat, "_Fstatistics.csv"), quote = FALSE)

# Fstat between locations
matFst <- genet.dist(genind, method = "WC84")
print(matFst)
datFst <- as.matrix(matFst)

boot.ppfst(dat=g_h,nboot=100,quant=c(0.025,0.975))
boot.ppfis(dat=g_h,nboot=100,quant=c(0.025,0.975))

write.csv(datFst, paste0("../results/fstat/", taxa, "_", missing.dat, "_LocationFst.csv"), quote = FALSE)

# Individual heterozygosity scores ======
mat <- as.matrix(genlight)
dat <- as.data.frame(mat)

hetero_loci <- row_count(dat, count = 1, var = "het", append = TRUE)
hetero_loci <- row_count(hetero_loci, count = NA, var = "missdat", append = TRUE)

for (i in 1:nrow(hetero_loci)) {
  hetero_loci$num_loci[i] <- genlight$n.loc - hetero_loci$missdat[i]
  hetero_loci$het_frac[i] <- hetero_loci$het[i]/hetero_loci$num_loci[i]
}

indiv_het <- dplyr::select(hetero_loci, het, het_frac, num_loci)
high_het <- indiv_het[indiv_het$het_frac > 0.05,]
low_het <- indiv_het[indiv_het$het_frac < 0.05,]

plot(indiv_het$het_frac)

write.csv(indiv_het, file = paste0("../results/fstat/", taxa, "_", missing.dat, "_indiv_het.csv"), quote = FALSE)

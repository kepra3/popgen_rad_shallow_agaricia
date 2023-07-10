# Packages ####
library("vcfR")
library("adegenet")
library("hierfstat")
library("dplyr")
library("sjmisc")
library("pegas")

# Setting variables
setwd("~/Dropbox/agaricia_project_2019/shalo_ag/gen_project/2 - Population Structure/2d - Popgen stats/")

# Cryptic taxa stats (without clones) =====
taxa <- "ah3"
missing.dat <- "20"

# Importing VCF and popfile =====
pop <- read.csv(paste0('pop_', taxa, ".txt"), header = FALSE, sep = "\t")
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
genind <- vcfR2genind(read.vcfR(paste0(taxa, "_1div_nc_", missing.dat, ".vcf")))
genind@pop <- pop$Loc
genlight <- vcfR2genlight(read.vcfR(paste0(taxa, "_1div_nc_", missing.dat, ".vcf")))

# Heterozygosity and diversity =====
sum <- summary(genind)

print(paste("Mean of observed Heterozygosity:", round(mean(sum$Hobs), 2)))
# aa1 0.22
# aa2 0.13
# ah1 0.08
# ah3 0.06
# al1 0.21
# al2 0.16
print(paste("Mean of expected Heterozygosity:", round(mean(sum$Hexp), 2)))
# aa1 0.23
# aa2 0.13
# ah3 0.29
# al1 0.21
# al2 0.17
plot(sum$Hexp, sum$Hobs)
abline(a = 0, b = 1, col = "red")


het.hw.test <- hw.test(genind, B=100)
length(het.hw.test[het.hw.test[,4] < 0.05,4])
mean(het.hw.test[,4])
hwe_names <- names(het.hw.test[het.hw.test[,4] < 0.05,4])
plot(sum$Hexp, sum$Hobs)
plot(sum$Hexp[hwe_names], sum$Hobs[hwe_names])

het_bartlett.test <- bartlett.test(list(sum$Hexp,sum$Hobs))
het_bartlett.test 
# equal for aa1, aa2, ah1, al2
# unequal ah2, ah3, al1

het_t.test <- t.test(sum$Hexp,sum$Hobs, pair = T, var.equal = T, alter = "greater")
het_t.test
# aa1
# t = 3.028, df = 486, p-value = 0.001296
# mean difference  0.0101724 
# aa2
# t = 8.4303, df = 1463, p-value < 2.2e-16
# mean difference 0.006920864 
# ah1
#t = 23.375, df = 808, p-value < 2.2e-16
#mean difference 0.1249934
# ah2
#t = 24.774, df = 601, p-value < 2.2e-16
#mean difference 0.1417464 
# ah3
# t = 33.903, df = 622, p-value < 2.2e-16
# mean difference  0.2345531 
# al1
# t = 0.8963, df = 578, p-value = 0.1852
# mean difference  0.003219509 


# Weir & Cockham's FST ======
print(wc(genind))
# aa1
#$FST
#[1] 0.0930458
#$FIS
#[1] -0.0004519145
# aa2
#$FST
#[1] 0.03888143
#$FIS
#[1] 0.02714732
#ah2
#$FST
#[1] 0.06441596
#$FIS
#[1] 0.5226727

ftab <- Fst(as.loci(genind))
ftab <- as.data.frame(ftab)
ftab_removed <- ftab[!is.na(ftab$Fis) & !is.na(ftab$Fit) &  !is.na(ftab$Fst),]
print(colMeans(ftab_removed))
# aa1
# Fit          Fst          Fis 
# 0.075695288  0.078431864 -0.007040724 
# aa2
#Fit        Fst        Fis 
#0.05123940 0.02592600 0.02624643


# Heirfstat ####
g_h <- genind2hierfstat(genind)

statistics <- as.data.frame(boot.vc(g_h[1], g_h[-1])$ci)
round(statistics$`Hobs`, 3)
round(statistics$`H-Total`, 3)
round(statistics$`F-pop/Total`, 3)
round(statistics$`F-Ind/Total`, 3)
round(statistics$`F-Ind/pop`, 3)

print(boot.vc(g_h[1], g_h[-1])$ci)
# aa1
#       H-Total F-pop/Total F-Ind/Total  H-pop F-Ind/pop   Hobs
#2.5%   0.2303      0.0778      0.0642 0.2093   -0.0263 0.2082
#50%    0.2455      0.0922      0.0915 0.2231   -0.0011 0.2232
#97.5%  0.2606      0.1074      0.1197 0.2368    0.0252 0.2380
# aa2
#       H-Total F-pop/Total F-Ind/Total  H-pop F-Ind/pop   Hobs
#2.5%   0.1286      0.0353      0.0538 0.1235    0.0162 0.1197
#50%    0.1367      0.0390      0.0652 0.1313    0.0274 0.1278
#97.5%  0.1450      0.0428      0.0774 0.1393    0.0388 0.1356
# ah1
#       H-Total F-pop/Total F-Ind/Total  H-pop F-Ind/pop   Hobs
#2.5%   0.1971      0.0695      0.6023 0.1812    0.5686 0.0677
#50%    0.2103      0.0813      0.6409 0.1931    0.6087 0.0756
#97.5%  0.2233      0.0939      0.6724 0.2050    0.6420 0.0846
# ah2
#       H-Total F-pop/Total F-Ind/Total  H-pop F-Ind/pop   Hobs
#2.5%   0.2729      0.0495      0.5174 0.2553    0.4856 0.1180
#50%    0.2874      0.0639      0.5534 0.2687    0.5230 0.1281
#97.5%  0.3015      0.0809      0.5848 0.2820    0.5548 0.1400
# al1
#     H-Total F-pop/Total F-Ind/Total  H-pop F-Ind/pop   Hobs
#2.5%   0.2008     -0.0124     -0.0005 0.2016    0.0034 0.1922
#50%    0.2142     -0.0045      0.0332 0.2153    0.0376 0.2070
#97.5%  0.2276      0.0037      0.0671 0.2292    0.0730 0.2230
# al2
#     H-Total F-pop/Total F-Ind/Total  H-pop F-Ind/pop   Hobs
#2.5%   0.1584     -0.0052      0.0217 0.1591    0.0238 0.1498
#50%    0.1703     -0.0022      0.0481 0.1707    0.0499 0.1623
#97.5%  0.1812      0.0007      0.0702 0.1817    0.0728 0.1742

boot.ppfst(dat=g_h,nboot=100,quant=c(0.025,0.975))
boot.ppfis(dat=g_h,nboot=100,quant=c(0.025,0.975))

basic.stats(g_h)


matFst <- genet.dist(genind, method = "WC84")
print(matFst)
datFst <- as.matrix(matFst)
write.csv(datFst, paste0(taxa, "_", missing.dat, "_matFst.csv"), quote = FALSE)

# Heterozygosity scores ======
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

write.csv(indiv_het, file = paste0(taxa, "_", missing.dat, "_indiv_het.csv"), quote = FALSE)


# plotting heterozygosity ####
ah1_50_indiv_het <- read.csv("~/Dropbox/agaricia_project_2019/shalo_ag/gen_project/2 - Population Structure/2d - Popgen stats/ah1_50_indiv_het.csv")
ah1_50_indiv_het <- ah1_50_indiv_het %>%
  separate(X, into = c("Sample", "Species", "Site"))
ah1_50_indiv_het$Site <- as.factor(ah1_50_indiv_het$Site)
ggplot(data = ah1_50_indiv_het, aes(x = Site, y = het_frac, colour = Site)) + geom_point()

ah2_50_indiv_het <- read.csv("~/Dropbox/agaricia_project_2019/shalo_ag/gen_project/2 - Population Structure/2d - Popgen stats/ah2_50_indiv_het.csv")
ah2_50_indiv_het <- ah2_50_indiv_het %>%
  separate(X, into = c("Sample", "Species", "Site"))
ah2_50_indiv_het$Site <- as.factor(ah2_50_indiv_het$Site)
ggplot(data = ah2_50_indiv_het, aes(x = Site, y = het_frac, colour = Site)) + geom_point()

ah3_50_indiv_het <- read.csv("~/Dropbox/agaricia_project_2019/shalo_ag/gen_project/2 - Population Structure/2d - Popgen stats/ah3_50_indiv_het.csv")
ah3_50_indiv_het <- ah3_50_indiv_het %>%
  separate(X, into = c("Sample", "Species", "Site"))
ah3_50_indiv_het$Site <- as.factor(ah3_50_indiv_het$Site)
ggplot(data = ah3_50_indiv_het, aes(x = Site, y = het_frac, colour = Site)) + geom_point()

ah3_20_indiv_het <- read.csv("~/Dropbox/agaricia_project_2019/shalo_ag/gen_project/2 - Population Structure/2d - Popgen stats/ah3_20_indiv_het.csv")
ah3_20_indiv_het <- ah3_20_indiv_het %>%
  separate(X, into = c("Sample", "Species", "Site"))
ah3_20_indiv_het$Site <- as.factor(ah3_20_indiv_het$Site)
ggplot(data = ah3_20_indiv_het, aes(x = Site, y = het_frac, colour = Site)) + geom_point()


ah_pop <- read.csv("~/Dropbox/agaricia_project_2019/shalo_ag/gen_project/data/hu_1div_nc_20_4.csv")

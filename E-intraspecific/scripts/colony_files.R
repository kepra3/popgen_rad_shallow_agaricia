# Convert 2 COLONY
library(vcfR)
library(adegenet)

args = commandArgs(TRUE)
vcf_name <- args[1]
# if want to do manual
#vcf_name <- "aa1-wp_1div_nc_20"
setwd("~/Dropbox/agaricia_project_2019/shalo_ag/gen_project/3 - Spatial/3c - Kinship/")

# Make sure vcf is in the data directory
genind <- vcfR2genind(read.vcfR(paste0("../../data/", vcf_name, ".vcf")))

df <- genind2df(genind)
print(dim(df))

loci_list <- colnames(df)

# Convert 0 to 01 and 1 to 02 and NA's to 0000
mat <- as.matrix(df)
mat <- apply(mat, FUN = gsub, MARGIN = 2, pattern = "0", replacement = "A")
mat <- apply(mat, FUN = gsub, MARGIN = 2, pattern = "1", replacement = "T")
mat <- apply(mat, FUN = gsub, MARGIN = 2, pattern = "A", replacement = "1")
mat <- apply(mat, FUN = gsub, MARGIN = 2, pattern = "T", replacement = "2")
mat <- apply(mat, FUN = gsub, MARGIN = 2, pattern = "11", replacement = "1 1")
mat <- apply(mat, FUN = gsub, MARGIN = 2, pattern = "22", replacement = "2 2")
mat <- apply(mat, FUN = gsub, MARGIN = 2, pattern = "12", replacement = "1 2")
mat <- apply(mat, FUN = gsub, MARGIN = 2, pattern = "21", replacement = "2 1")
mat[is.na(mat)] = "0 0"

write.table(mat, file = paste0("genotypes", "_", vcf_name, ".txt"), sep = "\t",
      quote = FALSE, row.names = TRUE, col.names = FALSE)

write.table(t(loci_list), file = paste0("loci", "_", vcf_name, ".txt"), sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = FALSE)


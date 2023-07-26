# Convert 2 COLONY
# R 4.2.0
library(vcfR) # vcfR 1.12.0 (http://dx.doi.org/10.1111/1755-0998.12549)
library(adegenet) # adegenet 2.1.7 (doi: 10.1093/bioinformatics/btr521)

args = commandArgs(TRUE)
vcf_name <- args[1]

# Make sure vcf is in the data directory
genind <- vcfR2genind(read.vcfR(paste0("../data/", vcf_name, ".vcf")))

df <- genind2df(genind)
print(paste(vcf_name, dim(df)))

loci_list <- colnames(df)

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

write.table(mat, file = paste0("../data/genotypes", "_", vcf_name, ".txt"), sep = "\t",
      quote = FALSE, row.names = TRUE, col.names = FALSE)

write.table(t(loci_list), file = paste0("../data/loci", "_", vcf_name, ".txt"), sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = FALSE)


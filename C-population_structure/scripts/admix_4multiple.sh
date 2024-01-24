#!/usr/bin/env bash
# Author: Katharine Prata
# Date created: 03/03/21
# Last edit: 28/06/23

# argument variables
VCF_NAME=$1
STAGE=$2
PARAMS=$3
MAX_K=$4
percent=$5
run=$6

# 1. thin dataset to be unlinked
# use Pim's random sample script
vcf_single_snp.py ../data/${VCF_NAME}_${STAGE}_${PARAMS}_${percent}.vcf > ${VCF_NAME}_${STAGE}_${PARAMS}_${percent}_temp-unlinked.vcf

# 2. make directories for each analysis
mkdir ../results/admix_runs/${percent}percent_${VCF_NAME}_${STAGE}_${PARAMS}_${run}

# 3. convert to plink
plink --vcf ${VCF_NAME}_${STAGE}_${PARAMS}_${percent}_temp-unlinked.vcf --const-fid --allow-extra-chr 0 --make-bed --out ../results/admix_runs/${percent}percent_${VCF_NAME}_${STAGE}_${PARAMS}_${run}/${VCF_NAME}_2bi_${STAGE}_${PARAMS}_${percent}
echo "---------------------------------- Converted to plink ----------------------------------"

# 4. order binary file by population
# Ensure you have a 0 column for family ID before sample name!
# plink --bfile ${percent}percent_${VCF_NAME}_${PARAMS}/${VCF_NAME}_2bi_${PARAMS}_$percent --const-fid --indiv-sort f samples_${VCF_NAME}_1d_${PARAMS}.txt --allow-extra-chr 0 --make-bed --out ${percent}percent_${VCF_NAME}_${PARAMS}/${VCF_NAME}_2bi_${PARAMS}_${percent}_ordered
# echo "---------------------------------- Ordered by population ----------------------------------"

# 5. do admixture for k
echo "---------------------------------- Starting admixture ----------------------------------"
start=1
end=${MAX_K}
for ((K=start; K<=end; K++))
do cd ../results/admix_runs/${percent}percent_${VCF_NAME}_${STAGE}_${PARAMS}_${run}
admixture ${VCF_NAME}_2bi_${STAGE}_${PARAMS}_${percent}.bed $K -j8 --cv=10 -C=100 | tee log.$K.out
grep -h CV log*.out
grep -h ^Log log*.out
echo "---------------------------------- Finished admixture ----------------------------------"
done

grep -h CV log*.out >> CV_error.txt
grep -h ^Log log*.out >> loglikelihood.txt
cd ../../../scripts
rm ${VCF_NAME}_${STAGE}_${PARAMS}_${percent}_temp-unlinked.vcf

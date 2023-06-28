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

# 1. make directories for each analysis
mkdir ../results/${percent}percent_${VCF_NAME}_${STAGE}_${PARAMS}

# 2. convert to plink
plink --vcf ../vcf/${VCF_NAME}_${STAGE}_${PARAMS}_${percent}.vcf --const-fid --allow-extra-chr 0 --make-bed --out ../results/${percent}percent_${VCF_NAME}_${STAGE}_${PARAMS}/${VCF_NAME}_2bi_${STAGE}_${PARAMS}_${percent}
echo "---------------------------------- Converted to plink ----------------------------------"

# 3. do admixture for k
echo "---------------------------------- Starting admixture ----------------------------------"
cd ../results/${percent}percent_${VCF_NAME}_${STAGE}_${PARAMS}
touch log_cv.${VCF_NAME}_${STAGE}_${PARAMS}.txt
start=1
end=${MAX_K}
for ((K=start; K<=end; K++))
do cd ../results/${percent}percent_${VCF_NAME}_${STAGE}_${PARAMS}
admixture ${VCF_NAME}_2bi_${STAGE}_${PARAMS}_${percent}.bed $K -j8 --cv=10 -C=100 | tee log.$K.out
if [ $K -lt 2 ]
then
tail -3 log.$K.out | head -2 >> log_cv.${VCF_NAME}_${STAGE}_${PARAMS}.txt
else
tail -$((5 + $K)) log.$K.out | head $((4 + $K)) >> log_cv.${VCF_NAME}_${STAGE}_${PARAMS}.txt
fi
echo "\n\n" >> log_cv.${VCF_NAME}_${STAGE}_${PARAMS}.txt
grep -h CV log*.out
grep -h ^Log log*.out
cd ..
echo "---------------------------------- Finished admixture ----------------------------------"
done

cd ../results/${percent}percent_${VCF_NAME}_${STAGE}_${PARAMS}
grep -h CV log*.out >> CV_error.txt
grep -h ^Log log*.out >> loglikelihood.txt
cd ..

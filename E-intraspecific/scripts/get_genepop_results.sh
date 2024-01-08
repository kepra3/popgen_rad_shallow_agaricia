#!/bin/sh

#  get_genepop_results.sh
#  
#
#  Created by Katharine Prata on 1/6/22.
#  
taxa=$1 # "AA2"
cat=$2 # "all"
scale=$3 # "within"
dim=$4 # "2D"

results_file="../results/ibd/${taxa}/${taxa}_${cat}.${scale}.${dim}.results"

touch "${results_file}_short.txt"

printf %s "slope,s.lowCI,s.highCI,intercept,i.lowCI,i.highCI,p.slope,p.IBD,other.p" >> "${results_file}_short.txt"

echo "" >> "${results_file}_short.txt"

grep -A 2 'ABC bootstrap results for SLOPE:' "${results_file}.txt" | tail -n 1 | awk '{printf $1}' >> "${results_file}_short.txt"

printf %s "," >> "${results_file}_short.txt"

grep -A 2 'ABC bootstrap results for SLOPE:' "${results_file}.txt" | tail -n 1 | awk '{printf $3}' >> "${results_file}_short.txt"

printf %s "," >> "${results_file}_short.txt"

grep -A 2 'ABC bootstrap results for SLOPE:' "${results_file}.txt" | tail -n 1 | awk '{printf $5}' >> "${results_file}_short.txt"

printf %s "," >> "${results_file}_short.txt"

grep -A 2 'ABC bootstrap results for INTERCEPT:' "${results_file}.txt" | tail -n 1  | awk '{printf $1}' >> "${results_file}_short.txt"

printf %s "," >> "${results_file}_short.txt"

grep -A 2 'ABC bootstrap results for INTERCEPT:' "${results_file}.txt" | tail -n 1  | awk '{printf $3}' >> "${results_file}_short.txt"

printf %s "," >> "${results_file}_short.txt"

grep -A 2 'ABC bootstrap results for INTERCEPT:' "${results_file}.txt" | tail -n 1  | awk '{printf $5}' >> "${results_file}_short.txt"

printf %s "," >> "${results_file}_short.txt"

grep -h 'Uni' "${results_file}.txt" | awk '{printf $8}' >> "${results_file}_short.txt"

printf %s "," >> "${results_file}_short.txt"

grep -A 1 'Test of isolation by distance' "${results_file}.txt" | tail -n 1  | awk '{printf $5}'  >> "${results_file}_short.txt"

printf %s "," >> "${results_file}_short.txt"

grep -A 1 'Other one tailed Pvalue' "${results_file}.txt" | tail -n 1  | awk '{printf $5}'  >> "${results_file}_short.txt"

echo "Finished!"


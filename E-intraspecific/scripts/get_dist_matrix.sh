#!/bin/sh

#  get_dist_matrix.sh
#  
#
#  Created by Katharine Prata on 30/4/22.
#  

taxa=$1 # "AA1"
depth=$2 # "all"
scale=$3 # "within"
indiv=$4 # number of individuals

count=`echo -n "${scale}" | wc -c`

cat ../results/ibd/${taxa}/${taxa}_${depth}.${scale}.results_copy.txt | head -n $((${indiv} + 2)) | tail -n $((${indiv} - 1)) | tr -s ' ' | tr -s ' ' '\t' | sed 's/^[ \t]*//'> "../results/ibd/${taxa}/${taxa}_${depth}.a_copy.txt"

if [ $count -ge 10 ]
then
last_lines=27
else
last_lines=25
fi

cat ../results/ibd/${taxa}/${taxa}_${depth}.${scale}.results_copy.txt | tail -n $((${indiv} + ${last_lines})) | head -n $((${indiv} - 1 ))  | tr -s ' ' | tr -s ' ' '\t' | sed 's/^[ \t]*//'> ../results/ibd/${taxa}/${taxa}_${depth}.lndist_copy.txt



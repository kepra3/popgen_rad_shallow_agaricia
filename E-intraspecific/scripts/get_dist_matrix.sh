#!/bin/sh

#  get_dist_matrix.sh
#  
#
#  Created by Katharine Prata on 30/4/22.
#  

taxa=$1 # "AA2"
depth=$2 # "20
scale=$3 # "between"
indiv=$4 # number of individuals

count=`echo -n "${scale}" | wc -c`

cat ../results/${taxa}/${taxa}_${depth}.${scale}.results.txt | head -n $((${indiv} + 2)) | tail -n $((${indiv} - 1)) | tr -s ' ' | tr -s ' ' '\t' | sed 's/^[ \t]*//'> "../results/${taxa}/${taxa}_${depth}.a.txt"

if [ $count -ge 10 ]
then
last_lines=27
else
last_lines=25
fi

cat ../results/${taxa}/${taxa}_${depth}.${scale}.results.txt | tail -n $((${indiv} + ${last_lines})) | head -n $((${indiv} - 1 ))  | tr -s ' ' | tr -s ' ' '\t' | sed 's/^[ \t]*//'> ../results/${taxa}/${taxa}_${depth}.lndist.txt



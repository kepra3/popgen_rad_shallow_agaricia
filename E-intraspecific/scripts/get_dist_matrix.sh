#!/bin/sh

#  get_dist_matrix.sh
#  
#
#  Created by Katharine Prata on 30/4/22.
#  

taxa=$1 # "AA1"
category=$2 # "all"
scale=$3 # all, within, between
dimension=$4 # 1D or 2D

indiv=`head -n 1 ../results/ibd/${taxa}/${taxa}_${category}.${scale}.${dimension}.results.txt | cut -d" " -f1`

cat ../results/ibd/${taxa}/${taxa}_${category}.${scale}.${dimension}.results.txt | head -n $((${indiv} + 2)) | tail -n $((${indiv} - 1)) | cut -d ' ' -f3- | tr -s ' ' '\t' > "../results/ibd/${taxa}/${taxa}_a_copy.txt"

last_lines=27

cat ../results/ibd/${taxa}/${taxa}_${category}.${scale}.${dimension}.results.txt | tail -n $((${indiv} + ${last_lines})) | head -n $((${indiv} - 1 )) | cut -d ' ' -f4- | tr -s ' ' '\t' > ../results/ibd/${taxa}/${taxa}_${dimension}_copy.txt

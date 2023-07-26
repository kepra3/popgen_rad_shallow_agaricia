#!/bin/sh

#  edit_dat_file.sh
#  
#
#  Created by Katharine Prata on 9/6/22.
#


data_name=$1
offspring_numb=$2
loci_numb=$3
prob=$4

quote_data_name="'${data_name}'"
quote_output_name="'${data_name}_${prob}'"
output_name="${data_name}_$prob"
genotypes_filename="../data/genotypes_${data_name}.txt"
loci_filename="../data/loci_${data_name}.txt"
nl=$'\n'

# Name line numbers and number of offspring and loci
sed -e 1"s/.*/${quote_data_name} &/" ../data/new.dat | sed -e 2"s/.*/${quote_output_name} &/" > ${output_name}_1.dat


NUMBER_OFFSPRING_LINE=$(grep -n '! Number of offspring' ${output_name}_1.dat | cut -d : -f 1)
NUMBER_LOCI_LINE=$(grep -n '! Number of loci' ${output_name}_1.dat | cut -d : -f 1)
MARKER_LINE=$(grep -n '!Marker names' ${output_name}_1.dat | cut -d : -f 1)

sed -e ${NUMBER_OFFSPRING_LINE}"s/.*/${offspring_numb} &/" ${output_name}_1.dat | sed -e ${NUMBER_LOCI_LINE}"s/.*/${loci_numb} &/" | sed -e ${MARKER_LINE}"s/.*/`cat ${loci_filename}` &/" > ${output_name}_2.dat

# Insert offspring and parent genotypes
sed -e "N;/\n.*\!prob. of dad/{r ${genotypes_filename}" -e '};P;D' ${output_name}_2.dat | sed -e "N;/\n.*!#known father/{r ${genotypes_filename}" -e '};P;D' > ${output_name}_3.dat

FATER_LINE=$(grep -n '!#known father' ${output_name}_3.dat | cut -d : -f 1)

sed -e ${FATER_LINE}"s/.*/\\${nl} &/" ${output_name}_3.dat >  ${output_name}_4.dat


sed -e "N;/\n.*!#known father/{r ${genotypes_filename}" -e '};P;D' ${output_name}_4.dat > ${output_name}_5.dat

# PARAMETERS FOR PROB AND NUMBER OF CAND
PROB_OF_PARENTS_LINE=$(grep -n '!prob. of dad/mum' ${output_name}_5.dat | cut -d : -f 1)
sed -e ${PROB_OF_PARENTS_LINE}"s/.*/${prob} ${prob}&/" ${output_name}_5.dat > ${output_name}_6.dat

NUMBERS_OF_CAND_LINE=$(grep -n '!numbers of candiadte' ${output_name}_6.dat | cut -d : -f 1)
sed -e ${NUMBERS_OF_CAND_LINE}"s/.*/${offspring_numb} ${offspring_numb}&/" ${output_name}_6.dat > ../results/kinship/${output_name}.dat


rm ${output_name}_*

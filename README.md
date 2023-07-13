# Data accessibility for Prata et al (2023)

TODO:

- [x] PCA
- [x] Admixture
  - [x] log/CV error
- [ ] Combination
- [x] Spatial coords
- [x] FST
- [x] RDA
- [ ] IBD
- [ ] Kinship
- [ ] Clone distances
- [ ] Change names of all files
- [ ] Add all software versions and requirements

## raw files

Raw fastq files and initial vcf output from ipyrad are stored on the University of Queensland eSpace repository (hyperlink)

TODO: add ipyrad settings file to RDM & here

Ipyrad settings:

```

```



# popgen_rad_shallow_agaricia

A git hub repository for analysis datasets, scripts and results

## A - Filtered vcf file

**Files structure:**

|--blast_files/

​			|--sym_loci_to_remove_all-aga_nogaps_min4_n.txt [match e from Cladocopium genome]

​			|--sym_loci_to_remove_all-aga_SymRAD16_n.txt [match e from RAD symbiont extracts todo: hyperlink to pim github]

​			|--RAD_other_loci_to_remove.txt [match e from NCBI nucleotide database]

|--high-miss-indiv_all-aga_1ci.txt [removed individuals with high missing data]

|--pcadapt_outliers/

​		|--data/

​		|--results/

​				|--ac_pcadapt_outliers.txt

​				|--hu_pcadapt_outliers.txt

​				|--lm_pcadapt_outliers.txt

​		|--scripts/

​				|--pcadapt_shallow.R

**Notes**

Steps are all written out instead of all code commands as the commands used are relatviely simple and repetitive.

**Steps**

1. Activate conda environment to access genetics softwares.

   ```bash
   $ conda activate radkat
   ```

1. The initial vcf [found on RDM] was filtered for minimum and maximum depth (> 3x mean depth), minimum allele count, maximum missing SNPs.

```bash
$ vcftools --vcf all-aga_min4_renamed.vcf --min-meanDP 5 --max-meanDP 1046 --mac 3 --max-missing 0.5 --min-alleles 2 --max-alleles 2 --recode --stdout > all-aga_1ci.vcf
```

2. Symbiont, other contamination (see `blast_files/`) and high missing data individuals (see `high-miss-indiv_all-aga_1ci.txt`) were removed, then final filters were applied, 5, 10, 20% maximum missing data, & monomorphic SNPs removed.
3. Initial PCA, NJ-tree and Admixture analyses were conducted in order to separate into species datasets for outlier removal (see *C - Population structure*). VCF files found within *C - Population structure* already have outliers removed, there were no significant differences among the structure results for these datsets. Individuals that assigned to each of the species were subset from the initial vcf and filtering, using the same steps as above, was repeated.
4. For most analyses outliers were removed and datasets were subset to one SNP per RAD contig. Outlier SNPs were discovered per species dataset using pcadapt (see `pcadapt_outliers/`).

## B - Clones

**File structure**

|--data/

​		|--ac_1d_wc_20_gdmatrix.tsv

​		|--clone_matches_ac_1d_wc_20.csv

​		|--clone_matches_hu_1d_wc_20.csv

​		|--clone_matches_lm_1d_wc_20.csv

​		|--hu_1d_wc_20_gdmatrix.tsv

​		|--lm_1d-pure_wc_20_gdmatrix.tsv

​		|--pop_ac_1d_wc_20.txt

​		|--pop_hu_1d_wc_20.txt

​		|--pop_lm_1d-pure_wc_20.txt

|--scripts

​		|--make_clone_groups.py

|--results/

​		|--chosen_individuals_ac.csv

​		|--chosen_individuals_hu.csv

​		|--chosen_individuals_lm.csv

​		|--clone_groups_ac.csv

​		|--clone_groups_hu.csv

​		|--clone_groups_lm.csv

​		|--distribution_pairwise_sum_ac_1d_wc.pdf

​		|--distribution_pairwise_sum_hu_1d_wc.pdf

​		|--distribution_pairwise_sum_lm_1d-pure_wc.pdf

​		|--individuals2keep_ac.txt

​		|--individuals2keep_hu.txt

​		|--individuals2keep_lm.txt

**Notes**

The script `vcf_find_clones.py` from www.github.com/pimbongaerts/radseq makes the `clone_matches_*` files.

**Steps**

These datasets were not filtered for unlinked and neutral SNPs.

1. Each separate species dataset was tranformed into a genetic distance matrix (Hamming's distance)
2. Run make_clone_groups.py to assess similarity distribution threshold, based upon breaks in distribution and similarity to technical replicates
3. Retain one individual per clone groups with the highest number of SNPs to create no clones datasets ('nc')

**Code**

```bash
$ conda activate radkat
$ make_clone_groups.py ac_1d_wc_20
$ make_clone_groups.py hu_1d_wc_20
$ make_clone_groups.py lm_1d-pure_wc_20
```

## C - Population structure

**File structure**

C-population_structure/

​	|--data/

​			|--ac_1div_nc_20_4.csv

​			|--all-aga_1div_nc-wnr_20_4.csv

​			|--hu_1div_nc_20_4.csv

​			|--lm_1div_nc_20_2.csv

​			|--lm_1div_nc-wnr_20_2.csv

​			|--pop_lm_1div_nc-wnr.txt

​			|--pop_lm_1div_nc.txt

​			|--ac_1div_nc_20.vcf

​			|--all-aga_1div_nc_wnr.vcf 

​			|--hu_1div_nc_20.vcf

​			|--lm_1div_nc_20.vcf

​			|--lm_1div_nc-wnr_20.vcf

​	|--results/

​			|--admix_plots/



​			|--admix_runs/

​			|--pca/

​					|--ac_stats/

​							|--ac_1div_nc_20_depth_facet_pca-x12.pdf

​							|--ac_2div_nc_20_depth_facet_pca-x34.pdf

​							|--ac_1div_nc_20_eig_plot.pdf

​							|--ac_1div_nc_20_eig.csv

​							|--ac_1div_nc_20_pcscores4.csv

​					|--all-aga_stats/

​					|--hu_stats/

​					|--lm_stats/

​	|--scripts/

​			|--admix_4multiple.sh

​			|--admix_plots2.R

​			|--admix_unlinked.sh

​			|--basic_pca.R

​			|--compile_logs.R

​			|--func-k.R

​			|--Qvalues.R

**Notes**

These analysis datasets do not include clones (removed during B - Clones), have 1 SNP per locus and have outliers removed from pcadapt.

Need to use `vcf_single_snp.py` from Pim Bongaerts github for creating unlinked datasets when using admix_4multiple.sh. Running admixture will mean results will slightly differ from manuscript results as runs are unseeded and random draws come from `vcf_single_snp.py`.

Haven't not included `all-aga_1diii_nc-wnr_20.vcf` due to being >100Mb (larger than the size limit for github).

**Steps**

1. Separate and filter species specific datasets

2. PCA

   -run using unlinked and neutral data and with depth category for individual species and species category for all-aga dataset.

3. Admixture

   -run admix_4multiple to create CV error and log-likelihood

   -make plots for CV error, log-likelihood and Qvalue taxa thresholds

   -run admix_unlinked to run on analysis dataset (*i.e.*, the sampled dataset used for all other analyses)

   -make admixture plots

4. PCA again

   -run with admixture setting to see concordance and create input files for combination plots

5. Combination plots (with NJ-tree)

**Code**

```bash
$ conda activate radkat
# TODO: show filtered final results files
```

### PCA

```bash
$ Rscript basic_pca.R ac_1div_nc_20 4 depth
$ Rscript basic_pca.R hu_1div_nc_20 4 depth
$ Rscript basic_pca.R lm_1div_nc-wnr_20 4 depth
$ Rscript basic_pca.R lm_1div_nc_20 2 depth
$ Rscript basic_pca.R all-aga_1div_nc-wnr_20 6 species
```

### Admixture

```bash
$ ./admix_unlinked.sh ac 1div nc 10 20
$ Rscript admix_plots2.R ac 1div nc 20 4 no
$ ./admix_unlinked.sh hu 1div nc 10 20
$ Rscript admix_plots2.R hu 1div nc 20 4 no
$ ./admix_unlinked.sh all-aga 1div nc-wnr 10 20
$ Rscript admix_plots2.R all-aga 1div nc-wnr 20 6 no
$ ./admix_unlinked.sh lm 1div nc-wnr 10 20
$ Rscript admix_plots2.R lm 1div nc-wnr 20 4 no
$ ./admix_unlinked.sh lm 1div nc 10 20
$ Rscript admix_plots2.R lm 1div nc 20 4 no
### Run Admixture on multiple datasets per species 
$ for run in 2 3 4 5 6 7 8 9 10; do ./admix_4multiple.sh ac 1diii nc 10 20 $run; done
$ for run in 2 3 4 5 6 7 8 9 10; do ./admix_4multiple.sh hu 1diii nc 10 20 $run; done
$ for run in 2 3 4 5 6 7 8 9 10; do ./admix_4multiple.sh lm 1diii nc 10 20 $run; done
$ for run in 2 3 4 5 6 7 8 9 10; do ./admix_4multiple.sh lm 1diii nc-wnr 10 20 $run; done
# $ for run in 2 3 4 5 6 7 8 9 10; do ./admix_4multiple.sh all-aga 1diii nc-wnr 10 20 $run; done
### CV error and log-likelihood plots
$ Rscript complie_logs.R
### Q plots
$ Rscript Qvalues.R all-aga_2bi_1div_nc-wnr_20 3
$ Rscript Qvalues.R ac_2bi_1div_nc_20 2
$ Rscript Qvalues.R hu_2bi_1div_nc_20 3
$ Rscript Qvalues.R lm_2bi_1div_nc_20 2
```

### Combination plots

```bash
$
```



## D - Spatial coordinates 

**File structure:**

|--data/

​		|--

|--scripts/

​		|--

|--results/



**Notes:**



**Steps:**

1. Colonies were annotated on point clouds using CloudCompare as three points, one centroid and two on the longest edge of colony.
2. Two transformations were performed on coordinates for each analysis:

   i) For Redundancy analyses, points annotated were rotated to align with 'real world up', thus the z-axis reflected changes in depth. Then plots were oriented according to geographical arrangments.

   ii) For isolation-by-distance, points annotated were rotated to onto the 2D XY-plane, thus any deviation in depth along the slope was reduced to 0. Then plots were oriented accordint to geographical arrangements.

**Code:**

```bash
$ conda activate open3d
# For redundancy analysis
# in data directory
$ cd data
$ json_file=(`ls *.json | tr '\n' ' '`)
$ text_file=(`ls *.txt | tr '\n' ' '`)
# in script directory
$ cd ../scripts
$ for ((i=0;i<${#json_file[@]};i++)); 
do python rotate_annotations_depthoriented.py ${json_file[$i]} ${text_file[$i]};
done
$ rotate_parallel_depthoriented.py
# For isolation-by-distance
$ for ((i=0;i<${#json_file[@]};i++)); 
do python rotate_annotations_2D.py ${json_file[$i]} ${text_file[$i]};
done
$ python rotate_annotations_2D.py cur_sna_20m_20200303_subsets.json cur_sna_20m_20200303_decvis_02.txt
$ python rotate_parallel_2D.py
$ Rscript set_distances.R all_annotations_X_HORIZ_parallel
$ Rscript set_distancces.R all_annotations_X_DEPTH_parallel
```

TODO: naming of files and inconsistency of scripts...

## E - Intraspecific analyses

**File structure:**

|--data/

​		|--

|--scripts/

​		|--

|--results/

​		|--

Notes:



**Steps:**

1. Remove mislabels and outgroups
2. Make separate taxa files for Fstatistics and Kinship
3. Run redundacy analysis
4. Run isolation-by-distance analysis
5. F-statistics
6. Kinship

**Code:**

```bash
$ conda activate radkat
# Remove mislabels with vcftools, adding a MAF filter as previously variant sites can become nonvariable after removing individuals.
$ vcftools --vcf ac_1div_nc_20.vcf --remove ac_mislabels.txt --min-alleles 2 --max-alleles 2 --maf 0.000000001 --recode --stdout > ac_3b_nc_20.vcf
#After filtering, kept 327 out of 335 Individuals
#After filtering, kept 1604 out of a possible 1606 Sites
$ vcftools --vcf hu_1div_nc_20.vcf --remove hu_mislabels.txt --min-alleles 2 --max-alleles 2 --maf 0.000000001 --recode --stdout > hu_3b_nc_20.vcf
#AAfter filtering, kept 120 out of 121 Individuals
#After filtering, kept 1280 out of a possible 1282 Sites
# For lamarcki need to also remove nextrad samples that were used for cryptic lineage assignment from Prata et al., (2022)
$ vcftools --vcf ../data/lm_1div_nc-wnr_20.vcf --remove lm_mislabels.txt --min-alleles 2 --max-alleles 2 --recode --stdout > lm_3b_nc_20_x.vcf
#After filtering, kept 85 out of 98 Individuals
#Outputting VCF file...
#After filtering, kept 940 out of a possible 940 Sites
#Run Time = 0.00 seconds
$ vcftools --vcf lm_3b_nc_20_x.vcf --remove ../data/nextrad_lm-gr.txt --maf 0.000000001 --min-alleles 2 --max-alleles 2 --recode --stdout > lm_3b_nc_20.vcf
#After filtering, kept 79 out of 85 Individuals
#Outputting VCF file...
#After filtering, kept 929 out of a possible 940 Sites
#Run Time = 0.00 seconds
$ rm lm_3b_nc_20_x.vcf
```

**Making separate vcf files for each taxon, for Kinship and Fstatistics**

```bash
# All individuals for f-statistics
# A. agaricites
$ vcftools --vcf ac_1div_nc_20.vcf --keep aa1.names.txt --min-alleles 2 --max-alleles 2 --max-missing 0.8 --maf 0.000000001 --stdout --recode > aa1_1div_nc_20.vcf
#After filtering, kept 35 out of 335 Individuals
#After filtering, kept 465 out of a possible 1606 Sites
$ vcftools --vcf ac_1div_nc_20.vcf --keep aa2.names.txt --min-alleles 2 --max-alleles 2 --max-missing 0.8 --maf 0.000000001 --stdout --recode > aa2_1div_nc_20.vcf
#After filtering, kept 300 out of 335 Individuals
#After filtering, kept 1461 out of a possible 1606 Sites

# A. humilis
$ vcftools --vcf hu_1div_nc_20.vcf --keep ah1.names.txt --min-alleles 2 --max-alleles 2 --max-missing 0.8 --maf 0.000000001 --stdout --recode > ah1_1div_nc_20.vcf
#After filtering, kept 61 out of 121 Individuals
#After filtering, kept 803 out of a possible 1282 Sites

$ vcftools --vcf hu_1div_nc_20.vcf --keep ah2.names.txt --min-alleles 2 --max-alleles 2 --max-missing 0.8 --maf 0.000000001 --stdout --recode > ah2_1div_nc_20.vcf
#After filtering, kept 29 out of 121 Individuals
#After filtering, kept 565 out of a possible 1282 Sites

$ vcftools --vcf hu_1div_nc_20.vcf --keep ah3.names.txt --min-alleles 2 --max-alleles 2 --max-missing 0.8 --maf 0.000000001 --stdout --recode > ah3_1div_nc_20.vcf
#After filtering, kept 15 out of 121 Individuals
#After filtering, kept 580 out of a possible 1282 Sites

# A. lamarcki
$ vcftools --vcf lm_1div_nc_20.vcf --keep al1.names.txt --min-alleles 2 --max-alleles 2 --max-missing 0.8 --maf 0.000000001 --stdout --recode > al1_1div_nc_20.vcf
#After filtering, kept 28 out of 92 Individuals
#After filtering, kept 544 out of a possible 945 Sites

$ vcftools --vcf lm_1div_nc_20.vcf --keep al2.names.txt --min-alleles 2 --max-alleles 2 --max-missing 0.8 --maf 0.000000001 --stdout --recode > al2_1div_nc_20.vcf
#After filtering, kept 62 out of 92 Individuals
#After filtering, kept 734 out of a possible 945 Sites

# Spatial individuals for kinship # TODO: finish separating files
$ vcftools --vcf ac_3b_nc_20.vcf --keep AA1_names.txt --min-alleles 2 --max-alleles 2 --max-missing 0.5 --maf 0.000000001 --stdout --recode > aa1_3b_nc_50.vcf
#After filtering, kept 30 out of 327 Individuals
#After filtering, kept 482 out of a possible 1604 Sites
```

### Redundancy analysis

```bash
$ for taxa in AA1 AA2 AH1 AH2 AH3 AL1 AL2;
do Rscript rda.R $taxa;
done
```

### Isolation-by-distance

```bash
# Running genepop across all locations
$ for taxa in AA1 AH1 AH2 AH3 AL1 AL2;
do Rscript genepop.R $taxa all all; # TODO: need to alter vcf names within script, ac_3b_nc_20
done
# Running genepop across all locations, for within 'locations', e.g., all depths at one spatial location 
$ for taxa in AA1 AH1 AH2 AH3 AL1 AL2;
do Rscript genepop.R $taxa all within;
done
# Summarise results
$ for taxa in AA1 AH1 AH2 AH3 AL1 AL2;
do ./get_genepop_results.sh $taxa all within;
done
$ for taxa in AA1 AH1 AH2 AH3 AL1 AL2;
do ./get_genepop_results.sh $taxa all all;
done
# Make plots
$ for taxa in AA1 AH1 AH2 AH3 AL1 AL2;
do Rscript ibd_plots_new.R $taxa all all;
done
Rscript ibd_plots_new.R  AA1 all within
```

### F-statistics

```bash
$ for taxa in aa1 aa2 ah1 ah2 ah3 al1 al2;
do Rscript popgen_stats-taxa.R $taxa 20;
done
```

### Kinship

```bash
$ Rscript colony_files.R aa1_1div_nc_50
$ for param in 0.5 0.2 0.1 0.05;
do ./edit_dat_file.sh aa1_1div_nc_50 35 487 $param
$ ./colony2s.out IFN:aa1_1div_nc_50_0.2.dat
$ ./colony2s.out IFN:aa1_1div_nc_50_0.5.dat
$ ./colony2s.out IFN:aa1_1div_nc_50_0.05.dat
$ ./colony2s.out IFN:aa1_1div_nc_50_0.1.dat
# all other taxa run the same except AA2 where datafile was split per location for computation time speed up.
```

### Clone distances

```

```



## Appendix

### Conda environments

Details for different conda environments

### Script descriptions



### Filename descriptions


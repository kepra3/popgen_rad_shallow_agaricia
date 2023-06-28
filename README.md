# Data accessibility for Prata et al (2023)

TODO:

- [x] place on github
- [x] clean make_clone_groups.py
  - [x] make sure everything runs
  - [x] annotate everything
  - [ ] check inconsistencies with newly generated files and original files
  - [x] add versions for each package
- [ ] clean PCA script
  - [ ] make sure everything runs
  - [ ] annotate everything
  - [ ] add versions for each package
- [ ] clean Admixture script(s)
  - [ ] make sure everything runs
  - [ ] annotate everything
  - [ ] add versions for each package
- [ ] clean combination scripts
  - [ ] make sure everything runs
  - [ ] annotate everything
  - [ ] add versions for each package
- [ ] clean RDA scripts
  - [ ] make sure everything runs
  - [ ] annotate everything
  - [ ] add versions for each package
- [ ] clean IBD scripts
  - [ ] make sure everything runs
  - [ ] annotate everything
  - [ ] add versions for each package
- [ ] clean kinship scripts
  - [ ] make sure everything runs
  - [ ] annotate everything
  - [ ] add versions for each package
- [ ] clean FST scripts
  - [ ] make sure everything runs
  - [ ] annotate everything
  - [ ] add versions for each package
- [ ] clean rotate annotation scripts
  - [ ] make sure everything runs
  - [ ] annotate everything
  - [ ] add versions for each package



- script process, copy over files from previous directory, then edit to make result files



## raw files

Raw fastq files and initial vcf output from ipyrad are stored on the University of Queensland eSpace repository (hyperlink)

TODO: add ipyrad settings file to RDM

- can upload large files special way (LFS)

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

​		|--scripts/

​				|--pcadapt_shallow.R

**Steps**

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

**Steps**

These datasets were not filtered for unlinked and neutral SNPs.

1. Each separate species dataset was tranformed into a genetic distance matrix (Hamming's distance)
2. Run make_clone_groups.py to assess similarity distribution threshold, based upon breaks in distribution and similarity to technical replicates
3. Retain one individual per clone groups with the highest number of SNPs to create no clones datasets ('nc')

**Notes**

This was run using a conda environment radkat see appendix for details

The script `vcf_find_clones.py` from www.github.com/pimbongaerts/radseq makes the `clone_matches_*` files

**Code**

```bash
$ make_clone_groups.py "ac_1d_wc_20" # make files for A. agaricites
$ make_clone_groups.py "hu_1d_wc_20" # make files for A. humilis
$ make_clone_groups.py "lm_1d-pure_wc_20" # make files for A. lamarcki
```

## C - Population structure

**File structure**

population_structure/

​	|--metadata/

​			|--ac_1div_nc_20_4.csv [*A. agaricites* metadata]

​			|--all-aga_1div_nc-wnr_20_4.csv [all species with NextRAD samples* metadata]

​			|--hu_1div_nc_20_4.csv [*A. humilis* metadata]

​			|--lm_1div_nc_20_2.csv [*A. lamarcki* metadata]

​			|--lm_1div_nc-wnr_20_2.csv [*A. lamarcki* with NextRAD samples* metadata]

​			|--pop_lm_1div_nc-wnr.txt

​			|--pop_lm_1div_nc.txt

​	|--scripts/

​			|--basic_pca.R

​			|--func-k.R

​	|--vcf/

​			|--ac_1div_nc_20.vcf [*A. agaricites* vcf]

​			|--all-aga_1div_nc_wnr.vcf [all species with NexRAD samples* vcf]

​			|--hu_1div_nc_20.vcf [*A. humilis* vcf]

​			|--lm_1div_nc_20.vcf [*A. lamarcki* vcf]

​			|--lm_1div_nc-wnr_20.vcf [*A. lamarcki* with NextRAD samples* vcf]

\* NextRAD samples include X *A. grahamae* and X *A. lamarcki* samples resequence with current data and provided from Prata et al. (2022) (todo: add hyperlink)

These analysis datasets do not include clones, have 1 SNP per locus and have outliers removed from pcadapt.

**Steps**

1. PCA
2. Admixture
3. Combination plot

**Note**



**Code**

```bash
## PCA
$ Rscript basic_pca.R ac_1div_nc_20 4 depth
$ Rscript basic_pca.R hu_1div_nc_20 4 depth
$ Rscript basic_pca.R lm_1div_nc-wnr_20 4 depth
$ Rscript basic_pca.R lm_1div_nc_20 2 depth
$ Rscript basic_pca.R all-aga_1div_nc-wnr_20 6 species
## Admixture

## Combination plot
```



## D - Spatial coordinates 

**File structure:**

|--data/

​		|--

|--scripts/

​		|--

|--results/



What to do?

-- mention that 3d models will be made available??

-- provide annotation files in data/

-- scripts for transformations in scripts/

-- results fro transformation metdata/

**Steps**

1. Colonies were annotated on point clouds using CloudCompare as three points, one centroid and two on the longest edge of colony.
2. Two transformations 



## E - Intraspecific analyses

**File structure:**

|--data/

​		|--

|--scripts/

​		|--

|--results/



(todo: clean FST, KINSHIP, RDA, IBD)

**Steps**

For IBD and Kinship vcf files had to be converted into



## Appendix

Details for different conda environments

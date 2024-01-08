TODO:

- [ ] Fix get_dist_matrix.sh in E-intraspecific

# Data accessibility for Prata et al (202X)

This repository contains all scripts and datasets required to run all analyses performed within Prata et al (202X) Some dominant reef-builders disperse only metres per generation. Large files such as raw sequence data and intial vcf files are not present within the repository but links to the sources to access these files are provided.

Blocks of code are provided in README on how to perform analyses for each section. For each block of code to run activate the appropriate conda environment or software versions (see *Appendix*) and change to the scripts directory within the appropriate heading directory. Files will either be accessed in data or results and will write files into data or results.

## raw files and de novo assembly

Raw fastq files and initial vcf output from ipyrad are stored on the University of Queensland eSpace repository and sequences will be uploaded to NCBI SRA database.

*De novo* assembly of RAD contigs was conducted using ipyrad (v.09.67) and performed on HPCs provided by the University of Queensland and California Academy Sciences.

Ipyrad settings:

```
------- ipyrad params file (v.0.9.67)-------------------------------------------
all_aga                        ## [0] [assembly_name]: Assembly name. Used to name output directories for assembly steps
/gpfs1/scratch/30days/uqkprat2/pyrad ## [1] [project_dir]: Project dir (made in curdir if not present)
                               ## [2] [raw_fastq_path]: Location of raw non-demultiplexed fastq files
                               ## [3] [barcodes_path]: Location of barcodes file
/gpfs1/scracth/30days/uqkprat2/stacks/KP01_14_set/*fq.gz                               ## [4] [sorted_fastq_path]: Location of demultiplexed/sorted fastq files
denovo                         ## [5] [assembly_method]: Assembly method (denovo, reference)
                               ## [6] [reference_sequence]: Location of reference sequence file
pairddrad                            ## [7] [datatype]: Datatype (see docs): rad, gbs, ddrad, etc.
TGCAG,CGG                         ## [8] [restriction_overhang]: Restriction overhang (cut1,) or (cut1, cut2)
5                              ## [9] [max_low_qual_bases]: Max low quality base calls (Q<20) in a read
33                             ## [10] [phred_Qscore_offset]: phred Q score offset (33 is default and very standard)
10                              ## [11] [mindepth_statistical]: Min depth for statistical base calling
10                              ## [12] [mindepth_majrule]: Min depth for majority-rule base calling
10000                          ## [13] [maxdepth]: Max cluster depth within samples
0.85                           ## [14] [clust_threshold]: Clustering threshold for de novo assembly
0                              ## [15] [max_barcode_mismatch]: Max number of allowable mismatches in barcodes
2                              ## [16] [filter_adapters]: Filter for adapters/primers (1 or 2=stricter)
35                             ## [17] [filter_min_trim_len]: Min length of reads after adapter trim
2                              ## [18] [max_alleles_consens]: Max alleles per site in consensus sequences
0.05                           ## [19] [max_Ns_consens]: Max N's (uncalled bases) in consensus
0.05                           ## [20] [max_Hs_consens]: Max Hs (heterozygotes) in consensus
4                              ## [21] [min_samples_locus]: Min # samples per locus for output
0.2                            ## [22] [max_SNPs_locus]: Max # SNPs per locus
8                              ## [23] [max_Indels_locus]: Max # of indels per locus
0.5                            ## [24] [max_shared_Hs_locus]: Max # heterozygous sites per locus
0, 0, 0, 0                     ## [25] [trim_reads]: Trim raw read edges (R1>, <R1, R2>, <R2) (see docs)
0, 0, 0, 0                     ## [26] [trim_loci]: Trim locus edges (see docs) (R1>, <R1, R2>, <R2)
p, s, v                        ## [27] [output_formats]: Output formats (see docs)
                               ## [28] [pop_assign_file]: Path to population assignment file
                               ## [29] [reference_as_filter]: Reads mapped to this reference are removed in step 
```

## A - Filtered vcf file

**Notes**

Some steps are written out instead of all code commands due to utilising databases and other large files.

**Steps**

1. Activate conda environment to access genetics softwares.

   ```bash
   $ conda activate radkat
   ```

1. The initial VCF [on UQ eSpace]() was filtered for minimum and maximum depth (>3x mean depth), minimum allele count, maximum missing SNPs (50%).

```bash
$ vcftools --vcf all-aga_min4_renamed.vcf --min-meanDP 5 --max-meanDP 1046 --mac 3 --max-missing 0.5 --min-alleles 2 --max-alleles 2 --recode --stdout > all-aga_1ci.vcf
```

2. Symbiont and other contamination was found using the `blastn` tool (see `blast_files/` for results). The loci file from ipyrad output was converted into a fasta file using `pyrad2fasta.py`. For Symbiodiniacae the *Cladocopium* spp. genome fasta file was downloaded from https://www.ncbi.nlm.nih.gov/assembly/GCA_003297045.1 and a RAD reference from www.github.com/pimbongaerts/bermuda-rad, reads that matched an e-value of 1e^-15^ were removed. For the other symbionts and contamination, the web-based software [Galaxy](https://usegalaxy.org.au/) was used and the fasta file was matched with the `nt 17-Apr-2014` database for non-Scleractinia sequences using an e-value of 0.001 which were then removed.
3. High missing data individuals (>50% missing data, see `high-miss-indiv_all-aga_1ci.txt`) were removed, then final filters were applied, 5, 10, 20% maximum missing data, & monomorphic SNPs removed.
4. Initial PCA, Admixture, and genetic distance analyses were conducted in order to separate into species datasets for outlier removal (see *C - Population structure*). VCF files found within *C - Population structure* already have outliers removed, there were no significant differences among the structure results for these datsets. Individuals that assigned to each of the species were subset from the initial VCF and filtering, using the same steps as above, was repeated.
5. For all analyses (apart from clone analyses) outliers were removed and datasets were subset to one SNP per RAD contig. Outlier SNPs were discovered per species dataset using pcadapt.

**Code:**

```bash
$ for taxa in ac hu lm;
		do Rscript pcadapt.R ${taxa};
		done
```

### Base datasets for analyses

| Dataset name               | # individuals | # SNPs  |                           Details                            | Use                                           | Locations                                                    |
| -------------------------- | ------------- | ------- | :----------------------------------------------------------: | --------------------------------------------- | ------------------------------------------------------------ |
| all-aga_1d_wc-wnr_20.vcf   | 767           | 15659   | All *Agaricia* species samples, linked, no outlier removal, with clones and NextRAD samples and 20% max missing data | Clones                                        | -                                                            |
| all-aga_1div_nc-wnr_20.vcf | 557           | 751     | All *Agaricia* species samples, unlinked, neutral, no clones with NextRAD samples and 20% max missing data | Population structure                          | C - population structure/data                                |
| ac_1d_wc_20.vcf            | 512           | 22274   | *A. agaricites* samples, linked, no outlier removal, with clones and 20% max missing data | Clones                                        | Distance matrix stored in B - Clones/data                    |
| ac_1div_nc_20.vcf          | 335           | 1606    | *A. agaricites* samples, unlinked, neutral, no clones and 20% max missing data | Population structure & Intraspecific analyses | C - population strucutre/data, D - intraspeific analyses/data |
| hu_1d_wc_20.vcf            | 142           | 25380   | *A. humilis* samples, linked, no outlier removal, with clones and 20% max missing data | Clones                                        | Distance matrix stored in B - Clones/data                    |
| hu_1div_nc_20.vcf          | 121           | 1282    | *A. humilis* samples, unlinked, neutral, no clones and 20% max missing data | Population structure & Intraspecific analyses | C - population strucutre/data, D - intraspeific analyses/data |
| lm_1d-pure_wc_20.vcf       | 105           | *18541* | *A. lamarcki* samples (including *A. lamarcki* NextRAD samples), linked, no outlier removal, with clones and 20% max missing data | Clones                                        | Distance matrix stored in B - Clones/data                    |
| lm_1div_nc_20.vcf          | 92            | 941     | *A. lamarcki* samples, unlinked, neutral, no clones and 20% max missing data | Population structure & Intraspecific analyses | C - population strucutre/data, D - intraspeific analyses/data |
| lm_1div_nc-wnr_20.vcf      | 98            | 940     | *A. lamarcki* samples, unlinked, neutral, no clones with NextRAD samples, and 20% max missing data | Population structure & Intraspecific analyses | C - population strucutre/data, D - intraspeific analyses/data |

For intraspecific analyses these were further filtered see *E - Intraspecific analyses* section for more details.

## B - Clones

**Notes**

The script `vcf_find_clones.py` from www.github.com/pimbongaerts/radseq makes the `clone_matches_*` files.

**Steps**

These datasets were not filtered for unlinked and neutral SNPs.

1. Each separate species dataset was tranformed into a genetic distance matrix (Hamming's distance)
2. Run `make_clone_groups.py` to assess similarity distribution threshold, based upon breaks in distribution and similarity to technical replicates
3. Retain one individual per clone groups with the highest number of SNPs to create no clones datasets ('nc')

**Code**

```bash
$ conda activate radkat
$ make_clone_groups.py ac_1d_wc_20
$ make_clone_groups.py hu_1d_wc_20
$ make_clone_groups.py lm_1d-pure_wc_20
```

## C - Population structure

**Notes**

These analysis datasets do not include clones (removed during B - Clones), have 1 SNP per locus and have outliers removed from pcadapt.

Need to use `vcf_single_snp.py` from www.github.com/pimbongaerts/radseq for creating unlinked datasets when using `admix_4multiple.sh`. Running admixture will mean results will slightly differ from manuscript results as runs are unseeded and random draws come from `vcf_single_snp.py`.

Haven't not included `all-aga_1diii_nc-wnr_20.vcf` due to being >100Mb.

For creating the combination plots (*e.g.*, Admixture and NJ-tree summary plots) need to use `vcf_gdmatrix.py` & `gdmatrix2tree.py` from www.github.com/pimbongaerts/radseq.

**Steps**

1. Separate and filter species specific datasets (see *A - Filtered vcf* file section)

2. PCA

   -run `basic_pca.R` with depth category for individual species and species category for all-aga dataset (*e.g.*, points are either coloured by location and depth or species identity)

3. Admixture

   -run `admix_4multiple.sh` to create multiple replicate runs of randomly sampled SNPs

   -make plots for CV error, log-likelihood and Qvalue taxa thresholds using `compile_logs.R` and `Qvalues.R`	

   -run `admix_unlinked.sh` to run on analysis dataset (*i.e.*, the sampled dataset used for all other analyses)

   -make admixture plots

4. PCA again

   -run with admixture category to see concordance and create input files for combination plots

5. Make combination plots with Admixture and neighbour-joining tree.

**Code**

```bash
$ conda activate radkat
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
$ Rscript compile_logs.R
### Q plots
$ Rscript Qvalues.R all-aga_2bi_1div_nc-wnr_20 3
$ Rscript Qvalues.R ac_2bi_1div_nc_20 2
$ Rscript Qvalues.R hu_2bi_1div_nc_20 3
$ Rscript Qvalues.R lm_2bi_1div_nc_20 2
```

### Combination plots

```bash
# PCA with Admixture
$ Rscript basic_pca.R ac_1div_nc_20 4 admixture
$ Rscript basic_pca.R hu_1div_nc_20 4 admixture
$ Rscript basic_pca.R lm_1div_nc-wnr_20 4 admixture
$ Rscript basic_pca.R lm_1div_nc_20 2 admixture
# Make Neighbour Joining trees
$ for taxa in ac hu lm;
		do vcf_gdmatrix.py ../data/${taxa}_1div_nc_20.vcf ../data/pop_${taxa}_1div_nc.txt > ../results/${taxa}_1div_nc_20_gdmatrix.tsv; 
		done
$ for taxa in ac hu lm;
		do gdmatrix2tree.py ../results/${taxa}_1div_nc_20_gdmatrix.tsv ../results/${taxa}_1div_nc_20_tree.nex;
		done
# Run combination plot script
$ for taxa in ac hu lm;
		do Rscript combining_popstruc.R ${taxa} 1div nc Taxa no;
		done
```

## D - Spatial coordinates 

**Steps**

1. Colonies were annotated on point clouds using CloudCompare as three points, one centroid and two on the longest edge of colony.
2. Two transformations were performed on coordinates for each analysis:

   i) For Redundancy analyses, points annotated were rotated to align with 'real world up', thus the z-axis reflected changes in depth. Then plots were oriented according to geographical arrangments.

   ii) For isolation-by-distance, points annotated were rotated to onto the 2D XY-plane, thus any deviation in depth along the slope was reduced to 0. Then plots were oriented accordint to geographical arrangements.

**Code**

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

## E - Intraspecific analyses

**Steps**

1. Remove mislabels and outgroups for spatial analyses
2. Make separate taxa files for Fstatistics and Kinship
3. Redundancy analysis
4. Isolation-by-distance
5. F-statistics
6. Kinship (requires HPC to run quickly)

**Code**

```bash
$ conda activate radkat
# Remove mislabels with vcftools, adding a MAF filter as previously variant sites can become nonvariable after removing individuals.
$ vcftools --vcf ac_1div_nc_20.vcf --remove ac_mislabels.txt --min-alleles 2 --max-alleles 2 --maf 0.000000001 --recode --stdout > ac_3b_nc_20.vcf
#After filtering, kept 327 out of 335 Individuals
#After filtering, kept 1604 out of a possible 1606 Sites
$ vcftools --vcf hu_1div_nc_20.vcf --remove hu_mislabels.txt --min-alleles 2 --max-alleles 2 --maf 0.000000001 --recode --stdout > hu_3b_nc_20.vcf
#AAfter filtering, kept 120 out of 121 Individuals
#After filtering, kept 1280 out of a possible 1282 Sites
# For lamarcki need to also remove nextrad samples that were used for cryptic lineage assignment
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
# A. agaricites
$ vcftools --vcf ac_1div_nc_20.vcf --keep aa1.names.txt --min-alleles 2 --max-alleles 2 --max-missing 0.8 --maf 0.000000001 --stdout --recode > aa1_1div_nc_20.vcf
#After filtering, kept 35 out of 335 Individuals
#After filtering, kept 465 out of a possible 1606 Sites
$ vcftools --vcf ac_1div_nc_20.vcf --keep aa2.names.txt --min-alleles 2 --max-alleles 2 --max-missing 0.8 --maf 0.000000001 --stdout --recode > aa2_1div_nc_20.vcf
#After filtering, kept 300 out of 335 Individuals
#After filtering, kept 1461 out of a possible 1606 Sites

## Subsetting A. agaricities 2 for Kinship
$ for location in WP CA SB SQ;
		do grep ${location} aa2.names.txt > aa2-${location}.names.txt;
		vcftools --vcf ac_1div_nc_20.vcf --keep aa2-${location}.names.txt --min-alleles 2 --max-alleles 2 --max-missing 0.2 --maf 0.000000001 --stdout --recode > aa2-${location}_1div_nc_20.vcf;
		tail -n 4;
		done
### WP
#After filtering, kept 144 out of 335 Individuals
#After filtering, kept 945 out of a possible 1606 Sites
### CA
#After filtering, kept 53 out of 335 Individuals
#After filtering, kept 887 out of a possible 1606 Sites
#### SB
#After filtering, kept 71 out of 335 Individuals
#After filtering, kept 1005 out of a possible 1606 Sites
#### SQ
#After filtering, kept 30 out of 335 Individuals
#After filtering, kept 1052 out of a possible 1606 Sites

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
```

### Redundancy analysis

```bash
$ for taxa in AA1 AA2 AH1 AH2 AH3 AL1 AL2;
		do Rscript rda.R $taxa;
		done
```

### Isolation-by-distance

```bash
# Rousset and genepop
# NOTE: the AA2 analysis was run on HPC due to memory issues
## Running genepop across all locations
$ for taxa in AA1 AA2	AH1 AH2 AH3 AL1 AL2;
		do Rscript genepop.R $taxa all all;
		done

# Running genepop across all locations, for within 'locations', e.g., all depths at one spatial location 
$ for taxa in AA1 AA2 AH1 AH2 AH3 AL1 AL2;
		do Rscript quick_genepop.R $taxa all within;
		done

# Running genepop across all locations, for between 'locations', e.g., min scale 100 
$ for taxa in AA1 AA2 AH1 AH2 AH3 AL1 AL2;
		do Rscript quick_genepop.R $taxa all between;
		done

# create distance matrices for â and ln(geodist)
$ for taxa in AA1 AA2 AH1 AH2 AH3 AL1 Al2;
		do  ./get_dist_matrix.sh ${taxa} all within 2D;
		done
		
# create distance matrices for â and geodist (?)
$ for taxa in AA1 AA2 AH1 AH2 AH3 AL1 Al2;
		do ./get_dist_matrix.sh $taxa all between 1D;
		done
# Note: may need to open with excel to fix up, TODO: edit script so all lines that start with tab or space the tab or space is removed.

# Summarise results
$ for taxa in AA1 AA2 AH1 AH2 AH3 AL1 AL2;
		do ./get_genepop_results.sh $taxa all between 1D;
		done
$ for taxa in AA1 AA2 AH1 AH2 AH3 AL1 AL2;
		do ./get_genepop_results.sh $taxa all within 2D;
		done

# Make plots
$ for taxa in AA1 AA2 AH1 AH2 AH3 AL1 AL2;
		do Rscript ibd_plots.R $taxa all all;
		done
$ Rscript ibd_plots.R AA1 all within
```

**Using loiselle's kinship and SPAGedI v1.5**

```bash
# Loiselle and spagedi
# NOTE: the aa2 analysis was run on HPC due to memory issues
## Convert to spagedi format
$ for taxa in aa1 aa2 ah1 ah2 ah3 al1 al2;
		do Rscript vcf2spagedi.R $taxa;
		done
## Running spagedi and summarising results
$ for taxa in aa1 ah1 ah2 ah3;
		for i in 1 2;
			do spagedi < cmds/${taxa}.spagedi_cmds${i}.txt;
			mv ../data/${taxa}.spagedi-results${i}.txt ../results/ibd/${taxa}/;
			Rscript get_spagedi_results.R ${taxa} ${i};
		done 
$ for i in 1 2;
			Rscript get_spagedi_results.R aa2 ${i};
			done
$ for taxa in al1 al2;
		do spagedi < cmds/${taxa}.spagedi_cmds1.txt;
			mv ../data/${taxa}.spagedi-results1.txt ../results/ibd/${taxa}/;
			Rscript get_spagedi_results.R ${taxa} 1;
		done
```

### F-statistics

```bash
$ for taxa in aa1 aa2 ah1 ah2 ah3 al1 al2;
		do Rscript popgen_stats-taxa.R $taxa 20;
		done
```

### Kinship

Analyses were done using parallel processing on HPC, official results for each datasets can be found in `results/kinship/`. The [COLONY](https://www.zsl.org/about-zsl/resources/software/colony) software (v2.0.6.8) needs to be downloaded to use `./colony2s.out`.

**Analysis summary:**

| Taxon                                    | # individuals | # SNP loci | # Kin | Prob | *Ne* | FIS  |
| ---------------------------------------- | ------------- | ---------- | ----- | ---- | ---- | ---- |
| aa1-wp - *A. agarictes* 1 at West Point  | 20            | 361        | 7     | 0.2  | 99   | 0    |
| aa1-sq - *A. agaricites* 1 at Seaquarium | 7             | 314        | 0     | 0.2  | -    | 0    |
| aa2-wp - *A. agaricties* 2 at West Point | 144           | 938        |       |      |      |      |
| aa2-ca - *A. agaricites* 2 at Cas Abao   | 53            | 856        |       |      |      |      |
| aa2-sb - *A. agaricites* 2 at Snake Bay  | 71            | 977        |       |      |      |      |
| aa2-sq - *A. agaricites* 2 at Seaquarium | 30            | 1044       |       |      |      |      |
| ah1-wp - *A. humilis* 1 at West Point    | 29            | 644        | 0     | 0.2  | -    | 0.57 |
| ah1-ca - *A. humilis* 1 at Cas Abao      | 23            | 581        | 2     | 0.2  | 316  | 0.58 |
| ah2-wp - *A. humilis* 2 at West Point    | 17            | 466        |       |      |      |      |
| ah2-sb - *A. humilis* 2 at Snake Bay     | 12            | 486        |       |      |      |      |
| ah3-wp - *A. humilis* 3 at West Point    | 11            | 526        |       |      |      |      |
| al1 - *A. lamarck*i 1                    | 28            | 544        |       |      |      |      |
| al2 - *A. lamarcki* 2                    | 62            | 734        |       |      |      |      |

Kinship general settings file:

```
! 0 - not updating allele frequency
! 1 - Monoecious species
! 1 - Inbreeding
! 0 - Diploid species
! 0 / 0 - Polygamy for both MF
! 0 - No clone inference
! 0 1.0 1.0 - No sibship prior; mean parternal maternal size - not necessary if sibship prior is 0
! 0 - Unknown pop freq
! 10 - Number of runs
! 4 - Length of run
! 0 - Monitor method always choose 0 when not using windows GUI
! 10000 - Monitor interval (ignore)
! 0 - non-windows version
! 1 - Full likelihood
! 3 - High precision
```

Code:

```bash
# Prepare genotypes and loci files for COLONY input

$ for dataset in aa1 ah1 ah2 ah3 al1 al2 aa2-WP aa2-CA aa2-SB aa2-SQ;
		do Rscript colony_files.R ${dataset}_1div_nc_20;
		done
# EXAMPLE for one dataset
# Create COLONY input files
$ for param in 0.5 0.2 0.1 0.05;
		do ./edit_dat_file.sh aa1_1div_nc_20 35 465 $param;
		done
# Running analyses
$ for param in 0.5 0.1 0.05;
		do ./colony2s.out IFN:aa1_1div_nc_20_${param}.dat;
		done
```

### Clone distances

```bash
# clone distances
$ Rscript clone_spatial_plots.R
# all genotype distances
$ for taxa in AA1 AA2 AH1 AH2 AH3 AL1 AL2
		do Rscript genotype_spatial_plots.R $taxa;
		done
```

## Appendix

### Python packages

Details for different conda environments.

`radkat`

```bash
# packages in environment at /Users/kprata/anaconda3/envs/radkat:
#
# Name                    Version                   Build  Channel
admixture                 1.3.0                         0    bioconda
appnope                   0.1.0                 py37_1000    conda-forge
attrs                     19.3.0                     py_0    conda-forge
backcall                  0.1.0                      py_0    conda-forge
bcftools                  1.9                  h16e57c4_7    bioconda
biopython                 1.74             py37h01d97ff_0    conda-forge
blas                      2.12                   openblas    conda-forge
blast                     2.12.0               h0370960_3    bioconda
bleach                    3.1.0                      py_0    conda-forge
boost                     1.68.0          py37h9888f84_1001    conda-forge
boost-cpp                 1.68.0            h6f8c590_1000    conda-forge
bzip2                     1.0.8                h0b31af3_2    conda-forge
c-ares                    1.18.1               h0d85af4_0    conda-forge
ca-certificates           2021.10.8            h033912b_0    conda-forge
certifi                   2021.10.8        py37hf985489_1    conda-forge
curl                      7.82.0               h9dce1e3_0    conda-forge
cython                    0.29.14          py37h4a8c4bd_0    conda-forge
dadi                      2.0.3                    pypi_0    pypi
dbus                      1.13.6               h2f22bb5_0    conda-forge
decorator                 4.4.1                      py_0    conda-forge
defusedxml                0.6.0                      py_0    conda-forge
entrez-direct             16.2                 h193322a_0    bioconda
entrypoints               0.3                   py37_1000    conda-forge
expat                     2.2.9                h4a8c4bd_2    conda-forge
gettext                   0.19.8.1          hd1a6beb_1008    conda-forge
glib                      2.70.2               hcf210ce_4    conda-forge
glib-tools                2.70.2               hcf210ce_4    conda-forge
gmp                       6.2.0                h4a8c4bd_1    conda-forge
gnutls                    3.6.5             h53004b3_1002    conda-forge
gsl                       2.5                  ha2d443c_1    conda-forge
htslib                    1.9                  h356306b_9    bioconda
icu                       58.2              h0a44026_1000    conda-forge
importlib_metadata        1.5.0                    py37_0    conda-forge
inflect                   4.1.0                    py37_0    conda-forge
ipykernel                 5.1.4            py37h5ca1d4c_0    conda-forge
ipython                   7.8.0            py37h5ca1d4c_0    conda-forge
ipython_genutils          0.2.0                      py_1    conda-forge
ipywidgets                7.5.1                      py_0    conda-forge
jaraco.itertools          5.0.0                      py_0    conda-forge
jedi                      0.16.0                   py37_0    conda-forge
jinja2                    2.11.1                     py_0    conda-forge
jpeg                      9c                h1de35cc_1001    conda-forge
jsonschema                3.2.0                    py37_0    conda-forge
jupyter                   1.0.0                      py_2    conda-forge
jupyter_client            5.3.4                    py37_1    conda-forge
jupyter_console           6.0.0                      py_0    conda-forge
jupyter_core              4.6.1                    py37_0    conda-forge
krb5                      1.19.3               hb98e516_0    conda-forge
libblas                   3.8.0               12_openblas    conda-forge
libcblas                  3.8.0               12_openblas    conda-forge
libcurl                   7.82.0               h9dce1e3_0    conda-forge
libcxx                    13.0.1               hc203e6f_0    conda-forge
libdeflate                1.3                  h01d97ff_0    conda-forge
libedit                   3.1.20191231         h0678c8f_2    conda-forge
libev                     4.33                 haf1e3a3_1    conda-forge
libffi                    3.4.2                h0d85af4_5    conda-forge
libgfortran               3.0.1                         0    conda-forge
libglib                   2.70.2               hf1fb8c0_4    conda-forge
libiconv                  1.16                 haf1e3a3_0    conda-forge
liblapack                 3.8.0               12_openblas    conda-forge
liblapacke                3.8.0               12_openblas    conda-forge
libnghttp2                1.47.0               hca56917_0    conda-forge
libopenblas               0.3.7                hd44dcd8_1    conda-forge
libpng                    1.6.37               h2573ce8_0    conda-forge
libsodium                 1.0.17               h01d97ff_0    conda-forge
libssh2                   1.10.0               hd3787cc_2    conda-forge
libzlib                   1.2.11            h9173be1_1013    conda-forge
llvm-openmp               13.0.1               hcb1a161_1    conda-forge
markupsafe                1.1.1            py37h0b31af3_0    conda-forge
mistune                   0.8.4           py37h0b31af3_1000    conda-forge
more-itertools            8.2.0                      py_0    conda-forge
mpi                       1.0                     openmpi    conda-forge
nbconvert                 5.6.1                    py37_0    conda-forge
nbformat                  5.0.4                      py_0    conda-forge
ncurses                   6.3                  he49afe7_0    conda-forge
nettle                    3.4.1             h3efe00b_1002    conda-forge
notebook                  6.0.1                    py37_0    conda-forge
numpy                     1.18.0.dev0+463acda          pypi_0    pypi
openblas                  0.3.7                hd44dcd8_1    conda-forge
openjdk                   11.0.1            hbbe82c9_1018    conda-forge
openmpi                   3.1.4                ha90c164_0    conda-forge
openssl                   3.0.0                h0d85af4_2    conda-forge
pandoc                    2.9.1.1                       0    conda-forge
pandocfilters             1.4.2                      py_1    conda-forge
parso                     0.6.1                      py_0    conda-forge
pcre                      8.45                 he49afe7_0    conda-forge
perl                      5.32.1          2_h0d85af4_perl5    conda-forge
perl-archive-tar          2.40            pl5321hdfd78af_0    bioconda
perl-carp                 1.50            pl5321hd8ed1ab_0    conda-forge
perl-common-sense         3.75            pl5321hdfd78af_0    bioconda
perl-compress-raw-bzip2   2.101           pl5321h9722bc1_1    bioconda
perl-compress-raw-zlib    2.101           pl5321h9722bc1_2    bioconda
perl-encode               3.16            pl5321ha5712d3_1    bioconda
perl-exporter             5.74            pl5321hd8ed1ab_0    conda-forge
perl-exporter-tiny        1.002002        pl5321hdfd78af_0    bioconda
perl-extutils-makemaker   7.64            pl5321hd8ed1ab_0    conda-forge
perl-io-compress          2.102           pl5321h9722bc1_1    bioconda
perl-io-zlib              1.11            pl5321hdfd78af_0    bioconda
perl-json                 4.05            pl5321hdfd78af_0    bioconda
perl-json-xs              2.34            pl5321hcd10b59_5    bioconda
perl-list-moreutils       0.430           pl5321hdfd78af_0    bioconda
perl-list-moreutils-xs    0.430           pl5321ha5712d3_1    bioconda
perl-parent               0.238           pl5321hd8ed1ab_0    conda-forge
perl-pathtools            3.75            pl5321ha5712d3_3    bioconda
perl-scalar-list-utils    1.62            pl5321ha5712d3_0    bioconda
perl-types-serialiser     1.01            pl5321hdfd78af_0    bioconda
perl-xsloader             0.24            pl5321hd8ed1ab_0    conda-forge
pexpect                   4.8.0                    py37_0    conda-forge
pickleshare               0.7.5                 py37_1000    conda-forge
pip                       20.0.2                     py_2    conda-forge
plink                     1.90b6.21            hb4d813b_1    bioconda
prometheus_client         0.7.1                      py_0    conda-forge
prompt_toolkit            2.0.10                     py_0    conda-forge
ptyprocess                0.6.0                   py_1001    conda-forge
pygments                  2.5.2                      py_0    conda-forge
pyqt                      5.9.2            py37h2a560b1_4    conda-forge
pyrsistent                0.15.7           py37h0b31af3_0    conda-forge
python                    3.7.12          hf3644f1_100_cpython    conda-forge
python-dateutil           2.8.1                      py_0    conda-forge
python_abi                3.7                     1_cp37m    conda-forge
pyvcf                     0.6.8                 py37_1000    conda-forge
pyzmq                     18.1.1           py37h4bf09a9_0    conda-forge
qt                        5.9.7                h93ee506_2    conda-forge
qtconsole                 4.6.0                      py_0    conda-forge
readline                  8.1                  h05e3726_0    conda-forge
send2trash                1.5.0                      py_0    conda-forge
setuptools                45.1.0                   py37_0    conda-forge
sip                       4.19.8          py37h0a44026_1000    conda-forge
six                       1.14.0                   py37_0    conda-forge
snpsift                   4.3.1t                        1    bioconda
sqlite                    3.37.1               hb516253_0    conda-forge
structure                 2.3.4                h470a237_1    bioconda
terminado                 0.8.3                    py37_0    conda-forge
testpath                  0.4.4                      py_0    conda-forge
tk                        8.6.12               h5dbffcc_0    conda-forge
tornado                   6.0.3            py37h0b31af3_0    conda-forge
traitlets                 4.3.3                    py37_0    conda-forge
vcftools                  0.1.16               h5c9b4e4_3    bioconda
wcwidth                   0.1.8                      py_0    conda-forge
webencodings              0.5.1                      py_1    conda-forge
wheel                     0.34.2                     py_1    conda-forge
widgetsnbextension        3.5.1                    py37_0    conda-forge
xz                        5.2.5                haf1e3a3_1    conda-forge
zeromq                    4.3.2                h6de7cb9_2    conda-forge
zipp                      2.1.0                      py_0    conda-forge
zlib                      1.2.11            h9173be1_1013    conda-forge
```

`open3d`

```
# packages in environment at /Users/kprata/anaconda3/envs/open3d:
#
# Name                    Version                   Build  Channel
addict                    2.4.0                    pypi_0    pypi
alphashape                1.3.1              pyh44b312d_0    conda-forge
aom                       3.3.0                h96cf925_1    conda-forge
attrs                     21.4.0             pyhd8ed1ab_0    conda-forge
blosc                     1.21.0               he49afe7_0    conda-forge
boost-cpp                 1.74.0               hdbf7018_7    conda-forge
brotli                    1.0.9                h0d85af4_6    conda-forge
brotli-bin                1.0.9                h0d85af4_6    conda-forge
bzip2                     1.0.8                h0d85af4_4    conda-forge
c-ares                    1.18.1               h0d85af4_0    conda-forge
ca-certificates           2021.10.8            h033912b_0    conda-forge
cairo                     1.16.0            he01c77b_1009    conda-forge
certifi                   2021.10.8        py39h6e9494a_2    conda-forge
cfitsio                   4.0.0                hb20e66c_0    conda-forge
click                     8.0.3            py39h6e9494a_1    conda-forge
click-log                 0.3.2              pyh9f0ad1d_0    conda-forge
click-plugins             1.1.1                      py_0    conda-forge
cligj                     0.7.2              pyhd8ed1ab_1    conda-forge
curl                      7.81.0               h97da3c1_0    conda-forge
cycler                    0.11.0             pyhd8ed1ab_0    conda-forge
descartes                 1.1.0                      py_4    conda-forge
expat                     2.4.4                he49afe7_0    conda-forge
ffmpeg                    5.0.1                h9220da4_0    conda-forge
fiona                     1.8.21           py39haa9df5e_0    conda-forge
font-ttf-dejavu-sans-mono 2.37                 hab24e00_0    conda-forge
font-ttf-inconsolata      3.000                h77eed37_0    conda-forge
font-ttf-source-code-pro  2.038                h77eed37_0    conda-forge
font-ttf-ubuntu           0.83                 hab24e00_0    conda-forge
fontconfig                2.13.94              h10f422b_0    conda-forge
fonts-conda-ecosystem     1                             0    conda-forge
fonts-conda-forge         1                             0    conda-forge
fonttools                 4.28.4                   pypi_0    pypi
freetype                  2.10.4               h4cff582_1    conda-forge
freexl                    1.0.6                h0d85af4_0    conda-forge
fribidi                   1.0.10               hbcb3906_0    conda-forge
future                    0.18.2           py39h6e9494a_5    conda-forge
gdal                      3.4.1            py39h478c985_3    conda-forge
geopandas                 0.9.0              pyhd8ed1ab_0    conda-forge
geos                      3.10.2               he49afe7_0    conda-forge
geotiff                   1.7.0                h0ca5f94_6    conda-forge
gettext                   0.19.8.1          hd1a6beb_1008    conda-forge
giflib                    5.2.1                hbcb3906_2    conda-forge
gmp                       6.2.1                h2e338ed_0    conda-forge
gnutls                    3.6.13               h756fd2b_1    conda-forge
hdf4                      4.2.15               hefd3b78_3    conda-forge
hdf5                      1.12.1          nompi_hd9e8a45_103    conda-forge
icu                       69.1                 he49afe7_0    conda-forge
jbig                      2.1               h0d85af4_2003    conda-forge
joblib                    1.1.0                    pypi_0    pypi
jpeg                      9e                   h0d85af4_0    conda-forge
json-c                    0.15                 hcb556a6_0    conda-forge
kealib                    1.4.14               ha22a8b1_3    conda-forge
kiwisolver                1.3.2            py39hf018cea_1    conda-forge
krb5                      1.19.2               h289aae4_3    conda-forge
lame                      3.100             h35c211d_1001    conda-forge
lcms2                     2.12                 h577c468_0    conda-forge
lerc                      3.0                  he49afe7_0    conda-forge
libblas                   3.9.0           13_osx64_openblas    conda-forge
libbrotlicommon           1.0.9                h0d85af4_6    conda-forge
libbrotlidec              1.0.9                h0d85af4_6    conda-forge
libbrotlienc              1.0.9                h0d85af4_6    conda-forge
libcblas                  3.9.0           13_osx64_openblas    conda-forge
libcurl                   7.81.0               h97da3c1_0    conda-forge
libcxx                    12.0.1               habf9029_1    conda-forge
libdap4                   3.20.6               h3e144a0_2    conda-forge
libdeflate                1.10                 h0d85af4_0    conda-forge
libedit                   3.1.20191231         h0678c8f_2    conda-forge
libev                     4.33                 haf1e3a3_1    conda-forge
libffi                    3.4.2                h0d85af4_5    conda-forge
libgdal                   3.4.1                h467bfbe_3    conda-forge
libgfortran               5.0.0           9_3_0_h6c81a4c_23    conda-forge
libgfortran5              9.3.0               h6c81a4c_23    conda-forge
libglib                   2.70.2               hf1fb8c0_4    conda-forge
libiconv                  1.16                 haf1e3a3_0    conda-forge
libimagequant             2.17.0               h9bde063_1    conda-forge
libkml                    1.3.0             h8fd9edb_1014    conda-forge
liblapack                 3.9.0           13_osx64_openblas    conda-forge
libnetcdf                 4.8.1           nompi_h6609ca0_101    conda-forge
libnghttp2                1.46.0               hfd382f3_0    conda-forge
libopenblas               0.3.18          openmp_h3351f45_0    conda-forge
libpng                    1.6.37               h7cec526_2    conda-forge
libpq                     14.2                 h8ce7ee7_0    conda-forge
librttopo                 1.1.0                hec60dd8_9    conda-forge
libspatialindex           1.9.3                he49afe7_4    conda-forge
libspatialite             5.0.1               h88e7940_14    conda-forge
libssh2                   1.10.0               hd3787cc_2    conda-forge
libtiff                   4.3.0                h17f2ce3_3    conda-forge
libvpx                    1.11.0               he49afe7_3    conda-forge
libwebp                   1.2.2                h28dabe5_0    conda-forge
libwebp-base              1.2.2                h0d85af4_1    conda-forge
libxcb                    1.13              h0d85af4_1004    conda-forge
libxml2                   2.9.12               h7e28ab6_1    conda-forge
libzip                    1.8.0                h7e5727d_1    conda-forge
libzlib                   1.2.11            h9173be1_1013    conda-forge
llvm-openmp               13.0.1               hda6cdc1_0    conda-forge
lz4-c                     1.9.3                he49afe7_1    conda-forge
matplotlib-base           3.5.1            py39hb07454d_0    conda-forge
munch                     2.5.0                      py_0    conda-forge
munkres                   1.1.4              pyh9f0ad1d_0    conda-forge
ncurses                   6.2                  h2e338ed_4    conda-forge
nettle                    3.6                  hedd7734_0    conda-forge
networkx                  2.6.3              pyhd8ed1ab_1    conda-forge
nspr                      4.32                 hcd9eead_1    conda-forge
nss                       3.74                 h31e2bf1_0    conda-forge
numpy                     1.21.4                   pypi_0    pypi
open3d                    0.14.1                   pypi_0    pypi
openh264                  2.1.1                hfd3ada9_0    conda-forge
openjpeg                  2.4.0                h6e7aa92_1    conda-forge
openssl                   3.0.2                h6c3fc93_1    conda-forge
packaging                 21.3               pyhd8ed1ab_0    conda-forge
pandas                    1.3.5                    pypi_0    pypi
pcre                      8.45                 he49afe7_0    conda-forge
pillow                    8.4.0                    pypi_0    pypi
pip                       21.3.1             pyhd8ed1ab_0    conda-forge
pixman                    0.40.0               hbcb3906_0    conda-forge
poppler                   22.01.0              h9573804_0    conda-forge
poppler-data              0.4.11               hd8ed1ab_0    conda-forge
postgresql                14.2                 h0fd25fa_0    conda-forge
proj                      8.2.1                h1512c50_0    conda-forge
pthread-stubs             0.4               hc929b4f_1001    conda-forge
pyglet                    1.5.16           py39h6e9494a_1    conda-forge
pyparsing                 3.0.6                    pypi_0    pypi
pyproj                    3.3.0            py39h074cefc_1    conda-forge
python                    3.9.7           h38b4d05_3_cpython    conda-forge
python-dateutil           2.8.2              pyhd8ed1ab_0    conda-forge
python_abi                3.9                      2_cp39    conda-forge
pytz                      2021.3             pyhd8ed1ab_0    conda-forge
pyyaml                    6.0                      pypi_0    pypi
readline                  8.1                  h05e3726_0    conda-forge
rtree                     0.9.7            py39h7d0d40a_3    conda-forge
scikit-learn              1.0.1                    pypi_0    pypi
scipy                     1.7.3                    pypi_0    pypi
setuptools                59.6.0           py39h6e9494a_0    conda-forge
shapely                   1.8.0            py39hbfbc381_5    conda-forge
six                       1.16.0             pyh6c4a22f_0    conda-forge
sqlite                    3.37.0               h23a322b_0    conda-forge
svt-av1                   0.9.1                h96cf925_0    conda-forge
threadpoolctl             3.0.0                    pypi_0    pypi
tiledb                    2.6.2                h2a16ea5_1    conda-forge
tk                        8.6.11               h5dbffcc_1    conda-forge
tqdm                      4.62.3                   pypi_0    pypi
trimesh                   3.10.0             pyh6c4a22f_0    conda-forge
tzcode                    2021e                h0d85af4_0    conda-forge
tzdata                    2021e                he74cb21_0    conda-forge
unicodedata2              14.0.0           py39h89e85a6_0    conda-forge
wheel                     0.37.0             pyhd8ed1ab_1    conda-forge
x264                      1!161.3030           h0d85af4_1    conda-forge
x265                      3.5                  hbb4e6a2_3    conda-forge
xerces-c                  3.2.3                h6564042_4    conda-forge
xorg-libxau               1.0.9                h35c211d_0    conda-forge
xorg-libxdmcp             1.1.3                h35c211d_0    conda-forge
xz                        5.2.5                haf1e3a3_1    conda-forge
zlib                      1.2.11            h9173be1_1013    conda-forge
zstd                      1.5.2                h582d3a0_0    conda-forge
```

###  R packages

| Package      | Version | DOI or package link                                          |
| ------------ | ------- | ------------------------------------------------------------ |
| adegenet     | 2.1.7   | 10.1093/bioinformatics/btn129,  10.1093/bioinformatics/btr521 |
| ape          | 5.6-2   | [https://CRAN.R-project.org/package=ape](https://cran.r-project.org/package=ape) |
| corrplot     | 0.92    | https://github.com/taiyun/corrplot                           |
| dplyr        | 1.0.9   | [https://CRAN.R-project.org/package=dplyr](https://cran.r-project.org/package=dplyr) |
| genepop      | 1.1.7   | http://dx.doi.org/10.1111/j.1471-8286.2007.01931.x           |
| ggplot2      | 3.4.0   | [https://ggplot2.tidyverse.org](https://ggplot2.tidyverse.org/) |
| ggplotify    | 0.1.0   | [https://CRAN.R-project.org/package=ggplotify](https://cran.r-project.org/package=ggplotify) |
| ggstance     | 0.3.5   | [https://CRAN.R-project.org/package=ggstance](https://cran.r-project.org/package=ggstance) |
| ggtree       | 3.4.0   | 10.1111/2041-210X.12628                                      |
| hierfstat    | 0.5-11  | [https://CRAN.R-project.org/package=hierfstat](https://cran.r-project.org/package=hierfstat) |
| lattice      | 0.20-45 | [http://lmdvr.r-forge.r-project.org](http://lmdvr.r-forge.r-project.org/) |
| pcadapt      | 4.3.3   | https://doi.org/10.1093/molbev/msaa053                       |
| pegas        | 1.1     | [https://CRAN.R-project.org/package=pegas](https://cran.r-project.org/package=pegas) |
| permute      | 0.9-7   | [https://CRAN.R-project.org/package=permute](https://cran.r-project.org/package=permute) |
| phytools     | 1.0-3   | 10.1111/j.2041-210X.2011.00169.x                             |
| qvalue       | 2.28.0  | http://github.com/jdstorey/qvalue                            |
| RColorBrewer | 1.1-3   | [https://CRAN.R-project.org/package=RColorBrewer](https://cran.r-project.org/package=RColorBrewer) |
| readxl       | 1.4.0   | [https://CRAN.R-project.org/package=readxl](https://cran.r-project.org/package=readxl) |
| reshape2     | 1.4.4   | http://www.jstatsoft.org/v21/i12/                            |
| sjmisc       | 2.8.9   | 10.21105/joss.00754                                          |
| spaa         | 0.2.2   | [https://CRAN.R-project.org/package=spaa](https://cran.r-project.org/package=spaa) |
| stringr      | 1.4.0   | [https://CRAN.R-project.org/package=stringr](https://cran.r-project.org/package=stringr) |
| tidyr        | 1.2.0   | [https://CRAN.R-project.org/package=tidyr](https://cran.r-project.org/package=tidyr) |
| tidyverse    | 1.3.1   | 10.21105/joss.01686                                          |
| treeio       | 1.20.0  | 0.1093/molbev/msz240                                         |
| vcfR         | 1.12.0  | [http://dx.doi.org/10.1111/1755-0998.12549 ](http://dx.doi.org/10.1111/1755-0998.12549) |
| vegan        | 2.6-2   | [https://CRAN.R-project.org/package=vegan](https://cran.r-project.org/package=vegan) |


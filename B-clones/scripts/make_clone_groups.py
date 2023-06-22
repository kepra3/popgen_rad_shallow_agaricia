import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as ssd
from sklearn.cluster import AgglomerativeClustering
import argparse

# updated: 16/6/23
# Information:
# script runs from script directory
# takes one argument
# e.g, ac_1d_wc_20, lm_1d-pure_wc, hu_1d_wc

# Arguments
#parser = argparse.ArgumentParser(prog="Make numbered clone groups", usage="[options]")
#parser.add_argument("name")
#args = parser.parse_args()
#name = args.name
name = "ac_1d_wc"

# Conditions
if name == 'lm_1d-pure_wc':
    colour = "purple"
    clone_threshold = 0.98
    taxa_threshold = 0.946
    y_height = 3000
elif name == 'ac_1d_wc':
    colour = "green"
    clone_threshold = 0.98
    taxa_threshold = 0.945
    y_height = 70000
elif name == 'hu_1d_wc':
    colour = "brown"
    clone_threshold = 0.99
    taxa_threshold = 0.972
    taxa_threshold_2 = 0.963
    y_height = 3000
else:
    clone_threshold = 0.99
    taxa_threshold = 0.94
    y_height = 3000
    colour = "red"

# import data set (pairwise distance matrix)
gd_mat = pd.read_csv('../data/{}_20_gdmatrix.tsv'.format(name), index_col=0, sep='\t')
euc_mat = gd_mat.stack().reset_index().rename(columns={'level_0': 'Obj1', 'level_1': 'Obj2', 0: 'Dist'})

plt.hist(1 - euc_mat["Dist"], bins=50, color=colour, range=[0.9, 1])
plt.xlabel('Percentage of genetic similarity between shared SNPs (%)')
plt.ylabel('Number of pairwise comparisons')
plt.vlines(x=clone_threshold, ymax=y_height, ymin=0, linestyles="dashed", colors="red")
plt.text(x=clone_threshold - 0.003, y=y_height / 2, size=10, s="Clonal threshold", rotation='vertical')
plt.vlines(x=taxa_threshold, ymax=y_height, ymin=0, linestyles="dashed", colors="blue")
plt.text(x=taxa_threshold - 0.003, y=y_height / 2, size=10, s="Taxa threshold", rotation='vertical')
if name != 'hu_1d_wc':
    print('')
else:
    plt.vlines(x=taxa_threshold_2, ymax=y_height, ymin=0, linestyles="dashed", colors="green")
    plt.text(x=taxa_threshold_2 - 0.003, y=y_height / 2, size=10, s="Taxa threshold - 2", rotation='vertical')
# plt.title('Distribution of pairwise genetic similarity')


plt.savefig('../results/distribution_pairwise_sim_{}.pdf'.format(name), dpi=400)

gd_mat_cond = ssd.squareform(gd_mat)
gd_link = sch.linkage(gd_mat_cond)
gd_dendro = sch.dendrogram(gd_link)

# Grouping clusters by number of clusters
hc = AgglomerativeClustering(n_clusters=365, affinity = 'euclidean', linkage = 'ward')
y_hc = hc.fit_predict(gd_link)
print(y_hc)

# 5% clonal threshold
# hu - 0.005
# ac - 0.0025 clonal threshold not high enough
# lm - 0.005
groups = pd.DataFrame(sch.fcluster(gd_link, 0.005, criterion='distance'))

# Sort clone file with individual popfile
all_indiv = pd.read_csv('../data/pop_{}.txt'.format(name), sep='\t', header=None)
all_indiv['Groups'] = groups
all_indiv.columns = ['Sample', 'Pop', 'Groups']
all_indiv.to_csv('./results/clone_groups_{}.csv'.format(name), index=False)

# retain individuals from clone groups with highest number of sites
clones = pd.read_csv('../results/clone_matches_{}.csv'.format(name))

all_indiv['Sites'] = 0
# loop through whole data set for sample 1
for i in range(len(all_indiv['Sample'].values)):
    label = all_indiv['Sample'].values[i]
    ssites = clones['ind1_snps'].values[clones['# ind1'].values == label]
    print(label, ssites, i)
    if len(ssites) == 0:
        all_indiv['Sites'][i] = None
    elif len(ssites) == 1:
        all_indiv['Sites'][i] = ssites
    elif len(ssites) > 1:
        all_indiv['Sites'][i] = max(ssites)

# loop through whole data set for sample 2
for i in range(len(all_indiv['Sample'].values)):
    label = all_indiv['Sample'].values[i]
    ssites = clones['ind2_snps'].values[clones['ind2'].values == label]
    print(label, ssites, i)
    if len(ssites) == 0:
        pass
    elif len(ssites) == 1:
        all_indiv['Sites'][i] = ssites
    elif len(ssites) > 1:
        all_indiv['Sites'][i] = max(ssites)
# individuals with site values should be within clone groups

# For all groups
chosen_individuals = pd.DataFrame()
for group_int in all_indiv['Groups'].drop_duplicates():
    group = all_indiv[all_indiv['Groups'] == group_int]
    group = group[group['Sites'] == max(group['Sites'])]
    chosen_individuals = chosen_individuals.append(group, ignore_index=True)
# yay works :D

chosen_individuals.to_csv('../data/chosen_individuals_{}.csv'.format(name), index=False)
individuals2keep = chosen_individuals['Sample']

individuals2keep.to_csv('../data/individuals2keep_{}.txt'.format(name), index=False)

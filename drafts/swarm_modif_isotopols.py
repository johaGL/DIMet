# Credits : Macha Nikolski

import csv
import numpy as np
import pandas as pd
import argparse
import itertools
import scipy.stats as sta
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
import math
from sklearn.decomposition import PCA
from scipy.cluster import hierarchy as hac
from scipy.cluster.hierarchy import fcluster

from locale import *

parser = argparse.ArgumentParser()
# read the file
parser.add_argument("--file", "-f", help='Input file', required=True)
parser.add_argument("--impute", action=argparse.BooleanOptionalAction, default=True, help='Impute or not impute')
parser.add_argument("--xlab", help="isotopologue m+0 | isotopologue m+1 | ...")
parser.add_argument("--measure", help="abundance | isotopol_abund | mean_enrichment | ...")
parser.add_argument("--outdir", help="...")

# python ./foo.py --file meanEnrich_myelin_cell_4h.tsv

args = parser.parse_args()
data = args.file
impute = args.impute

print("Reading cells data...")
Metabolites = pd.read_csv(data, sep="\t", decimal=".", index_col=0)  # header=0)


print("Read data:", Metabolites.shape)

column_groups = defaultdict(list)
for sample in Metabolites.columns.values:
    type = sample.split('_')[0]
    column_groups[type].append(sample)

print(column_groups)


def compute_PCA(df_w_labels, df, file_name, targets, colors, select_on='cat'):
    print("Computing PCA")
    print(df_w_labels)
    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(df)
    principalDf = pd.DataFrame(data=principalComponents, columns=['principal component 1', 'principal component 2'])
    print("   Explained variance", pca.explained_variance_ratio_)

    finalDf = pd.concat([principalDf, df_w_labels[select_on].rename('labels')], axis=1)
    print(finalDf)

    plt.figure(figsize=(8, 8))
    ax = plt.axes()
    ax.set_xlabel('PC 1', fontsize=15)
    ax.set_ylabel('PC 2', fontsize=15)
    ax.set_title('Explained variance ' + '%.2f' % pca.explained_variance_ratio_[0] + ' and ' + '%.2f' %
                 pca.explained_variance_ratio_[1], fontsize=20)

    for target, color in zip(targets, colors):
        print(target, color)
        indicesToKeep = finalDf[finalDf['labels'] == target].index
        ax.scatter(finalDf.loc[indicesToKeep, 'principal component 1']
                   , finalDf.loc[indicesToKeep, 'principal component 2']
                   , alpha=0.8
                   , c=color
                   , s=35)
        ax.legend(targets)
        ax.grid()

    plt.savefig(file_name)
    plt.close()


def swarmplot(expr, data, myxlab, fig_name):
    print("\n*** Drawing swarmplot, splitting by group")
    measure_data = expr.copy()
    group_names = column_groups.keys()

    typeofmeasure =  args.measure# TODO:  transform into an option !

    measure_data = measure_data.assign(group=['' for x in np.arange(len(measure_data.index))])
    measure_data.loc[expr['cat'].str.contains('Vehicle'), 'group'] = 'Vehicle'
    measure_data.loc[expr['cat'].str.contains('Pellets'), 'group'] = 'Pellets'
    measure_data.loc[expr['cat'].str.contains('Myelin'), 'group'] = 'Myelin'
    print(measure_data)

    melted = measure_data.melt(id_vars=["group", "cat"], var_name='metabolite', value_name=typeofmeasure)
    print(melted)

    sample_colors = dict(
        [('Vehicle_cell_T4h-1', 'royalblue'), ('Vehicle_cell_T4h-2', 'blue'), ('Vehicle_cell_T4h-3', 'dodgerblue'),
         ('Pellets_cell_T4h-1', 'green'), ('Pellets_cell_T4h-2', 'darkgreen'), ('Pellets_cell_T4h-3', 'mediumseagreen'),
         ('Myelin_cell_T4h-1', 'orangered'), ('Myelin_cell_T4h-2', 'mediumvioletred'),
         ('Myelin_cell_T4h-3', 'crimson')])

    fig, ax = plt.subplots(figsize=(24, 6))
    melted["metabolite"] = melted["metabolite"].apply(lambda x: x[:10])  # shorten the metabolite names
    sns.swarmplot(x="metabolite", y=typeofmeasure, s=5, hue="cat", palette=sample_colors, data=melted)

    #plt.ylim(-0.5, 25)
    # rotate x-axis labels
    plt.xticks(rotation=45, fontsize=9)

    plt.xlabel(myxlab)
    fig.suptitle(myxlab)

    # Put a legend to the right of the current axis
    plt.legend(bbox_to_anchor=(0, 1.02, 1, 0.2), loc="lower left", borderaxespad=0, ncol=4)
    fig.tight_layout()

    fig.savefig(fig_name)
    plt.close()


# reduce the data
print("\n#######\n")
reduced = pd.DataFrame(columns=Metabolites.columns.values)
count = 0


for m in Metabolites.index.values:
    measures = Metabolites.loc[m].copy()
    print(measures)
    total = sum(measures.values)  # np.count_nonzero(np.isnan(replicate_measures))
    if total == 0:  # skip metabolites with no measure values
        continue

    if impute:
        #        mask = (measures == 0.0)
        num_zeros = (measures == 0).sum()
        if num_zeros > 0:  # impute
            print("Imputing values for", m)
            for grp in column_groups.keys():
                group_measures = measures[column_groups[grp]]
                for sample in column_groups[grp]:
                    if measures[sample] == 0.0:
                        measures.loc[sample] = np.nanmean(group_measures)

    reduced.loc[m, :] = measures / measures.std(ddof=0)

samples = ['Vehicle_cell_T4h-1', 'Vehicle_cell_T4h-2', 'Vehicle_cell_T4h-3',
           'Myelin_cell_T4h-1', 'Myelin_cell_T4h-2', 'Myelin_cell_T4h-3']
measure_data = reduced[samples].copy().T
values_no_labels = reduced.T.values
reduced.loc['cat', :] = reduced.columns.values

all_colors = ['royalblue', 'blue', 'dodgerblue'
              'orangered', 'mediumvioletred', 'crimson']

#compute_PCA(reduced.T.reset_index(), values_no_labels, "all_PCA.png", targets=Metabolites.columns.values,
#            colors=all_colors)

swarmplot(reduced.T.reset_index(drop=True), measure_data.reset_index(), args.xlab, f"{args.outdir}swarmplot_{args.xlab}.png" )

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Metabologram module
johagl 2023
"""
import os
import argparse
import yaml
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize
from matplotlib.colors import LinearSegmentedColormap
import matplotlib as mpolib
import functions_general as fg


def metabologram_args():
    parser = argparse.ArgumentParser(prog="python -m DIMet.src.abundance_bars",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('config', type=str,
                        help="configuration file in absolute path")

    return parser


def bars_args():
    parser = argparse.ArgumentParser(prog="python -m DIMet.src.abundance_bars",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('config', type=str,
                        help="configuration file in absolute path")

    parser.add_argument('--palette', action=argparse.BooleanOptionalAction, default="pastel",
                        help="qualitative or categorical palette name as in Seaborn library (Python)")

    return parser


def read_config(metabologram_config):
    try:
        with open(metabologram_config, "r") as f:
            confgramD = yaml.load(f, Loader=yaml.Loader)
    except Exception as err:
        print("Error when opening metabologram configuration file in {metabologram_dir}")
        print(err)
    return confgramD


def give_path_dico(pathways_files):
    #def openfile(filename):
    def get_as_dico(path_file):
        pathsdf = pd.read_csv(path_file, sep="\t", index_col=None)
        pathsdf = pathsdf.fillna("")
        preD = pathsdf.to_dict(orient='list')
        pathD = dict()
        for k in preD.keys():
            tmpl = [i for i in preD[k] if i != ""]
            deduplicated = set(tmpl)
            pathD[k] = list(deduplicated)
        return pathD

    genes_dico = get_as_dico(pathways_files['genes'])
    metabo_dico = get_as_dico(pathways_files['metabolites'])
    assert set(genes_dico.keys()) == set(metabo_dico.keys()), \
      "Error, pathway names in both pathway files are not identical"
    finalD = genes_dico.copy()
    for k in metabo_dico.keys():
        finalD[k] += metabo_dico[k]

    return genes_dico, metabo_dico, finalD


def deg_dam_corresp( DEG_tables,  DAM_tables, titles):
    assert ( set(DEG_tables.keys()) == set(DAM_tables.keys()) ) and \
    (set(titles.keys()) == set(DEG_tables.keys())),  "error with the numbering of DEG_tables DAM_tables and titles"
    comparisondico = dict()
    for k in titles.keys():
        comparisondico[k] = { 'gene' : DEG_tables[k],
                              'metabolite' : DAM_tables[k],
                              'title': titles[k] }
    return comparisondico


def pile_dfs_byconds(comparisondico, dir_table, columns_pick, subkey_comparisondico):
    df_out = pd.DataFrame({'name': [], 'log2FC': []})
    for comparison in comparisondico.keys():
        file_here = comparisondico[comparison][subkey_comparisondico]
        df = pd.read_csv(dir_table + file_here, sep='\t')
        df = df[columns_pick]
        df.columns = ['name', 'log2FC']
        df.log2FC = df.log2FC.round(2)
        df = df.assign(comparison=comparison, typemol=subkey_comparisondico)
        df_out = pd.concat([df_out, df], axis=0)
    return df_out


def filter_by_path_dico(df, dico):
    created_elems_list = set()
    for k in dico.keys():
        created_elems_list = created_elems_list.union(set(dico[k]))
    df = df.loc[df.name.isin(list(created_elems_list)),:]
    return df



def get_custom_color_palette_hash(lowcol, midcol, highcol):
    """
    courtesy from :
    https://richardhildebrand.wordpress.com/2019/09/18/create-a-custom-color-palette-with-matplotlib-and-seaborn/
    """
    colorlist = [lowcol, midcol, highcol]
    return LinearSegmentedColormap.from_list("", colorlist, N=256)


def rgbas2hex(rgbas_):
    colorsout = []
    for tup in rgbas_:
        tmp = mpolib.colors.to_hex(tup)
        colorsout.append(tmp)
    return colorsout


def values2rgbas(myvalues, mycmap, vmin, vmax, center):
    if center == 0:
        # Normalize data before giving colors, because map interval is [0,1] by matplotlib
        # https://stackoverflow.com/questions/25408393/getting-individual-colors-from-a-color-map-in-matplotlib
        norm = Normalize(vmin=vmin, vmax=vmax)
        rgba_tuples = mycmap(norm(myvalues))
        return rgba_tuples
    else :
        print("only center == 0 is handled here")


def inner_pie_colors(inner_dico, mycmap,  gabs, mabs):
    metaboval = inner_dico['metabo_mean_val']
    geneval = inner_dico['gene_mean_val']
    metabocolors_ = rgbas2hex(values2rgbas([metaboval], mycmap, -mabs, mabs, center=0))
    genecolors_ = rgbas2hex(values2rgbas([geneval], mycmap, -gabs, gabs, center=0))
    return {'metab_color' : metabocolors_[0], 'gene_color': genecolors_[0]}


def metabologram_run(confidic, dimensions_pdf):

    confD = confidic
    genes_dico, metabo_dico, pathdico = give_path_dico(confD['pathways_files'])

    comparisondico = deg_dam_corresp(confD['DEG_tables'],
                                     confD['DAM_tables'],
                                     confD['titles'])

    DEG_full = pile_dfs_byconds(comparisondico, confD['dir_deg'],
                                confD['columns_deg'], "gene")
    DAM_full = pile_dfs_byconds(comparisondico, confD['dir_dam'],
                                confD['columns_dam'], "metabolite")

    DEG_full = filter_by_path_dico(DEG_full, genes_dico)
    DAM_full = filter_by_path_dico(DAM_full, metabo_dico)

    mycmap = get_custom_color_palette_hash('#0070C0', 'white', '#D30000')

    mabs = max(abs(DAM_full['log2FC'])) # caution
    gabs = max(abs(DEG_full['log2FC'])) # caution

    DEG_full['mycolors'] = rgbas2hex(values2rgbas(DEG_full['log2FC'].to_numpy(), mycmap, -gabs, gabs, center=0))
    DAM_full['mycolors'] = rgbas2hex(values2rgbas(DAM_full['log2FC'].to_numpy(), mycmap, -mabs, mabs, center=0))

    gathered = pd.concat([DEG_full, DAM_full], axis=0)

    gathered['typemol'] = pd.Categorical(gathered['typemol'], categories=['metabolite', 'gene'])
    gathered = gathered.sort_values(by=['typemol', 'name'], ascending=[True, False])

    gathered["elem_tag"] = gathered["name"]
    add_complex_tags = False # for data deep viz set to True
    if add_complex_tags:
        log2FCstrli = [str(i) for i in gathered["log2FC"]]
        gathered["elem_tag"] = gathered["name"].str.cat(log2FCstrli, sep=": ")
    print("shaped data for metabologram")

    # also see :https://proplot.readthedocs.io/en/latest/why.html
    ###################################### complicated grid
    # PLOT GRID :
    # as many columns as comparisons,
    # as many rows as paths  + add supplementary row(s) for bars
    nbpaths = len(pathdico)

    nbcompars = len(comparisondico)  # this has to be deduced from tables files

    tf = True

    if nbcompars == 1:
        supplerows = 3
    elif nbcompars == 2:
        supplerows = 2
    else:
        supplerows = 1

    if dimensions_pdf is None:
        dimensions_pdf = tuple(nbpaths * 7, supplerows * 9)
    #sns.set_style({'font.family': 'serif', 'font.serif': ['Meiryo']})
    fig, axes = plt.subplots(nrows=nbpaths + supplerows,
                             ncols=nbcompars, figsize=dimensions_pdf)
    fig.subplots_adjust(bottom=0, top=0.9, left=0, right=1,
                        wspace=0.2, hspace=0.4)
    # prepare subsetters as indexesdico for axes.flat usage
    indexesdico = dict()
    indexer = 0
    for i in pathdico.keys():
        for j in comparisondico.keys():
            indexesdico[indexer] = {'path': i, 'comparison': j, 'title': comparisondico[j]['title'] }
            indexer += 1
    indexer = 0
    for ax in axes.flat[:len(indexesdico)]:
        ax = axes.flat[indexer]
        print()
        path_elems_here = pathdico[indexesdico[indexer]['path']]
        gatheredsub = gathered.loc[gathered['name'].isin(path_elems_here), :]  # this will mess up the order
        compari_here = indexesdico[indexer]['comparison']
        title_here = indexesdico[indexer]['title']
        #ax.set_title(f"{indexesdico[indexer]['path']}\n {compari_here}\n")
        ax.set_title(f"{indexesdico[indexer]['path']}\n {title_here}\n")
        gatheredsub = gatheredsub.loc[gatheredsub['comparison'] == compari_here, :]

        ##################
        #  donut
        ##################
        gatheredsub['circportion'] = ''
        genecircportion = 50/gatheredsub.loc[gatheredsub.typemol == "gene",:].shape[0]
        metabocircportion = 50/gatheredsub.loc[gatheredsub.typemol == "metabolite",:].shape[0]
        gatheredsub.loc[gatheredsub.typemol == "gene", "circportion"] = genecircportion
        gatheredsub.loc[gatheredsub.typemol == "metabolite", "circportion"] = metabocircportion

        sizes_list = gatheredsub["circportion"]
        annots = gatheredsub["elem_tag"]
        mappedcolors_list = gatheredsub["mycolors"]

        if tf == False:
            ax.pie(sizes_list,
                   colors=mappedcolors_list,
                   wedgeprops={'width': 1, 'edgecolor': 'black', 'linewidth': 0.8},
                   radius=1,
                   startangle=90)
        else:  # add metabolites to the plot
            ax.pie(sizes_list,
                   colors=mappedcolors_list,
                   wedgeprops={'width': 1, 'edgecolor': 'black', 'linewidth': 0.8},
                   radius=1,
                   startangle=90,

                   labels=annots,  ## this one yiels the  labels annotated in the plot
                   textprops={'fontsize': 8}
                   )
            ## white circles for artist patches
        my_circle2 = plt.Circle((0, 0), radius=0.47, edgecolor="black", linewidth=1.6)
        my_circle = plt.Circle((0, 0), radius=0.465, color="white")
        ax.add_patch(my_circle2)
        ax.add_patch(my_circle)
        inner_dico = {'metabo_mean_val': gatheredsub.loc[gatheredsub.typemol == 'metabolite', 'log2FC'].mean(),
                      'gene_mean_val': gatheredsub.loc[gatheredsub.typemol == 'gene', 'log2FC'].mean()}
        inner_colorsD = inner_pie_colors(inner_dico, mycmap, gabs, mabs)
        # internal pie
        ax.pie([50, 50],
               colors=[inner_colorsD['metab_color'], inner_colorsD['gene_color']],
               wedgeprops={'width': 0.41, 'edgecolor': 'black', 'linewidth': 0.7},
               radius=0.41,
               startangle=90,
               labels=np.array([inner_dico['metabo_mean_val'], inner_dico['gene_mean_val']]).round(1),
               labeldistance=0.2)
        ax.axis('equal')
        ax.legend('', frameon=False)  # https://www.statology.org/remove-legend-matplotlib/
        # ax.tight_layout()
        # ####
        # end donut
        ###
        indexer += 1
    # end for

    # fill last three panels with color bar key

    # do "fake" separated heatmaps to take the colorbar key separately for metabolites and for genes
    sns.heatmap([[]], ax=axes.flat[-3], cmap=mycmap, center=0, cbar=True,
                annot=False,
                square=True,
                vmin=-mabs, vmax=mabs, cbar_kws={'shrink': 0.9, 'aspect': 10,
                                                 'label': 'Metabolite',
                                                 'drawedges': False})
    # axes.flat[-2].text(-0.3, 0.7, "metabolite", rotation=90)

    sns.heatmap([[]], ax=axes.flat[-2], cmap=mycmap, center=0, cbar=True,
                annot=False,
                square=True,
                vmin=-gabs, vmax=gabs, cbar_kws={'shrink': 0.9, 'aspect': 10,
                                                 'label': 'Transcript',
                                                 'drawedges': False})
    # axes.flat[-1].text(-0.3, 0.7, "gene", rotation=90)
    out_file = os.path.expanduser(confD['metabologram_out_dir'])
    plt.savefig( out_file + "metabologram_plot.pdf")


if __name__ == "__main__":
    print("\nMetabologram\n")

    parser = metabologram_args()
    args = parser.parse_args()
    configfile = os.path.expanduser(args.config)
    confidic = fg.open_config_file(configfile)
    dimensions_pdf = (15, 20) # TODO transform into option from metabologram_config
    metabologram_run( confidic, dimensions_pdf)



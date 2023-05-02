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
import warnings


def metabologram_args():
    parser = argparse.ArgumentParser(prog="python -m DIMet.src.metabologram",
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

def path_to_metabologram_plots(out_plot_dir):
    out_plot_dir = os.path.expanduser(os.path.join(out_plot_dir, "metabologram/"))
    fg.detect_and_create_dir(out_plot_dir)
    return(out_plot_dir)

def write_metabologram_plot(fig, dir, fname):
    fname = os.path.join(dir, fname)
    fig.savefig(fname=fname)



def metabologram_run(confidic, dimensions_pdf,
                     edgecolors=('#cecece','#8d8d8d'), linewidths=(1.6,1.2)):
    if type(edgecolors) is str:
        edgecolors = (edgecolors,)*2 
    if type(linewidths) is float or type(linewidths) is int:
        linewidths = (linewidths,)*2 
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

    if dimensions_pdf is None:
        dimensions_pdf = tuple(7, 5)
    #sns.set_style({'f~/Documents/projects/DIMet/gb-ldh-TD/analysis01/results/plots/metabologramont.family': 'serif', 'font.serif': ['Meiryo']})
    fig = plt.figure(figsize=dimensions_pdf)
    fig.tight_layout()
    # prepare subsetters as indexesdico for axes.flat usage
    indexesdico = dict()
    indexer = 0
    for i in pathdico.keys():
        for j in comparisondico.keys():
            indexesdico[indexer] = {'path': i, 'comparison': j, 'title': comparisondico[j]['title'] }
            indexer += 1
    
    out_path = path_to_metabologram_plots(confD['metabologram_out_dir'])
    for indexer, subdict in enumerate(indexesdico.values()):
        fig = plt.figure(figsize=dimensions_pdf)
        fig.tight_layout()
        print(f"Plotting metabologram {indexer+1} out of {len(indexesdico)}", end="\r")
        path_elems_here = pathdico[subdict['path']]
        gatheredsub = gathered.loc[gathered['name'].isin(path_elems_here), :]  # this will mess up the order
        compari_here = indexesdico[indexer]['comparison']
        title_here = subdict['title']
        #ax.set_title(f"{indexesdico[indexer]['path']}\n {compari_here}\n")
        fig.suptitle(f"{subdict['path']}\n {title_here}\n")
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

        plt.pie(sizes_list,
                colors=mappedcolors_list,
                wedgeprops={'width': 1, 'edgecolor': edgecolors[0], 'linewidth': linewidths[0]},
                radius=1,
                startangle=90,
                labels=annots if tf else None,  ## this one yiels the  labels annotated in the plot
                textprops={'fontsize': 8} if tf else None)
        ## white circles for artist patches
        plot_circles= [
            plt.Circle((0, 0), radius=0.47, edgecolor="black", linewidth=1.6),
            plt.Circle((0, 0), radius=0.465, color="white")
        ]
        fig.patches.extend(plot_circles)
        inner_dico = {'metabo_mean_val': gatheredsub.loc[gatheredsub.typemol == 'metabolite', 'log2FC'].mean(),
                      'gene_mean_val': gatheredsub.loc[gatheredsub.typemol == 'gene', 'log2FC'].mean()}
        inner_colorsD = inner_pie_colors(inner_dico, mycmap, gabs, mabs)
        # internal pie
        plt.pie([50, 50],
                colors=[inner_colorsD['metab_color'], inner_colorsD['gene_color']],
                wedgeprops={'width': 0.41, 'edgecolor': edgecolors[1], 'linewidth': linewidths[1]},
                radius=0.41,
                startangle=90,
                labels=np.array([inner_dico['metabo_mean_val'], inner_dico['gene_mean_val']]).round(1),
                labeldistance=0.2)
        # ax.axis('equal')
        # ax.legend('', frameon=False)  # https://www.statology.org/remove-legend-matplotlib/
        # ax.tight_layout()
        # ####
        # end donut
        ###
        # save pdf
        fname=f'{subdict["path"]}_comparison{subdict["comparison"]}.pdf'
        write_metabologram_plot(fig,out_path,fname)
    # end for
    fig, axes = plt.subplots(ncols=2, nrows=1, figsize=dimensions_pdf)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # Remove warning internal heatmap (df.shape)
        vs=(mabs, gabs) ; labels=('Metabolite', 'Transcript')
        for ax, v, label in zip(axes, vs, labels):
            sns.heatmap([[]], ax=ax, cmap=mycmap, center=0, cbar=True,
                        annot=False, yticklabels=False,
                        square=True,
                        vmin=-v, vmax=v, cbar_kws={'shrink': 0.9, 'aspect': 10,
                                                        'label': label,
                                                        'drawedges': False})
    write_metabologram_plot(fig, out_path, 'legend.pdf')

    print("\nDone plotting!")




if __name__ == "__main__":
    print("\nMetabologram\n")

    parser = metabologram_args()
    args = parser.parse_args()
    configfile = os.path.expanduser(args.config)
    confidic = fg.open_config_file(configfile)
    out_path = os.path.expanduser(confidic['out_path'])
    dimensions_pdf = (7, 5) # TODO transform into option from metabologram_config
    metabologram_run( confidic, dimensions_pdf)



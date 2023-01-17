import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize
from matplotlib.colors import LinearSegmentedColormap
import matplotlib as mpolib



# from txt files make pathways dico (admits mixed genes and metabolites by pathway),
# same gene/met can be shared by two or more pathways, this is ok
pathways_files = ['pathway1.txt', 'pathTWO.txt']
pathdico = dict()
for fi in pathways_files:
    tmp = pd.read_csv(fi, header=None)
    path_name = fi.replace(".txt", "")
    pathdico[path_name] = tmp[0].to_list()

print(pathdico) # {'pathway1': ['GENE1', 'GENE2', ... 'metabo1', 'metabo2' ...], ...}

dimensions_pdf = (14,10)
# comparisons in yaml file : provide each comparison 'title', followed by [DEG table , DAM table ]
comparisondico = {'fake1' : ['genesfake.tsv' ,'metabosfake.tsv'],
                  'rigolo' : ['genesrigolo.tsv' ,'metabosrigolo.tsv'],}


DEG_full = pd.DataFrame({'name': [], 'log2FC': []})
DAM_full = pd.DataFrame({'name': [], 'log2FC': []})

print('\nVerifying input data for metabologram')
for comparison in comparisondico.keys():
    print("")
    DEGtable = comparisondico[comparison][0]
    DAMtable = comparisondico[comparison][1]
    print('comparison : ', comparison)
    print('DEG table :', DEGtable)
    print('DAM table :', DAMtable)
    genesdf = pd.read_csv(DEGtable, sep='\t')
    metsdf =  pd.read_csv(DAMtable, sep='\t')
    metsdf = metsdf[['mets', 'log2FC']]
    metsdf.columns = ['name', 'log2FC']
    genesdf.columns = ['name', 'log2FC']
    metsdf = metsdf.assign(comparison = comparison, typemol = "metabolite")
    genesdf = genesdf.assign(comparison = comparison, typemol = "gene")
    DEG_full = pd.concat([DEG_full, genesdf], axis = 0)
    DAM_full = pd.concat([DAM_full, metsdf], axis = 0)

print(DEG_full)
print(DAM_full)


###################################################### color functions


def get_custom_color_palette_hash(lowcol, midcol, highcol):
    """
    courtesy from :
    https://richardhildebrand.wordpress.com/2019/09/18/create-a-custom-color-palette-with-matplotlib-and-seaborn/
    """
    colorlist = [lowcol, midcol, highcol]
    return LinearSegmentedColormap.from_list("", colorlist, N=256)

mycmap = get_custom_color_palette_hash('#0070C0', 'white', '#D30000')

# max absolute value both genes and metabo
mabs = max(abs(DAM_full['log2FC']))
gabs = max(abs(DEG_full['log2FC']))


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
    metaboval = inner_dico['metabo_sum_val']
    geneval = inner_dico['gene_sum_val']
    metabocolors_ = rgbas2hex(values2rgbas([metaboval], mycmap, -mabs, mabs, center=0))
    genecolors_ = rgbas2hex(values2rgbas([geneval], mycmap, -gabs, gabs, center=0))
    return {'metab_color' : metabocolors_[0], 'gene_color': genecolors_[0]}


DEG_full['mycolors'] = rgbas2hex(values2rgbas(DEG_full['log2FC'].to_numpy(), mycmap, -gabs, gabs, center=0))
DAM_full['mycolors'] = rgbas2hex(values2rgbas(DAM_full['log2FC'].to_numpy(), mycmap, -mabs, mabs, center=0))


DAM_full['much'] = 50/(DAM_full.shape[0])
DEG_full['much'] = 50/(DEG_full.shape[0])
gathered = pd.concat([DAM_full,DEG_full], axis = 0)

gathered['typemol'] = pd.Categorical(gathered['typemol'], categories = ['metabolite', 'gene'])
gathered = gathered.sort_values(by='typemol')

log2FCstrli = [str(i) for i in gathered["log2FC"]]
gathered["blah"] =  gathered["name"].str.cat(log2FCstrli, sep = ": ")


print("preparation ok")
# also see :https://proplot.readthedocs.io/en/latest/why.html


###################################### complicated grid
# PLOT GRID :
# as many columns as comparisons,
# as many rows as paths  + add supplementary row(s) for bars
nbpaths = len(pathdico)

nbcompars = len(comparisondico) # this has to be deduced from tables files

tf = True

if  nbcompars == 1:
    supplerows = 2
else:
    supplerows = 1


if dimensions_pdf is None:
    dimensions_pdf = tuple(nbpaths*7, supplerows*9)
sns.set_style({'font.family' : 'serif', 'font.serif':['Meiryo']})
fig, axes = plt.subplots(nrows = nbpaths + supplerows,
                         ncols = nbcompars, figsize=dimensions_pdf)
fig.subplots_adjust(bottom=0, top=0.9, left=0, right=1,
                    wspace=0.2, hspace=0.4)
# prepare subsetters as indexesdico for axes.flat usage
indexesdico = dict()
indexer = 0
for i in pathdico.keys():
    for j in comparisondico.keys():
        indexesdico[indexer] = {'path': i, 'comparison': j}
        indexer += 1
indexer = 0
for ax in axes.flat[:len(indexesdico)]:
    print(indexer)
    ax = axes.flat[indexer]
    path_elems_here = pathdico[ indexesdico[indexer]['path'] ]
    gatheredsub = gathered.loc[gathered['name'].isin(path_elems_here),: ]
    compari_here = indexesdico[indexer]['comparison']
    ax.set_title(f"{indexesdico[indexer]['path'] }\n {compari_here}\n")
    gatheredsub = gatheredsub.loc[gatheredsub['comparison'] == compari_here,: ]
    print(gatheredsub)

    ##################
    #  donut
    ##################
    log2FC = gatheredsub["log2FC"]
    sizefoo = gatheredsub["much"]
    annots = gatheredsub["blah"]
    mappedcolorshaha = gatheredsub["mycolors"]

    if tf == False:
        ax.pie(sizefoo,
               colors=mappedcolorshaha,
               wedgeprops={'width': 1, 'edgecolor': 'black', 'linewidth': 0.8},
               radius=1,
               startangle=90)
    else:  # add metabolites to the plot
        ax.pie(sizefoo,
               colors=mappedcolorshaha,
               wedgeprops={'width': 1, 'edgecolor': 'black', 'linewidth': 0.8},
               radius=1,
               startangle=90,

               labels=annots,  ## this one yiels the  labels annotated in the plot
               textprops={'fontsize': 12}
               )
        ## white circles for artist patches
    my_circle2 = plt.Circle((0, 0), radius=0.47, edgecolor="black", linewidth=1.6)
    my_circle = plt.Circle((0, 0), radius=0.465, color="white")
    ax.add_patch(my_circle2)
    ax.add_patch(my_circle)
    inner_dico = { 'metabo_sum_val' : gatheredsub.loc[gatheredsub.typemol == 'metabolite', 'log2FC'].sum(),
                    'gene_sum_val' : gatheredsub.loc[gatheredsub.typemol == 'gene', 'log2FC'].sum() }
    inner_colorsD = inner_pie_colors(inner_dico, mycmap,  gabs, mabs)
    innerlabelsorder = [ inner_colorsD['metab_color'], inner_colorsD['gene_color'] ]
    print(inner_dico)
    # internal pie
    ax.pie([50, 50],
           colors=[ inner_colorsD['metab_color'], inner_colorsD['gene_color'] ],
           wedgeprops={'width': 0.41, 'edgecolor': 'black', 'linewidth': 0.7},
           radius=0.41,
           startangle=90,
           labels = np.array([inner_dico['metabo_sum_val'], inner_dico['gene_sum_val']]).round(1),
            labeldistance = 0.2)
    ax.axis('equal')
    ax.legend('', frameon=False)  # https://www.statology.org/remove-legend-matplotlib/
    # ax.tight_layout()
    # ####
    # end donut
    ###
    indexer += 1
# end for

# fill last two panels with color bar key


# do "fake" separated heatmaps to take the colorbar key separately for metabolites and for genes
sns.heatmap([[]], ax=axes.flat[-2], cmap=mycmap, center=0, cbar=True,
            annot=False,
            square=True,
            vmin=-mabs, vmax=mabs, cbar_kws={'shrink': 0.9, 'aspect': 10,
                                             'label': 'metabolite',
                                             'drawedges': False})
# axes.flat[-2].text(-0.3, 0.7, "metabolite", rotation=90)

sns.heatmap([[]], ax=axes.flat[-1], cmap=mycmap, center=0, cbar=True,
            annot=False,
            square=True,
            vmin=-gabs, vmax=gabs, cbar_kws={'shrink': 0.9, 'aspect': 10,
                                             'label': 'gene',
                                             'drawedges': False})
# axes.flat[-1].text(-0.3, 0.7, "gene", rotation=90)

plt.savefig("MYDONUTS.pdf")

# thanks also to:
#https://stackoverflow.com/questions/49199164/increasing-pie-chart-size-with-matplotlib-radius-parameter-appears-to-do-nothin

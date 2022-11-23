#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PCA plots

Created on  Nov 15  2022
@author: johanna
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import seaborn as sns
from scipy import stats

from matplotlib.patches import Ellipse


# calculate Ellipses , many thanks to :
# https://rayblick.gitbooks.io/my-python-scrapbook/content/analysis/plotting/scatterplot_ellipse.html
# use make_ellipse OR eigsorted, but not both


def make_ellipse(mean, cov,  level=0.95, color=None):
    """Support function for scatter_ellipse."""
    # ex 4 : https://www.programcreek.com/python/example/61396/matplotlib.patches.Ellipse
    v, w = np.linalg.eigh(cov)
    u = w[0] / np.linalg.norm(w[0])
    angle = np.arctan(u[1]/u[0])
    angle = 180 * angle / np.pi # convert to degrees
    v = 2 * np.sqrt(v * stats.chi2.ppf(level, 2)) #get size corresponding to level
    ell = Ellipse(mean[:2], v[0], v[1], 180 + angle, facecolor='none',
                  edgecolor="gray", alpha=0.4,
                  #ls='dashed',  #for debugging
                  lw=1.5)
    return ell

def eigsorted(cov):
    """
    used for calculating ellipses
    many thanks to :
https://rayblick.gitbooks.io/my-python-scrapbook/content/analysis/plotting/scatterplot_ellipse.html
    """
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    return vals[order], vecs[:,order]

def massage_datadf_4pca(df, metadatasub):
    df = df[metadatasub['sample']]

    advtests = False
    if advtests == False:
        df = df[~(df == 0).all(axis=1)] # drop zero rows
    elif advtests == True:
        print("Danger: advanced tests  set to true")
        # if row zero all, replace by nan, then replace by minimal of each col
        df = df.where(~(df==0).all(axis=1)) # zero row became nan row
        df = df.fillna(df[df>0].min()) # this is column wise !

    df = df.div(df.std(axis=1, ddof=0), axis=0)  # reduce rows
    df = df.dropna(axis=0)
    return df

def calcPCAand2Dplot(mymat, metadata, col1, col2, pointlabels, title, odir, *args):
    """
    *args :  options for advanced PCA plotting:
          args[0] is how many linear combinations to calculate, must be <= nb variables
          args[1] is column name for ellipses

    """
    col_ellipses = None
    desdim = 2

    try:
        desdim = args[0]
        col_ellipses = args[1]
    except IndexError as indexerr:
        # print("ellipses args not set == > ", indexerr)
        pass # args[1] was not set but we can continue
    except Exception as e:
        print("unknown error! ", e)


    X = np.transpose(np.array(mymat))
    pca = PCA(n_components=desdim)
    pc = pca.fit_transform(X)

    pc_df = pd.DataFrame(data=pc,
                         columns=['PC' + str(i) for i in range(1, desdim + 1)])
    pc_df = pc_df.assign(sample = mymat.columns)

    pc_df = pd.merge(pc_df, metadata, on= "sample")

    dfvare = pd.DataFrame({
        # 'var': pca.explained_variance_ratio_,
        'var': pca.explained_variance_ratio_ * 100,
        'PC': ['PC' + str(i) for i in range(1, desdim + 1)]})

    # barplot
    plt.figure()
    sns.barplot(x='PC', y="var", data=dfvare, color="cadetblue")
    plt.title("Principal components: explained variance")
    plt.ylabel("variance (%)")
    plt.savefig(odir+"variance_pca_"+title+".pdf")
    plt.close()


    # scatterplot

    # sns.set_style("whitegrid")
    fig,ax = plt.subplots()
    sns.scatterplot(x="PC1", y="PC2",
                    ax = ax,
                    data=pc_df,
                    hue=col1,
                    style=col2,
                    legend=True,
                    s=80, zorder=3)
    ax.axhline(0, ls="--", color="gray",zorder = 1)
    ax.axvline(0, ls="--", color="gray", zorder = 1)

    yesnolabel = "no"
    if pointlabels != "":
        yesnolabel = "yes"
        for l in range(0, pc_df.shape[0]):
            ax.text(pc_df.PC1[l] + 0.5, pc_df.PC2[l]+0.1, pc_df[pointlabels].tolist()[l],
                                    horizontalalignment='center', size='x-small')
    # end if
    plt.xlabel(f"{dfvare.iloc[0, :]['PC']} {round(dfvare.iloc[0, :]['var'], 2) } %")
    plt.ylabel(f'{dfvare.iloc[1, :]["PC"]} {round(dfvare.iloc[1, :]["var"], 2) } %')
    plt.title(title)

    # ellipses, if true
    if col_ellipses is not None:
        myellipsesnames = pc_df[col_ellipses].unique()
        for lab in myellipsesnames:
            xdata = pc_df.loc[pc_df[col_ellipses] == lab, 'PC1']
            ydata = pc_df.loc[pc_df[col_ellipses] == lab, 'PC2']

            # get values to build the ellipse
            cov = np.cov(xdata, ydata)
            #ell = make_ellipse((np.mean(xdata), np.mean(ydata)),cov,
            #                  level= 0.95, color=None)
            vals, vecs = eigsorted(cov)
            theta = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
            print(vals, vecs, theta)
            w, h = 2 * 2 * np.sqrt(vals)

            # create the ellipse
            ell = Ellipse(xy=(np.mean(xdata), np.mean(ydata)),
                          width=w, height=h,
                          angle=theta,
                          edgecolor="lightgray",
                  ls='dashed', facecolor='none')

            # ell.set_facecolor() # reference the colour for each factor
            ax.add_artist(ell)

    # end  if ellipses
    plt.savefig(f"{odir}pca_{title}_label{yesnolabel}.pdf", format="pdf")
    return pc_df, dfvare

if __name__ == "__main__":


    doiristest = True  # for debug
    if doiristest == True:
        iris = sns.load_dataset("iris")
        sns.relplot(data=iris, x="sepal_width", y="petal_width",
                    hue="species")
        iris = iris.assign(sample=[str(i) for i in iris.index])
        fakemeta = iris[["sample", "species"]]
        iris = iris.drop(columns=["sample", "species"])
        fakedf = iris.T  # variables rows, samples columns
        fakedf = fakedf.div(fakedf.std(axis=1, ddof=0), axis=0)

        fakedf.columns = [str(i) for i in iris.index]
        calcPCAand2Dplot(fakedf, fakemeta, "species", "species", "",
                         "irisTest", "/home/johanna/myelin_metabo/", 3, "species")


    # metadata = pd.read_csv("data/metadata_myelin.csv", header = 0)
    # tmpdir = "tmp/"
    # compartments_names = {"cells" : "cell",
    #                       "supernatant" : "sn"}
    # suffix = "myelin"
    # odir = "results/plots/pca/"
    # tps = ["T0h", "T4h", "T24h"]
    #
    # domytest = True
    # if domytest == True:
    #     for k in compartments_names.keys():
    #         co = compartments_names[k]
    #         df = pd.read_csv(f"{tmpdir}meanEnrich_{suffix}_{co}.tsv", sep='\t', header = 0, index_col = 0)
    #         metadatasub = metadata.loc[metadata['short_comp'] == co,:]
    #         df = df[metadatasub['sample']]
    #         # reduce by std each metabolite :
    #         #df = df + 1
    #         df = df.div(df.std(axis=1, ddof=0), axis=0) # reduce rows
    #         #df = df.div(df.sum(axis = 0 )) # normalize columns
    #         df = df.dropna(axis = 0 )
    #
    #         #calcPCAand2Dplot(df, metadatasub, col1 = "timepoint" , col2="condition",
    #         #                 pointlabels="",
    #         #                 title = f'{suffix}-{k}', odir=odir, 6)
    #         calcPCAand2Dplot(df, metadatasub, "timepoint" , "condition",
    #                         "",  f'{suffix}-{k}', odir, 6)
    #         calcPCAand2Dplot(df, metadatasub, "timepoint", "condition",
    #                          "sample_descrip", f'{suffix}-{k}', odir, 6)
    #
    #     print("now parsing time-points to yield separate time-point dataframe for plot")
    #
    #     for k in compartments_names.keys():
    #         co = compartments_names[k]
    #         for tp in tps:
    #             df = pd.read_csv(f"{tmpdir}meanEnrich_{suffix}_{co}.tsv", sep='\t', header = 0, index_col = 0)
    #             metadatasub = metadata.loc[(metadata['short_comp'] == co) & (metadata['timepoint'] == tp),:]
    #             df = df[metadatasub['sample']]
    #             print(df.shape)
    #             print(metadatasub.shape)
    #             # reduce by std each metabolite :
    #             #df = df + 1
    #             df = df.div(df.std(axis=1, ddof=0), axis=0) # reduce rows
    #             #df = df.div(df.sum(axis = 0 )) # normalize columns
    #             df = df.dropna(axis = 0 )
    #             # mymat= ,  df=,  col1=, col2=, pointlabels=, title=, odir=, *args= desdim, col_ellipses ...
    #             calcPCAand2Dplot(df, metadatasub,  "condition" , "condition",
    #                              "sample_descrip",
    #                             f'{suffix}-{k}-{tp}', odir, 6)






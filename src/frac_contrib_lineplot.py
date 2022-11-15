#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  9 11:20:09 2022

@author: johanna
"""
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns



def yieldfraccountrib(dicos, tablePicked, co):
    """
    description TODO
    """
    assert (
        "contribution" in tablePicked.lower()
    ), "this funcion only works \
        with fractionnal contribution table or similar"
    dfcompartment = dicos[co][tablePicked].T

    metabolites = dfcompartment.columns
    dfcompartment["sample"] = dfcompartment.index
    dfcompartment = pd.merge(dfcompartment, dicos[co]["metadata"], on="sample")

    # empty dataframe to fill
    df4plot = pd.DataFrame(
        columns=["timenum", "condition", "metabolite", "Fractional Contribution (%)"]
    )

    for z in range(len(metabolites)):
        subdf = dfcompartment.loc[:, [metabolites[z], "timenum", "condition"]]
        subdf["metabolite"] = metabolites[z]
        subdf["Fractional Contribution (%)"] = subdf[metabolites[z]] * 100
        subdf = subdf.drop(columns=[metabolites[z]])
        # df4plot = df4plot.append(subdf, ignore_index=True)
        df4plot = pd.concat((df4plot, subdf), ignore_index=True)
    return df4plot


def nestedDi2list(anestedD):
    elems_ = []
    for p in anestedD.values():
        for q in p:
            elems_.append(q)
    return elems_


def df2musddf(one_m):
    """
    input: dataframe by one metabolite:
      "timenum", "condition", "metabolite" "Fractional Contribution (%)"
    returns :
        dataframe by one metabolite, metabolite column deleted:
                condition  timenum    mean    sd
    108        Control     0  0.000893  0.002611
    111        Control     1  0.236453  0.023246
    ...
    123        Control    24  0.854101  0.055241
    126  L-Cycloserine     0  0.010083  0.003465
    129  L-Cycloserine     1  0.259570  0.008602
    ...
    141  L-Cycloserine    24  0.815613  0.050756
    """
    m_s = one_m[["condition", "timenum"]].drop_duplicates()
    m_s["mean"] = 0
    m_s["sd"] = 0
    for kc in m_s["condition"]:
        for hou in one_m["timenum"]:
            ss = one_m.loc[(one_m["timenum"] == hou) & (one_m["condition"] == kc), :]
            mymean = np.mean(ss["Fractional Contribution (%)"])
            mysd = np.std(ss["Fractional Contribution (%)"])
            m_s.loc[(m_s["condition"] == kc) & (m_s["timenum"] == hou), "mean"] = mymean
            m_s.loc[(m_s["condition"] == kc) & (m_s["timenum"] == hou), "sd"] = mysd
    return m_s


def complextimetracer(co, df4plot, grmetsD, mycolorsD, outfile):
    themets = nestedDi2list(grmetsD)
    hoursticks = df4plot['timenum'].unique()
    somem = df4plot.loc[df4plot["metabolite"].isin(themets)]
    m_s = pd.DataFrame(columns=["condition", "timenum", "mean", "sd", "metabolite"])
    for k in set(somem["metabolite"]):
        one_m = somem[somem["metabolite"] == k]
        m_s1 = df2musddf(one_m)
        m_s1["metabolite"] = k
        m_s = pd.concat([m_s, m_s1])  # mean and sd ( bio replicates)
    figziz = 6.2 * len(grmetsD)
    fig, axs = plt.subplots(2, len(grmetsD), sharey=False, figsize=(figziz, 11))
    def doemptyrow(axs, nbcolumns):
        for z in range(nbcolumns): # do empty row
            axs[0,z].set_axis_off()
        return axs

    axs = doemptyrow(axs, len(grmetsD))

    for z in range(len(grmetsD)):
        sns.lineplot(
            ax=axs[1,z],
            x="timenum",
            y="Fractional Contribution (%)",
            hue="metabolite",
            style="condition",
            err_style=None,
            alpha=0.9,
            linewidth = 4.5,
            palette=mycolorsD,
            data=somem.loc[somem["metabolite"].isin(grmetsD[z])],
            legend=True,  ## !!!!!!!!!!!!!!!!attention here   <=====
        )
        axs[1,z].set_xticks([int(i) for i in hoursticks])
        m_s1 = m_s.loc[m_s["metabolite"].isin(grmetsD[z])]
        axs[1,z].scatter(
            m_s1["timenum"], m_s1["mean"], s=22, facecolors="none", edgecolors="black"
        )
        axs[1,z].errorbar(
            m_s1["timenum"],
            m_s1["mean"],
            yerr=m_s1["sd"],
            fmt="none",
            capsize=3,
            ecolor="black",
            zorder=3,
        )
        axs[1,z].set(ylabel=None),
        axs[1,z].legend(loc="upper center", bbox_to_anchor= (0.5, 2), frameon=False)
        # ha , la = axs[z].get_legend_handles_labels()
        # handless.append(ha)
        # labels.append(la)
    fig.suptitle(co)
    plt.subplots_adjust(bottom=0.1, right=0.8, hspace=0.1, wspace=0.4)

    fig.text(
        0.1, 0.3, "Fractional Contribution (%)", va="center", rotation="vertical"
    )
    # loc="upper left"
    # fig.legend(handless, labels, bbox_to_anchor=(0.65,1.05))
    #fig.legend(handles=handless, labels=labels)
    fig.savefig(outfile, format="pdf")
    return 0


def savefraccontriplots(
    datadi, names_compartments, metadata, tableFC, namesuffix, gbycompD, coloreachmetab
):
    for co in names_compartments.values():
        adf = pd.read_csv(
            datadi + tableFC + "_" + namesuffix + "_" + co + ".tsv",
            sep="\t",
            index_col=0,
        )
        dicos = dict()
        dicos[co] = {}
        dicos[co]["metadata"] = metadata.loc[metadata.short_comp == co]
        dicos[co][tableFC] = adf[dicos[co]["metadata"]["sample"]]

        df4plot = yieldfraccountrib(dicos, tableFC, co)
        grmetsD = gbycompD[co]

        sns.set_style(
            {"font.family": "sans-serif", "font.sans-serif": "Liberation Sans"}
        )
        plt.rcParams.update({"font.size": 22})
        # https://stackoverflow.com/questions/53137983/define-custom-seaborn-color-palette

        odiric = "results/plots/fc/"
        if not os.path.exists(odiric):
            os.makedirs(odiric)
        ofile = "{}fc_plot_{}.pdf".format(odiric, co)
        complextimetracer(co, df4plot, grmetsD, coloreachmetab, ofile)

        #### uninteresting metabolites : quick ugly print :
        pickedmets_ = nestedDi2list(grmetsD)
        thiscompartment_mets = set(df4plot["metabolite"])

        trashmets = set(thiscompartment_mets) - set(pickedmets_)
        ppp = df4plot.loc[df4plot["metabolite"].isin(trashmets)]

        plt.rcParams.update({"font.size": 10})

        g = sns.FacetGrid(
            ppp, row="condition", col="metabolite", sharey=False, margin_titles=True
        )
        g.map(
            sns.lineplot,
            "timenum",
            "Fractional Contribution (%)",
            marker="o",
            err_style="bars",
            ci=95,
            color="black",
            alpha=0.4,
            size=5,
        )
        g.savefig(odiric + "fc_remaining_" + co + ".pdf", format="pdf")
    return 0

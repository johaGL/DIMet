#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 28 09:33:01 2022
This script exploits species and compartment specific tables at 'data/abufromperc/'
(resulting from abund_frompercentages)

@author: johanna
"""

import os

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from .fun_fm import *


def fi2_smx_compartment(filename, indexsmx, indexcompartment):
    """example output: 'm+3' , 'cell'"""
    elems = filename.replace(".tsv", "").split("_")
    smx = elems[indexsmx]
    compartment = elems[indexcompartment]
    return smx, compartment


def stackallabundace(abu_sel, metada_sel):
    """
    all this function is inspired from yieldfraccontrib (frac_contrib_v4)
    """
    dfcompartment = abu_sel.T
    metabolites = dfcompartment.columns
    dfcompartment["sample"] = dfcompartment.index
    dfcompartment = pd.merge(dfcompartment, metada_sel, on="sample")
    dafull = pd.DataFrame(columns=["timepoint", "condition", "metabolite", "abundance"])
    for z in range(len(metabolites)):
        subdf = dfcompartment.loc[:, [metabolites[z], "timepoint", "condition"]]
        subdf["metabolite"] = metabolites[z]
        subdf["abundance"] = subdf[metabolites[z]] * 100
        subdf = subdf.drop(columns=[metabolites[z]])
        dafull = pd.concat(
            [dafull, subdf], ignore_index=True
        )  # TODO : replace df4plot.append by pd.concat([]) elsewhere
    return dafull


def mean_sd_D(newdf, ametadasel):
    """
    returns dico (oD) with keys: metabolite,
    mu|condition_timepoint ,
    sd|condition_timepoint

    All in same order as metabolites
    example : sd|L-Cycloserine_T0
    """
    oD = dict()
    condis = set(ametadasel["condition"])
    times = set(ametadasel["timepoint"])
    for cd in condis:
        for t in times:
            sasas = ametadasel.loc[
                (ametadasel["condition"] == cd) & (ametadasel["timepoint"] == t)
            ]

            tmp = newdf[sasas["sample"]]
            means_here = []
            sd_here = []
            for i, row in tmp.iterrows():  # iterates metabolites
                means_here.append(np.mean(row))
                sd_here.append(np.std(row))
            oD["metabolite"] = tmp.index
            oD[f"sd|{cd}_{t}"] = sd_here
            oD[f"mu|{cd}_{t}"] = means_here
    return oD


def tmpstack(oD):
    oudf = pd.DataFrame(
        data={"metabolite": [], "condition": [], "timepoint": [], "mean": [], "sd": []}
    )
    mets = oD["metabolite"]
    for k in list(oD.keys())[1:]:
        cd, t = k.split("|")[1].split("_")  # this separator "_" by default
        tmp = pd.DataFrame(
            data={
                "metabolite": mets,
                "condition": [cd for i in range(len(mets))],
                "timepoint": [t for i in range(len(mets))],
                "mean": [x for x in oD[f"mu|{cd}_{t}"]],
                "sd": [x for x in oD[f"sd|{cd}_{t}"]],
            }
        )
        oudf = pd.concat([oudf, tmp])
    oudf = oudf.drop_duplicates()
    return oudf


def printtestpluslegend(zoo, selectedmets, CO, SMX, plotwidth, odirbars):
    ili = selectedmets
    zooi = zoo.loc[zoo["metabolite"].isin(ili), :]
    sns.set_style({"font.family": "sans-serif", "font.sans-serif": "Liberation Sans"})
    plt.rcParams.update({"font.size": 21})
    YLABE = "Abundance"
    fig, axs = plt.subplots(1, len(ili), sharey=False, figsize=(plotwidth, 7))
    #  constrained_layout=True)
    for il in range(len(ili)):
        ici = zooi.loc[zooi["metabolite"] == ili[il], :]
        ici = ici.reset_index()  # NEW
        sns.barplot(
            ax=axs[il],
            data=ici,
            x="condition",
            y="mean",
            hue="timepoint",
            palette="pastel",
            alpha=1,
            edgecolor="black",
        )
        axs[il].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        # axs[il].legend_.remove()
        axs[il].set(title=" " + ili[il] + "\n")
        axs[il].set(ylabel="")
        axs[il].set(xlabel="")
        sns.despine(ax=axs[il])
        axs[il].tick_params(axis="x", labelrotation=20)
        axs[il].set_ylim(bottom=0)  # set minimal val display : y axis : 0
        imaxmean = ici["mean"].idxmax()  # NEW
        maxva = ici.loc[imaxmean, "mean"]  # NEW
        sdofmaxva = ici.loc[imaxmean, "sd"]  # NEW
        axs[il].set_ylim(top=maxva + sdofmaxva)  # NEW
        isd = 0
        # note: as hue plotted by timepoint: sort_values to have good sd by bar
        icialter = ici.sort_values(by=["timepoint", "condition"])
        # icialter = ici
        for p in axs[il].patches:
            xp = p.get_x()
            wp = p.get_width()
            yp = p.get_height()
            axs[il].errorbar(
                xp + (wp / 2),
                yp,
                yerr=icialter["sd"].tolist()[isd],
                fmt="none",
                capsize=10,
                ecolor="black",
                alpha=0.9,
                zorder=1,
            )  # zorder -1 behind color bar
            isd += 1
        # end for p
    # end for il
    fig.text(-0.1, 0.5, YLABE, va="center", rotation="vertical", size=26)
    fig.suptitle(f"{CO} {SMX}".upper())
    plt.savefig(f"{odirbars}bars_legend.pdf", format="pdf")
    return 0


def printabundbarswithdots(piled_sel, selectedmets, CO, SMX, plotwidth, odirbars):
    li_ = selectedmets
    sns.set_style({"font.family": "sans-serif", "font.sans-serif": "Liberation Sans"})
    plt.rcParams.update({"font.size": 21})
    YLABE = "Abundance"
    fig, axs = plt.subplots(1, len(li_), sharey=False, figsize=(plotwidth, 5.5))
    #  constrained_layout=True)
    for il in range(len(li_)):
        herep = piled_sel.loc[piled_sel["metabolite"] == li_[il], :]
        herep = herep.reset_index()
        sns.barplot(
            ax=axs[il],
            data=herep,
            x="condition",
            y="abundance",
            hue="timepoint",
            palette="pastel",
            alpha=1,
            edgecolor="black",
            errcolor="black",
            errwidth=1.7,
            ci="sd",
            capsize=0.2,
        )
        sns.stripplot(
            ax=axs[il],
            data=herep,
            x="condition",
            y="abundance",
            hue="timepoint",
            palette="pastel",
            dodge=True,
            edgecolor="black",
            linewidth=1.5,
            alpha=1,
        )
        axs[il].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        axs[il].legend_.remove()
        axs[il].set(title=" " + li_[il] + "\n")
        axs[il].set(ylabel="")
        axs[il].set(xlabel="")
        sns.despine(ax=axs[il])
        axs[il].tick_params(axis="x", labelrotation=20)
        axs[il].set_ylim(bottom=0)  # set minimal val display : y axis : 0

    # end for il
    # fig.text(0.1, 0.5, YLABE, va = "center", rotation = 'vertical', size=26)
    # xlaloc = (0.025 * len(selectedmets)) - 0.15
    fig.text(0.04, 0.5, YLABE, va="center", rotation="vertical", size=26)
    fig.suptitle(f"{CO} {SMX}".upper())
    plt.subplots_adjust(top=0.76, bottom=0.2, wspace=0.6, hspace=0.5)
    # plt.tight_layout(pad = 0.01, w_pad = -2, h_pad=0.1)
    plt.savefig(f"{odirbars}bars_{CO}_{SMX}.pdf", format="pdf")
    return 0


###############
# https://stackoverflow.com/questions/42017049/how-to-add-error-bars-on-a-grouped-barplot-from-a-column
# give allways scientific notation :
# https://www.adamsmith.haus/python/answers/how-to-scale-an-axis-to-scientific-notation-in-a-matplotlib-plot-in-python

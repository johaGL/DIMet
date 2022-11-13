#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 28 09:33:01 2022
This script exploits species and compartment specific tables at 'data/abufromperc/'
(resulting from abund_frompercentages)

@author: johanna
"""


import matplotlib.pyplot as plt

import seaborn as sns

from .fun_fm import *

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



def printabundbarswithdots(piled_sel, selectedmets, CO, SMX, col1, col2, plotwidth, odirbars):
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
            x=col1,
            y="abundance",
            hue=col2,
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
            x=col1,
            y="abundance",
            hue=col2,
            palette="pastel",
            dodge=True,
            edgecolor="black",
            linewidth=1.5,
            alpha=1,
        )
        axs[il].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        axs[il].set(title=" " + li_[il] + "\n")
        axs[il].set(ylabel="")
        axs[il].set(xlabel="")
        sns.despine(ax=axs[il])
        axs[il].tick_params(axis="x", labelrotation=20)
        axs[il].set_ylim(bottom=0)  # set minimal val display : y axis : 0

    thehandles, thelabels = axs[-1].get_legend_handles_labels()
    for il in range(len(li_)):
        axs[il].legend_.remove()

    # end for il
    # fig.text(0.1, 0.5, YLABE, va = "center", rotation = 'vertical', size=26)
    # xlaloc = (0.025 * len(selectedmets)) - 0.15

    fig.text(0.1, 0.5, YLABE, va="center", rotation="vertical", size=26)
    fig.suptitle(f"{CO} {SMX}".upper())
    plt.subplots_adjust(top=0.76, bottom=0.2, wspace=0.3, hspace=1)
    plt.legend(handles=thehandles, labels=thelabels, loc='upper right',
               bbox_to_anchor=(-plotwidth/3, 1))
    # plt.tight_layout(pad = 0.01, w_pad = -2, h_pad=0.1)
    plt.savefig(f"{odirbars}bars_{CO}_{SMX}.pdf", format="pdf")
    plt.close()
    plt.figure()
    plt.legend(handles=thehandles, labels=thelabels, loc='upper right')
    plt.savefig(f"{odirbars}legend.pdf", format="pdf")

    return 0


###############
# https://stackoverflow.com/questions/42017049/how-to-add-error-bars-on-a-grouped-barplot-from-a-column
# give allways scientific notation :
# https://www.adamsmith.haus/python/answers/how-to-scale-an-axis-to-scientific-notation-in-a-matplotlib-plot-in-python

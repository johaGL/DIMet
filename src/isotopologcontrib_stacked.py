#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 18 11:14:34 2022

@author: johanna
"""

import os

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.ticker as mticker

def icontrib_2df4plot(dicos, tablePicked, co, levelshours_str):
    """
    Imitates de behaviour of a 'melt', but this function is more secure
    example:
      pd.merge output, where 1st column has the values we want:

        L-Phenylalanine_C13-label-2    timenum   condition ...  sample descr...
        0.01                           0          Control        xxxxxx

      is transformed into:
          timenum    condition    isotopolgFull        Isotopologue Contrib
          0           control      L-Phenylala...    0.01

    """
    dfcompartment = dicos[co][tablePicked].T

    metabolites = dfcompartment.columns
    dfcompartment["sample"] = dfcompartment.index
    dfcompartment = pd.merge(dfcompartment, dicos[co]["metadata"], on="sample")
    # empty dataframe to fill
    df4plot = pd.DataFrame(
        columns=[
            "timenum",
            "condition",
            "isotopolgFull",
            "Isotopologue Contribution (%)",
        ]
    )
    df4plot["timenum"] = pd.Categorical(df4plot["timenum"], levelshours_str)
    # iteratively pile up
    for z in range(len(metabolites)):  # corrected! error I had : -30)): now ok :)
        subdf = dfcompartment.loc[:, [metabolites[z], "timenum", "condition"]]
        subdf['timenum'] = subdf['timenum'].astype(str)
        subdf["isotopolgFull"] = metabolites[z]  # 1st colname as cell value, reps
        subdf["Isotopologue Contribution (%)"] = subdf[metabolites[z]] * 100
        subdf = subdf.drop(columns=[metabolites[z]])  # 1st col no longer needed
        df4plot = pd.concat([df4plot, subdf])
        del subdf

    return df4plot


def massageisotopologues(df4plot):
    """
    returns dataframe splitting metabolite and m+x into two separate columns
    and also correcting weird values
    """
    xu = {"name": [], "m+x": []}
    for ch in df4plot["isotopolgFull"]:
        elems = ch.split("_m+")
        xu["name"].append(elems[0])
        xu["m+x"].append("m+" + elems[1])
    df4plot["metabolite"] = xu["name"]
    df4plot["m+x"] = xu["m+x"]
    # dealing with weird values: bigger than 100 and less than 0 :
    df4plot.loc[
        df4plot["Isotopologue Contribution (%)"] > 100, "Isotopologue Contribution (%)"
    ] = 100

    df4plot.loc[
        df4plot["Isotopologue Contribution (%)"] < 0, "Isotopologue Contribution (%)"
    ] = 0

    return df4plot


def preparemeansreplicates(df4plot, selectedmets):
    """
    returns a dictionary of dataframes, keys are metabolites
    """
    dfcopy = df4plot.copy()
    dfcopy = dfcopy.groupby(["condition", "metabolite", "m+x", "timenum"]).mean()
    dfcopy = dfcopy.reset_index()

    dfs_Dico = dict()
    for i in selectedmets:
        tmp = dfcopy.loc[dfcopy["metabolite"] == i, ].reset_index()
        dfs_Dico[i] = tmp.sort_values(by="m+x", axis=0, ascending=True, inplace=False)
    return dfs_Dico

def addcombinedconditime(dfs_Dico, combined_tc_levels ):
    """
    add column string '{timepoint} : {condition}' to each dataframe

    """
    for metab in dfs_Dico.keys():
        #dfs_Dico['time_cat'] = dfs_Dico['timepoint'].astype(str).replace("T", "")
        #dfs_Dico['time_cat'] = dfs_Dico['time_cat'].str.replace("h", "")
        dfs_Dico[metab]["timeANDcondi"] = dfs_Dico[metab]["timenum"] + " : " + \
            dfs_Dico[metab]["condition"]
        dfs_Dico[metab]["timeANDcondi"] = pd.Categorical(dfs_Dico[metab]["timeANDcondi"],
                                                         combined_tc_levels)
    return dfs_Dico





def complexstacked(
    co, selectedmets, dfs_Dico, darkbarcolor, palsD,
        outfilename, figu_width, xlabyesno
):
    """plot highly custom, recommended that selectedmets <= 6 subplots"""
    ### set font style
    sns.set_style({"font.family": "sans-serif", "font.sans-serif": "Liberation Sans"})
    f, axs = plt.subplots(1, len(selectedmets), sharey=False, figsize=(figu_width, 4.8))
    plt.rcParams.update({"font.size": 20})

    for z in range(len(selectedmets)):
        # sns.set_style({ 'font.family': 'sans-serif',
        #                'font.sans-serif' : 'Liberation Sans'   })
        axs[z].set_title(selectedmets[z])
        sns.histplot(
            ax=axs[z],
            data=dfs_Dico[selectedmets[z]],
            x="timeANDcondi",
            # Use the value variable here to turn histogram counts into weighted
            # values.
            weights="Isotopologue Contribution (%)",
            hue="m+x",
            multiple="stack",
            palette=palsD,  # ['#0000c0',  '#410257' ] , # '#440154FF'],
            # Add  borders to the bars.
            edgecolor="black",
            # Shrink the bars a bit so they don't touch.
            shrink=0.85,
            alpha=1,
            legend=False,
        )
        #
        axs[z].tick_params(axis="x", labelrotation=90, labelsize=18)
        axs[z].tick_params(axis="y", length=11, labelsize=19 )
        axs[z].set_ylim([0,100])
        for bar in axs[z].patches:
            selcol = "black"
            # print(bar.get_facecolor())
            herergba = bar.get_facecolor() # returns rgba !
            if herergba in darkbarcolor:
                selcol = "white"
            thebarvalue = round(bar.get_height(), 1)
            if thebarvalue >= 100:
                thebarvalue = 100  # no decimals if 100
            if round(bar.get_height(), 1) > 4:
                axs[z].text(
                    # Put the text in the middle of each bar. get_x returns the start
                    # so we add half the width to get to the middle.
                    bar.get_x() + bar.get_width() / 2,
                    # Vertically, add the height of the bar to the start of the bar,
                    # along with the offset.
                    (bar.get_height() / 2) + (bar.get_y()) + 2,  #
                    # This is actual value we'll show.
                    thebarvalue,
                    # Center the labels and style them a bit.
                    ha="center",
                    color=selcol,
                    size=int((figu_width / len(selectedmets)) * 2) ,
                )  # end axs[z].text
            else:
                continue
            # end if round(...)
        # end for bar

        axs[z].set_ylabel("", size=20)
        axs[z].xaxis.set_tick_params(length=0) # no need of x ticks
        axs[z].set_xlabel("", size=13)
    # end for z

    [ax.invert_yaxis() for ax in axs] # invert y, step 1


    for ax in axs:
        ylabels = ax.get_yticks().tolist()
        ax.yaxis.set_major_locator(mticker.FixedLocator(ylabels))
        ax.set_yticklabels([100 - int(i) for i in ylabels])  # invert y , step2

    if xlabyesno == "no":
        for ax in axs:
            xlabelshere = ax.get_xticks()
            ax.set_xticklabels(["" for i in xlabelshere])

    f.subplots_adjust(hspace=0.5, wspace=0.25, top = 0.85, bottom=0.26, left=0.15, right=0.99)
    f.suptitle(f"compartment : {co.upper()}  \n", fontsize=20)
    f.text(0.03, 0.57, "Isotopologue Contribution (%)\n", va="center", rotation="vertical", size=20)
    f.savefig(outfilename, format="pdf")
    plt.close()
    return 0


def default_colors_stacked():
    """
    Returns colors defined by default
     important : darkbarcolor sets which rgba obliges text inside bar to be white color
    """
    darkbarcolor = [
        (0.0, 0.0, 0.7529411764705882, 1.0), # equivalent to #410257
        (0.2549019607843137, 0.00784313725490196, 0.3411764705882353, 1.0), # equivalent to #0000c0
    ]
    palsD = {
        "m+0": "#410257",
        "m+1": "#0000c0",
        "m+2": "#55a0fb",
        "m+3": "#41f0ae",
        "m+4": "#addead",
        "m+5": "#eec046",
        "m+6": "#ffa040",
    }
    return darkbarcolor, palsD

def add_joker_tolabs(condilevels, levelshours_str):
    condilevels += ["joker"]  # joker to add space among time categories
    combined_tc_levels = list()
    for x in levelshours_str:
        for y in condilevels:
            if y == "joker":
                combined_tc_levels.append(str(x))
            else:
                combined_tc_levels.append(f'{x} : {y}')
    return condilevels, combined_tc_levels

def simplelabs(condilevels, levelshours_str):
    combined_tc_levels = list()
    for x in levelshours_str:
        for y in condilevels:
            combined_tc_levels.append(f'{x} : {y}')
    return condilevels, combined_tc_levels

def saveisotopologcontriplot(
    datadi,
    tablePicked,
    names_compartments,
    levelstimepoints_,
    namesuffix,
    metadata,
    selbycompD,
    darkbarcolor,
    palsD,
    condilevels,
):
    levelshours_str = [str(i) for i in sorted(metadata['timenum'].unique())]

    # condilevels, combined_tc_levels = add_joker_tolabs(condilevels, levelshours_str)
    condilevels, combined_tc_levels = simplelabs(condilevels, levelshours_str)

    for co in names_compartments.values():  #
        print(co)

        adf = pd.read_csv(
            datadi + tablePicked + "_" + namesuffix + "_" + co + ".tsv",
            sep="\t",
            index_col=0,
        )
        # note that pandas automatically transform any 99.9% in decimal 0.999

        dicos = dict()
        dicos[co] = {}
        dicos[co]["metadata"] = metadata.loc[metadata.short_comp == co]
        dicos[co][tablePicked] = adf[dicos[co]["metadata"]["sample"]]

        # call complicated functions

        df4plot = icontrib_2df4plot(dicos, tablePicked, co, levelshours_str)
        df4plot = massageisotopologues(df4plot)
        ####
        # conditions to plot in desired order :
        ####
        odiric = "results/plots/ic/"
        if not os.path.exists(odiric):
            os.makedirs(odiric)
        metscustomgroups = selbycompD[co]


        for j in range(len(metscustomgroups)):
            selectedmets = metscustomgroups[j]
            print(selectedmets)
            outfname = "{}ic_{}_group{}.pdf".format(odiric, co,  j)
            print(outfname)
            dfs_Dico = preparemeansreplicates( df4plot,  selectedmets )

            dfs_Dico = addcombinedconditime(dfs_Dico, combined_tc_levels)
            dfs_Dico.keys()  # just the metabolites subframes, one co
            figu_width = 4.3 * len(selectedmets)  # note, change width
            complexstacked(
                co, selectedmets, dfs_Dico, darkbarcolor, palsD, outfname,
                figu_width, xlabyesno="yes"
            )
            plt.close()
            # new :
            outfnameNoXlab = "{}ic_{}_group{}_noxlab.pdf".format(odiric, co,  j)
            complexstacked(
                co, selectedmets, dfs_Dico, darkbarcolor, palsD, outfnameNoXlab,
                figu_width, xlabyesno="no"
            )

        # legend alone
        plt.figure()
        myhandless = []
        for c in palsD.keys():
            paobj = mpatches.Patch(facecolor=palsD[c], label=c, edgecolor="black")
            myhandless.append(paobj)
        plt.legend(handles=myhandless)
        plt.savefig(f"{odiric}ic_legend.pdf", format="pdf")

    return 0

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 20:38:00 2022

@author: johanna
"""

import os

import numpy as np
import pandas as pd


def extrulist2dico(extrulist):
    extruD = dict()
    for i, r in extrulist.iterrows():
        metaboliteok = r["metabolite"].replace('"', "")
        if r["compartment"] not in extruD.keys():
            extruD[r["compartment"]] = [metaboliteok]
        else:
            extruD[r["compartment"]].append(metaboliteok)
    return extruD


def detectifinbadlist(al_, amet):
    amet = amet.replace("_PARENT", "")
    amet = amet.split("_C13")[0]

    if amet in al_:
        return 0
    else:
        return 1


def save_new_dfs(
    datadi, names_compartments, filename, metadata, extrulist_fi, diroutput
):
    extrulist = pd.read_csv(datadi + extrulist_fi, sep=",")
    extruD = extrulist2dico(extrulist)
    pa_i = pd.read_csv(datadi + filename, sep="\t", index_col=0)

    for compartment in names_compartments.values():
        f_o = filename.split(".")[0] + "_" + compartment + ".tsv"

        # using metadata, to set apart (pxpx) compartment specific samples
        samplestokeep = metadata.loc[metadata["short_comp"] == compartment, "sample"]
        pxpx = pa_i[samplestokeep]

        # create column to mark rows to delete
        # pxpx['todelete'] = np.nan
        pxpx = pxpx.assign(todelete=np.nan)
        for i, r in pxpx.iterrows():
            # print(i) # i is the i metabolite ....
            hereboolnum = detectifinbadlist(extruD[compartment], i)
            # print(hereboolnum)
            pxpx.loc[i, "todelete"] = hereboolnum
            # = hereboolnum # 0 or 1
        # end for iterrows
        # pxpx.todelete.value_counts() # visualize "table" on 'todelete' colummn
        pxpx_o = pxpx[pxpx["todelete"] != 0]
        pxpx_o = pxpx_o.drop(columns=["todelete"])
        pxpx_o.to_csv(diroutput + f_o, sep="\t")
    # end for
    return 0

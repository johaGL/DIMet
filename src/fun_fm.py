#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 31 12:24:01 2022

@author: johanna
"""

import os

import numpy as np
import pandas as pd

def autodetect_isotop_nomenclature(datadi, tableIC, namesuffix) :
    icdf = pd.read_csv(datadi + tableIC + "_" + namesuffix +".tsv", sep="\t", index_col=0)
    isotopololist_ = icdf.index.tolist()
    vibg1 = [i for i in isotopololist_ if  "_C13-label-" in i]
    vibg2 = [i for i in isotopololist_ if "_PARENT" in i]
    if (len(vibg1) + len(vibg2)) == len(isotopololist_):
        return "VIB"
    else:
        styuni = [ i for i in isotopololist_ if "_label" in i]
        if len(styuni) == len(isotopololist_):
            return "generic"
        else:
            return "style not recognized"


def yieldrowdataB(newdf):
    xu = {"metabolite": [], "m+x": [], "isotopolgFull": []}
    for ch in newdf.index:
        elems = ch.split("_m+")
        xu["metabolite"].append(elems[0])
        xu["m+x"].append("m+{}".format(elems[-1].split("-")[-1]))
        xu["isotopolgFull"].append(ch)
    rowdata = pd.DataFrame.from_dict(xu)
    return rowdata

def prepare4contrast(idf, ametadata, newcateg, contrast):
    """
    newcateg : example :  ['condition', 'timepoint' ]
    contrast : example : ["treatment_t12h", "control_t12h" ]
    creates column "newcol" : aa_bb <as in newcateg> suitable to contrast
    """
    cc = ametadata.copy()
    l_ = (ametadata[newcateg[0]] + "_" + ametadata[newcateg[1]]).tolist()
    cc["newcol"] = l_
    metas = cc.loc[cc["newcol"].isin(contrast), :]
    newdf = idf[metas["sample"]]
    return newdf, metas


def calcgeommean(avector, eps):
    vech = np.array(avector)
    vech[vech == 0] = eps  # replace any zeroes
    return np.exp(np.mean(np.log(vech)))


def splitrowbynewcol(row, metas):
    """
    Returns : miniD
    example : Control_T24h : [0.0, 0.0, 0.0] , L-Cycloserine_T24h : [0.0, 0.0, 0.0]
    """
    newcoluniq = list(set(metas["newcol"]))
    miniD = dict()
    for t in newcoluniq:
        # print(metas.loc[metas["newcol"] == t,:])
        koo = metas.loc[metas["newcol"] == t, :]
        selsams = koo["sample"]
        miniD[t] = row[selsams].tolist()
    return miniD


def jitterzero(avector):
    """only for internal tests, not a validated method"""
    if sum(avector) == 0:  # all vector is zero
        # avector = [7e-15,  4e-15,  9e-15, 5e-15] # if many small values
        avector = [1500000, 1600000, 1700000]  # if huge values
    return avector


def a12(lst1, lst2, rev=True):
    """
    Non-parametric hypothesis testing using Vargha and Delaney's A12 statistic:
    how often is x in lst1 greater than y in lst2?
    == > it gives a size effect, good to highlight potentially real effects <==
    """
    more = same = 0.0
    for x in lst1:
        for y in lst2:
            if x == y:
                same += 1
            elif rev and x > y:
                more += 1
            elif not rev and x < y:
                more += 1
    return (more + 0.5 * same) / (len(lst1) * len(lst2))

def detect_fraccontrib_missing(the_tsv_files_):
    fraccontribfile = ""
    for fi in the_tsv_files_:
        if "frac" in fi.lower():
            fraccontribfile = fi
    if fraccontribfile != "":
        return True
    return False


##################
# def yieldrowdata_old(newdf):
#     xu = {"metabolite": [], "m+x": [], "isotopolgFull": []}
#     for ch in newdf.index:
#         if "_C13-label-" in ch:
#             elems = ch.split("_C13-label-")
#             xu["metabolite"].append(elems[0])
#             xu["m+x"].append("m+{}".format(elems[-1].split("-")[-1]))
#             xu["isotopolgFull"].append(ch)
#         elif "_PARENT" in ch:
#             elems = ch.split("_PARENT")
#             xu["metabolite"].append(elems[0])
#             xu["m+x"].append("m+0")
#             xu["isotopolgFull"].append(ch)
#     rowdata = pd.DataFrame.from_dict(xu)
#     return rowdata
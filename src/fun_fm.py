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

def yieldrowdata(newdf):
    xu = {"metabolite": [], "m+x": [], "isotopolgFull": []}
    for ch in newdf.index:
        if "_C13-label-" in ch:
            elems = ch.split("_C13-label-")
            xu["metabolite"].append(elems[0])
            xu["m+x"].append("m+{}".format(elems[-1].split("-")[-1]))
            xu["isotopolgFull"].append(ch)
        elif "_PARENT" in ch:
            elems = ch.split("_PARENT")
            xu["metabolite"].append(elems[0])
            xu["m+x"].append("m+0")
            xu["isotopolgFull"].append(ch)
    rowdata = pd.DataFrame.from_dict(xu)
    return rowdata

def yieldrowdataB(newdf):
    xu = {"metabolite": [], "m+x": [], "isotopolgFull": []}
    for ch in newdf.index:
        elems = ch.split("_m+")
        xu["metabolite"].append(elems[0])
        xu["m+x"].append("m+{}".format(elems[-1].split("-")[-1]))
        xu["isotopolgFull"].append(ch)
    rowdata = pd.DataFrame.from_dict(xu)
    print(rowdata)
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


def fi2_smx_compartment(filename, indexsmx, indexcompartment):
    """example output: 'm+3' , 'cell'"""
    elems = filename.replace(".tsv", "").split("_")
    smx = elems[indexsmx]
    compartment = elems[indexcompartment]
    return smx, compartment


def reducebysd_2(avector, eps):
    """
    avector : list
    eps : value to replace zero
    returns :  numpy array
    """
    Y = np.array(avector)
    Y = np.where(Y == 0, eps, Y)
    stdevi = np.std(Y)
    if stdevi == 0:
        stdevi = eps
    # try:
    #     output = Y / np.std(Y)
    # except ZeroDivisionError():
    #     output = Y / 1e-05
    return Y / stdevi


def compute_reduction(df, ddof):
    """from proteomiX"""
    res = df.copy()
    for protein in df.index.values:
        # get array with abundances values
        protein_values = np.array(
            df.iloc[protein].map(lambda x: locale.atof(x) if type(x) == str else x)
        )

        # return array with each value divided by standard deviation of the whole array
        reduced_abundances = protein_values / np.nanstd(protein_values, ddof=ddof)

        # replace values in result df
        res.loc[protein] = reduced_abundances
    return res


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

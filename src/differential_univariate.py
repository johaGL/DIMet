#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 14:32:27 2022

@author: johanna
"""
import numpy as np
import scipy.stats
import statsmodels.stats.multitest as ssm
from .functions_diffmode import outStat_df, compute_padj_version2
from .fun_fm import *



def outFC_df(newdf, metas, contrast):
    """
    Calculates fold change
    returns : pandas dataframe
    """
    mets = []
    FC_geommu = []
    metaboliteshere = newdf.index
    for i in metaboliteshere:
        mets.append(i)
        row = newdf.loc[i,:]  # remember row is a seeries, colnames pass to index

        columnsInterest = metas.loc[metas["newcol"] == contrast[0], "sample"]
        columnsBaseline = metas.loc[metas["newcol"] == contrast[1], "sample"]

        vInterest = row[columnsInterest].to_numpy()
        vBaseline = row[columnsBaseline].to_numpy()

        # Calculate Fold Change
        #geommInterest = calcgeommean(vInterest, eps)
        #geommBaseline = calcgeommean(vBaseline, eps)

        val_interest = compute_gmean_nonan(vInterest)
        val_baseline = compute_gmean_nonan(vBaseline)
        if val_baseline == 0:
            ratioval = val_interest / 1e-10
        else:
            ratioval = val_interest / val_baseline

        FC_geommu.append(ratioval)


    fctab = pd.DataFrame(
        data={"metabolite": mets, "FC_geommu": FC_geommu, "log2FC": np.log2(FC_geommu)}
    )

    return fctab





def rundiffer(datadi, tablePicked, namesuffix, metadata, newcateg, contrast,
                   whichtest,  co, outdiffdirs, autochoice):
    """
    runs functions above,
    saves DAM (Differentially abundant metabolites/isotopologues)
    into tables
    """
    eps = 1
    if autochoice == "TOTAL":
        abun = pd.read_csv(
            datadi + tablePicked + "_" + namesuffix + "_" + co + ".tsv",
            sep="\t",
            index_col=0,
        )
        outkey = "TOTAL"
    else:
        abun = pd.read_csv(datadi + tablePicked, sep="\t", index_col=0)
        outkey = tablePicked.split("_")[2]  # the species m+x as saved
    print("\nprocessing : ", outkey, co)
    #tokeep = [i for i in abun.index if (i not in technical_toexclude)]
    #abun = abun.loc[tokeep, :]
    metadahere = metadata.loc[metadata.short_comp == co]
    selecols = metadahere["sample"]

    newdf_tot, metas_tot = prepare4contrast(abun, metadahere, newcateg, contrast)

    df = newdf_tot.copy()
    mets = df.index
    df.index = range(len(mets))  # to make a 1st column numeric positions, compute_reduction accepted
    ddof = 0
    tmp = compute_reduction(df, ddof)
    tmp.index = newdf_tot.index
    newdf_tot_red = tmp

    prediffresult = outStat_df(newdf_tot_red, metas_tot, contrast, whichtest)

    prediffresult = compute_padj_version2(prediffresult)

    FCresult = outFC_df(newdf_tot_red, metas_tot, contrast)

    DIFFRESULT = pd.concat(
        [prediffresult.set_index("metabolite"), FCresult.set_index("metabolite")],
        axis=1,
        join="inner").reset_index()

    ocols = ["metabolite", "stat", "p-value", "padj", "FC_geommu", "log2FC"]
    OUTPUT = DIFFRESULT[ocols]

    if outkey.startswith("m+") or (outkey =="mktot") :
        name_iso_long = [f"{m}_{outkey}" for m in OUTPUT.metabolite.tolist()]
        OUTPUT["metabolite"] = name_iso_long

    contrastword = "_".join(contrast)

    extended_dir = [ i for i in outdiffdirs if "extended" in i]
    signi_dir = [i for i in outdiffdirs if "significant" in i]

    ofi = f"{extended_dir[0]}{co}_{outkey}_{contrastword}_{whichtest}.tsv"
    OUTPUT.to_csv(ofi, sep="\t")

    sigoutput = OUTPUT.loc[OUTPUT["padj"] <= 0.05,:]
    if sigoutput.shape[0] > 0:
        ofisig = f"{signi_dir[0]}{co}_{outkey}_{contrastword}_{whichtest}_sig.tsv"
        sigoutput.to_csv(ofisig, sep="\t")
    return 0

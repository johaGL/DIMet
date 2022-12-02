#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 14:32:27 2022

@author: johanna
"""


import scipy.stats
import statsmodels.stats.multitest as ssm
import locale
from .fun_fm import *



def outFC_df(newdf, metas, contrast, eps):
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
        geommInterest = calcgeommean(vInterest, eps)
        geommBaseline = calcgeommean(vBaseline, eps)
        try:
            FC_geommu.append(geommInterest / geommBaseline)
        except ZeroDivisionError():
            FC_geommu.append(geommInterest / eps)

    fctab = pd.DataFrame(
        data={"mets": mets, "FC_geommu": FC_geommu, "log2FC": np.log2(FC_geommu)}
    )

    return fctab


def outStat_df(newdf, metas, contrast, whichtest):
    """
    Parameters
    ----------
    newdf : pandas
        1st element output of prepare4contrast.
    metas : pandas
        2nd element output of prepare4contrast..
    contrast : a list
    Returns
    -------
    DIFFRESULT : pandas
        THE PAIR-WISE DIFFERENTIAL ANALYSIS RESULTS.
    """
    mets = []
    stare = []
    pval = []
    FC_geommu = []
    metaboliteshere = newdf.index
    for i in metaboliteshere:
        mets.append(i)
        row = newdf.loc[
            i,
        ]  # remember row is a seeries, colnames pass to index

        columnsInterest = metas.loc[metas["newcol"] == contrast[0], "sample"]
        columnsBaseline = metas.loc[metas["newcol"] == contrast[1], "sample"]

        vInterest = row[columnsInterest].to_numpy()
        vBaseline = row[columnsBaseline].to_numpy()

        if whichtest == "MW":
            # vInterest = jitterzero(vInterest)
            # vBaseline = jitterzero(vBaseline)
            # Calculate Mann–Whitney U test (a.k.a Wilcoxon rank-sum test,
            # or Wilcoxon–Mann–Whitney test, or Mann–Whitney–Wilcoxon (MWW/MWU), )
            # DO NOT : use_continuity False AND method "auto" at the same time.
            # because "auto" will set continuity depending on ties and sample size.
            # If ties in the data  and method "exact" (i.e use_continuity False) pvalues cannot be calculated
            # check scipy doc
            usta, p = scipy.stats.mannwhitneyu(
                vInterest,
                vBaseline,
                # method = "auto",
                use_continuity=False,
                alternative="less",
            )
            usta2, p2 = scipy.stats.mannwhitneyu(
                vInterest,
                vBaseline,
                # method = "auto",
                use_continuity=False,
                alternative="greater",
            )
            usta3, p3 = scipy.stats.mannwhitneyu(
                vInterest,
                vBaseline,
                # method = "auto",
                use_continuity=False,
                alternative="two-sided",
            )

            # best (smaller pvalue) among all tailed tests
            pretups = [(usta, p), (usta2, p2), (usta3, p3)]
            tups = []
            for t in pretups:  # make list of tuples with no-nan pvalues
                if not np.isnan(t[1]):
                    tups.append(t)

            if len(tups) == 0:  # if all pvalues are nan assign two sided result
                tups = [(usta3, p3)]

            stap_tup = min(tups, key=lambda x: x[1])  # if any nan, will always pick nan as min
            stare.append(stap_tup[0])
            pval.append(stap_tup[1])

        elif whichtest == "Tt":
            tstav, pvaltv = scipy.stats.ttest_ind(vInterest, vBaseline, alternative="two-sided")
            stare.append(tstav)
            pval.append(pvaltv)
        # end if
    # end for
    prediffr = pd.DataFrame(data={"mets": mets, "stat": stare, "pval": pval})
    return prediffr


def detect_and_create_dir(namenesteddir):
    if not os.path.exists(namenesteddir):
        os.makedirs(namenesteddir)


def rundiffer(datadi, tablePicked, namesuffix, metadata, newcateg, contrast,
                   whichtest,  co, outdiffdirs, autochoice):
    """
    runs functions above,
    saves DAM (Differentially abundant metabolites/isotopologues)
    into tables
    """
    reduce = True
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
    if reduce:
        df = newdf_tot.copy()
        mets = df.index
        df.index = range(len(mets))  # to make a 1st column numeric positions, compute_reduction accepted
        ddof = 0
        tmp = compute_reduction(df, ddof)
        tmp.index = newdf_tot.index
        newdf_tot_red = tmp

    prediffresult = outStat_df(newdf_tot_red, metas_tot, contrast, whichtest)

    cpdiff = prediffresult.copy()
    cpdiff["pval"] = cpdiff[["pval"]].fillna(1)  # inspired from R documentation in p.adjust
    (sgs, corrP, _, _) = ssm.multipletests(cpdiff["pval"], method="fdr_bh")  # Benjamini-Hochberg
    cpdiff["padj"] = corrP
    truepadj = []
    for v, w in zip(prediffresult["pval"], cpdiff["padj"]):
        if np.isnan(v):
            truepadj.append(v)
        else:
            truepadj.append(w)

    prediffresult["padj"] = truepadj

    FCresult = outFC_df(newdf_tot, metas_tot, contrast, eps)

    DIFFRESULT = pd.concat(
        [prediffresult.set_index("mets"), FCresult.set_index("mets")],
        axis=1,
        join="inner").reset_index()

    ocols = ["mets", "stat", "pval", "padj", "FC_geommu", "log2FC"]
    OUTPUT = DIFFRESULT[ocols]
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

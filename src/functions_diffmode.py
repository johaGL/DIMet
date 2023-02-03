#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 14:32:27 2022

@author: johanna
"""


import pandas as pd
import statsmodels.stats.multitest as smm
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from .distrib_fit_fromProteomix import compute_z_score, find_best_distribution, compute_p_value
from .fun_fm import prepare4contrast, give_reduced_df, \
    give_coefvar_new, give_geommeans_new, give_ratios_df
import scipy.stats
import statsmodels.stats.multitest as ssm




def compute_overlap(df: pd.DataFrame, group1, group2, overlap_method: str) -> pd.DataFrame:
    # Credits: Claire Lescoat, Macha Nikolski, Benjamin Dartigues, Cedric Usureau, Aurélien Barré, Hayssam Soueidan
    for i in df.index.values:
        group1_values = np.array(group1.iloc[i])
        group2_values = np.array(group2.iloc[i])

        if overlap_method == "symmetric":
            df.loc[i, 'distance'] = overlap_symmetric(group1_values, group2_values)
        else:
            df.loc[i, 'distance'] = overlap_asymmetric(group1_values, group2_values)

    return df


def overlap_symmetric(x: np.array, y: np.array) -> int:
    # Credits: Claire Lescoat, Macha Nikolski, Benjamin Dartigues, Cedric Usureau, Aurélien Barré, Hayssam Soueidan
    a = [np.nanmin(x), np.nanmax(x)]
    b = [np.nanmin(y), np.nanmax(y)]

    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)

    overlap = np.nanmax([a[0], b[0]]) - np.nanmin([a[1], b[1]])
    return overlap


def overlap_asymmetric(x: np.array, y: np.array) -> int:
    # Credits: Claire Lescoat, Macha Nikolski, Benjamin Dartigues, Cedric Usureau, Aurélien Barré, Hayssam Soueidan
    # x is the reference group
    overlap = np.nanmin(y) - np.nanmax(x)
    return overlap



def compute_p_adjusted(df: pd.DataFrame, correction_method: str) -> pd.DataFrame:
    rej, pval_corr = smm.multipletests(df['pvalue'].values, alpha=float('0.05'), method=correction_method)[:2]
    df['padj'] = pval_corr
    return df



def divide_groups(df4c, metad4c, selected_contrast):
        """split into two df"""
        sc = selected_contrast
        sam0 = metad4c.loc[metad4c['newcol'] == sc[0], "sample"] # interest
        sam1 = metad4c.loc[metad4c['newcol'] == sc[1], "sample"] # control
        group_interest = df4c[sam0]
        group_control = df4c[sam1]
        group_interest.index = range(group_interest.shape[0])
        group_control.index = range(group_control.shape[0])

        return group_interest, group_control


def calc_reduction( df, metad4c, selected_contrast):
    def renaming_original_col_sams(df):
        newcols = ["copy_" + i for i in df.columns]
        df.columns = newcols
        # origdf.index = metabolites
        return df

    ddof = 0  # for compute reduction
    df4c = df[metad4c['sample']]

    df4c = give_reduced_df(df4c, ddof)

    df_orig_vals = renaming_original_col_sams(df[metad4c['sample']])

    df4c = pd.merge(df_orig_vals, df4c, left_index=True, right_index=True)

    return df4c

def calc_ratios(df4c, metad4c, selected_contrast):


    # df4c = pd.merge(df4c, cv_df, left_index=True, right_index=True)

    c_interest = selected_contrast[0]
    c_control = selected_contrast[1]
    # geometric means
    df4c, geominterest, geomcontrol = give_geommeans_new(df4c,
                                                         metad4c, 'newcol', c_interest, c_control)


    df4c = give_ratios_df(df4c, geominterest, geomcontrol)
    return df4c



def distance_or_overlap(df4c, metad4c, selected_contrast):
    # distance (syn: overlap)
    groupinterest, groupcontrol = divide_groups(df4c, metad4c, selected_contrast)
    rownames = df4c.index
    tmp_df = df4c.copy()
    tmp_df.index = range(len(rownames))
    tmp_df = compute_overlap(tmp_df, groupcontrol, groupinterest, "symmetric")
    tmp_df.columns = [*tmp_df.columns[:-1], "distance"]
    tmp_df.index = rownames
    df4c = tmp_df.copy()
    return df4c


def compute_span_incomparison(df, metadata, contrast):
    expected_samples = metadata.loc[metadata['newcol'].isin(contrast), "sample"]
    selcols_df = df[expected_samples].copy()
    for i in df.index.values:
        df.loc[i, 'span_allsamples'] = max(selcols_df.loc[i,:]) - min(selcols_df.loc[i,:])
    return df

def split_byalert_df(df):
    good_df = df.loc[df['alert'] == '', :]
    bad_df = df.loc[df['alert'] != '', :]

    return good_df, bad_df


def auto_detect_tailway(good_df, best_distribution, args_param):
    min_pval_ = list()
    for tail_way in ["two-sided", "right-tailed"]:
        tmp = compute_p_value(good_df, tail_way, best_distribution, args_param)

        min_pval_.append(tuple([tail_way, tmp["pvalue"].min()]))

    return min(min_pval_, key=lambda x: x[1])[0]


def compute_log2FC(df4c):
    for i, r in df4c.iterrows():
        df4c.loc[i,'log2FC'] = np.log(df4c.loc[i,'ratio'].to_numpy() , 2)

    return df4c


def outStat_df(newdf, metas, contrast, whichtest):
    """
    Parameters
    ----------
    newdf : pandas
        reduced dataframe
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

        vInterest = vInterest[~np.isnan(vInterest)]
        vBaseline = vBaseline[~np.isnan(vBaseline)]


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
            if (len(vInterest) >= 2) and (len(vBaseline) >= 2): # new: 2 or more samples required
                tstav, pvaltv = scipy.stats.ttest_ind(vInterest, vBaseline, alternative="two-sided")
            else:
                tstav = np.nan
                pvaltv = np.nan
            stare.append(tstav)
            pval.append(pvaltv)
        # end if
    # end for
    prediffr = pd.DataFrame(data={"metabolite": mets, "stat": stare, "p-value": pval})
    return prediffr


def compute_padj_version2(df):
    df["p-value"] = df[["p-value"]].fillna(1)  # inspired from R documentation in p.adjust
    (sgs, corrP, _, _) = ssm.multipletests(df["p-value"], method="fdr_bh")  # Benjamini-Hochberg
    df["padj"] = corrP
    truepadj = []
    for v, w in zip(df["p-value"], df["padj"]):
        if np.isnan(v):
            truepadj.append(v)
        else:
            truepadj.append(w)
    df["padj"] = truepadj

    return df


if __name__ == "__main__":

    myconff = "configs/cmyelin.yml"
    import yaml
    os.chdir(os.path.expanduser("~/myelin_metabo/"))
    with open(myconff, "r") as f:
        confidic = yaml.load(f, Loader=yaml.Loader)

    namesuffix = confidic["namesuffix"]
    datadi = confidic["datadi"]
    extrudf_fi = confidic["extrulist_fi"]
    names_compartments = confidic["names_compartments"]
    metadata_fi = confidic["metadata_fi"]
    levelstimepoints_ = confidic["levelstime"]

    whichtest = confidic["whichtest"]

    metadata = pd.read_csv(datadi + metadata_fi, index_col=False)

    tableIC = confidic["name_isotopologue_contribs"].split(".")[0]  # no extension

    tableAbund = confidic["name_abundances"].split(".")[0]  # no extension
    max_m_species = confidic["max_m_species"]

    dirtmpdata = "tmp/"
    abunda_species_4diff_dir = dirtmpdata + "abufromperc/"

    ### functions specific for differential abundance by disfit

    spefiles = [i for i in os.listdir(abunda_species_4diff_dir)]

    newcateg = confidic["newcateg"]
    contrasts_ = confidic["contrasts"]
    autochoice = "TOTAL"

    # outdiffdir = "results/tables/"
    # if not os.path.exists(outdiffdir):
    #     os.makedirs(outdiffdir)
    # outputsubdirs = ["m+" + str(i) + "/" for i in range(max_m_species + 1)]
    # outputsubdirs.append("totmk/")
    # outputsubdirs.append("TOTAL/")
    # alloutdirs = list()
    # for exte_sig in ["extended/", "significant/"]:
    #     for subdir_spec in outputsubdirs:
    #         x = outdiffdir + exte_sig + subdir_spec
    #         alloutdirs.append(x)
    #         if not os.path.exists(x):
    #             os.makedirs(x)

    print(newcateg)

    # print a table with reduded values by time
    # red_cv_g_bytime_dffull(names_compartments, dirtmpdata, tableAbund,
    #                   namesuffix, metadata, levelstimepoints_)

    # new loop for plotting?


    #use contrasts for selecting the data and calculate all
    selected_contrast = ['Myelin_T0h', 'Vehicle_T0h']
    #assert selected_contrast[1] == selected_contrast[1], f"error, geomean table wont be find, check {newcateg}"

    for co in names_compartments.values():
        if autochoice == "TOTAL":
            filehere = f"{dirtmpdata}{tableAbund}_{namesuffix}_{co}.tsv"
        else:
            filehere = f"{dirtmpdata}IFORGOTHERE{autochoice}_{co}.tsv"

        ratiosdf = calcs_before_fit(filehere, co,
                          metadata, newcateg, selected_contrast)


        strcontrast = "_vs_".join(selected_contrast)
        out_histo_file = f"results/plots/distrib-{strcontrast}-{co}.pdf"
        print(out_histo_file)

        print("fitting to distributions to find the best ... ")
        # **
        ratiosdf.to_csv(f"{autochoice}{strcontrast}_{co}_prep.tsv",
                           header=True, sep='\t')


        # **
        # best_distribution, args_param = find_best_distribution(ratiosdf,
        #                                                        out_histogram_distribution= out_histo_file)
        # #argsided = "right-tailed" #"two-sided" # or "rigth-tailed"
        # #ratiosdf = compute_p_value(ratiosdf, argsided, best_distribution, args_param )
        # #ratiosdf.to_csv("test_right.csv",header=True)
        #
        # argsided = "two-sided"  # "two-sided" # or "rigth-tailed"
        # ratiosdf = compute_p_value(ratiosdf, argsided, best_distribution, args_param)
        #
        # ratiosdf = compute_p_adjusted(ratiosdf, "fdr_bh")
        #

        # ratiosdf.to_csv(f"{outdiffdir}{autochoice}{strcontrast}_{co}_twosided.tsv",
        #                   header=True, sep='\t')



# end
# ----- older tests
# for co in names_compartments.values():
    # for t in ['T0h']:
    #     outfi_geomean = f"{dirtmpdata}abund_reduced_geomean_{t}_{co}.tsv"
    #     df_t_red_geomean = pd.read_csv(outfi_geomean, index_col=0)
    #     c_interest = selected_contrast[0] + "_geomean"
    #     c_control = selected_contrast[1] +"_geomean"
    #     ratiosdf = give_ratios_df(df_t_red_geomean, c_interest, c_control)
    # ...




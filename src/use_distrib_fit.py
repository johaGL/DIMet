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




def detect_and_create_dir(namenesteddir): # TODO: deduplicate, same in other scripts
    if not os.path.exists(namenesteddir):
        os.makedirs(namenesteddir)


def compute_overlap(df: pd.DataFrame, group1, group2, overlap_method: str) -> pd.DataFrame:
    # Credits: Claire Lescoat, Macha Nikolski, Benjamin Dartigues, Cedric Usureau, Aurélien Barré, Hayssam Soueidan
    for i in df.index.values:
        group1_values = np.array(group1.iloc[i])
        group2_values = np.array(group2.iloc[i])

        if overlap_method == "symmetric":
            df.loc[i, 'score_overlap'] = overlap_symmetric(group1_values, group2_values)
        else:
            df.loc[i, 'score_overlap'] = overlap_asymmetric(group1_values, group2_values)

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



def plot_overlap_hist(df_overls, fileout):
    values_sym = df_overls["score_overlap_symmetric"]
    a = pd.DataFrame({'value' : values_sym,
                      'type_overlap': ["symm" for i in range(len(values_sym))] })
    vasym = df_overls["s_o_Assymetric"]
    b = pd.DataFrame({'value': vasym,
                      'type_overlap': ["assym" for i in range(len(vasym))]})
    dfplotov = pd.concat([a,b], ignore_index=True, axis=0)

    #fig, ax = plt.subplots()
    with sns.axes_style("darkgrid"):
        sns.displot(data=dfplotov, x = 'value', hue='type_overlap',
                       kde=False)
        plt.savefig(fileout)
    plt.close()
    return 0


def calcs_red_to_ratios( df4c, co,
                      metad4c, newcateg, selected_contrast):

    ddof = 0  # for compute reduction
    df4c = df4c[metad4c['sample']]

    df4c_red = give_reduced_df(df4c, ddof)

    cv_df = give_coefvar_new(df4c_red, metad4c, 'newcol')

    df4c_red_cv = pd.merge(df4c_red, cv_df, left_index=True, right_index=True)

    groupinterest, groupcontrol = divide_groups(df4c_red_cv, metad4c, selected_contrast)

    # overlap # TODO : pick symmetric or assymetric
    rownames = df4c_red_cv.index
    df4c_red_cv.index = range(len(rownames))
    o_sym = compute_overlap(df4c_red_cv, groupcontrol, groupinterest, "symmetric")
    o_sym.columns = [*o_sym.columns[:-1], "score_overlap_symmetric"]
    o_sym = compute_overlap(o_sym,  groupcontrol, groupinterest , "assymmetric")
    o_sym.columns = [*o_sym.columns[:-1], "s_o_Assymetric"]
    o_sym.index = rownames

    df4c_red_cv_o = o_sym.copy()

    c_interest = selected_contrast[0]
    c_control = selected_contrast[1]

    # geometric means
    df4c_red_cv_o_g, geominterest, geomcontrol = give_geommeans_new(df4c_red_cv_o,
                                         metad4c, 'newcol', c_interest, c_control)
    #df4c_red_cv_o_g = pd.merge(df4c_red_cv_o, geomdf, left_index=True, right_index=True)


    # ratio (adds this column to df in construction)
    ratiosdf = give_ratios_df(df4c_red_cv_o_g, geominterest, geomcontrol)

    return ratiosdf



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




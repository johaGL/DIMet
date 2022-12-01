#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 14:32:27 2022

@author: johanna
"""


import locale
from scipy.stats.mstats import gmean
from distrib_fit_fromProteomix import compute_z_score, find_best_distribution, compute_p_value
from fun_fm import prepare4contrast

def compute_reduction(df, ddof):
    """
    modified, original from ProteomiX
    johaGL 2022: if all row is zeroes, set same protein_values
    """
    res = df.copy()
    for protein in df.index.values:
        # get array with abundances values
        protein_values = np.array(
            df.iloc[protein].map(lambda x: locale.atof(x) if type(x) == str else x) )
        # return array with each value divided by standard deviation of the whole array
        if np.nanstd(protein_values, ddof=ddof) == 0:
            reduced_abundances = protein_values  # because all row is zeroes
        else:
            reduced_abundances = protein_values / np.nanstd(protein_values, ddof=ddof)

        # replace values in result df
        res.loc[protein] = reduced_abundances
    return res


def compute_cv(reduced_abund):
    reduced_abu_np = reduced_abund.to_numpy().astype('float64')
    if np.nanmean(reduced_abu_np) != 0:
        return np.nanstd(reduced_abu_np) / np.nanmean(reduced_abu_np)
    elif np.nanmean(reduced_abu_np) != 0 and np.nanstd(reduced_abu_np) == 0:
        return 0
    else:
        return np.nan



def give_coefvar_bycond(df_red, red_meta, t):
    condilist = red_meta["condition"].unique()
    tmpdico = dict()
    for condi in condilist:
        samplesthiscondi = red_meta.loc[red_meta['condition']==condi, "sample"]
        subdf = df_red[samplesthiscondi]
        subdf = subdf.assign(CV=subdf.apply(compute_cv, axis=1))
        #print(subdf.head())
        tmpdico[f"{condi}_{t}_CV"] = subdf.CV.tolist()

    dfout = pd.DataFrame.from_dict(tmpdico)
    dfout.index = df_red.index
    return dfout

def give_coefvar_new(df_red, red_meta, newcol : str):
    groups_ = red_meta[newcol].unique()

    tmpdico = dict()
    for group in groups_:
        samplesthisgroup = red_meta.loc[red_meta[newcol] == group, "sample"]
        subdf = df_red[samplesthisgroup]
        subdf = subdf.assign(CV=subdf.apply(compute_cv, axis=1))
        #print(subdf.head())
        tmpdico[f"{group}_CV"] = subdf.CV.tolist()

    dfout = pd.DataFrame.from_dict(tmpdico)
    dfout.index = df_red.index
    return dfout

def give_geommeans_new(df_red, red_meta, newcol : str):
    groups_ = red_meta[newcol].unique()
    tmpdico = dict()
    for group in groups_:
        samplesthisgroup = red_meta.loc[red_meta[newcol] == group, "sample"]
        subdf = df_red[samplesthisgroup]
        subdf = subdf.assign(geomean=subdf.apply(gmean, axis=1))
        # print(subdf.head())
        tmpdico[f"{group}_gm"] = subdf.geomean.tolist()

    dfout = pd.DataFrame.from_dict(tmpdico)
    dfout.index = df_red.index
    return dfout


def give_geommeans(df_red, red_meta, t ):
    condilist = red_meta["condition"].unique()
    tmpdico = dict()
    for condi in condilist:
        samplesthiscondi = red_meta.loc[red_meta['condition'] == condi, "sample"]
        subdf = df_red[samplesthiscondi]
        subdf = subdf.assign(geomean=subdf.apply(gmean, axis=1))
        # print(subdf.head())
        tmpdico[f"{condi}_{t}_geomean"] = subdf.geomean.tolist()

    dfout = pd.DataFrame.from_dict(tmpdico)
    dfout.index = df_red.index

    return dfout

def give_ratios_df(df, conditionInterest, conditionControl):
    df = df.assign(ratio=df[conditionInterest] / df[conditionControl])
    return df


def detect_and_create_dir(namenesteddir): # TODO: deduplicate, same in other scripts
    if not os.path.exists(namenesteddir):
        os.makedirs(namenesteddir)


def test_bytime_dffull(names_compartments, dirtmpdata, tableAbund,
                       namesuffix, metadata, levelstimepoints_):
    # calculate and save : reduced data, coefficient of variation, splitting by timepoint, here only T0h test
    ddof = 0
    for co in names_compartments.values():
        df = pd.read_csv(f"{dirtmpdata}{tableAbund}_{namesuffix}_{co}.tsv", sep='\t', header=0, index_col=0)

        metada_sel = metadata.loc[metadata['short_comp']==co, :]
        #get reduced rows , cv and geommeans,
        for t in ['T0h']: # for t in levelstimepoints_
            print(t)
            samples_t = metada_sel.loc[metada_sel['timepoint'] == t, "sample"]
            samples_t = sorted(samples_t)
            df_t = df[samples_t]
            rownames = df_t.index
            df_t.index = range(len(rownames))#  index must be numeric because compute reduction accepted
            df_t_red = compute_reduction(df_t, ddof)
            df_t_red.index = rownames
            #outfilereduced = f"{dirtmpdata}abund_reduced_{t}_{co}.tsv"
            # df_t_red.to_csv(outfilereduced, header=True, sep='\t')
            # add coefficient of variation, by condition
            red_meta = metada_sel.loc[metada_sel["sample"].isin(df_t_red.columns) ,:]


            df_cv = give_coefvar_bycond(df_t_red, red_meta, t )
            df_t_red_cv = pd.merge(df_t_red, df_cv , left_index=True, right_index=True)
            #outfi_coefvar = f"{dirtmpdata}abund_reduced_coefvar_{t}_{co}.tsv"
            #df_t_red_cv.to_csv(outfi_coefvar)

            # save intervals overlap # TODO?... but here we have 3 conditions (in this timepoint)

            # save geometric means table, and the ratio
            df_t_red_geomean = give_geommeans(df_t_red, red_meta, t)
            print(df_t_red_geomean.head())
            dfo1 = pd.merge(df_t_red_cv, df_t_red_geomean, left_index=True, right_index=True)
            outfi_geomean = f"{dirtmpdata}abund_reduced_geomean_{t}_{co}.tsv"
            dfo1.to_csv(outfi_geomean, header=True)
    return 0

if __name__ == "__main__":

    myconff = "configs/cmyelin.yml"
    import yaml
    import os
    import pandas as pd
    import numpy as np
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

    ddof = 0 # for compute reduction


    outdiffdir = "results/tables/"
    if not os.path.exists(outdiffdir):
        os.makedirs(outdiffdir)
    outputsubdirs = ["m+" + str(i) + "/" for i in range(max_m_species + 1)]
    outputsubdirs.append("totmk/")
    outputsubdirs.append("TOTAL/")
    alloutdirs = list()
    # for exte_sig in ["extended/", "significant/"]:
    #     for subdir_spec in outputsubdirs:
    #         x = outdiffdir + exte_sig + subdir_spec
    #         alloutdirs.append(x)
    #         if not os.path.exists(x):
    #             os.makedirs(x)

    print(newcateg)
    #test_bytime_dffull(names_compartments, dirtmpdata, tableAbund,
    #                   namesuffix, metadata, levelstimepoints_)

    # new loop for plotting
    # ...


    #use contrasts for selecting the data and calculate all
    selected_contrast = ['Myelin_T0h', 'Vehicle_T0h']
    #assert selected_contrast[1] == selected_contrast[1], f"error, geomean table wont be find, check {newcateg}"
    for co in names_compartments.values():
        # for t in ['T0h']:
        #     outfi_geomean = f"{dirtmpdata}abund_reduced_geomean_{t}_{co}.tsv"
        #     df_t_red_geomean = pd.read_csv(outfi_geomean, index_col=0)
        #     c_interest = selected_contrast[0] + "_geomean"
        #     c_control = selected_contrast[1] +"_geomean"
        #     ratiosdf = give_ratios_df(df_t_red_geomean, c_interest, c_control)

        df = pd.read_csv(f"{dirtmpdata}{tableAbund}_{namesuffix}_{co}.tsv", sep='\t', header=0, index_col=0)

        metada_sel = metadata.loc[metadata['short_comp'] == co, :]

        df4c, metad4c = prepare4contrast(df, metada_sel, newcateg, selected_contrast)

        # sort them by condition
        metad4c = metad4c.sort_values("condition")
        df4c = df4c[metad4c['sample']]

        c_interest = selected_contrast[0]
        c_control = selected_contrast[1]

        rownames = df4c.index
        df4c.index = range(len(rownames))  # index must be numeric because compute reduction accepted
        df4c_red = compute_reduction(df4c, ddof)
        df4c_red.index = rownames
        print(metad4c.columns)
        cv_df = give_coefvar_new(df4c_red, metad4c, 'newcol')
        df4c_red_cv = pd.merge(df4c_red, cv_df, left_index=True, right_index=True)

        #overlap # TODO

        geomdf = give_geommeans_new(df4c_red_cv, metad4c, 'newcol')
        df4c_red_cv_g = pd.merge(df4c_red_cv, geomdf, left_index=True, right_index=True)

        interest_gm = c_interest + "_gm"
        control_gm = c_control + "_gm"
        ratiosdf = give_ratios_df(df4c_red_cv_g, interest_gm, control_gm )

        ratiosdf = compute_z_score(ratiosdf)
        print(ratiosdf.head())
        ratiosdf.to_csv( f"xx_{co}",header=True)

        print("fitting to distributions to find the best ... ")
        import sys
        sys.exit()
        strcontrast = "_vs_".join(selected_contrast)
        out_histo_file = f"results/plots/distrib-{strcontrast}-{co}.pdf"
        best_distribution, args_param = find_best_distribution(ratiosdf,
                                                               out_histogram_distribution= out_histo_file)
        argsided = "right-tailed" #"two-sided" # or "rigth-tailed"
        ratiosdf = compute_p_value(ratiosdf, argsided, best_distribution, args_param )
        ratiosdf.to_csv("test_right.csv",header=True)

        argsided = "two-sided"  # "two-sided" # or "rigth-tailed"
        ratiosdf = compute_p_value(ratiosdf, argsided, best_distribution, args_param)
        ratiosdf.to_csv("test_twosided.csv", header=True)







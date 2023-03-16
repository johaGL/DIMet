#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 14:32:27 2022

@author: johanna
"""
import os
import sys
import argparse
import pandas as pd
import numpy as np
import scipy.stats
import statsmodels.stats.multitest as ssm
import matplotlib.pyplot as plt
import seaborn as sns
sys.path.append(os.path.dirname(__file__))
import functions_general as fg
from distrib_fit_fromProteomix import compute_z_score, find_best_distribution, compute_p_value


    #prepare4contrast, give_reduced_df,
   # give_coefvar_new, give_geommeans_new, give_ratios_df



def diff_args():
    parser = argparse.ArgumentParser(prog="python -m DIMet.src.differential_analysis",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('config', type=str,
                        help="configuration file in absolute path")

    parser.add_argument('--abundances', action=argparse.BooleanOptionalAction, default=True)

    parser.add_argument('--meanEnrichment_or_fractionalContribC',
                        action=argparse.BooleanOptionalAction, default=True)

    parser.add_argument('--isotopologues', action=argparse.BooleanOptionalAction, default=True)

    parser.add_argument("--qualityDistanceOverSpan", default=-0.3, type=float,
                        help="For each metabolite, for samples (x, y)  the\
                            distance is calculated, and span is max(x U y) - min(x U y).  A 'distance/span' \
                            inferior to this value excludes the metabolite from testing (pvalue=NaN).")


    return parser


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


def compute_padj_version2(df, correction_alpha, correction_method):
    tmp = df.copy()
    tmp["pvalue"] = tmp[["pvalue"]].fillna(1)  # inspired from R documentation in p.adjust

    (sgs, corrP, _, _) = ssm.multipletests(tmp["pvalue"], alpha=float(correction_alpha),
                                           method=correction_method)
    df["padj"] = corrP
    truepadj = []
    for v, w in zip(df["pvalue"], df["padj"]):
        if np.isnan(v):
            truepadj.append(v)
        else:
            truepadj.append(w)
    df["padj"] = truepadj

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


def calc_reduction(df, metad4c, selected_contrast):
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
    c_interest = selected_contrast[0]
    c_control = selected_contrast[1]
    # geometric means
    df4c, geominterest, geomcontrol = give_geommeans_new(df4c, metad4c,
                                                         'newcol', c_interest, c_control)

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


def auto_detect_tailway(good_df, best_distribution, args_param):
    min_pval_ = list()
    for tail_way in ["two-sided", "right-tailed"]:
        tmp = compute_p_value(good_df, tail_way, best_distribution, args_param)

        min_pval_.append(tuple([tail_way, tmp["pvalue"].min()]))

    return min(min_pval_, key=lambda x: x[1])[0]


def steps_fitting_method(ratiosdf, out_histo_file):
    ratiosdf = compute_z_score(ratiosdf)
    best_distribution, args_param = find_best_distribution(ratiosdf,
                                 out_histogram_distribution=out_histo_file)
    autoset_tailway = auto_detect_tailway(ratiosdf, best_distribution, args_param)
    print("auto, best pvalues calculated :", autoset_tailway)
    ratiosdf = compute_p_value(ratiosdf, autoset_tailway, best_distribution, args_param)

    return ratiosdf


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

        vInterest = np.array(row[columnsInterest], dtype=float)
        vBaseline = np.array(row[columnsBaseline], dtype=float)

        vInterest = vInterest[~np.isnan(vInterest)]
        vBaseline = vBaseline[~np.isnan(vBaseline)]

        if (len(vInterest) >= 2) and (len(vBaseline) >= 2):  #  2 or more samples required

            if whichtest == "MW":
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
                stat_result = stap_tup[0]
                pval_result = stap_tup[1]

            elif whichtest == "Tt":
                stat_result, pval_result = scipy.stats.ttest_ind(vInterest, vBaseline, alternative="two-sided")

            elif whichtest == "K":
                stat_result, pval_result = scipy.stats.kruskal(vInterest, vBaseline)


            stare.append(stat_result)
            pval.append(pval_result)
        else:
            stare.append(np.nan)
            pval.append(np.nan)

        # end if

    # end for
    prediffr = pd.DataFrame(data={"metabolite": mets, "stat": stare, "pvalue": pval})
    return prediffr


def flag_has_replicates(ratiosdf):
    bool_results = list()
    tuples_list = ratiosdf['count_nan_samples'].tolist()
    for tup in tuples_list:
        group_x = tup[0].split("/")
        group_y = tup[1].split("/")
        x_usable = int(group_x[1]) - int(group_x[0])
        y_usable = int(group_y[1]) - int(group_y[0])
        if x_usable <= 1 or y_usable <= 1 : # if any side only has one replicate
            bool_results.append(0)
        else:
            bool_results.append(1)
    return bool_results


def divide_before_stats(ratiosdf, quality_dist_span):
    ratiosdf['has_replicates'] = flag_has_replicates(ratiosdf)
    m = 0
    try:
        quality_dist_span = float(quality_dist_span)
        good_df = ratiosdf.loc[(ratiosdf['distance/span'] >= quality_dist_span) &
                        (ratiosdf['has_replicates'] == 1) ]

        undesired_mets = set(ratiosdf.index) - set(good_df.index)
        bad_df = ratiosdf.loc[list(undesired_mets)]
        good_df = good_df.drop(columns=['has_replicates'])
        bad_df = bad_df.drop(columns=['has_replicates'])
    except Exception as e:
        print(e)
        print("Unknown error in divide_before_stats, check qualityDistanceOverSpan arg")

    return good_df, bad_df


def complete_columns_for_bad(bad_df, ratiosdf):
    columns_missing =  set(ratiosdf.columns) - set(bad_df.columns)
    for col in list(columns_missing):
        bad_df[col] = np.nan
    return bad_df


def filter_diff_results(ratiosdf, padj_cutoff, log2FC_cutoff, distance_span):
    ratiosdf['abslfc'] = ratiosdf['log2FC'].abs()
    ratiosdf = ratiosdf.loc[( ratiosdf['padj'] <= padj_cutoff ) &
                        ( ratiosdf['abslfc'] >= log2FC_cutoff ) &
                        ( ratiosdf['distance/span'] >= distance_span ), :]
    ratiosdf = ratiosdf.drop(columns=['abslfc'])
    return ratiosdf



def save_each_df(good_df, bad_df, outdiffdir,
                 co, autochoice, strcontrast):
    rn = f"{outdiffdir}/extended/{autochoice}"
    good_o = f"{rn}/{co}_{autochoice}_{strcontrast}_good.tsv"
    bad_o = f"{rn}/{co}_{autochoice}_{strcontrast}_bad.tsv"
    good_df.to_csv(good_o, sep='\t', header=True)
    bad_df.to_csv(bad_o, sep='\t', header=True)
    return "saved to results"

# --

def whatever_the_table(df, metadatadf, out_file_name, confidic, args):
    quality_dist_span = args.qualityDistanceOverSpan
    print(quality_dist_span)
    huge_df = 2
    return huge_df


def blah_abund(clean_tables_path, table_prefix, metadatadf, out_diff_abun, confidic, args):
    test_abu = confidic['statistical_test']['abundances']
    suffix = confidic['suffix']
    compartments = metadatadf['short_comp'].unique().tolist()
    # dynamically open the file based on prefix, compartment and suffix:
    for co in compartments:
        meta_sub = metadatadf.loc[metadatadf['short_comp'] == co, :]
        fn = f'{clean_tables_path}{table_prefix}--{co}--{suffix}.tsv'
        measurements = pd.read_csv(fn, sep='\t', header=0, index_col=0)
        df, meta_sub = fg.prepare4contrast(measurements, meta_sub,
                                           confidic['grouping'], confidic['contrasts'])
    out_file_name = 44
    huge_df = whatever_the_table(df, meta_sub, out_file_name, confidic, args )
    return 0



if __name__ == "__main__":
    print("\n  -*- searching for Differentially Abundant-or-Marked Metabolites (DAM) -*-\n")
    parser = diff_args()
    args = parser.parse_args()
    configfile = os.path.expanduser(args.config)
    confidic = fg.open_config_file(configfile)
    fg.auto_check_validity_configuration_file(confidic)
    out_path = os.path.expanduser(confidic['out_path'])
    meta_path = os.path.expanduser(confidic['metadata_path'])
    clean_tables_path = out_path + "results/prepared_tables/"

    # tpd : tables prefixes dictionary
    tpd = fg.clean_tables_names2dict(f'{clean_tables_path}TABLESNAMES.csv')
    metadatadf = fg.open_metadata(meta_path)


    print(confidic['contrasts'])
    # abu
    abund_tab_prefix = tpd['name_abundance']
    out_diff_abun = out_path + "results/differential_analysis/abundance/"
    fg.detect_and_create_dir(out_diff_abun)
    blah_abund(clean_tables_path, abund_tab_prefix, metadatadf, out_diff_abun, confidic, args)


    # ME or FC
    fraccon_tab_prefix = tpd['name_meanE_or_fracContrib']

    # ii
    isos_abs_tab_prefix = tpd['name_isotopologue_abs']
    if isos_abs_tab_prefix == "None":  # is string because prepare.py did it
        print("fun in prop")
        # run here the function on proportions
        isos_prop_tab_prefix = tpd['name_isotopologue_prop']
    else:
        print("fun on absolutes")
        # run here the function on absolute vals

        # for each of these make a separate function :
    #1 -abundance : allows all tests
    #2 -ME or FC : allows all BUT not T-test
    #3 -Isotopologues:
        # -- absolute ( if available) allows all, otherwise:
        #  - - isotopologues percentages  allows all BUT not T-test


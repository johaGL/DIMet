#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 14:32:27 2022

@author: johanna
"""
import os
import argparse
import pandas as pd
import numpy as np
import scipy.stats
import statsmodels.stats.multitest as ssm
import functions_general as fg
from distrib_fit_fromProteomix import compute_z_score, find_best_distribution, compute_p_value


def diff_args():
    parser = argparse.ArgumentParser(prog="python -m DIMet.src.differential_analysis",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('config', type=str,
                        help="configuration file in absolute path")

    # Before reduction and gmean, way to replace zeros:
    parser.add_argument("--abundance_replace_zero_with", default="min", type=str, help="min | min/n | VALUE")

    parser.add_argument( "--meanEnrichOrFracContrib_replace_zero_with", metavar="MEorFC_replace_zero_with",
                         default="min", type=str,
                         help="min | min/n | VALUE")

    parser.add_argument("--isotopologueProp_replace_zero_with", default="min", type=str,
                        help="min | min/n | VALUE")

    parser.add_argument("--isotopologueAbs_replace_zero_with", default="min", type=str,
                        help="min | min/n | VALUE")

    # by default include all the types of measurements:

    parser.add_argument('--abundances', action=argparse.BooleanOptionalAction, default=True)

    parser.add_argument('--meanEnrich_or_fracContrib',
                        action=argparse.BooleanOptionalAction, default=True)

    parser.add_argument('--isotopologues', action=argparse.BooleanOptionalAction, default=True)


    parser.add_argument("--qualityDistanceOverSpan", default=-0.3, type=float,
                        help="For each metabolite, for samples (x, y)  the\
                            distance is calculated, and span is max(x U y) - min(x U y).  A 'distance/span' \
                            inferior to this value excludes the metabolite from testing (pvalue=NaN).")

    parser.add_argument("--multitest_correction", default='fdr_bh',
                        help="see https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html")
    return parser


def validate_zero_repl_arg(zero_repl_arg: str) -> None:
    isvalid = False
    if zero_repl_arg == "min":
        isvalid = True
    elif zero_repl_arg.startswith("min/"):
        try:
            float(str(zero_repl_arg.split("/")[1]))
            isvalid = True
        except ValueError:
            pass
    else:
        try:
            float(args.zero_repl_arg)
            isvalid = True
        except ValueError:
            pass

    assert isvalid, "Your value for zero replacement is not valid"


def arg_repl_zero2value(argum_zero_rep: str, df: pd.DataFrame) -> float:
    if argum_zero_rep.startswith("min"):
        try:
            divisor_min = float(argum_zero_rep.split("/")[1])
        except IndexError:
            divisor_min = 1
        value_out = df[df > 0].min().min() / divisor_min
    else:
        value_out = float(args.zero_repl_arg)
    return value_out


def flag_has_replicates(ratiosdf: pd.DataFrame):
    bool_results = list()
    tuples_list = ratiosdf['count_nan_samples'].tolist()
    for tup in tuples_list:
        group_x = tup[0].split("/")
        group_y = tup[1].split("/")
        x_usable = int(group_x[1]) - int(group_x[0])
        y_usable = int(group_y[1]) - int(group_y[0])
        if x_usable <= 1 or y_usable <= 1:  # if any side only has one replicate
            bool_results.append(0)
        else:
            bool_results.append(1)
    return bool_results


def compute_span_incomparison(df: pd.DataFrame, metadata: pd.DataFrame, contrast: list):
    expected_samples = metadata.loc[metadata['newcol'].isin(contrast), 'name_to_plot']
    selcols_df = df[expected_samples].copy()
    for i in df.index.values:
        values_this_comparison = selcols_df.loc[i,:].to_numpy()
        df.loc[i, 'span_allsamples'] = values_this_comparison.max() - values_this_comparison.min()

    return df


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


def calc_reduction(df, metad4c):
    def renaming_original_col_sams(df):
        newcols = ["input_" + i for i in df.columns]
        df.columns = newcols

        return df

    ddof = 0  # for compute reduction
    df4c = df[metad4c['name_to_plot']]

    df4c = fg.give_reduced_df(df4c, ddof)

    df_orig_vals = renaming_original_col_sams(df[metad4c['name_to_plot']])

    df4c = pd.merge(df_orig_vals, df4c, left_index=True, right_index=True)

    return df4c


def calc_ratios(df4c, metad4c, selected_contrast):
    c_interest = selected_contrast[0]  # columns names interest
    c_control = selected_contrast[1]  # columns names control

    df4c, col_g_interest, col_g_control = fg.give_geommeans_new(df4c, metad4c,
                                                         'newcol', c_interest, c_control)
    df4c = fg.give_ratios_df(df4c, col_g_interest, col_g_control)

    return df4c


def divide_groups(df4c, metad4c, selected_contrast):
    """split into two df"""
    sc = selected_contrast
    sam0 = metad4c.loc[metad4c['newcol'] == sc[0], 'name_to_plot']  # interest
    sam1 = metad4c.loc[metad4c['newcol'] == sc[1], 'name_to_plot']  # control
    group_interest = df4c[sam0]
    group_control = df4c[sam1]
    group_interest.index = range(group_interest.shape[0])
    group_control.index = range(group_control.shape[0])

    return group_interest, group_control


def distance_or_overlap(df4c, metad4c, selected_contrast):
    """ calculate distance between groups (synonym: overlap) """
    groupinterest, groupcontrol = divide_groups(df4c, metad4c, selected_contrast)
    rownames = df4c.index
    tmp_df = df4c.copy()
    tmp_df.index = range(len(rownames))
    tmp_df = compute_overlap(tmp_df, groupcontrol, groupinterest, "symmetric")
    tmp_df.columns = [*tmp_df.columns[:-1], "distance"]
    tmp_df.index = rownames
    df4c = tmp_df.copy()

    return df4c


def separate_before_stats(ratiosdf, quality_dist_span):
    ratiosdf['has_replicates'] = flag_has_replicates(ratiosdf)
    try:
        quality_dist_span = float(quality_dist_span)
        good_df = ratiosdf.loc[(ratiosdf['distance/span'] >= quality_dist_span) &
                        (ratiosdf['has_replicates'] == 1) ] # has_replicates: 1 true

        undesired_mets = set(ratiosdf.index) - set(good_df.index)
        bad_df = ratiosdf.loc[list(undesired_mets)]
        good_df = good_df.drop(columns=['has_replicates'])
        bad_df = bad_df.drop(columns=['has_replicates'])
    except Exception as e:
        print(e)
        print("Error in separate_before_stats, check qualityDistanceOverSpan arg")

    return good_df, bad_df


def compute_mann_whitney_allH0(vInterest, vBaseline):
    # Calculate Mann–Whitney U test (a.k.a Wilcoxon rank-sum test,
    # or Wilcoxon–Mann–Whitney test, or Mann–Whitney–Wilcoxon (MWW/MWU), )
    # DO NOT : use_continuity False AND method "auto" at the same time.
    # because "auto" will set continuity depending on ties and sample size.
    # If ties in the data  and method "exact" (i.e use_continuity False) pvalues cannot be calculated
    # check scipy doc
    usta, p = scipy.stats.mannwhitneyu(vInterest,  vBaseline,  # method = "auto",
                                       use_continuity=False, alternative="less")

    usta2, p2 = scipy.stats.mannwhitneyu(vInterest,  vBaseline,  # method = "auto",
                                         use_continuity=False, alternative="greater")

    usta3, p3 = scipy.stats.mannwhitneyu(vInterest,  vBaseline,  # method = "auto",
                                         use_continuity=False, alternative="two-sided")

    # best (smaller pvalue) among all tailed tests
    pretups = [(usta, p), (usta2, p2), (usta3, p3)]
    tups = []
    for t in pretups:  # make list of tuples with no-nan pvalues
        if not np.isnan(t[1]):
            tups.append(t)

    if len(tups) == 0:  # if all pvalues are nan assign two sided result
        tups = [(usta3, p3)]

    stap_tup = min(tups, key=lambda x: x[1]) # nan already excluded
    stat_result = stap_tup[0]
    pval_result = stap_tup[1]

    return stat_result, pval_result


def compute_ranksums_allH0(vInterest: np.array, vBaseline: np.array):
    # The Wilcoxon rank-sum test tests the null hypothesis that two sets of measurements are drawn from the same distribution
    # ‘two-sided’: one of the distributions (underlying x or y) is stochastically greater than the other.
    # ‘less’: the distribution underlying x is stochastically less than the distribution underlying y.
    #  ‘greater’: the distribution underlying x is stochastically greater than the distribution underlying y.
    vInterest = vInterest[~np.isnan(vInterest)]
    vBaseline = vBaseline[~np.isnan(vBaseline)]
    sta, p = scipy.stats.ranksums(vInterest,  vBaseline,  alternative="less")
    sta2, p2 = scipy.stats.ranksums(vInterest,  vBaseline,  alternative="greater")
    sta3, p3 = scipy.stats.ranksums(vInterest,  vBaseline, alternative="two-sided")

    # best (smaller pvalue) among all tailed tests
    pretups = [(sta, p), (sta2, p2), (sta3, p3)]
    tups = []
    for t in pretups:  # make list of tuples with no-nan pvalues
        if not np.isnan(t[1]):
            tups.append(t)

    if len(tups) == 0:  # if all pvalues are nan assign two sided result
        tups = [(sta3, p3)]

    stap_tup = min(tups, key=lambda x: x[1]) # nan already excluded
    stat_result = stap_tup[0]
    pval_result = stap_tup[1]

    return stat_result, pval_result


def compute_wilcoxon_allH0(vInterest: np.array, vBaseline: np.array):
    #  Wilcoxon signed-rank test
    vInterest = vInterest[~np.isnan(vInterest)]
    vBaseline = vBaseline[~np.isnan(vBaseline)]
    sta, p = scipy.stats.wilcoxon(vInterest,  vBaseline,  alternative="less")
    sta2, p2 = scipy.stats.wilcoxon(vInterest,  vBaseline,  alternative="greater")
    sta3, p3 = scipy.stats.wilcoxon(vInterest,  vBaseline, alternative="two-sided")

    # best (smaller pvalue) among all tailed tests
    pretups = [(sta, p), (sta2, p2), (sta3, p3)]
    tups = []
    for t in pretups:  # make list of tuples with no-nan pvalues
        if not np.isnan(t[1]):
            tups.append(t)

    if len(tups) == 0:  # if all pvalues are nan assign two sided result
        tups = [(sta3, p3)]

    stap_tup = min(tups, key=lambda x: x[1]) # nan already excluded
    stat_result = stap_tup[0]
    pval_result = stap_tup[1]

    return stat_result, pval_result


def compute_brunnermunzel_allH0(vInterest: np.array, vBaseline: np.array):
    vInterest = vInterest[~np.isnan(vInterest)]
    vBaseline = vBaseline[~np.isnan(vBaseline)]
    sta, p = scipy.stats.brunnermunzel(vInterest,  vBaseline,  alternative="less")
    sta2, p2 = scipy.stats.brunnermunzel(vInterest,  vBaseline,  alternative="greater")
    sta3, p3 = scipy.stats.brunnermunzel(vInterest,  vBaseline, alternative="two-sided")

    # best (smaller pvalue) among all tailed tests
    pretups = [(sta, p), (sta2, p2), (sta3, p3)]
    tups = []
    for t in pretups:  # make list of tuples with no-nan pvalues
        if not np.isnan(t[1]):
            tups.append(t)

    if len(tups) == 0:  # if all pvalues are nan assign two sided result
        tups = [(sta3, p3)]

    stap_tup = min(tups, key=lambda x: x[1]) # nan already excluded
    stat_result = stap_tup[0]
    pval_result = stap_tup[1]

    return stat_result, pval_result


def statistic_absolute_geommean_diff(b_values: np.array, a_values: np.array):
    m_b = fg.compute_gmean_nonan(b_values)
    m_a = fg.compute_gmean_nonan(a_values)
    #denom = m_a + m_b
    #diff_normalized = abs((m_b - m_a) / denom)
    diff_absolute = abs(m_b - m_a)
    return diff_absolute


def run_statistical_test(redu_df, metas, contrast, whichtest):
    """
    Parameters
    ----------
    redu_df : pandas, reduced dataframe
    metas : pandas, 2nd element output of prepare4contrast.
    contrast : a list
    Returns
    -------
    DIFFRESULT : pandas
        THE PAIR-WISE DIFFERENTIAL ANALYSIS RESULTS.
    """
    mets = []
    stare = []
    pval = []
    metaboliteshere = redu_df.index
    for i in metaboliteshere:
        mets.append(i)
        row = redu_df.loc[i, :]  #  row is a series, colnames pass to index

        columnsInterest = metas.loc[metas["newcol"] == contrast[0], 'name_to_plot']
        columnsBaseline = metas.loc[metas["newcol"] == contrast[1], 'name_to_plot']

        vInterest = np.array(row[columnsInterest], dtype=float)
        vBaseline = np.array(row[columnsBaseline], dtype=float)

        vInterest = vInterest[~np.isnan(vInterest)]  # exclude nan elements
        vBaseline = vBaseline[~np.isnan(vBaseline)]  # exclude nan elements

        if (len(vInterest) >= 2) and (len(vBaseline) >= 2):  

            if whichtest == "MW":
                stat_result, pval_result = compute_mann_whitney_allH0(vInterest, vBaseline)

            elif whichtest == "Tt":
                stat_result, pval_result = scipy.stats.ttest_ind(vInterest, vBaseline, alternative="two-sided")

            elif whichtest == "KW":
                stat_result, pval_result = scipy.stats.kruskal(vInterest, vBaseline)

            elif whichtest == "ranksum":
                stat_result, pval_result = compute_ranksums_allH0(vInterest, vBaseline)

            elif whichtest == "Wcox": # signed-rank test: one sample (independence), or two paired or related samples
                stat_result, pval_result = compute_wilcoxon_allH0(vInterest, vBaseline)

            elif whichtest == "BrMu":
                stat_result, pval_result = compute_brunnermunzel_allH0(vInterest, vBaseline)

            elif whichtest == "prm-scipy":
                # test statistic is absolute geommean differences, so "greater" satisfy
                prm_res = scipy.stats.permutation_test((vInterest, vBaseline),
                                                       statistic=statistic_absolute_geommean_diff,
                                                       permutation_type='independent',
                                                       vectorized=False,
                                                       n_resamples=9999,
                                                       batch=None,
                                                       alternative='greater')
                stat_result, pval_result = prm_res.statistic, prm_res.pvalue


            stare.append(stat_result)
            pval.append(pval_result)

        else:
            stare.append(np.nan)
            pval.append(np.nan)

        # end if
    # end for
    prediffr = pd.DataFrame(data={"metabolite": mets, "stat": stare, "pvalue": pval})
    return prediffr


def steps_fitting_method(ratiosdf, out_histo_file):

    def auto_detect_tailway(good_df, best_distribution, args_param):
        min_pval_ = list()
        for tail_way in ["two-sided", "right-tailed"]:
            tmp = compute_p_value(good_df, tail_way, best_distribution, args_param)

            min_pval_.append(tuple([tail_way, tmp["pvalue"].min()]))

        return min(min_pval_, key=lambda x: x[1])[0]


    ratiosdf = compute_z_score(ratiosdf)
    best_distribution, args_param = find_best_distribution(ratiosdf,
                                 out_histogram_distribution=out_histo_file)
    autoset_tailway = auto_detect_tailway(ratiosdf, best_distribution, args_param)
    print("auto, best pvalues calculated :", autoset_tailway)
    ratiosdf = compute_p_value(ratiosdf, autoset_tailway, best_distribution, args_param)

    return ratiosdf


def compute_p_adjusted(df: pd.DataFrame, correction_method: str) -> pd.DataFrame:
    rej, pval_corr = ssm.multipletests(df['pvalue'].values, alpha=float('0.05'), method=correction_method)[:2]
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


def complete_columns_for_bad(bad_df, ratiosdf):
    columns_missing =  set(ratiosdf.columns) - set(bad_df.columns)
    for col in list(columns_missing):
        bad_df[col] = np.nan
    return bad_df


def filter_diff_results(ratiosdf, padj_cutoff, log2FC_abs_cutoff):
    ratiosdf['abslfc'] = ratiosdf['log2FC'].abs()
    ratiosdf = ratiosdf.loc[( ratiosdf['padj'] <= padj_cutoff ) &
                        ( ratiosdf['abslfc'] >= log2FC_abs_cutoff ), :]
    ratiosdf = ratiosdf.sort_values(['padj', 'pvalue', 'distance/span'], ascending=[True, True, False])
    ratiosdf = ratiosdf.drop(columns=['abslfc'])

    return ratiosdf


def reorder_columns_diff_end(df):
    standard_cols = [
        'count_nan_samples',
        'distance',
        'span_allsamples',
        'distance/span',
        'stat',
        'pvalue',
        'padj',
        'log2FC',
        'compartment']  # geommean_ratio not to be reordered, so not here

    desired_order = [
        'log2FC',
        'stat',
        'pvalue',
        'padj',
        'distance/span',
        'distance',
        'span_allsamples',
        'count_nan_samples',
        'compartment' ]

    standard_df = df[standard_cols]
    df = df.drop(columns=standard_cols)
    # reorder the standard part
    standard_df = standard_df[desired_order]
    # re-join them, indexes are the metabolites
    df = pd.merge(standard_df, df, left_index=True, right_index=True, how='left')
    return df


def run_differential_steps(measurements: pd.DataFrame, metadatadf: pd.DataFrame,
                       out_file_elements: dict, confidic: dict, whichtest:str,  args) -> None:

    out_dir = out_file_elements['odir']
    prefix = out_file_elements['prefix']
    co = out_file_elements['co']
    suffix = out_file_elements["suffix"]

    fg.detect_and_create_dir(f"{out_dir}/extended/")
    fg.detect_and_create_dir(f"{out_dir}/filtered/")

    for contrast in confidic['comparisons']:
        strcontrast = '_'.join(contrast)
        df4c, metad4c = fg.prepare4contrast(measurements, metadatadf,
                                       confidic['grouping'], contrast)
        df4c = df4c[(df4c.T != 0).any()]  # delete rows being zero everywhere
        # sort them by 'newcol' the column created by prepare4contrast
        metad4c = metad4c.sort_values("newcol")
        df4c = calc_reduction(df4c, metad4c)
        df4c = fg.countnan_samples(df4c, metad4c)  # adds nan_count_samples column
        df4c = distance_or_overlap(df4c, metad4c, contrast)
        df4c = compute_span_incomparison(df4c, metad4c, contrast)
        df4c['distance/span'] = df4c.distance.div(df4c.span_allsamples)
        ratiosdf = calc_ratios(df4c, metad4c, contrast)
        ratiosdf, df_bad = separate_before_stats(ratiosdf, args.qualityDistanceOverSpan)
        if whichtest == "disfit":
            out_histo_file = f"{out_dir}/extended/{prefix}--{co}--{suffix}-{strcontrast}_fitdist_plot.pdf"
            ratiosdf = steps_fitting_method(ratiosdf, out_histo_file)
            ratiosdf = compute_padj_version2(ratiosdf, 0.05, args.multitest_correction)

        else:
            extract_test_df = run_statistical_test(ratiosdf, metad4c, contrast, whichtest)
            extract_test_df = compute_padj_version2(extract_test_df, 0.05, args.multitest_correction)
            extract_test_df.set_index("metabolite", inplace=True)
            ratiosdf = pd.merge(ratiosdf, extract_test_df, left_index=True, right_index=True)

        ratiosdf["log2FC"] = np.log2(ratiosdf['geommean_ratio'])
        df_bad = complete_columns_for_bad(df_bad, ratiosdf)
        if df_bad.shape[0] >= 1:
            ratiosdf = pd.concat([ratiosdf, df_bad])

        ratiosdf["compartment"] = co
        ratiosdf = reorder_columns_diff_end(ratiosdf)
        ratiosdf = ratiosdf.sort_values(['padj', 'distance/span'], ascending=[True, False])
        ratiosdf.to_csv(f"{out_dir}/extended/{prefix}--{co}--{suffix}-{strcontrast}-{whichtest}.tsv",
                        index_label="metabolite", header=True, sep='\t')
        #filtered by thresholds :
        filtered_df = filter_diff_results(ratiosdf,
                                          confidic['thresholds']['padj'],
                                          confidic['thresholds']['absolute_log2FC'])
        filfi = f"{out_dir}/filtered/{prefix}--{co}--{suffix}-{strcontrast}-{whichtest}_filter.tsv"
        filtered_df.to_csv(filfi, index_label="metabolite", header=True, sep='\t')


def wrapper_for_abund(clean_tables_path, table_prefix, metadatadf,  confidic, args) -> None:
    out_diff_abun = out_path + "results/differential_analysis/abundance/"
    fg.detect_and_create_dir(out_diff_abun)
    whichtest = confidic['statistical_test']['abundances']
    suffix = confidic['suffix']
    compartments = metadatadf['short_comp'].unique().tolist()
    # dynamically open the file based on prefix, compartment and suffix:
    for co in compartments:
        meta_co = metadatadf.loc[metadatadf['short_comp'] == co, :]
        fn = f'{clean_tables_path}{table_prefix}--{co}--{suffix}.tsv'
        measurements = pd.read_csv(fn, sep='\t', header=0, index_col=0)  # compartment specific
        val_instead_zero = arg_repl_zero2value(args.abundance_replace_zero_with, measurements)
        measurements = measurements.replace(to_replace=0, value=val_instead_zero)
        out_file_elems = {'odir': out_diff_abun, 'prefix': table_prefix, 'co': co, "suffix":suffix}
        run_differential_steps(measurements, meta_co,
                                     out_file_elems,
                                     confidic, whichtest, args)


def wrapper_for_mefc(clean_tables_path, table_prefix, metadatadf,  confidic, args) -> None:
    out_diff = out_path + "results/differential_analysis/meanE_fracContr/"
    fg.detect_and_create_dir(out_diff)
    whichtest = confidic['statistical_test']['meanE_or_fracContrib']
    suffix = confidic['suffix']
    compartments = metadatadf['short_comp'].unique().tolist()
    # dynamically open the file based on prefix, compartment and suffix:
    for co in compartments:
        meta_co = metadatadf.loc[metadatadf['short_comp'] == co, :]
        fn = f'{clean_tables_path}{table_prefix}--{co}--{suffix}.tsv'
        measurements = pd.read_csv(fn, sep='\t', header=0, index_col=0) # compartment specific
        val_instead_zero = arg_repl_zero2value(args.meanEnrichOrFracContrib_replace_zero_with,
                                               measurements)
        measurements = measurements.replace(to_replace=0, value=val_instead_zero)
        out_file_elems = {'odir': out_diff, 'prefix': table_prefix, 'co': co, "suffix":suffix}
        run_differential_steps(measurements, meta_co,
                                     out_file_elems,
                                     confidic, whichtest, args)


def wrapper_for_isoAbsol(clean_tables_path, table_prefix, metadatadf, confidic, args) -> None:
    out_diff = out_path + "results/differential_analysis/isotopol_abs/"
    fg.detect_and_create_dir(out_diff)
    whichtest = confidic['statistical_test']['isotopologue_abs']
    suffix = confidic['suffix']
    compartments = metadatadf['short_comp'].unique().tolist()
    # dynamically open the file based on prefix, compartment and suffix:
    for co in compartments:
        meta_co = metadatadf.loc[metadatadf['short_comp'] == co, :]
        fn = f'{clean_tables_path}{table_prefix}--{co}--{suffix}.tsv'
        measurements = pd.read_csv(fn, sep='\t', header=0, index_col=0)  # compartment specific
        val_instead_zero = arg_repl_zero2value(args.isotopologueAbs_replace_zero_with,
                                               measurements)
        measurements = measurements.replace(to_replace=0, value=val_instead_zero)
        out_file_elems = {'odir': out_diff, 'prefix': table_prefix, 'co': co, "suffix": suffix}
        run_differential_steps(measurements, meta_co,
                           out_file_elems,
                           confidic, whichtest, args)


def wrapper_for_isoProp(clean_tables_path, table_prefix, metadatadf, confidic, args) -> None:
    out_diff = out_path + "results/differential_analysis/isotopol_prop/"
    fg.detect_and_create_dir(out_diff)
    whichtest = confidic['statistical_test']['isotopologue_prop']
    suffix = confidic['suffix']
    compartments = metadatadf['short_comp'].unique().tolist()
    # dynamically open the file based on prefix, compartment and suffix:
    for co in compartments:
        meta_co = metadatadf.loc[metadatadf['short_comp'] == co, :]
        fn = f'{clean_tables_path}{table_prefix}--{co}--{suffix}.tsv'
        measurements = pd.read_csv(fn, sep='\t', header=0, index_col=0)  # compartment specific
        val_instead_zero = arg_repl_zero2value(args.isotopologueProp_replace_zero_with,
                                               measurements)
    
        measurements = measurements.replace(to_replace=0, value=val_instead_zero)
        out_file_elems = {'odir': out_diff, 'prefix': table_prefix, 'co': co, "suffix": suffix}
        run_differential_steps(measurements, meta_co,
                               out_file_elems,
                               confidic, whichtest, args)


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

    tables_prefixes_dico = fg.clean_tables_names2dict(f'{clean_tables_path}TABLESNAMES.csv')
    metadatadf = fg.open_metadata(meta_path)

    # 1- abund
    if args.abundances:
        print("processing abundances")
        validate_zero_repl_arg(args.abundance_replace_zero_with)
        abund_tab_prefix = tables_prefixes_dico['name_abundance']
        wrapper_for_abund(clean_tables_path, abund_tab_prefix, metadatadf,  confidic, args)

    # 2- ME or FC
    if args.meanEnrich_or_fracContrib:
        print("processing mean enrichment or fractional contributions")
        validate_zero_repl_arg(args.meanEnrichOrFracContrib_replace_zero_with)
        fraccon_tab_prefix = tables_prefixes_dico['name_meanE_or_fracContrib']
        wrapper_for_mefc(clean_tables_path, fraccon_tab_prefix, metadatadf,  confidic, args)

    # 3- isotopologues
    if args.isotopologues:
        isos_abs_tab_prefix = tables_prefixes_dico['name_isotopologue_abs']
        isos_prop_tab_prefix = tables_prefixes_dico['name_isotopologue_prop']
        if (isos_abs_tab_prefix is not np.nan) and \
            (isos_abs_tab_prefix != "None") and \
            (isos_abs_tab_prefix is not None):
            print("processing absolute isotopologues")
            validate_zero_repl_arg(args.isotopologueAbs_replace_zero_with)
            wrapper_for_isoAbsol(clean_tables_path, isos_abs_tab_prefix,
                                 metadatadf,  confidic, args)
        else:
            print("processing isotopologues (values given as proportions)")
            validate_zero_repl_arg(args.isotopologueProp_replace_zero_with)
            wrapper_for_isoProp(clean_tables_path, isos_prop_tab_prefix, metadatadf,  confidic, args)

    print("end")




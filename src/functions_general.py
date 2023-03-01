#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 31 12:24:01 2022

@author: johanna
"""

import os
import yaml
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
import warnings


def wdir_configpaths_validate(wdir, config) -> None:
    aproval = [True,True, True]
    if not os.path.isdir(wdir):
        print(f"not a directory: {wdir}")
        aproval[0] = False
    if not os.path.isfile(config):
        print(f"not a configuration file: {config}")
        aproval[1] = False
    if not os.path.exists(wdir + "data/"):
        print(f"data/ folder is missing in {wdir}")
    if not (aproval[0] and aproval[1] and aproval[2]):
        raise ValueError("\nDid you inverted the order of directory and config file?")


def open_config_file(confifile):
    try:
        with open(confifile, "r") as f:
            confidic = yaml.load(f, Loader=yaml.Loader)
    except yaml.YAMLError as yam_err:
        print(yam_err)
        confidic = None
    except Exception as e:
        print(e)
        confidic = None

    if confidic is None:
        raise ValueError("\nimpossible to read configuration file")

    return confidic


def detect_and_create_dir(namenesteddir):
    if not os.path.exists(namenesteddir):
        os.makedirs(namenesteddir)


def fullynumeric(mystring):
    try:
        float(mystring)
        return True
    except Exception as e:
        return False


def open_metadata(workingdir, confidic):
    try:
        metadata = pd.read_csv(workingdir + "data/" + confidic['metadata_file'])
        return metadata
    except Exception as e:
        print(e)
        print('problem with opening metadata file')
        metadata = None
    if metadata is None:
        raise ValueError("\nproblem opening configuration file")


def verify_metadata_sample_not_duplicated(metadata_df) -> None:
    def yield_repeated_elems(mylist):
        occur_dic = dict(map(lambda x: (x, list(mylist).count(x)), mylist)) # credits: w3resource.com
        repeated_elems = list()
        for k in occur_dic.keys():
            if occur_dic[k] > 1:
                repeated_elems.append(k)
        return repeated_elems

    sample_duplicated = yield_repeated_elems(list(metadata_df['sample']))
    if len(sample_duplicated) > 0:
        txt_errors = f"-> duplicated sample names: {sample_duplicated}\n"
        raise ValueError(f"Error, found these conflicts in your metadata:\n{txt_errors}")




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



def a12(lst1, lst2, rev=True):
    """
    Non-parametric hypothesis testing using Vargha and Delaney's A12 statistic:
    how often is x in lst1 greater than y in lst2?
    == > it gives a size effect, good to highlight potentially real effects <==
    """
    # credits : Macha Nikolski
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
        # return array with each value divided by standard deviation, row-wise
        if sum(protein_values) == 0:
            reduced_abundances = protein_values  # because all row is zeroes
        else:
            reduced_abundances = protein_values / np.nanstd(protein_values, ddof=ddof)

        # replace values in result df
        res.loc[protein] = reduced_abundances
    return res


def give_reduced_df( df, ddof ):
    rownames = df.index
    df.index = range(len(rownames))  # index must be numeric because compute reduction accepted
    df_red = compute_reduction(df, ddof)  # reduce
    df_red.index = rownames
    return df_red


def compute_cv(reduced_abund):
    reduced_abu_np = reduced_abund.to_numpy().astype('float64')
    if np.nanmean(reduced_abu_np) != 0:
        return np.nanstd(reduced_abu_np) / np.nanmean(reduced_abu_np)
    elif np.nanmean(reduced_abu_np) == 0 and np.nanstd(reduced_abu_np) == 0:
        return 0
    else:
        return np.nan


def give_coefvar_new(df_red, red_meta, newcol : str):
    print("give cv")

    groups_ = red_meta[newcol].unique()
    tmpdico = dict()
    for group in groups_:
        samplesthisgroup = red_meta.loc[red_meta[newcol] == group, "sample"]
        subdf = df_red[samplesthisgroup]
        subdf = subdf.assign(CV=subdf.apply(compute_cv, axis=1))
        tmpdico[f"CV_{group}"] = subdf.CV.tolist()

    dfout = pd.DataFrame.from_dict(tmpdico)
    dfout.index = df_red.index
    return dfout

def compute_gmean_nonan(anarray):
    anarray = np.array(anarray, dtype=float)
    anarray = anarray[~np.isnan(anarray)]
    if sum(anarray) == 0:  # replicates all zero
        outval = 0
    else:
        outval = stats.gmean(anarray)
    return outval

def give_geommeans_new(df_red, metad4c, newcol : str , c_interest, c_control):
    """
    output: df, str, str
    """
    print("give GEO")

    sams_interest = metad4c.loc[metad4c['newcol'] == c_interest, "sample"]
    sams_control = metad4c.loc[metad4c['newcol'] == c_control, "sample"]
    dfout = df_red.copy()
    geomcol_interest = "gm_"+c_interest
    geomcol_control = "gm_" + c_control
    dfout[geomcol_interest] = [np.nan for i in range(dfout.shape[0])]
    dfout[geomcol_control] = [np.nan for i in range(dfout.shape[0])]

    for i, row in df_red.iterrows():
        metabolite = i
        vec_interest = np.array(row[sams_interest])  #[ sams_interest]
        vec_control = np.array(row[sams_control])

        val_interest = compute_gmean_nonan(vec_interest)
        val_control = compute_gmean_nonan(vec_control)

        dfout.loc[i, geomcol_interest] = val_interest
        dfout.loc[i, geomcol_control] = val_control

    return dfout, geomcol_interest, geomcol_control



def give_ratios_df(df1, geomInterest, geomControl):
    print("give RATIO")
    df = df1.copy()
    df = df.assign(ratio=[np.nan for i in range(df.shape[0])])
    for i, row in df1.iterrows():
        intere = row[geomInterest]
        contr = row[geomControl]
        if contr == 0 :
            df.loc[i, "ratio"] = intere / 1e-10
        else:
            df.loc[i, "ratio"] = intere / contr

    return df


def countnan_samples(df, metad4c):
    vecout = []
    grs = metad4c['newcol'].unique()
    gr1 = metad4c.loc[metad4c['newcol']  == grs[0], "sample"]
    gr2 = metad4c.loc[metad4c['newcol']  == grs[1], "sample"]

    for i, row in df.iterrows():
        vec1 = row[gr1].tolist()
        vec2 = row[gr2].tolist()
        val1 = np.sum(np.isnan(vec1))
        val2 = np.sum(np.isnan(vec2))
        vecout.append(tuple([str(val1)+'/'+str(len(vec1)),
                             str(val2)+'/'+str(len(vec2))]))

    df['count_nan_samples'] = vecout#[str(tup) for tup in vecout]
    return df

def add_alerts(df, metad4c):
    # deprecated
    df['alert'] = ''
    df.loc[df["distance"] < 0, "alert"] = "overlap"
    alert_reps = list()
    for i in df['count_nan_samples'].tolist():
        if i[0] >= 2 or i[1] >= 2:
            alert_reps.append("no replicates")
        else:
            alert_reps.append('')

    df['foo'] = alert_reps
    df.loc[df['foo'] != '', "alert"] = "no replicates"
    df = df.drop(columns = ['foo'])
    return df


def calcgeommean(avector, eps):
    # TODO: old, used in differential_univariate, work pending
    vech = np.array(avector)
    vech[vech == 0] = eps  # replace any zeroes
    return np.exp(np.mean(np.log(vech)))



def plot_overlap_hist(df_overls, colname_symetric, colname_assymetric, fileout):
    import seaborn as sns
    import matplotlib as plt
    """just for debugging or other tests"""
    values_sym = df_overls[colname_symetric]
    a = pd.DataFrame({'value' : values_sym,
                      'type_overlap': ["symm" for i in range(len(values_sym))] })
    vasym = df_overls[colname_assymetric]
    b = pd.DataFrame({'value': vasym,
                      'type_overlap': ["assym" for i in range(len(vasym))]})
    dfplotov = pd.concat([a,b], ignore_index=True, axis=0)

    with sns.axes_style("darkgrid"):
        sns.displot(data=dfplotov, x = 'value', hue='type_overlap',
                       kde=False)
        plt.savefig(fileout)
    plt.close()
    return 0


# from here, functions for isotopologue preview (mode None)

def add_metabolite_column(df):
    theindex = df.index
    themetabolites = [i.split("_m+")[0] for i in theindex]
    df = df.assign(metabolite=themetabolites)

    return df


def add_isotopologue_type_column(df):
    theindex = df.index
    preisotopologue_type = [i.split("_m+")[1] for i in theindex]
    theisotopologue_type = [int(i) for i in preisotopologue_type]
    df = df.assign(isotopologue_type=theisotopologue_type)

    return df


def save_heatmap_sums_isos(thesums, figuretitle, outputfigure):
    fig, ax = plt.subplots( figsize=(9,10))
    sns.heatmap(thesums,
                annot=True, fmt=".1f", cmap="crest",
                square = True,
                annot_kws = {
                    'fontsize' : 6
                },
                ax= ax)
    plt.xticks(rotation=90)
    plt.title(figuretitle)
    plt.savefig(outputfigure)
    plt.close()

    return 0


def givelevels(melted):
    another = melted.copy()
    another = another.groupby('metabolite').min()
    another = another.sort_values(by='value', ascending=False)
    levelsmetabolites = another.index
    tmp = melted['metabolite']
    melted['metabolite'] = pd.Categorical(tmp, categories=levelsmetabolites)

    return melted


def table_minimalbymet(melted, fileout):
    another = melted.copy()
    another = another.groupby('metabolite').min()
    another = another.sort_values(by='value', ascending=False)
    another.to_csv(fileout, sep='\t', header=True)


def save_rawisos_plot(dfmelt, figuretitle, outputfigure):
    fig, ax = plt.subplots(1, 1, figsize=(16, 10))
    sns.stripplot(ax=ax, data=dfmelt, x="value", y="metabolite", jitter=False,
                  hue="isotopologue_type", size=4, palette="tab20")
    plt.axvline(x=0,
                ymin=0,
                ymax=1,
                linestyle="--", color="gray")
    plt.axvline(x=1,
                ymin=0,
                ymax=1,
                linestyle="--", color="gray")
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
    plt.title(figuretitle)
    plt.xlabel("fraction")
    plt.savefig(outputfigure)
    plt.close()
    return 0

# end functions for isotopologue preview
# END
# dico with nb occur: https://www.w3resource.com/python-exercises/lambda/python-lambda-exercise-49.php
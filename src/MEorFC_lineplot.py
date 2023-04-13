#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  9 11:20:09 2022

@author: johanna
"""

import os
import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import functions_general as fg


def lineplot_args():
    parser = argparse.ArgumentParser(prog="python -m DIMet.src.abundance_bars",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('config', type=str,
                        help="configuration file in absolute path")

    parser.add_argument('--color_metabolites', type=str,  default="auto_multi_color",
                        help="auto_multi_color | PATH_TO_COLOR_TABLE\n\
                        \nIf you add a PATH_TO_COLOR_TABLE, must be a .csv \
                        with header: metabolite,color")

    parser.add_argument("--xaxis_title", type=str, default='time',
                        help="the x axis lab title, for example, \
                        'time', 'hours', 'minutes', 'seconds'")

    parser.add_argument("--plot_the_rest_of_metabolites",
                        action=argparse.BooleanOptionalAction, default=False,
                        help="plot a simple facetgrid with the metabolites that weren't\
                        in your config (can be slow if many metabolites)")

    parser.add_argument('--alpha', type=int, default=1)

    return parser


def yieldfraccountrib(df_co, metada_co, co):

    dfcompartment = df_co.T

    metabolites = dfcompartment.columns
    dfcompartment['name_to_plot'] = dfcompartment.index
    dfcompartment = pd.merge(dfcompartment, metada_co, on='name_to_plot')

    # empty dataframe to fill
    df4plot = pd.DataFrame(
        columns=["timenum", "condition", "metabolite", "Fractional Contribution (%)"]
    )

    for z in range(len(metabolites)):
        subdf = dfcompartment.loc[:, [metabolites[z], "timenum", "condition"]]
        subdf["metabolite"] = metabolites[z]
        subdf["Fractional Contribution (%)"] = subdf[metabolites[z]] * 100
        subdf = subdf.drop(columns=[metabolites[z]])
        # df4plot = df4plot.append(subdf, ignore_index=True)
        df4plot = pd.concat((df4plot, subdf), ignore_index=True)
    return df4plot


def nestedDi2list(anestedD):
    elems_ = []
    for p in anestedD.values():
        for q in p:
            elems_.append(q)
    return elems_


def df2musddf(one_m):
    """
    input: dataframe by one metabolite:
      "timenum", "condition", "metabolite" "Fractional Contribution (%)"
    returns :
$        dataframe by one metabolite, metabolite column deleted:
                condition  timenum    mean    sd
    108        Control     0  0.000893  0.002611
    111        Control     1  0.236453  0.023246
    ...
    123        Control    24  0.854101  0.055241
    126  L-Cycloserine     0  0.010083  0.003465
    129  L-Cycloserine     1  0.259570  0.008602
    ...
    141  L-Cycloserine    24  0.815613  0.050756
    """
    m_s = one_m[["condition", "timenum"]].drop_duplicates()
    m_s["mean"] = 0
    m_s["sd"] = 0
    for kc in m_s["condition"]:
        for hou in one_m["timenum"]:
            ss = one_m.loc[(one_m["timenum"] == hou) & (one_m["condition"] == kc), :]
            mymean = np.mean(ss["Fractional Contribution (%)"])
            mysd = np.std(ss["Fractional Contribution (%)"])
            m_s.loc[(m_s["condition"] == kc) & (m_s["timenum"] == hou), "mean"] = mymean
            m_s.loc[(m_s["condition"] == kc) & (m_s["timenum"] == hou), "sd"] = mysd
    return m_s


def complextimetracer_plot(co, df4plot, grmetsD, mycolorsD, outfile, args):
    themets = nestedDi2list(grmetsD)
    hoursticks = df4plot['timenum'].unique()
    somem = df4plot.loc[df4plot["metabolite"].isin(themets), :].copy()
    m_s = pd.DataFrame(columns=["condition", "timenum", "mean", "sd", "metabolite"])
    for k in set(somem["metabolite"]):
        one_m = somem[somem["metabolite"] == k]
        m_s1 = df2musddf(one_m)
        m_s1["metabolite"] = k
        m_s = pd.concat([m_s, m_s1])  # mean and sd ( bio replicates)
    figziz = 6.2 * len(grmetsD)
    fig, axs = plt.subplots(2, len(grmetsD), sharey=False, figsize=(figziz, 11))

    def doemptyrow(axs, nbcolumns):
        for z in range(0, nbcolumns): # do empty row
            axs[0, z].set_axis_off()
        return axs

    axs = doemptyrow(axs, len(grmetsD))

    for z in range(len(grmetsD)):
        sns.lineplot(
            ax=axs[1, z],
            x="timenum",
            y="Fractional Contribution (%)",
            hue="metabolite",
            style="condition",
            err_style=None,
            alpha=args.alpha,
            linewidth=4.5,
            palette=mycolorsD,
            zorder=1,
            data=somem.loc[somem["metabolite"].isin(grmetsD[z])],
            legend=True,
        )
        axs[1,z].set_xticks([int(i) for i in hoursticks])
        m_s1 = m_s.loc[m_s["metabolite"].isin(grmetsD[z])]
        axs[1,z].scatter(
            m_s1["timenum"], m_s1["mean"], s=22, facecolors="none", edgecolors="black"
        )
        axs[1,z].errorbar(
            m_s1["timenum"],
            m_s1["mean"],
            yerr=m_s1["sd"],
            fmt="none",
            capsize=3,
            ecolor="black",
            zorder=2
        )
        axs[1, z].set(ylabel=None),
        axs[1, z].set(xlabel=args.xaxis_title)
        axs[1, z].legend(loc="upper center", bbox_to_anchor= (0.5, 2), frameon=False)

    fig.suptitle(co)
    plt.subplots_adjust(bottom=0.1, right=0.8, hspace=0.1, wspace=0.4)

    fig.text(
        0.1, 0.3, "Fractional Contribution (%)", va="center", rotation="vertical"
    )
    fig.savefig(outfile, format="pdf")
    return 0


def give_colors_by_arg(argument_color_metabolites, metabos_nested_dico):

    handycolors = ["rosybrown", "lightcoral", "brown", "firebrick", "tomato", "coral", "sienna",
                  "darkorange",   "peru", "darkgoldenrod", "gold", "darkkhaki", "olive",
                  "yellowgreen", "limegreen", "green", "lightseagreen",  "mediumturquoise",
                  "darkcyan", "teal", "cadetblue", "slategrey", "stellblue", "navy",
                  "darkslateblue", "blueviolet",
                  "darkochid", "purple", "mediumvioletred", "crimson"]

    dicoout = dict()

    if argument_color_metabolites =="auto_multi_color":
        allmets = set()
        for co in metabos_nested_dico.keys():
            for k in metabos_nested_dico[co].keys():
                allmets.update(set(metabos_nested_dico[co][k]))
        mets_l = sorted(list(allmets))
        if len(mets_l) > 12:
            for i in range(len(mets_l)):
                dicoout[mets_l[i]] = handycolors[i]
        else:
            palettecols = sns.color_palette("Paired", 12)
            for i in range(len(mets_l)):
                dicoout[mets_l[i]] = palettecols[i]

    else: # argument_color_metabolites is a csv file
        try:
            df = pd.read_csv(argument_color_metabolites, header=True)
            for i, row in df.iterrows():
                metabolite = df.iloc[i, 0] # first column metabolite
                color = df.iloc[i, 1] # second column is color
                dicoout[metabolite] = color
        except Exception as e:
            print(e, "\n could not assign color to ")
            dicoout =  None

    return dicoout


def undesired_metabs_fast_plot(df4plot, trashmets, co, out_plot_dir):
    ppp = df4plot.loc[df4plot["metabolite"].isin(trashmets)]

    plt.rcParams.update({"font.size": 10})

    g = sns.FacetGrid(
        ppp, row="condition", col="metabolite", sharey=False, margin_titles=True
    )
    g.map(
        sns.lineplot,
        "timenum",
        "Fractional Contribution (%)",
        size=5,
        alpha=1,
        # error bars config :
        marker="o",
        err_style="bars",
        errorbar=('ci', 95),
        color="darkgray"
    )
    g.savefig(out_plot_dir + "undesired_" + co + ".pdf", format="pdf")


def savefraccontriplots(table_prefix, metadatadf, out_plot_dir, confidic, args):
    try:
        metabos_nested_dico = confidic["groups_toplot_frac_contribs"]
    except KeyError:
        try:
            metabos_simple_dico = confidic["metabolites_to_plot"]
            single_metab_by_facet_dico = dict()
            for co in metabos_simple_dico:
                single_metab_by_facet_dico[co] = dict()
                i = 0
                for metabo in confidic["metabolites_to_plot"][co]:
                    single_metab_by_facet_dico[co][i] = [metabo]
                    i += 1
            metabos_nested_dico = single_metab_by_facet_dico
        except KeyError:
            print("No metabolites for plotting in your config file")
    # end try

    condilevels = confidic["conditions"]  # <= locate where it is used

    out_path = os.path.expanduser(confidic['out_path'])
    suffix = confidic['suffix']
    compartments = metadatadf['short_comp'].unique().tolist()

    coloreachmetab = give_colors_by_arg(args.color_metabolites, metabos_nested_dico)

    # dynamically open the file based on prefix, compartment and suffix:
    for co in compartments:
        metada_co = metadatadf.loc[metadatadf['short_comp'] == co, :]
        fn = f'{out_path}results/prepared_tables/{table_prefix}--{co}--{suffix}.tsv'
        adf = pd.read_csv(fn, sep='\t', header=0, index_col=0)

        df4plot = yieldfraccountrib(adf, metada_co , co)
        df4plot["condition"] = pd.Categorical(df4plot["condition"], condilevels)
        grmetsD = metabos_nested_dico[co]

        sns.set_style(
            {"font.family": "sans-serif", "font.sans-serif": "Liberation Sans"}
        )
        plt.rcParams.update({"font.size": 22})
        # https://stackoverflow.com/questions/53137983/define-custom-seaborn-color-palette

        ofile = "{}MEorFC_line_plot_{}.pdf".format(out_plot_dir, co)
        complextimetracer_plot(co, df4plot, grmetsD, coloreachmetab, ofile, args)

         #### uninteresting metabolites, quick ugly print :
        pickedmets_ = nestedDi2list(grmetsD)
        thiscompartment_mets = set(df4plot["metabolite"])

        trashmets = set(thiscompartment_mets) - set(pickedmets_)
        if args.plot_the_rest_of_metabolites:
            undesired_metabs_fast_plot(df4plot, trashmets, co, out_plot_dir)

    return 0


if __name__ == "__main__":
    parser = lineplot_args()
    args = parser.parse_args()
    configfile = os.path.expanduser(args.config)
    confidic = fg.open_config_file(configfile)
    fg.auto_check_validity_configuration_file(confidic)
    out_path = os.path.expanduser(confidic['out_path'])
    meta_path = os.path.expanduser(confidic['metadata_path'])
    clean_tables_path = out_path + "results/prepared_tables/"

    # tpd : tables prefixes dictionary
    tpd = fg.clean_tables_names2dict(f'{out_path}results/prepared_tables/TABLESNAMES.csv')
    metadatadf = fg.open_metadata(meta_path)

    table_prefix = tpd['name_meanE_or_fracContrib']
    out_plot_dir = out_path + "results/plots/lineplots_MEorFC/"
    fg.detect_and_create_dir(out_plot_dir)

    print(" Mean Enrichment -or Fractional contributions- plots \n")
    savefraccontriplots(table_prefix, metadatadf, out_plot_dir, confidic, args)





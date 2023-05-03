#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 18 11:14:34 2022

@author: johanna
"""

import os
import argparse
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import matplotlib.ticker as mticker
import functions_general as fg


def stacked_args():
    parser = argparse.ArgumentParser(
        prog="python -m DIMet.src.abundance_bars",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('config', type=str,
                        help="configuration file in absolute path")

    parser.add_argument("--separated_plots_by_condition",
                        action=argparse.BooleanOptionalAction,
                        default=False,
                        help="print one pdf plot file by condition")

    parser.add_argument("--max_nb_carbons_possible", type=int,
                        default=24,
                        help="number of carbons, defines colors")

    parser.add_argument("--plots_height", type=float,
                        default=4.8,
                        help="height for all the subplots")

    parser.add_argument("--sharey",
                        action=argparse.BooleanOptionalAction,
                        default=False,
                        help="share y axis across subplots")

    parser.add_argument("--x_ticks_text_size", type=float,
                        default=18)

    parser.add_argument("--y_ticks_text_size", type=float,
                        default=19)

    return parser


def isotopol_prop_2df4plot(df_co, metada_co, levelshours_str):
    """
    Imitates de behaviour of a 'melt', but this function is more secure
    example:
      pd.merge output, where 1st column has the values we want:

        L-Phenylalanine_m+2    timenum   condition ...  sample descr...
        0.01                           0          Control        xxxxxx

      is transformed into:
          timenum    condition    isotopolgFull        Isotopologue
          0           control      L-Phenylala...    0.01

    """
    dfcompartment = df_co.T

    metabolites = dfcompartment.columns
    dfcompartment['name_to_plot'] = dfcompartment.index
    dfcompartment = pd.merge(dfcompartment, metada_co, on='name_to_plot')
    # empty dataframe to fill
    df4plot = pd.DataFrame(
        columns=[
            "timenum",
            "condition",
            "isotopolgFull",
            "Isotopologue Contribution (%)",
        ]
    )
    df4plot["timenum"] = pd.Categorical(df4plot["timenum"], levelshours_str)
    # iteratively pile up
    for z in range(len(metabolites)):
        subdf = dfcompartment.loc[:, [metabolites[z], "timenum", "condition"]]
        subdf['timenum'] = subdf['timenum'].astype(str)
        # 1st colname as cell value, reps
        subdf["isotopolgFull"] = metabolites[z]
        subdf["Isotopologue Contribution (%)"] = subdf[metabolites[z]] * 100
        # 1st col no longer needed
        subdf = subdf.drop(columns=[metabolites[z]])
        df4plot = pd.concat([df4plot, subdf])
        del subdf

    return df4plot


def massageisotopologues(df4plot):
    """
    returns dataframe splitting metabolite and m+x into two separate columns
    and also correcting weird values
    """
    xu = {"name": [], "m+x": []}
    for ch in df4plot["isotopolgFull"]:
        elems = ch.split("_m+")
        xu["name"].append(elems[0])
        xu["m+x"].append("m+" + elems[1])
    df4plot["metabolite"] = xu["name"]
    df4plot["m+x"] = xu["m+x"]
    # dealing with weird values: bigger than 100 and less than 0 :
    df4plot.loc[
        df4plot["Isotopologue Contribution (%)"] > 100,
        "Isotopologue Contribution (%)"
    ] = 100

    df4plot.loc[
        df4plot["Isotopologue Contribution (%)"] < 0,
        "Isotopologue Contribution (%)"
    ] = 0

    return df4plot


def preparemeansreplicates(df4plot, selectedmets):
    """
    returns a dictionary of dataframes, keys are metabolites
    """
    dfcopy = df4plot.copy()
    dfcopy = dfcopy.groupby(
        ["condition", "metabolite", "m+x", "timenum"])\
        .mean("Isotopologue Contribution %")
    dfcopy = dfcopy.reset_index()

    dfs_Dico = dict()
    for i in selectedmets:
        tmp = dfcopy.loc[dfcopy["metabolite"] == i, ].reset_index()
        # set m+x as numeric to avoid any bad reordering of stacked m+x
        tmp["m+x"] = tmp["m+x"].str.split("m+", regex=False).str[1]
        tmp["m+x"] = tmp["m+x"].astype(int)

        dfs_Dico[i] = tmp
    return dfs_Dico


def addcombinedconditime(dfs_Dico, combined_tc_levels):
    """
    add column 'timeANDcondi' to each  metabolite dataframe in Dico
    """
    for metab in dfs_Dico.keys():
        dfs_Dico[metab]["timeANDcondi"] = \
            dfs_Dico[metab]["timenum"] + " : " + dfs_Dico[metab]["condition"]

        dfs_Dico[metab]["timeANDcondi"] = pd.Categorical(
            dfs_Dico[metab]["timeANDcondi"],
            combined_tc_levels)
    return dfs_Dico


def addcategoricaltime(dfs_Dico, levelstime_str):
    """
    existing column 'timenum' as categorical, with specific order.
    Each metabolite dataframe in Dico is processed.
    """
    for metab in dfs_Dico.keys():
        dfs_Dico[metab]["timenum"] = pd.Categorical(
            dfs_Dico[metab]["timenum"],
            levelstime_str)
    return dfs_Dico


def give_colors_carbons(nb_of_carbons):
    color_d = dict()
    color_d[0] = "lightgray"  # m+0
    # set colors m+1 to m+8 from Spectral palette,
    # with custom spaced selected colors (validated)
    spectralPal = sns.color_palette("Spectral", 30)
    color_d[1] = spectralPal[29]
    color_d[2] = spectralPal[26]
    color_d[3] = spectralPal[21]
    color_d[4] = spectralPal[17]
    color_d[5] = spectralPal[10]
    color_d[6] = spectralPal[6]
    color_d[7] = spectralPal[2]
    color_d[8] = spectralPal[0]
    # rest of the colors from tab20b palette
    added_pal = sns.color_palette("tab20b", 20)
    i = 9
    j = 19
    while i <= nb_of_carbons:
        color_d[i] = added_pal[j]
        j = j - 1
        i += 1

    return color_d


def complexstacked(co, selectedmets, dfs_Dico, outfilename, args,
                   figu_width, xlabyesno,  wspace_stacks, numbers_size,
                   x_to_plot, x_ticks_text_tilt):
    """plot highly custom, recommended that selectedmets <= 6 subplots"""
    palsautoD = give_colors_carbons(args.max_nb_carbons_possible)

    sns.set_style({"font.family": "sans-serif",
                   "font.sans-serif": "Liberation Sans"})
    f, axs = plt.subplots(1, len(selectedmets), sharey=args.sharey,
                          figsize=(figu_width, args.plots_height))
    plt.rcParams.update({"font.size": 20})

    for z in range(len(selectedmets)):

        axs[z].set_title(selectedmets[z])
        sns.histplot(
            ax=axs[z],
            data=dfs_Dico[selectedmets[z]],
            x=x_to_plot,
            # Use the value variable here to turn histogram counts into
            # weighted values.
            weights="Isotopologue Contribution (%)",
            hue="m+x",
            multiple="stack",
            palette=palsautoD,
            # Add  borders to the bars.
            edgecolor="black",
            # Shrink the bars a bit so they don't touch.
            shrink=0.85,
            alpha=1,
            legend=False,
        )
        #
        axs[z].tick_params(axis="x",
                           labelrotation=x_ticks_text_tilt,
                           labelsize=args.x_ticks_text_size)
        axs[z].tick_params(axis="y", length=3,
                           labelsize=19)
        axs[z].set_ylim([0, 100])

        for bar in axs[z].patches:
            # assign stacked bars text color
            thebarvalue = round(bar.get_height(), 1)
            if thebarvalue >= 100:
                thebarvalue = 100  # no decimals if 100
            if round(bar.get_height(), 1) >= 4:
                axs[z].text(
                    # Put the text in the middle of each bar. get_x returns
                    # the start so we add half the width to get to the middle.
                    bar.get_x() + bar.get_width() / 2,
                    # Vertically, add the height of the bar to the start of
                    # the bar, along with the offset.
                    (bar.get_height() / 2) + (bar.get_y()) + 2,  #
                    # This is actual value we'll show.
                    thebarvalue,
                    # Center the labels and style them a bit.
                    ha="center",
                    size=numbers_size,
                )  # end axs[z].text
            else:
                continue
            # end if round
        # end for bar

        axs[z].set_ylabel("", size=20)
        axs[z].xaxis.set_tick_params(length=0)  # no need of x ticks
        axs[z].set_xlabel("", size=13)
    # end for z

    [ax.invert_yaxis() for ax in axs]  # invert y, step 1

    for ax in axs:
        ylabels = ax.get_yticks().tolist()
        ax.yaxis.set_major_locator(mticker.FixedLocator(ylabels))
        ax.set_yticklabels([100 - int(i) for i in ylabels])  # invert y, step2

    if xlabyesno == "no":
        for ax in axs:
            xlabelshere = ax.get_xticks()
            ax.set_xticklabels(["" for i in xlabelshere])

    f.subplots_adjust(hspace=0.5, wspace=wspace_stacks, top=0.85,
                      bottom=0.26, left=0.15, right=0.99)
    # f.suptitle(f"compartment {co.upper()}  \n", fontsize=20)
    f.text(0.03, 0.57, "Isotopologue Contribution (%)\n", va="center",
           rotation="vertical", size=20)
    f.savefig(outfilename,
              bbox_inches="tight", format="pdf")
    plt.close()
    return 0  # end the new version foo


def add_xemptyspace_tolabs(condilevels, levelshours_str):
    condilevels += ["xemptyspace"]  # to add space among time categories
    combined_tc_levels = list()
    for x in levelshours_str:
        for y in condilevels:
            if y == "xemptyspace":
                combined_tc_levels.append(str(x))
            else:
                combined_tc_levels.append(f'{x} : {y}')
    return condilevels, combined_tc_levels


def time_plus_condi_labs(condilevels, levelshours_str):
    combined_tc_levels = list()
    for x in levelshours_str:
        for y in condilevels:
            combined_tc_levels.append(f'{x} : {y}')
    return combined_tc_levels


def givelabelstopalsD(palsautoD):
    tmp = dict()
    for k in palsautoD.keys():
        tmp["m+"+str(k)] = palsautoD[k]
    return tmp


def save_isotopol_stacked_plot(table_prefix, metadatadf,
                               out_plot_dir, confidic, args):
    out_path = os.path.expanduser(confidic['out_path'])
    suffix = confidic['suffix']
    compartments = metadatadf['short_comp'].unique().tolist()

    condilevels = confidic["conditions"]  # <= locate where it is used
    width_each_stack = float(confidic["width_each_stack"])

    wspace_stacks = float(confidic["wspace_stacks"])
    numbers_size = int(confidic["numbers_size"])

    levelshours_str = [str(i) for i in sorted(metadatadf['timenum'].unique())]

    # dynamically open the file based on prefix, compartment and suffix:
    for co in compartments:
        metada_co = metadatadf.loc[metadatadf['short_comp'] == co, :]
        the_folder = f'{out_path}results/prepared_tables/'
        fn = f'{the_folder}{table_prefix}--{co}--{suffix}.tsv'
        adf = pd.read_csv(fn, sep='\t', header=0, index_col=0)

        # note that pandas automatically transform any 99.9% in decimal 0.999
        df4plot = isotopol_prop_2df4plot(adf, metada_co, levelshours_str)

        df4plot = massageisotopologues(df4plot)

        selectedmets = confidic['metabolites_to_plot'][co]

        dfs_Dico = preparemeansreplicates(df4plot,  selectedmets)

        # adapt width to nb of metabolites
        figu_width = width_each_stack * len(selectedmets)

        if args.separated_plots_by_condition:
            for condition in condilevels:
                outfname = "{}isotopologues_stack_{}--{}.pdf".\
                    format(out_plot_dir, condition, co)
                metada_this_condi = metada_co.loc[
                    metada_co['condition'] == condition, :]
                df_this_condi = adf[metada_this_condi['name_to_plot']]
                df4plot = isotopol_prop_2df4plot(df_this_condi,
                                                 metada_this_condi,
                                                 levelshours_str)
                df4plot = massageisotopologues(df4plot)
                dfs_Dico = preparemeansreplicates(df4plot,  selectedmets)
                dfs_Dico = addcategoricaltime(dfs_Dico, levelshours_str)
                complexstacked(
                    co, selectedmets, dfs_Dico, outfname, args,
                    figu_width, xlabyesno="yes", wspace_stacks=wspace_stacks,
                    numbers_size=numbers_size, x_to_plot="timenum",
                    x_ticks_text_tilt=0
                )
                plt.close()
        else:
            # condilevels, combined_tc_levels = add_xemptyspace_tolabs(
            #                               condilevels, levelshours_str)
            combined_tc_levels = time_plus_condi_labs(condilevels,
                                                      levelshours_str)

            dfs_Dico = addcombinedconditime(dfs_Dico, combined_tc_levels)

            outfname = "{}isotopologues_stack--{}.pdf".format(out_plot_dir,
                                                              co)
            complexstacked(
                co, selectedmets, dfs_Dico, outfname, args,
                figu_width, xlabyesno="yes",  wspace_stacks=wspace_stacks,
                numbers_size=numbers_size, x_to_plot="timeANDcondi",
                x_ticks_text_tilt=90
            )
            plt.close()
            # new :
            outfnameNoXlab = "{}isotopologues_stack--{}_noxlab.pdf".\
                format(out_plot_dir, co)
            complexstacked(
                co, selectedmets, dfs_Dico,  outfnameNoXlab, args,
                figu_width, xlabyesno="no",  wspace_stacks=wspace_stacks,
                numbers_size=numbers_size, x_to_plot="timeANDcondi",
                x_ticks_text_tilt=90
            )

        # legend alone
        plt.figure(figsize=(4, args.max_nb_carbons_possible * 0.6))
        palsautoD = give_colors_carbons(args.max_nb_carbons_possible)
        palsautoD_labeled = givelabelstopalsD(palsautoD)
        myhandless = []
        for c in palsautoD_labeled.keys():
            paobj = mpatches.Patch(facecolor=palsautoD_labeled[c],
                                   label=c, edgecolor="black")
            myhandless.append(paobj)
        plt.legend(handles=myhandless, labelspacing=0.01)
        plt.axis("off")
        plt.savefig(f"{out_plot_dir}legend_isotopologues_stackedbars.pdf",
                    format="pdf")

    return 0


if __name__ == "__main__":
    parser = stacked_args()
    args = parser.parse_args()
    configfile = os.path.expanduser(args.config)
    confidic = fg.open_config_file(configfile)
    fg.auto_check_validity_configuration_file(confidic)
    confidic = fg.remove_extensions_names_measures(confidic)

    out_path = os.path.expanduser(confidic['out_path'])
    meta_path = os.path.expanduser(confidic['metadata_path'])
    clean_tables_path = out_path + "results/prepared_tables/"

    metadatadf = fg.open_metadata(meta_path)

    # isotopologues proportions
    isotopol_prop_tab_prefix = confidic['name_isotopologue_prop']

    out_plot_dir = out_path + "results/plots/stacked_Isotopologue_prop/"
    fg.detect_and_create_dir(out_plot_dir)
    save_isotopol_stacked_plot(isotopol_prop_tab_prefix, metadatadf,
                               out_plot_dir, confidic, args)
# END


############################################################################
## old functions keept in this file because reproduce the figures of very first versions
############################################################################
############################################################################
def custom_colors_stacked():
    """
    Returns colors defined by default
     important : darkbarcolor sets which rgba obliges text inside bar to be white color
    """
    darkbarcolor = [
        (0.0, 0.0, 0.7529411764705882, 1.0),  # equivalent to #410257
        (0.2549019607843137, 0.00784313725490196, 0.3411764705882353, 1.0),  # equivalent to #0000c0
    ]
    palsD = {
        "m+0": "#410257",
        "m+1": "#0000c0",
        "m+2": "#55a0fb",
        "m+3": "#41f0ae",
        "m+4": "#addead",
        "m+5": "#eec046",
        "m+6": "#ffa040",
    }
    return darkbarcolor, palsD


def complexstacked_old(
        co, selectedmets, dfs_Dico, darkbarcolor, palsD,
        outfilename, figu_width, xlabyesno
):
    """plot highly custom, recommended that selectedmets <= 6 subplots"""
    ### set font style
    sns.set_style({"font.family": "sans-serif", "font.sans-serif": "Liberation Sans"})
    f, axs = plt.subplots(1, len(selectedmets), sharey=False, figsize=(figu_width, 4.8))
    plt.rcParams.update({"font.size": 20})

    for z in range(len(selectedmets)):
        # sns.set_style({ 'font.family': 'sans-serif',
        #                'font.sans-serif' : 'Liberation Sans'   })
        axs[z].set_title(selectedmets[z])
        sns.histplot(
            ax=axs[z],
            data=dfs_Dico[selectedmets[z]],
            x="timeANDcondi",
            # Use the value variable here to turn histogram counts into weighted
            # values.
            weights="Isotopologue Contribution (%)",
            hue="m+x",
            multiple="stack",
            palette=palsD,  # ['#0000c0',  '#410257' ] , # '#440154FF'],
            # Add  borders to the bars.
            edgecolor="black",
            # Shrink the bars a bit so they don't touch.
            shrink=0.85,
            alpha=1,
            legend=False,
        )
        #
        axs[z].tick_params(axis="x", labelrotation=90, labelsize=18)
        axs[z].tick_params(axis="y", length=11, labelsize=19)
        axs[z].set_ylim([0, 100])
        for bar in axs[z].patches:
            # assign stacked bars text color
            selcol = "black"
            # print(bar.get_facecolor())
            herergba = bar.get_facecolor()  # returns rgba !
            if herergba in darkbarcolor:
                selcol = "white"
            thebarvalue = round(bar.get_height(), 1)
            if thebarvalue >= 100:
                thebarvalue = 100  # no decimals if 100
            if round(bar.get_height(), 1) > 4:
                axs[z].text(
                    # Put the text in the middle of each bar. get_x returns the start
                    # so we add half the width to get to the middle.
                    bar.get_x() + bar.get_width() / 2,
                    # Vertically, add the height of the bar to the start of the bar,
                    # along with the offset.
                    (bar.get_height() / 2) + (bar.get_y()) + 2,  #
                    # This is actual value we'll show.
                    thebarvalue,
                    # Center the labels and style them a bit.
                    ha="center",
                    color=selcol,
                    size=int((figu_width / len(selectedmets)) * 2),
                )  # end axs[z].text
            else:
                continue
            # end if round(...)
        # end for bar

        axs[z].set_ylabel("", size=20)
        axs[z].xaxis.set_tick_params(length=0)  # no need of x ticks
        axs[z].set_xlabel("", size=13)
    # end for z

    [ax.invert_yaxis() for ax in axs]  # invert y, step 1

    for ax in axs:
        ylabels = ax.get_yticks().tolist()
        ax.yaxis.set_major_locator(mticker.FixedLocator(ylabels))
        ax.set_yticklabels([100 - int(i) for i in ylabels])  # invert y , step2

    if xlabyesno == "no":
        for ax in axs:
            xlabelshere = ax.get_xticks()
            ax.set_xticklabels(["" for i in xlabelshere])

    f.subplots_adjust(hspace=0.5, wspace=0.25, top=0.85, bottom=0.26, left=0.15, right=0.99)
    f.suptitle(f"compartment : {co.upper()}  \n", fontsize=20)
    f.text(0.03, 0.57, "Isotopologue Contribution (%)\n", va="center", rotation="vertical", size=20)
    f.savefig(outfilename, format="pdf")
    plt.close()
    return 0

    def saveisotopologcontriplot_old(
            datadi,
            tablePicked,
            names_compartments,
            namesuffix,
            metadata,
            selbycompD,
            darkbarcolor,
            palsD,
            condilevels,
    ):
        levelshours_str = [str(i) for i in sorted(metadata['timenum'].unique())]

    # condilevels, combined_tc_levels = add_joker_tolabs(condilevels, levelshours_str)
    condilevels, combined_tc_levels = simplelabs(condilevels, levelshours_str)

    for co in names_compartments.values():  #
        print(co)

        adf = pd.read_csv(
            datadi + tablePicked + "_" + namesuffix + "_" + co + ".tsv",
            sep="\t",
            index_col=0,
        )
        # note that pandas automatically transform any 99.9% in decimal 0.999

        dicos = dict()
        dicos[co] = {}
        dicos[co]["metadata"] = metadata.loc[metadata.short_comp == co]
        dicos[co][tablePicked] = adf[dicos[co]["metadata"]["sample"]]

        # call complicated functions

        df4plot = icontrib_2df4plot(dicos, tablePicked, co, levelshours_str)
        df4plot = massageisotopologues(df4plot)
        ####
        # conditions to plot in desired order :
        ####
        odiric = "results/plots/ic/"
        if not os.path.exists(odiric):
            os.makedirs(odiric)
        metscustomgroups = selbycompD[co]

        for j in range(len(metscustomgroups)):
            selectedmets = metscustomgroups[j]
            print(selectedmets)
            outfname = "{}ic_{}_group{}.pdf".format(odiric, co, j)
            print(outfname)
            dfs_Dico = preparemeansreplicates(df4plot, selectedmets)

            dfs_Dico = addcombinedconditime(dfs_Dico, combined_tc_levels)
            dfs_Dico.keys()  # just the metabolites subframes, one co
            figu_width = 4.3 * len(selectedmets)  # note, change width
            complexstacked_old(
                co, selectedmets, dfs_Dico, darkbarcolor, palsD, outfname,
                figu_width, xlabyesno="yes"
            )
            plt.close()
            # new :
            outfnameNoXlab = "{}ic_{}_group{}_noxlab.pdf".format(odiric, co, j)
            complexstacked_old(
                co, selectedmets, dfs_Dico, darkbarcolor, palsD, outfnameNoXlab,
                figu_width, xlabyesno="no"
            )

        # legend alone
        plt.figure()
        myhandless = []
        for c in palsD.keys():
            paobj = mpatches.Patch(facecolor=palsD[c], label=c, edgecolor="black")
            myhandless.append(paobj)
        plt.legend(handles=myhandless)
        plt.savefig(f"{odiric}ic_legend.pdf", format="pdf")

    return 0
############################################################################
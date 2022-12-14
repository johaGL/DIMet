"""
needed for DIMet/__main__.py to be able to import .py files in this location (DIMet/src/)
"""
import argparse
import os
import shutil

import numpy as np
import yaml
from .pca_fun import massage_datadf_4pca, calcPCAand2Dplot
from .differential_univariate import *
from .abund_frompercentages import calculate_meanEnri, split_mspecies_files
from .prep_steps import *
from .abundances_bars import *
from .frac_contrib_lineplot import *
from .isotopologcontrib_stacked import *
from .use_distrib_fit import *
from .fun_fm import countnan_samples, add_alerts


def dimet_message():
    return "DIMet: *D*ifferential *I*sotopically-labeled targeted *Met*abolomics\n"


parser = argparse.ArgumentParser()
parser.add_argument("--mywdir", help="working directory (with data and results subfolders)")
parser.add_argument("--config", help = "path to configuration yaml file")
parser.add_argument("--mode", help = "prepare | PCA | diffabund | abundplots \
                            | timeseries_isotopologues |  timeseries_fractional")
parser.add_argument("--stomp_values", help = "[Y/N]. Optional. To use with mode prepare.  \
                                  Stomps isotopologues' proportions to max 1.0 and min 0.0. Default : Y")
parser.add_argument("--proportion_cutoff", help = "Optional. To use with mode prepare. For isotopologues, \
                this cutoff is the min accepted value for values different than zero.  \
                Values (proportions) under cutoff become zero. Default : 0.001")

args = parser.parse_args()

if args.proportion_cutoff is None:
    PROPORTIONCUTOFF = 0.001
else:
    PROPORTIONCUTOFF = float(args.proportion_cutoff)


confifile = os.path.expanduser(args.config)
with open(confifile, "r") as f:
    confidic = yaml.load(f, Loader=yaml.Loader)

namesuffix = confidic["namesuffix"]
datadi = confidic["datadi"]
extrudf_fi = confidic["extrulist_fi"]
names_compartments = confidic["names_compartments"]
metadata_fi = confidic["metadata_fi"]
levelstimepoints_ = confidic["levelstime"]

whichtest = confidic["whichtest"]

tableIC = confidic["name_isotopologue_contribs"].split(".")[0] # no extension


tableAbund = confidic["name_abundances"].split(".")[0] # no extension
max_m_species = confidic["max_m_species"]

# set working directory as current
os.chdir(os.path.expanduser(args.mywdir))

metadata = pd.read_csv(datadi + metadata_fi, index_col=False)

isotopolog_style = autodetect_isotop_nomenclature(datadi, tableIC, namesuffix)

allfi = os.listdir(datadi)
dirtmpdata = "tmp/"
abunda_species_4diff_dir = dirtmpdata + "abufromperc/"


if args.mode == "prepare":
    print("\nPreparing dataset for analysis\n")

    print(" [any temporary files (tmp) are being deleted by default]")
    if os.path.exists(dirtmpdata):
        shutil.rmtree(dirtmpdata) # clear at each run

    os.makedirs(dirtmpdata)

    tsvfi = [i for i in allfi if ".tsv" in i]
    print("Your .tsv files in data folder: ", tsvfi, "\n")

    # using list of metabolites to exclude, compartment aware:
    print("using list of undesired metabolites to drop off from data\t")
    extrudf = pd.read_csv(datadi + extrudf_fi, sep=",")

    stomp_values = args.stomp_values
    if stomp_values == 'N':
        print("are you sure not stomping values ? \
        \naberrant proportion values (negative and superior to 1)  will remain,"
              "and good quality analysis is not guaranteed")
    for filename in tsvfi:
        save_new_dfsB(datadi, names_compartments, filename, metadata, extrudf,
                      dirtmpdata, isotopolog_style, stomp_values, PROPORTIONCUTOFF)


    # NOTE : for abundances bars and Differential,
    # compulsory step: calculate isotopologues abundances from IC percentages
    calculate_meanEnri(
        dirtmpdata,
        tableAbund,
        tableIC,
        metadata,
        names_compartments,
        namesuffix,
        dirtmpdata
    )

    split_mspecies_files(dirtmpdata, names_compartments, namesuffix,
                   abunda_species_4diff_dir)

    if detect_fraccontrib_missing(tsvfi) is False:
        print("Warning !: you do not have fractional contributions file")

    print("\nsplited (by compartment) and clean files in tmp/ ready for analysis\n")


if args.mode == "PCA":
    #picked_for_pca = "meanEnrich"  # TODO: allow to pick Abundance or meanEnrich
    picked_for_pca = tableAbund
    odirpca = "results/plots/pca/"
    if not os.path.exists(odirpca):
        os.makedirs(odirpca)
    print(f"\n plotting pca(s) to: {odirpca}\n")
    for k in names_compartments.keys():
        co = names_compartments[k]
        file4pca = f"{dirtmpdata}{picked_for_pca}_{namesuffix}_{co}.tsv"
        df = pd.read_csv(file4pca, sep='\t', header=0, index_col=0)
        metadatasub = metadata.loc[metadata['short_comp'] == co, :]
        dfa = massage_datadf_4pca(df, metadatasub)
        pc_df, dfvare = calcPCAand2Dplot(dfa, metadatasub, "timepoint", "condition",
                         "", f'{picked_for_pca}-{namesuffix}-{k}', odirpca, 6)
        pc_df, dfvare = calcPCAand2Dplot(dfa, metadatasub, "timepoint", "condition",
                         "sample_descrip", f'{picked_for_pca}-{namesuffix}-{k}', odirpca, 6)
        for tp in levelstimepoints_:
            metadatasub = metadata.loc[(metadata['short_comp'] == co) & (metadata['timepoint'] == tp), :]
            dfb = massage_datadf_4pca(df, metadatasub)
            pc_df, dfvare = calcPCAand2Dplot(dfb, metadatasub, "condition", "condition",
                             "sample_descrip", f'{picked_for_pca}-{namesuffix}-{k}-{tp}', odirpca, 6)


if args.mode == "diffabund" and whichtest != "disfit":
    print("\n testing for Differentially Abundant Metabolites [or Isotopologues] : DAM\n")
    spefiles = [i for i in os.listdir(abunda_species_4diff_dir)]


    newcateg = confidic["newcateg"]  # see yml in example/configs/
    contrasts_ = confidic["contrasts"]

    outdiffdir = "results/tables/"
    if not os.path.exists(outdiffdir):
        os.makedirs(outdiffdir)
    outputsubdirs = ["m+"+str(i)+"/" for i in range(max_m_species+1)]
    outputsubdirs.append("totmk/")
    outputsubdirs.append("TOTAL/")
    alloutdirs = list()
    for exte_sig in ["extended/", "significant/"]:
        for subdir_spec in outputsubdirs:
            x = outdiffdir + exte_sig + subdir_spec
            alloutdirs.append(x)
            if not os.path.exists(x):
                os.makedirs(x) # each m+x output directory

    outdirs_total_abund_res_ = [d for d in alloutdirs if "TOTAL" in d]
    for contrast in contrasts_:
        print("\n    comparison ==>", contrast[0] ,"vs",contrast[1] )
        for co in names_compartments.values():
            rundiffer(
                dirtmpdata,
                tableAbund,
                namesuffix,
                metadata,
                newcateg,
                contrast,
                whichtest,
                co,
                outdirs_total_abund_res_,
                "TOTAL",
            )
    
            tableabuspecies_co_ = [i for i in spefiles if co in i]
            # any "m+x" where x > max_m_species, must be excluded
            donotuse = [ k for k in tableabuspecies_co_ if "m+" in k.split("_")[2]
                        and int(k.split("_")[2].split("+")[-1]) > max_m_species ]
            tabusp_tmp_ = set(tableabuspecies_co_) - set(donotuse)
            tableabuspecies_co_good_ = list(tabusp_tmp_)
            for tabusp in tableabuspecies_co_good_:
                outkey = tabusp.split("_")[2]  # the species m+x as saved
                outdiffdirs = [d for d in alloutdirs if outkey in d]
                rundiffer(
                    abunda_species_4diff_dir,
                    tabusp,
                    namesuffix,
                    metadata,
                    newcateg,
                    contrast,
                    whichtest,
                    co,
                    outdiffdirs,
                    outkey
                )
                # end for tabusp
            # end for co
        # end for contrast
    print("\nended differential analysis")
    # end if args.mode == "diffabund"


if args.mode == "abundplots":
    odirbars = "results/plots/abundbars/"
    if not os.path.exists(odirbars):
        os.makedirs(odirbars)
    xticks_text_ = confidic["xticks_text"]
    axisx_labeltilt = confidic["axisx_labeltilt"]
    time_sel = confidic["time_sel"]  # locate where it is used
    selectedmetsD = confidic["selectedmets_forbars"]  # locate where it is used
    condilevels = confidic["conditions"]  # <= locate where it is used
    vizorder = confidic["axisx_barcolor"]
    col1 = vizorder[0]
    col2 = vizorder[1]
    # in a first time print the TOTAL abundances, selectedmets_forbars
    for CO in names_compartments.values():
        file_total_co_ = [i for i in os.listdir(dirtmpdata) if tableAbund in i and CO in i]
        print(file_total_co_)

        abutab = pd.read_csv(dirtmpdata + file_total_co_[0], sep="\t", index_col=0)
        metada_sel = metadata.loc[metadata["sample"].isin(abutab.columns), :]

        # metadata and abundances time of interest
        metada_sel = metada_sel.loc[metadata["timepoint"].isin(time_sel), :]
        abu_sel = abutab[metada_sel["sample"]]

        # total piled-up data:
        piled_sel = stackallabundace(abu_sel, metada_sel)
        piled_sel["condition"] = pd.Categorical(piled_sel["condition"], condilevels)
        piled_sel["timepoint"] = pd.Categorical(piled_sel["timepoint"], time_sel)

        plotwidth = 4 * len(selectedmetsD[CO])
        print(f"sending to plot file  :  {selectedmetsD[CO]}")
        printabundbarswithdots(piled_sel, selectedmetsD[CO], CO, "TOTAL",
                               col1, col2, plotwidth, odirbars, xticks_text_ , axisx_labeltilt)


if args.mode == "timeseries_fractional":
    print(" Fractional contributions plots \n")
    tableFC = confidic["name_fractional_contribs"].split(".")[0] # no extension
    gbycompD = confidic["groups_toplot_frac_contribs"]

    def yieldcolorsbymet(): # TODO: get out colors from here, set them in yaml
        # coloreachmetab dictionary: contains many metabolites, just defines colors.
        coloreachmetab = {
            "L-Lactic_acid": "blue",
            "Citric_acid": "green",
            "Oxoglutaric_acid": "#903ee0",
            "Succinic_acid": "#019085",
            "Fumaric_acid": "#01668a",
            "L-Malic_acid": "#afc98f",
            "L-Alanine": "blueviolet",
            "Pyruvic_acid": "#33d3e3",
            "L-Glutamic_acid": "#643561",
            "L-Glutamine": "#950065",
            "L-Aspartic_acid": "#9d6d00",
            "L-Asparagine": "#ff7c7c",  # greenfluo : '#11dc79'
            "Stearic_acid" : "gray",
            "Palmitic_acid" : "orange",
            "Coenzyme_A" : "gray",
            "Acetyl-CoA" : "orange"
        }
        return coloreachmetab
    coloreachmetab = yieldcolorsbymet()

    savefraccontriplots(dirtmpdata, names_compartments,
                        metadata, tableFC, namesuffix,
                        gbycompD, coloreachmetab)


if args.mode == "timeseries_isotopologues":
    print(" Isotopologue's Contributions plots \n")

    condilevels = confidic["conditions"]  # <= locate where it is used

    #darkbarcolor, palsD = custom_colors_stacked()
    selbycompD = confidic["groups_toplot_isotopol_contribs"]
    # saveisotopologcontriplot_old(dirtmpdata, tableIC, names_compartments,
    #                           namesuffix, metadata, selbycompD,
    #                          darkbarcolor, palsD, condilevels )
    for co in selbycompD.keys():
        #mets_byco = get_metabolites(tableIC) # TODO: make this function
        for group in selbycompD[co].keys():
            pass
            #print([met for met in selbycompD[co][group]])
            #notfound = set([met for met in group]) - set(mets_byco)
    saveisotopologcontriplot(  dirtmpdata,
    tableIC,
    names_compartments,
    namesuffix,
    metadata,
    selbycompD,
    condilevels )


def save_each_df(good_df, bad_df, outdiffdir,
                 co, autochoice, strcontrast):
    rn = f"{outdiffdir}/extended/{autochoice}"
    good_o = f"{rn}/{co}_{autochoice}_{strcontrast}_good.tsv"
    bad_o = f"{rn}/{co}_{autochoice}_{strcontrast}_bad.tsv"
    good_df.to_csv(good_o, sep='\t', header=True)
    bad_df.to_csv(bad_o, sep='\t', header=True)
    return "saved to results"



if args.mode == "diffabund" and whichtest == "disfit":
    print("\nDistribution fitting (of ratios)\n")
    newcateg = confidic["newcateg"]
    contrasts_ = confidic["contrasts"]
    outdiffdir = "results/tables/"

    if not os.path.exists(outdiffdir):
        os.makedirs(outdiffdir)

    if not os.path.exists(f'{dirtmpdata}/preDiff/'):
        os.makedirs(f'{dirtmpdata}/preDiff/')

    spefiles = [i for i in os.listdir(abunda_species_4diff_dir)]

    if not os.path.exists(f'{outdiffdir}extended/TOTAL/'):
        os.makedirs(f'{outdiffdir}extended/TOTAL/')

    for contrast in contrasts_:
        strcontrast = "_".join(contrast)
        print("\n    comparison ==>", contrast[0] ,"vs",contrast[1] )


        for co in names_compartments.values():
            autochoice = "TOTAL"
            filehere = f"{dirtmpdata}{tableAbund}_{namesuffix}_{co}.tsv"

            df = pd.read_csv(filehere, sep='\t', header=0, index_col=0)

            metada_sel = metadata.loc[metadata['short_comp'] == co, :]

            df4c, metad4c = prepare4contrast(df, metada_sel, newcateg, contrast)

            # delete rows being zero everywhere in this TOTAL df
            df4c = df4c[(df4c.T != 0).any()]

            # sort them by 'newcol' the column created by prepare4contrast
            metad4c = metad4c.sort_values("newcol")
            minimal_val = df4c[df4c > 0].min().min()
            print(minimal_val)
            # replace zeroes by that minimal val
            df4c[df4c==0] = minimal_val
            ratiosdf = calcs_red_to_ratios(df4c,
                                        metad4c,
                                        contrast )


            ratiosdf = add_alerts(ratiosdf, metad4c)

            ratiosdf["compartment"] = co

            # print("fitting to distributions to find the best ... ")
            # **
            pre_out = f"{co}_{autochoice}_{strcontrast}_prep_fit.tsv"
            ratiosdf.to_csv(f"{dirtmpdata}/preDiff/{pre_out}",
                                index_label="metabolite", header=True, sep='\t')


            good_df, bad_df = split_byalert_df(ratiosdf)
            # good_df = compute_z_score(good_df)
            # fiting et blahblah
            save_each_df(good_df, bad_df, outdiffdir,
                             co, autochoice, strcontrast)

            # out_histo_file = f"{outdiffdir}/extended/{autochoice}/\
            #                     {co}_{autochoice}_{strcontrast}_fitdist_plot.pdf"
            # best_distribution, args_param = find_best_distribution(good_df,
            #                           out_histogram_distribution= out_histo_file)
            # argsided = "two-sided" #"two-sided" # or "rigth-tailed"
            # good_df = compute_p_value(good_df, argsided, best_distribution, args_param)
            # good_df = compute_p_adjusted(good_df, "fdr_bh")
            #
            # final_total_diff = good_df.copy()
            #
            # fout = f"{autochoice}/{co}_{autochoice}_{strcontrast}_fitted.tsv"
            # final_total_diff.to_csv(f"{outdiffdir}extended/{fout}",
            #                   index_label="metabolite", header=True, sep='\t')



            # --------------------- for isotopologues ---------------------------:

            tableabuspecies_co_ = [i for i in spefiles if co in i]
            # any "m+x" where x > max_m_species, must be excluded
            donotuse = [k for k in tableabuspecies_co_ if "m+" in k.split("_")[2]
                        and int(k.split("_")[2].split("+")[-1]) > max_m_species]
            tabusp_tmp_ = set(tableabuspecies_co_) - set(donotuse)
            tableabuspecies_co_good_ = list(tabusp_tmp_)

            isos_togetherD = dict()
            for tabusp in tableabuspecies_co_good_:
                autochoice = tabusp.split("_")[2]  # the species m+x as saved
                filehere = f"{abunda_species_4diff_dir}{tabusp}"
                print(filehere)

                df = pd.read_csv(filehere, sep='\t', header=0, index_col=0)

                metada_sel = metadata.loc[metadata['short_comp'] == co, :]

                df4c, metad4c = prepare4contrast(df, metada_sel, newcateg, contrast)
                df4c = df4c[(df4c.T != 0).any()]  # delete rows being zero all

                print(df4c.shape, " <<<<< ")
                # sort them by 'newcol' the column created by prepare4contrast
                metad4c = metad4c.sort_values("newcol")

                indexfull = [f'{m}_{autochoice}' for m in df4c.index]
                df4c.index = indexfull
                # ratiosdf = ratiosdf.assign(isotype=[autochoice for k in range(ratiosdf.shape[0])])
                # out_histo_file = f"results/plots/distrib-{strcontrast}-{autochoice}-{co}.pdf"
                # #fout = f"{autochoice}/{co}_{autochoice}_{strcontrast}_fitted.tsv"
                # #ratiosdf.to_csv(f"{outdiffdir}extended/{fout}", header=True, sep='\t')

                df4c = df4c.assign(isotype=[autochoice for k in range(df4c.shape[0])])

                if autochoice.startswith("m+"):
                    isos_togetherD[f"{autochoice}-{co}"] = df4c
                elif autochoice == "totmk" :
                    if confidic['also_total_marked'] == "Yes":
                        isos_togetherD[f"{autochoice}-{co}"] = df4c  # ratiosdf


                    elif confidic['also_total_marked'] == "No":
                        pass # total marked is not of interest for user

            # new: pile the isotopologues dataframes
            isos_piledup = pd.concat(isos_togetherD.values(), axis=0)


            # this dataframe has a non numeric column 'isotype', so select samples :
            tmp = isos_piledup[metad4c['sample']]
            minimal_val = tmp[tmp > 0].min().min()

            tmp.replace(to_replace=0, value=minimal_val, inplace= True)


            ratiosdf2 = calcs_red_to_ratios(tmp,
                                           metad4c,
                                           contrast)

            col_isotype = isos_piledup['isotype'].tolist()
            ratiosdf2 = ratiosdf2.assign(isotype=col_isotype)
            print(ratiosdf2.shape, "!*")


            ratiosdf2 = add_alerts(ratiosdf2, metad4c)

            ratiosdf2["compartment"] = co



            plt.figure(figsize=(4,4))
            sns.histplot(x=ratiosdf2['ratio'])
            plt.title(f"all isotopologues ratios {co} {strcontrast}")
            plt.savefig(f"ratios_isotopologues-{co}-{strcontrast}.pdf")
            plt.close()
            ratiosdf2.to_csv(f"{dirtmpdata}/preDiff/{co}_m+x_{strcontrast}_prep_fit.tsv",
                                index_label="isotopologue", header=True, sep='\t')

            #
            # # only do fitting in those not having alerts
            good_df, bad_df = split_byalert_df(ratiosdf2)
            good_df = compute_z_score(good_df)
            #
            rn = f"{outdiffdir}/extended/"
            good_o = f"{rn}/{co}_m+x_{strcontrast}_good.tsv"
            bad_o = f"{rn}/{co}_m+x_{strcontrast}_bad.tsv"
            good_df.to_csv(good_o, sep='\t', header=True)
            bad_df.to_csv(bad_o, sep='\t', header=True)


            # # do fitting on the good_df if metabolites number >= ??
            out_histo_file = f"{outdiffdir}/extended/{co}_m+x_{strcontrast}_fitdist_plot.pdf"
            best_distribution, args_param = find_best_distribution(good_df,
                                                                    out_histogram_distribution=out_histo_file)
            argsided = "right-tailed"  # "two-sided" # or "rigth-tailed"
            good_df = compute_p_value(good_df, argsided, best_distribution, args_param)
            good_df = compute_p_adjusted(good_df, "fdr_bh")
            #
            print("!!")
            #final_total_diff = good_df.copy()

            fout = f"/{co}_m+x_{strcontrast}_fitted.tsv"
            good_df.to_csv(f"{outdiffdir}extended/{fout}",
                                   index_label="metabolite", header=True, sep='\t')

    dico_colnames = dict()
    for fi in os.listdir(f'{dirtmpdata}preDiff/'):
        print(fi)

    """
    vec = 
    crazy = pd.DataFrame({'foo' : [i for i in range(len(vec))],
                          'zscore' : vec})
    plt.figure(figsize=(4, 4))
    sns.histplot(x=crazy['ratio'])
    plt.title(f"all ratios {strcontrast}")
    plt.savefig(f"ratios-{strcontrast}.pdf")
    """









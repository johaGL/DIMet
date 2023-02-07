"""
needed for DIMet/__main__.py to be able to import .py files in this location (DIMet/src/)
@author: johanna
"""
import argparse
import shutil
import yaml
from .pca_fun import massage_datadf_4pca, calcPCAand2Dplot
from .differential_univariate import *
from .abund_frompercentages import calc_isos_absolute, split_mspecies_files
from .prep_steps import *
from .abundances_bars import *
from .frac_contrib_lineplot import *
from .isotopologcontrib_stacked import *
from .functions_diffmode import *
from .fun_fm import countnan_samples, add_alerts, add_metabolite_column,\
     add_isotopologue_type_column, save_heatmap_sums_isos, \
     givelevels, table_minimalbymet, save_rawisos_plot

from .metabologram import metabologram_run


def dimet_message():
    return "DIMet: *D*ifferential *I*sotopically-labeled targeted *Met*abolomics\n"


parser = argparse.ArgumentParser()
parser.add_argument("--mywdir", help="working directory (with data and results subfolders)",
                    required=True)

parser.add_argument("--config", help="path to configuration yaml file", required=True)

parser.add_argument("--isotopologue-preview", action=argparse.BooleanOptionalAction,
                    default=True, dest="isotopreview",
                    help="Initial visualization of isotopologues' values across samples. NOT to use with --mode")

parser.add_argument("--mode",  help="prepare | PCA | diffabund | abundplots \
                            | timeseries_isotopologues |  timeseries_fractional | metabologram")

parser.add_argument("--absolute-isotopologues-available", default="N",
                    dest="absolute_isos_avail", help="[Y/N]. To use with mode prepare.\
                      Y if you have  absolute values (not in percentage) for isotopologues. Default : N")

parser.add_argument("--isos_detect_cutoff", default=-0.05, type= float,
                    help="With mode prepare. When in percentages, metabolites whose isotopologues exhibit\
                     proportions less or equal this cutoff are removed (neg can occur). Default: -0.05")

parser.add_argument("--stomp_values", default="Y", help="[Y/N]. To use with mode prepare. \
                    Stomps isotopologues' proportions to max 1.0 and min 0.0. Default : Y")

parser.add_argument("--zero_to_nan", default=0, type=int,
                    help="Optional, with mode prepare. Only keep zero if entire group is zeros\
                    i.e. given 3 bio replicates:[0 22 0] change to [NaN 22 NaN].\
                    Will be performed in abundances, isotopologue and rawintensity. \
                     Default : 0 (not to perform) ")

parser.add_argument("--totalAbund_zero_replace", default="min", type=str,
                    help="To use with mode diffabund, treat zeroes before reduction and gmean. \n\
                         min : replaces zero by the minimal value of the entire dataframe (recommended)\n. \
                         Or a number to replace it with (eg: 1e-10, 1e-3, 3, 3000) if you know well your data. \
                         Default : min")

parser.add_argument("--isotoAbsolutes_zero_replace", default="min", type=str,
                    help="To use with mode diffabund, treat isotopologue zeroes before reduction and gmean.\
                         set to 'min' or 'min/n' where n can be integer or decimal (recommended);\
                         or a number to replace it with (eg: 1e-10, 1e-3, 3, 3000) if you know well your data. \
                         Default : min" )

parser.add_argument("--qualityDistanceOverSpan", default=-0.3, type=float,
                    help="With mode diffabund. For each metabolite, for samples (x, y)  the\
                        distance is calculated, and span is max(x U y) - min(x U y).  A 'distance/span' \
                        inferior to this value excludes the metabolite from testing (pvalue NaN).")

parser.add_argument("--padj", default=0.05, type=float, help="padj filter for diffabund results")
parser.add_argument("--log2FC", default=0.5, type=float, help="abs(log2FC) filter for diffabund results")
parser.add_argument("--DistanceOverSpan", default=0.15, type=float,
                    help="'distance/span' filter for diffabund results")

args = parser.parse_args()

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
metadata['timepoint'] = metadata['timepoint'].astype(str) # new
isotopolog_style = autodetect_isotop_nomenclature(datadi, tableIC, namesuffix)

allfi = os.listdir(datadi)
dirtmpdata = "tmp/"
abunda_species_4diff_dir = dirtmpdata + "abufromperc/"


if args.isotopreview and (args.mode is None):
    print("This option is dedicated to isotopologues in percentages when no absolute values given.",
          "Do not use for isotopologues absolute values")
    print("\nsending new files to isotop_preview/")
    if os.path.exists("isotop_preview/"):
        shutil.rmtree("isotop_preview/")  # clear at each run

    os.makedirs("isotop_preview/")

    extrudf = pd.read_csv(datadi + extrudf_fi, sep=",")

    fn = tableIC + "_" +  namesuffix + ".tsv"
    # set specia_zero_tonan to False as here we need original values unaltered
    save_new_dfsB(datadi, names_compartments, fn, metadata, extrudf, "isotop_preview/",
                  isotopolog_style, stomp_values='N', special_zero_tonan=False) # ok

    for long_compartment_name in names_compartments.keys():
        k = names_compartments[long_compartment_name]
        fnn = f"isotop_preview/{tableIC}_{namesuffix}_{k}.tsv"
        df = pd.read_csv(fnn, sep='\t', header=0, index_col=0)
        df = add_metabolite_column(df)
        df = add_isotopologue_type_column(df)

        thesums = df.groupby('metabolite').sum()
        thesums = thesums.drop(columns=['isotopologue_type'])
        thesums = thesums.round(3)
        ff = f"isotop_preview/sums_Iso_{k}.pdf"
        figuretitle = f"Sums of isotopologue proportions ({k}) "
        save_heatmap_sums_isos(thesums, figuretitle, ff)

        dfmelt = pd.melt(df, id_vars=['metabolite', 'isotopologue_type'])
        dfmelt = givelevels(dfmelt)
        table_minimalbymet(dfmelt, f"isotop_preview/minextremesIso_{k}.tsv")
        outputfigure = f"isotop_preview/allsampleIsos_{k}.pdf"
        figtitle = f"{long_compartment_name} compartment, all samples"
        save_rawisos_plot(dfmelt, figuretitle=figtitle, outputfigure=outputfigure)

    dfDict = dict()
    for k in names_compartments.values():
        df = pd.read_csv(f"isotop_preview/minextremesIso_{k}.tsv",
                         sep='\t', header=0, index_col=None)
        df['compartment'] = k
        df = df.drop(columns=['isotopologue_type', 'variable'])
        dfDict[k] = df

    joineddf = pd.concat(dfDict.values())
    joineddf.to_csv(f"isotop_preview/minextremesIso.csv",
                 header=True)


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

    # if any other list generated in the isotop_preview folder:
    try:
        extremeisos = pd.read_csv("isotop_preview/minextremesIso.csv",
                                  header=0, index_col=0)
        print("\nUser argument isos_detect_cutoff =", args.isos_detect_cutoff)
        extremeisos = extremeisos.loc[extremeisos['value'] <= args.isos_detect_cutoff, :]
        print("also the following metabolites will be excluded : ")
        extrtocat = extremeisos[['metabolite', 'compartment']]
        print(extrtocat)
        extrudf = pd.concat([extrudf, extrtocat], axis = 0)

    except Exception as e:
        pass

    print("\nWhole list of metabolites being excluded : ")
    extrudf.drop_duplicates()
    extrudf = extrudf.sort_values(by='compartment')
    print(extrudf)

    if args.absolute_isos_avail == 'Y':
        print("your table", tableIC, "is in absolute values")
        absolute_IC = tableIC
        percent_IC = "calcPercent" # TODO: add function to calc % (from given absolute) and saveit to data
        for filename in tsvfi :
            save_new_dfsB(datadi, names_compartments, filename, metadata, extrudf, dirtmpdata,
                          isotopolog_style, stomp_values=False,  special_zero_tonan=bool(args.zero_to_nan))

    elif args.absolute_isos_avail == 'N': # default
        if args.stomp_values == 'N':
            print("Warning : are you sure not stomping values ? \
            \naberrant proportion values (negative and superior to 1)  will remain,"
                  "and good quality analysis is not guaranteed")
        for filename in tsvfi:
            save_new_dfsB(datadi, names_compartments, filename, metadata, extrudf, dirtmpdata,
                          isotopolog_style, args.stomp_values,  bool(args.zero_to_nan))

        absolute_IC = calc_isos_absolute(dirtmpdata, tableAbund, tableIC, metadata,
                                         names_compartments, namesuffix, dirtmpdata)
        percent_IC = tableIC

    else:
        print("argument absolute-isotopologues-available not recognized")
    try:
        split_mspecies_files(dirtmpdata, names_compartments, namesuffix,
                        abunda_species_4diff_dir)
    except:
        print("[Note FOR DEV: something WRONG with split_mspecies_files]")
    if detect_fraccontrib_missing(tsvfi) is False:
        print("Warning !: you do not have fractional contributions file")

    detect_and_create_dir("results/")
    save_mini_report([tableAbund, absolute_IC, percent_IC],
                     namesuffix, names_compartments,  dirtmpdata)


    print("\nsplited (by compartment) and clean files in tmp/ ready for analysis\n")



if args.mode == "PCA":
    #picked_for_pca = "meanEnrich"  # TODO: allow to pick Abundance or meanEnrich, first fix meanEnrich
    picked_for_pca = tableAbund
    nbcomps= 6 # TODO make it an option from user config
    odirpca = "results/plots/pca/"
    detect_and_create_dir(odirpca)
    print("hahahah")
    advanced_test = False # for dev test

    print(f"\n plotting pca(s) to: {odirpca}\n")
    for k in names_compartments.keys():
        co = names_compartments[k]
        file4pca = f"{dirtmpdata}{picked_for_pca}_{namesuffix}_{co}.tsv"
        df = pd.read_csv(file4pca, sep='\t', header=0, index_col=0)
        metadatasub = metadata.loc[metadata['short_comp'] == co, :]
        dfa = massage_datadf_4pca(df, metadatasub, advanced_test)
        pc_df, dfvare = calcPCAand2Dplot(dfa, metadatasub, "timepoint", "condition",
                         "", f'{picked_for_pca}-{namesuffix}-{k}', odirpca, nbcomps)
        pc_df, dfvare = calcPCAand2Dplot(dfa, metadatasub, "timepoint", "condition",
                         "sample_descrip", f'{picked_for_pca}-{namesuffix}-{k}', odirpca, nbcomps)
        for tp in levelstimepoints_:
            metadatasub = metadata.loc[(metadata['short_comp'] == co) & (metadata['timepoint'] == tp), :]
            dfb = massage_datadf_4pca(df, metadatasub, advanced_test)
            pc_df, dfvare = calcPCAand2Dplot(dfb, metadatasub, "condition", "condition",
                             "sample_descrip", f'{picked_for_pca}-{namesuffix}-{k}-{tp}', odirpca, nbcomps)


if args.mode == "abundplots":
    odirbars = "results/plots/abundbars/"
    detect_and_create_dir(odirbars)
    xticks_text_ = confidic["xticks_text"]
    axisx_labeltilt = confidic["axisx_labeltilt"]
    time_sel = confidic["time_sel"]  # locate where it is used
    selectedmetsD = confidic["selectedmets_forbars"]  # locate where it is used
    condilevels = confidic["conditions"]  # <= locate where it is used
    vizorder = confidic["axisx_barcolor"]
    col1 = vizorder[0]
    col2 = vizorder[1]

    width_each_subfig = float(confidic["width_each_abu"])
    wspace_subfigs = float(confidic["wspace_abus"])


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

        plotwidth =  width_each_subfig * len(selectedmetsD[CO])
        print(f"sending to plot file  :  {selectedmetsD[CO]}")
        printabundbarswithdots(piled_sel, selectedmetsD[CO], CO, "TOTAL",
                               col1, col2, plotwidth, odirbars, xticks_text_ , axisx_labeltilt,
                               wspace_subfigs)


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
    width_each_stack = float(confidic["width_each_stack"])
    wspace_stacks = float(confidic["wspace_stacks"])
    numbers_size = confidic["numbers_size"]

    selbycompD = confidic["groups_toplot_isotopol_contribs"]

    if args.absolute_isos_avail == 'N':
        percent_IC = tableIC
    else:
        percent_IC = "calcPercent" # the suffix of the saved percentage (by prep mode)

    saveisotopologcontriplot(dirtmpdata,
                                percent_IC,
                                names_compartments,
                                namesuffix,
                                metadata,
                                selbycompD,
                                condilevels, width_each_stack, wspace_stacks, numbers_size )


def steps_fitting_method(ratiosdf, out_histo_file):
    ratiosdf = compute_z_score(ratiosdf)
    best_distribution, args_param = find_best_distribution(ratiosdf,
                                 out_histogram_distribution=out_histo_file)
    autoset_tailway = auto_detect_tailway(ratiosdf, best_distribution, args_param)
    print("auto, best pvalues calculated :", autoset_tailway)
    ratiosdf = compute_p_value(ratiosdf, autoset_tailway, best_distribution, args_param)

    return ratiosdf


if args.mode == "diffabund" :
    print("\n  -*- searching for Differentially Abundant Metabolites[or Isotopologues] (DAM) -*-\n")
    newcateg = confidic["newcateg"]
    contrasts_ = confidic["contrasts"]
    outdiffdir = "results/tables/"

    detect_and_create_dir(outdiffdir)

    # detect_and_create_dir(f'{dirtmpdata}/preDiff/')

    spefiles = [i for i in os.listdir(abunda_species_4diff_dir)]

    detect_and_create_dir(f'{outdiffdir}/extended/TOTAL/')
    detect_and_create_dir(f'{outdiffdir}/extended/isos/')

    detect_and_create_dir(f'{outdiffdir}/filtered/TOTAL/')
    detect_and_create_dir(f'{outdiffdir}/filtered/isos/')

    for contrast in contrasts_:
        strcontrast = "_".join(contrast)
        print("\ncomparison ==>", contrast[0] ,"vs",contrast[1] , '\n--------------')

        for co in names_compartments.values():
            print(". compartment:", co)
            autochoice = "TOTAL"
            print(".. ", autochoice)
            filehere = f"{dirtmpdata}{tableAbund}_{namesuffix}_{co}.tsv"

            df = pd.read_csv(filehere, sep='\t', header=0, index_col=0)

            metada_sel = metadata.loc[metadata['short_comp'] == co, :]

            df4c, metad4c = prepare4contrast(df, metada_sel, newcateg, contrast)

            # delete rows being zero everywhere in this TOTAL df
            df4c = df4c[(df4c.T != 0).any()]

            # sort them by 'newcol' the column created by prepare4contrast
            metad4c = metad4c.sort_values("newcol")
            if args.totalAbund_zero_replace == "min":
                minimal_val = df4c[df4c > 0].min().min()
                print("using minimal value to replace zeroes :", minimal_val)
            else :
                try:
                    minimal_val = float(args.totalAbund_zero)
                except:
                    print("Warning: --totalAbund_zero_replace argument not 'min' nor numeric, see help")
                    minimal_val = df4c[df4c > 0].min().min()

            df4c = df4c.replace(to_replace=0, value=minimal_val)
            df4c = calc_reduction(df4c, metad4c, contrast)
            df4c = countnan_samples(df4c, metad4c)  # adds nan_count_samples column
            df4c = distance_or_overlap(df4c, metad4c, contrast)
            df4c = compute_span_incomparison(df4c, metad4c, contrast)
            df4c['distance/span'] = df4c.distance.div(df4c.span_allsamples)
            ratiosdf = calc_ratios(df4c,  metad4c, contrast)
            ratiosdf, df_bad = divide_before_stats(ratiosdf, args.qualityDistanceOverSpan)  #HERE
            # specific test as defined by user
            if whichtest == "disfit" :
                out_histo_file = f"{outdiffdir}/extended/{autochoice}/" \
                                 f"{co}_{autochoice}_{strcontrast}_fitdist_plot.pdf"
                ratiosdf = steps_fitting_method(ratiosdf, out_histo_file)
                ratiosdf = compute_padj_version2(ratiosdf, 0.05, "fdr_bh")
            else:
                extract_test_df = outStat_df(ratiosdf, metad4c, contrast, whichtest)
                extract_test_df = compute_padj_version2(extract_test_df, 0.05, "fdr_bh")
                extract_test_df.set_index("metabolite", inplace=True)
                ratiosdf = pd.merge(ratiosdf, extract_test_df, left_index=True, right_index=True)

            ratiosdf["log2FC"] = np.log2(ratiosdf['ratio'])
            df_bad = complete_columns_for_bad(df_bad, ratiosdf) #HERE
            if df_bad.shape[0] >= 1:
                ratiosdf = pd.concat([ratiosdf, df_bad])  # HERE

            ratiosdf["compartment"] = co


            ratiosdf.to_csv(f"{outdiffdir}/extended/{autochoice}/{co}_{autochoice}_{strcontrast}_{whichtest}.tsv",
                            index_label="metabolite", header=True, sep='\t')
            #HERE :
            filtered_df = filter_diff_results(ratiosdf, args.padj, args.log2FC, args.DistanceOverSpan)
            filfi = f"{outdiffdir}/filtered/{autochoice}/{co}_{autochoice}_{strcontrast}_{whichtest}_filter.tsv"
            filtered_df.to_csv(filfi, index_label="metabolite", header=True, sep='\t' )
            del(ratiosdf, filfi)


            # --------------------- for isotopologues ---------------------------:

            print(".. Isotopologues")
            tableabuspecies_co_ = [i for i in spefiles if co in i]
            donotuse = [k for k in tableabuspecies_co_ if "m+" in k.split("_")[2]
                        and int(k.split("_")[2].split("+")[-1]) > max_m_species]
            tabusp_tmp_ = set(tableabuspecies_co_) - set(donotuse)
            tableabuspecies_co_good_ = list(tabusp_tmp_)

            isos_togetherD = dict()
            for tabusp in tableabuspecies_co_good_:
                autochoice = tabusp.split("_")[2]  # the species m+x as saved
                filehere = f"{abunda_species_4diff_dir}{tabusp}"

                df = pd.read_csv(filehere, sep='\t', header=0, index_col=0)

                metada_sel = metadata.loc[metadata['short_comp'] == co, :]

                df4c, metad4c = prepare4contrast(df, metada_sel, newcateg, contrast)
                df4c = df4c[(df4c.T != 0).any()]  # delete rows being zero all

                # sort them by 'newcol' the column created by prepare4contrast
                metad4c = metad4c.sort_values("newcol")

                indexfull = [f'{m}_{autochoice}' for m in df4c.index]
                df4c.index = indexfull

                df4c = df4c.assign(isotype=[autochoice for k in range(df4c.shape[0])])

                if autochoice.startswith("m+"):
                    isos_togetherD[f"{autochoice}-{co}"] = df4c
                elif autochoice == "totmk" :
                    if confidic['also_total_marked'] == "Yes":
                        isos_togetherD[f"{autochoice}-{co}"] = df4c
                    elif confidic['also_total_marked'] == "No":
                        pass # total marked is not of interest for user

            isos_piledup = pd.concat(isos_togetherD.values(), axis=0)

            tmp = isos_piledup[metad4c['sample']]

            if args.isotoAbsolutes_zero_replace == "min":
                mini_val_isos = tmp[tmp > 0 ].min().min()
            elif args.isotoAbsolutes_zero_replace.startswith("min/"):
                try:
                    denominator_str = str(args.isotoAbsolutes_zero_replace.split("/")[1])
                    mini_val_isos = tmp[tmp > 0].min().min() / float(denominator_str)
                except ValueError:
                    mini_val_isos = tmp[tmp > 0].min().min() / 2
                    print("Bad denominator --isotoAbsolutes_zero_replace, defaulting to ", 2)
            else:
                try:
                    mini_val_isos = float(args.isotoAbsolutes_zero_replace)
                except:
                    "--isotoAbsolutes_zero_replace argument value not recognized"
            # end if

            tmp = tmp.replace(to_replace=0, value=mini_val_isos)
            tmp = calc_reduction(tmp, metad4c, contrast)
            tmp = countnan_samples(tmp, metad4c)
            tmp = distance_or_overlap(tmp, metad4c, contrast)
            tmp = compute_span_incomparison(tmp, metad4c, contrast)
            tmp['distance/span'] = tmp.distance.div(tmp.span_allsamples)
            ratiosdf2 = calc_ratios(tmp, metad4c, contrast)
            ratiosdf2, df_bad = divide_before_stats(ratiosdf2, args.qualityDistanceOverSpan)  # HERE

            if whichtest == "disfit":
                out_histo_file = f"{outdiffdir}/extended/isos/{co}_m+x_{strcontrast}_fitdist_plot.pdf"
                ratiosdf2 = steps_fitting_method(ratiosdf2, out_histo_file)
                ratiosdf2 = compute_padj_version2(ratiosdf2, 0.05, "fdr_bh")
            else:
                extract_test_df = outStat_df(ratiosdf2, metad4c, contrast, whichtest)
                extract_test_df = compute_padj_version2(extract_test_df, 0.05, "fdr_bh")
                extract_test_df.set_index("metabolite", inplace=True)
                ratiosdf2 = pd.merge(ratiosdf2, extract_test_df, left_index=True, right_index=True)


            ratiosdf2["log2FC"] = np.log2(ratiosdf2['ratio'])
            df_bad = complete_columns_for_bad(df_bad, ratiosdf2)  # HERE
            if df_bad.shape[0] >= 1:
                ratiosdf2 = pd.concat([ratiosdf2, df_bad])  # HERE

            ratiosdf2["compartment"] = co
            col_isotype = [i.split("_m+")[1] for i in ratiosdf2.index]
            col_isotype = ["m+"+str(i) for i in col_isotype]
            ratiosdf2 = ratiosdf2.assign(isotype=col_isotype)


            ratiosdf2.to_csv(f"{outdiffdir}/extended/isos/{co}_isos_{strcontrast}_{whichtest}.tsv",
                            index_label="metabolite", header=True, sep='\t')
            # HERE :
            filtered_df = filter_diff_results(ratiosdf2, args.padj, args.log2FC, args.DistanceOverSpan)
            filfi = f"{outdiffdir}/filtered/isos/{co}_isos_{strcontrast}_{whichtest}_filter.tsv"
            filtered_df.to_csv(filfi, index_label="metabolite", header=True, sep='\t')

        # end for
    # end for
    print('\n')
# end if mode diffabund


if args.mode =="metabologram":
    print("\nMetabologram\n")
    metabologram_dir = confidic['metabologram_dir']
    metabologram_config = confidic['metabologram_config']
    dimensions_pdf = (15, 20) # TODO transform into option from metabologram_config
    metabologram_run(metabologram_dir, metabologram_config, dimensions_pdf)











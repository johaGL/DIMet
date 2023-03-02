import os
import sys
sys.path.append(os.path.dirname(__file__))
import argparse
import functions_general as fg
import shutil
import yaml
import pandas as pd
import numpy as np
import re
import warnings


def prep_args():
    parser = argparse.ArgumentParser(prog="python -m DIMet.src.prepare",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('wdir', type=str,
                        help="working directory, absolute path")
    parser.add_argument('config', type=str,
                        help="configuration file, absolute path")

    # for abundance
    parser.add_argument("--under_detection_limit_set_nan", action=argparse.BooleanOptionalAction, default=True,
                        help="On VIB results. Any abundance < LOD (Limit Of Detection) is set as NaN") # np.nan exactly

    parser.add_argument("--auto_drop_metabolite_LOD_based", action=argparse.BooleanOptionalAction, default=True,
                        help="On VIB results.  By compartment, a metabolite is automatically rejected \
                        if all abundances are under the detection limit. Has effect in all tables.")

    parser.add_argument("--subtract_blankavg", action=argparse.BooleanOptionalAction, default=True,
                        help="On VIB results. From samples' abundances, subtract the average of the blanks")

    parser.add_argument("--alternative_div_amount_material", action=argparse.BooleanOptionalAction, default=False,
                        help="On VIB results, when dividing abundances by amount of material, calculate\
                        (abundance_i/amountMaterial_i) * mean(amountMaterial) to stay in abundance units")

    parser.add_argument("--use_internal_standard", default=None, type=str,
                        help='Internal Standard for performing the division: abundances/internal_standard, \
                        example: --use_internal_standard Myristic_acid_d27. By default is not performed')

    # for isotopologues
    parser.add_argument("--isos_preview", action=argparse.BooleanOptionalAction, default=False,
                        help="Plot for isotopologue values overview")

    parser.add_argument("--isosprop_min_admitted", default=-0.05, type=float,
                        help="Metabolites whose isotopologues are less or equal this cutoff are removed" )

    parser.add_argument("--stomp_values", action=argparse.BooleanOptionalAction, default=True,
                        help="Stomps isotopologues' proportions to max 1.0 and min 0.0")

    # for all
    parser.add_argument("--remove_these_metabolites", default=None, type=str,
                        help="absolute path to the .csv file with columns: compartment, metabolite. This file contains \
                        metabolites to be completely excluded from all analysis (you know what you are doing)")
                        # all tables affected

    return parser


def vibexcel2frames_dic(workingdir: str, confidic: dict) -> dict:
    frames_dic = dict()
    vib_excel = workingdir + "data/" + confidic['input_file']
    xl = pd.ExcelFile(vib_excel)
    sheetsnames = xl.sheet_names
    list_config_tabs = [confidic['name_abundance'],
                        confidic['name_meanE_or_fracContrib'],
                        confidic['name_isotopologue_prop'],
                        confidic['name_isotopologue_abs']]
    list_config_tabs = [i for i in list_config_tabs if i is not None]

    def check_config_and_sheets_match(sheetsnames, list_config_tabs):
        name_notfound = set(list_config_tabs) - set(sheetsnames)
        message =  f"One or more name_ arguments in config file not matching \
        \nthe excel sheets names (vib results):  {name_notfound}. Check spelling!"
        if len(list(name_notfound)) > 0:
            print(message)
        assert len(list(name_notfound)) == 0, message

    check_config_and_sheets_match(sheetsnames, list_config_tabs)

    for i in list_config_tabs:
        tmp = pd.read_excel(vib_excel, sheet_name=i, engine='openpyxl',
                            header=0,  index_col=0)

        badcols = [i for i in list(tmp.columns) if i.startswith("Unnamed")]
        tmp = tmp.loc[:, ~tmp.columns.isin(badcols)]
        tmp.columns = tmp.columns.str.replace(" ", "_")
        tmp.index = tmp.index.str.replace(" ", "_")
        tmp = tmp.replace(" ", "_", regex=False)
        tmp = tmp.dropna(axis=0, how="all")
        frames_dic[i] = tmp

    return frames_dic


def pull_LOD_blanks_IS(abund_df) -> tuple[pd.Series, pd.DataFrame, pd.DataFrame, dict]:
    internal_st_precol = tuple()
    pathways_by_vib = list()
    for i in range(len(abund_df.columns)): # 00_Internal_Standard, ...., 01_Glycolysis ....
        if "internal_standard" in str(abund_df.columns[i].lower()):
            internal_st_precol = (i, abund_df.columns[i])
        elif re.search(".._", abund_df.columns[i]) and\
                fg.fullynumeric(abund_df.columns[i].split("_")[0]) : # often '01_Glycolysis' and so on
            pathways_by_vib.append((i, abund_df.columns[i]))

    icolIS = range(internal_st_precol[0]+1, pathways_by_vib[0][0])
    colIS = [abund_df.columns[i] for i in icolIS]
    internal_standards_df = abund_df[colIS]

    blanks_rows = [i for i in abund_df.index if (i.lower().startswith("blank") or
                                                 i.lower().startswith("mock"))] #synonyms: blank, mock
    blanks_df = abund_df.loc[blanks_rows, :]

    lod_values = abund_df.loc['LOD', :]
    # refine dfs
    elems_x_todrop = [internal_st_precol[1]] + [i[1] for i in pathways_by_vib] +\
                            list(internal_standards_df.columns)
    elems_y_todrop = ['LOD'] + blanks_rows
    lod_values = lod_values.loc[~lod_values.index.isin(elems_x_todrop)]
    internal_standards_df = internal_standards_df.loc[~internal_standards_df.index.isin(elems_y_todrop)]
    blanks_df = blanks_df.loc[:,~blanks_df.columns.isin(elems_y_todrop)]

    todrop_x_y = {'x': elems_x_todrop,
                     'y': elems_y_todrop} # this x and y just as vib originals (not transposed)

    return lod_values, blanks_df, internal_standards_df, todrop_x_y


def reshape_frames_dic_elems(frames_dic: dict, metadata: pd.DataFrame, todrop_x_y: dict):
    trans_dic = dict()
    for k in frames_dic.keys():
        df = frames_dic[k]
        df = df.loc[:, ~df.columns.isin(todrop_x_y["x"])]
        df = df.loc[~df.index.isin(todrop_x_y["y"]), :]
        compartments = metadata['short_comp'].unique().tolist()
        trans_dic[k] = dict()
        for co in compartments:
            metada_co = metadata.loc[metadata['short_comp'] == co, :]
            df_co = df.loc[metada_co['former_name'], :]
            trans_dic[k][co] = df_co

    frames_dic = trans_dic.copy()
    return frames_dic


def abund_under_lod_set_nan(confidic, frames_dic,
                       lod_values, under_detection_limit_set_nan) -> dict:
    if under_detection_limit_set_nan:
        for co in confidic['compartments']:
            abund_co = frames_dic[confidic['name_abundance']][co]
            abund_coT = abund_co.T
            for i, r in abund_coT.iterrows():
                abund_coT.loc[i, :].loc[abund_coT.loc[i, :] < lod_values[i]] = np.nan
            frames_dic[confidic['name_abundance']][co] = abund_coT.T

    return frames_dic


def drop__metabolites_not_transposed(frames_dic_orig, bad_metabolites_dic):
    frames_dic = frames_dic_orig.copy()
    for tab_name in frames_dic.keys():
        for co in frames_dic[tab_name].keys():
            tmpdf = frames_dic[tab_name][co]
            if not "isotopolog" in tab_name.lower():
                to_drop_now =  bad_metabolites_dic[co]
                tmpdf = tmpdf.drop(columns=to_drop_now)
                frames_dic[tab_name][co] = tmpdf
            else:
                to_drop_now_isos = list()
                for i in tmpdf.columns:
                    for j in bad_metabolites_dic[co]:
                        if i.startswith(j):
                            to_drop_now_isos.append(i)
                tmpdf = tmpdf.drop(columns=to_drop_now_isos)
                frames_dic[tab_name][co] = tmpdf
    return frames_dic


def auto_drop_metabolites_uLOD(confidic, frames_dic, lod_values,
                               auto_drop_metabolite_LOD_based: bool) -> dict:
    # affects all the datasets in frames_dic
    auto_bad_metabolites = dict()
    compartments = confidic['compartments']
    for k in compartments:
        auto_bad_metabolites[k] = list()

    if auto_drop_metabolite_LOD_based:
        for co in compartments:
            abund_co = frames_dic[confidic['name_abundance']][co]
            abund_coT = abund_co.T
            for i, r in abund_coT.iterrows():
                nb_nan = abund_coT.loc[i, :].isna().sum()
                nb_under_LOD = (abund_coT.loc[i, :] <= lod_values[i]).sum()
                if (nb_under_LOD >= r.size - 1) or (nb_nan >= r.size - 1):
                    auto_bad_metabolites[co].append(i)
        frames_dic = drop__metabolites_not_transposed(frames_dic, auto_bad_metabolites)

    return frames_dic


def abund_subtract_blankavg(frames_dic: dict, confidic: dict,
                             blanks_df: pd.Series, subtract_blankavg: bool):
    abund_dic = frames_dic[confidic['name_abundance']].copy()
    if subtract_blankavg:
        for compartment in abund_dic.keys():
            blanks_df_s = blanks_df[list(abund_dic[compartment].columns)]
            blanksAvg_s = blanks_df_s.mean(axis=0)
            abu_df_T = abund_dic[compartment].T
            tmp = abu_df_T.subtract(blanksAvg_s, axis = 'index')
            tmp[tmp < 0] = 0
            abund_dic[compartment] = tmp.T

        frames_dic[confidic['name_abundance']] = abund_dic

    return frames_dic


def abund_divideby_amount_material(frames_dic: dict, workingdir: str,
                                   amount_material: str, alternative_method: bool):
    if amount_material is not None:
        try:
            file = workingdir + "data/" + amount_material
            material_df = pd.read_csv(file, index_col=0)
        except FileNotFoundError as err_file:
            print(err_file)
        except Exception as e:
            print(e)

        assert material_df.shape[1] == 1, "amountMaterial table must have only 2 columns"
        assert (material_df.iloc[:, 0] <= 0).sum() == 0, "\
                    amountMaterial table must not contain zeros nor negative numbers"
        abund_dic = frames_dic[confidic['name_abundance']].copy()
        for compartment in abund_dic.keys():
            material_df_s = material_df.loc[list(abund_dic[compartment].index), :]
            if alternative_method:
                material_avg = material_df_s.iloc[:, 0].mean()
                material_avg_ser = pd.Series([float(material_avg) for i in range(material_df_s.shape[0])],
                                             index = material_df_s.index)
                tmp = abund_dic[compartment].div(material_df_s.iloc[:, 0], axis=0)
                tmp = tmp.mul(material_avg_ser, axis=0)
            else:
                tmp = abund_dic[compartment].div(material_df_s.iloc[:, 0], axis=0)

            frames_dic[confidic['name_abundance']][compartment] = tmp

    return frames_dic


def abund_divideby_internalStandard(frames_dic, confidic,
                                    internal_standards_df, use_internal_standard: [str, None]):
    # only useful for in vib preparation
    if use_internal_standard is None:
        return frames_dic
    else:
        picked_internal_standard = use_internal_standard
        assert picked_internal_standard in internal_standards_df.columns, "\
               Error, that internal standard is not present in the data"
        abund_dic = frames_dic[confidic['name_abundance']].copy()
        for compartment in abund_dic.keys():
            inte_sta_co = internal_standards_df.loc[abund_dic[compartment].index, :]
            inte_sta_co_is = inte_sta_co[picked_internal_standard]
            # replace zeros to avoid zero division, uses min of the pd series :
            inte_sta_co_is[inte_sta_co_is == 0] = inte_sta_co_is[inte_sta_co_is > 0].min()
            tmp = abund_dic[compartment].div(inte_sta_co_is, axis = 0)
            frames_dic[confidic['name_abundance']][compartment] = tmp
        return frames_dic


def transformmyisotopologues(isos_list, style):
    if "vib" in style.lower():
        outli = list()
        for ch in isos_list:
            if "_C13-label-" in ch:
                elems = ch.split("_C13-label-")
                metabolite = elems[0]
                species = "m+{}".format(elems[-1].split("-")[-1])
            elif "_PARENT" in ch:
                elems = ch.split("_PARENT")
                metabolite = elems[0]
                species = "m+0"
            else:
                metabolite = ch
                species = "m+?"
            outli.append(metabolite + "_" + species)
    elif "generic" in style.lower():
        try:
            outli = [i.replace("label", "m+") for i in isos_list]
        except Exception as e:
            print(e)
            print("not possible to change the isotopologues name style")
            outli = isos_list
    else:
        outli = isos_list
        raise ValueError("isotopologues style not vib nor generic")
    return outli


def save_isos_preview(dic_isos_prop, metadata, workingdir, output_plots_dir, the_boolean_arg):
    if the_boolean_arg:
        for k in metadata['short_comp'].unique():
            df = dic_isos_prop[k]
            sples_co = metadata.loc[metadata["short_comp"] == k, "former_name"]
            df = df.loc[sples_co,:]
            df = df.T # transpose
            df = fg.add_metabolite_column(df)
            df = fg.add_isotopologue_type_column(df)

            thesums = df.groupby('metabolite').sum()
            thesums = thesums.drop(columns=['isotopologue_type'])
            thesums = thesums.round(3)
            ff = f"{workingdir}{output_plots_dir}sums_Iso_{k}.pdf"
            figuretitle = f"Sums of isotopologue proportions ({k}) "
            fg.save_heatmap_sums_isos(thesums, figuretitle, ff)

            dfmelt = pd.melt(df, id_vars=['metabolite', 'isotopologue_type'])
            dfmelt = fg.givelevels(dfmelt)
            # fg.table_minimalbymet(dfmelt, f"isotop_preview/minextremesIso_{k}.tsv")
            outputfigure = f"{workingdir}{output_plots_dir}allsampleIsos_{k}.pdf"
            figtitle = f"{k} compartment, Isotopologues (proportions) across all samples"
            fg.save_rawisos_plot(dfmelt, figuretitle=figtitle, outputfigure=outputfigure)


def set_samples_names(frames_dic, metadata):
    # do smartest as possible: detect if dfs' rownames match to former_name or sample or none,
    # if match sample do nothing, if match former_name set sample as rownames, finally if none, error stop
    compartments = metadata['short_comp'].unique().tolist()

    for tab in frames_dic.keys():
        for co in compartments:
            metada_co = metadata.loc[metadata['short_comp'] == co, :]
            df = frames_dic[tab][co]
            df.reset_index(inplace=True)
            df.rename(columns={ df.columns[0]: "former_name" }, inplace = True)
            careful_samples_order = pd.merge(df.iloc[:, 0], metada_co[['sample', 'former_name']],
                                             how="left", on="former_name")
            df = df.assign(sample=careful_samples_order['sample'])
            df = df.set_index('sample')
            df = df.drop(columns=['former_name'])
            frames_dic[tab][co] = df

    return frames_dic


def drop_metabolites_infile(frames_dic, workingdir, exclude_list: [str, None]):
    if exclude_list is not None:
        # a file name to open
        try:
            file = workingdir + "data/" + exclude_list
            exclude_df = pd.read_csv(file,  header=0)
        except FileNotFoundError as err_file:
            print(err_file)
        except Exception as e:
            print(e)

        unwanted_metabolites = dict()
        for co in exclude_df["short_comp"].unique():
            mets_l = exclude_df.loc[exclude_df["short_comp"] == co, 'metabolite'].tolist()
            unwanted_metabolites[co] = mets_l

        frames_dic = drop__metabolites_not_transposed(frames_dic, unwanted_metabolites)

    return frames_dic


def propagate_nan(confidic, frames_dic): # TODO
    # propagates nan from abundance to isotopologues and fractional contributions
    for co in confidic['compartments']:
        abu_co = frames_dic[confidic['name_abundance']][co]
        frac_co = frames_dic[confidic['name_meanE_or_fracContrib']][co]
        frac_co.mask(abu_co.isnull())
        # pending the propagation to isotopologues, think !

    return frames_diccolumns


def prep_isotopologues_df(isos_df, isosprop_min_admitted: float, stomp_values: bool): # TODO
    # threshold for negatives
    # stomp values
    return 0




def do_vib_prep(workingdir, args, confidic, output_plots_dir):
    # the order of the steps is the one recommended by VIB
    frames_dic = vibexcel2frames_dic(workingdir,  confidic)
    metadata = fg.open_metadata(workingdir, confidic)
    fg.verify_metadata_sample_not_duplicated(metadata)
    abundance_df = frames_dic[confidic['name_abundance']]
    lod_values, blanks_df, internal_standards_df, bad_x_y = pull_LOD_blanks_IS(abundance_df)
    frames_dic = reshape_frames_dic_elems(frames_dic, metadata, bad_x_y)
    frames_dic = abund_under_lod_set_nan(confidic, frames_dic, lod_values,
                                         args.under_detection_limit_set_nan)
    frames_dic = auto_drop_metabolites_uLOD(confidic, frames_dic,  lod_values,
                                                  args.auto_drop_metabolite_LOD_based)
    frames_dic = abund_subtract_blankavg(frames_dic, confidic,
                                              blanks_df, args.subtract_blankavg)
    frames_dic = abund_divideby_amount_material(frames_dic, workingdir, confidic['amountMaterial'],
                                                args.alternative_div_amount_material)
    frames_dic = abund_divideby_internalStandard(frames_dic, confidic, internal_standards_df,
                                                 args.use_internal_standard)

    # transform isotopologues names to the easier "m+x" style:
    for tab in frames_dic.keys():
        if "isotopol" in tab.lower():
            for co in frames_dic[tab]:
                tmp = frames_dic[tab][co]
                new_col = transformmyisotopologues(tmp.columns, "vib")
                tmp.columns = new_col
                frames_dic[confidic['name_isotopologue_prop']][co] = tmp
    # end for

    save_isos_preview(frames_dic[confidic['name_isotopologue_prop']], metadata,
                      workingdir, output_plots_dir, args.isos_preview) # TODO: reactivate
    frames_dic = set_samples_names(frames_dic, metadata)
    frames_dic = drop_metabolites_infile(frames_dic, workingdir,  args.remove_these_metabolites)

    # propagate_nan
    # prep_isotopologues_df(isos_df, isosprop_min_admitted, stomp_values) **
    return frames_dic


def generictype_div_intern_stand(frames_dic, confidic, picked_internal_standard, ISdf):
    if picked_internal_standard is not None:
        frames_dic = abund_divideby_internalStandard(frames_dic, confidic, ISdf,
                                                     picked_internal_standard)
    return frames_dic


def compute_abund_from_absolute_isotopol(df, metabos_isos_df):
    df = df.T
    metabos_uniq = metabos_isos_df['metabolite'].unique()
    abundance = pd.DataFrame(index= metabos_uniq, columns=df.columns)
    for m in metabos_uniq:
        isos_here = metabos_isos_df.loc[metabos_isos_df['metabolite'] == m, 'isotopologue_name']
        sub_df = df.loc[isos_here, :]
        sub_df_sum = sub_df.sum(axis=0)
        abundance.loc[m, :] = sub_df_sum
    return abundance.T


def do_isocorOutput_prep(workingdir, args, confidic, output_plots_dir):

    def df2dic_comp(df, metadata): # TODO change name
        # splits df into dictionary of dataframes, each for one compartment
        out_dic = dict()
        for co in metadata['short_comp'].unique():
            metada_co = metadata.loc[metadata['short_comp'] == co, :]
            df_co = df.loc[metada_co['former_name'], :]
            out_dic[co] = df
        return out_dic

    def isocorOutput2frames_dic(isocor_out_df, metadata, confidic, internal_standard: [str, None]):
        df = isocor_out_df[['sample', 'metabolite', 'isotopologue', 'corrected_area',
                            'isotopologue_fraction', 'mean_enrichment']]
        isonames = df.metabolite.str.cat(df.isotopologue.astype(str), sep="_m+")
        df = df.assign(isotopologue_name=isonames)
        samples = isocor_out_df['sample'].unique()

        metabos_isos_df = df[['metabolite', 'isotopologue_name']]
        metabos_isos_df = metabos_isos_df.drop_duplicates()

        me_or_fc_melted = df[['sample', 'metabolite', 'mean_enrichment']]
        me_or_fc_melted = me_or_fc_melted.drop_duplicates()
        me_or_fc = me_or_fc_melted.pivot(index='sample', columns='metabolite')
        me_or_fc = me_or_fc['mean_enrichment'].reset_index()
        me_or_fc = me_or_fc.set_index('sample')

        isos_prop_melted = df[['sample', 'isotopologue_name', 'isotopologue_fraction']]
        isos_prop_melted = isos_prop_melted.drop_duplicates()
        isos_prop = isos_prop_melted.pivot(index='sample', columns='isotopologue_name')
        isos_prop = isos_prop['isotopologue_fraction'].reset_index()
        isos_prop = isos_prop.set_index('sample')

        isos_absolute_melted = df[['sample', 'isotopologue_name', 'corrected_area']]
        isos_absolute_melted = isos_absolute_melted.drop_duplicates()
        isos_absolute = isos_absolute_melted.pivot(index='sample', columns='isotopologue_name')
        isos_absolute = isos_absolute['corrected_area'].reset_index()
        isos_absolute = isos_absolute.set_index('sample')

        abundance = compute_abund_from_absolute_isotopol(isos_absolute, metabos_isos_df)
        if internal_standard is not None:
            intesta_abun_df = abundance[[internal_standard]]
        else:
            intesta_abun_df = None

        frames_dic = dict()

        frames_dic[confidic['name_meanE_or_fracContrib']] = df2dic_comp(me_or_fc, metadata)
        frames_dic[confidic['name_isotopologue_prop']] = df2dic_comp(isos_prop, metadata)
        frames_dic[confidic['name_isotopologue_abs']] = df2dic_comp(isos_absolute, metadata)
        frames_dic[confidic['name_abundance']] = df2dic_comp(abundance, metadata)

        return frames_dic, intesta_abun_df

    isocor_output_df = pd.read_excel(workingdir + "data/" + confidic['input_file'])
    metadata = fg.open_metadata(workingdir, confidic)
    fg.verify_metadata_sample_not_duplicated(metadata)
    frames_dic, intesta_abun_df = isocorOutput2frames_dic(isocor_output_df, metadata, confidic,
                                         args.use_internal_standard)

    save_isos_preview(frames_dic[confidic['name_isotopologue_prop']], metadata,
                      workingdir, output_plots_dir, args.isos_preview)

    frames_dic = abund_divideby_amount_material(frames_dic, workingdir, confidic['amountMaterial'],
                                                         args.alternative_div_amount_material)
    frames_dic = generictype_div_intern_stand(frames_dic, confidic, args.use_internal_standard, intesta_abun_df)

    frames_dic = set_samples_names(frames_dic, metadata)
    print(args.remove_these_metabolites)
    frames_dic = drop_metabolites_infile(frames_dic, workingdir,  args.remove_these_metabolites)
    # propagate_nan
    # prep_isotopologues_df(isos_df, isosprop_min_admitted, stomp_values)
    aa = 1
    return frames_dic

def do_generic_prep(): #  because user has a 'generic' set of csv tables, and does not have blanks nor LOD infos
    # inte_sta_direct_df = take it from abundance !
    # frames_dicto = funcitonharvesting the tables into a dictionary : first checks if metadata.sample and
                                       # rows in the tables match perfectly, if not, print and stop
    # frames_dic_trans = abund_divideby_amount_material(frames_dic_trans, workingdir, confidic['amountMaterial'],
    #                                                     args.alternative_div_amount_material)
    # frames_dic = generictype_div_intern_stand(frames_dic, confidic, args.use_internal_standard)
    # metadata
    # fg.verify_metadata_sample_not_duplicated(metadata)
    # frames_dicto = drop_metabolites_infile(df, exclude_list)
    # prep_isotopologues_df(isos_df, isosprop_min_admitted, stomp_values)
    return 0


def compute_isotopologues_proportions_from_absolute(frames_dic):
    return frames_dic


def compute_MEorFC_from_isotopologues_proportions(frames_dic):
    return frames_dic


def complete_missing_dfs(confidic, frames_dic, metabolites_isos_df):
    list_config_tabs = [confidic['name_abundance'],
     confidic['name_meanE_or_fracContrib'],
     confidic['name_isotopologue_prop'],
     confidic['name_isotopologue_abs']]

    existing = list(frames_dic.keys())
    key_name_suffixes = {'name_abundance': confidic['name_abundance'],
                         'name_meanE_or_fracContrib' : confidic['name_meanE_or_fracContrib'],
                         'name_isotopologue_prop' : confidic['name_isotopologue_prop'],
                         'name_isotopologue_abs': confidic['name_isotopologue_abs']}

    list_none_tabs = [i for i in list_config_tabs if i not in existing ]
    #ok_tabs = [i for i in dic_config_tabs.keys() if dic_config_tabs[i] is not None]
    ok_tabs = [i for i in list_config_tabs if i in existing]
    compartments = confidic['compartments']
    if 'name_abundance' in list_none_tabs:
        if 'name_isotopologue_abs' in ok_tabs:
            for co in compartments:
                df_co = frames_dic[confidic['name_isotopologue_abs']][co]
                tmp_co = compute_abund_from_absolute_isotopol(df_co, metabolites_isos_df)
                frames_dic["abundance_computed"][co] = tmp_co
                key_name_suffixes['name_abundance'] = "abundance_computed"
        elif 'name_isotopologue_abs' in list_none_tabs:
            print(" isotopologues' absolute values not available, impossible to get proportions")
    if 'name_isotopologue_prop' in list_none_tabs :
        if 'name_isotopologue_abs' in ok_tabs:
            for co in compartments:
                df_co = frames_dic[confidic['name_isotopologue_abs']][co]
                tmp_co = compute_isotopologues_proportions_from_absolute(df_co)
                frames_dic["isotopologues_props_computed"][co] = tmp_co
                key_name_suffixes['name_isotopologue_prop'] = "isotopologues_props_computed"
        elif 'name_isotopologue_abs' in list_none_tabs:
            print(" isotopologues' absolute values not available, impossible to get proportions")
    if 'name_meanE_or_fracContrib' in list_none_tabs:
        try:
            for co in compartments:
                df_co = frames_dic[confidic['name_isotopologue_prop']][co]
                tmp_co = compute_MEorFC_from_isotopologues_proportions(df_co)
                frames_dic["meanEnr_or_FracC_computed"][co] = tmp_co
                key_name_suffixes['name_meanE_or_fracContrib'] = "meanEnr_or_FracC_computed"
        except Exception as e:
            print("impossible to calculate: mean enrichment  or  fractional contribution\
                  isotopologues proportions not found")
            print(e)
    return frames_dic, key_name_suffixes

def perform_type_prep(workingdir, args, confidic):
    output_plots_dir = args.config.split("/")[-2] + "/results/plots/preview/"
    fg.detect_and_create_dir(workingdir + output_plots_dir)
    output_tabs_dir = args.config.split("/")[-2] + "/results/prepared_tables/"
    fg.detect_and_create_dir(workingdir + output_tabs_dir)

    if confidic['typeprep'] == 'IsoCor_output_prep':
        frames_dic  = do_isocorOutput_prep(workingdir, args, confidic,
                                           output_plots_dir)
    elif confidic['typeprep'] == 'vib_prep': # ?? 'VIBMEC_output_prep' !! yes # TODO change in config file
        frames_dic = do_vib_prep(workingdir, args, confidic,
                                 output_plots_dir)
    elif confidic['typeprep'] == 'generic_prep':
        frames_dic = do_generic_prep(workingdir, args, confidic,
                                     output_plots_dir)

    frames_dic, the_finale_blah = complete_missing_dfs(confidic, frames_dic, metabolites_isos_df)
    # for saving, transpose as finally desired: met x sam
    the_finale_blah = dict()
    txt_real_available_tables : ""
    for k in the_finale_blah:
        pass
        #txt_real_available_tables += f'{k},{}\n'
    with open(workingdir + output_tabs_dir + "nohope.csv") as f:
        f.write(txt_real_available_tables)
    # TODO : create a prefixes.txt file that contains the all tables, existing and completedmissing:  name_annnana, thePrefix\n .....


    # save all tables :
    return 0


if __name__ == "__main__":

    parser = prep_args()
    args = parser.parse_args()

    workingdir = os.path.expanduser(args.wdir)
    if not workingdir.endswith("/"):
        workingdir += "/"

    configfile = os.path.expanduser(args.config)
    fg.wdir_configpaths_validate(workingdir, configfile)
    confidic = fg.open_config_file(configfile)

    perform_type_prep(workingdir, args, confidic)
    print("\n,**Hey : ??at some point: replace toy1 by toy_VIBMEC,, toy2 by toy_IsoCor, ?? create toy_generic")







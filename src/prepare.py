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



def prep_args():
    parser = argparse.ArgumentParser(prog="python -m DIMet.src.prepare")

    parser.add_argument('wdir', type=str,
                        help="working directory, absolute path")
    parser.add_argument('config', type=str,
                        help="configuration file, also absolute path")

    # for abundance
    parser.add_argument("--under_detection_limit_set_nan", action=argparse.BooleanOptionalAction, default=True,
                        help="On VIB results. Any abundance < LOD (Limit Of Detection) is set as NaN") # np.nan exactly

    parser.add_argument("--auto_drop_metabolite_LOD_based", action=argparse.BooleanOptionalAction, default=True,
                        help="On VIB results.  By compartment, a metabolite is automatically rejected \
                        if all abundances are under the detection limit. Has effect in all tables.")

    parser.add_argument("--subtract_blankavg", action=argparse.BooleanOptionalAction, default=True,
                        help="On VIB results. From samples' abundances, subtract the average of the blanks")

    parser.add_argument("--amount_material_div_alternative", action=argparse.BooleanOptionalAction, default=False,
                        help="On VIB results, when dividing abundances by amount of material, \
                             stay in the abundance scale: (abundance_i/amountMaterial_i) * mean(amountMaterial)")

    parser.add_argument("--use_internal_standard", default=None, type=str,
                        help='Internal Standard for performing the division: abundances/internal_standard, \
                        example : "Myristic_acid_d27". By default is not performed')

    # for isotopologues
    parser.add_argument("--isosprop_min_admitted", default=-0.05, type=float,
                        help="Negative proportions can occur, for example, after El-Maven isotopologues correction.  \
                        Setting a cutoff, metabolites whose isotopologues are less or equal this cutoff are removed")

    parser.add_argument("--stomp_values", action=argparse.BooleanOptionalAction, default=True, \
                        help="Stomps isotopologues' proportions to max 1.0 and min 0.0")

    # for all
    parser.add_argument("--exclude_list", default=None, type=str,
                        help="path to .csv file with columns: compartment, metabolite. Metabolites to be excluded")
                        # all tables affected

    return parser


def excel2frames_dic(workingdir: str, confidic: dict) -> dict:
    frames_dic = dict()
    vib_excel = workingdir + "data/" + confidic['input_file']
    xl = pd.ExcelFile(vib_excel)
    sheetsnames = xl.sheet_names
    list_config_names = [confidic['name_abundance'],
                        confidic['name_meanE_or_fracContrib'],
                        confidic['name_isotopologue_prop'],
                        confidic['name_isotopologue_abs']]
    list_config_names = [i for i in list_config_names if i is not None]

    def check_config_and_sheets_match(sheetsnames, list_config_names):
        name_notfound = set(list_config_names) - set(sheetsnames)
        message =  f"One or more name_ arguments in config file not matching \
        \nthe excel sheets names (vib results):  {name_notfound}. Check spelling!"
        if len(list(name_notfound)) > 0:
            print(message)
        assert len(list(name_notfound)) == 0, message

    check_config_and_sheets_match(sheetsnames, list_config_names)

    for i in list_config_names:
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


def auto_drop_metabolites(confidic, frames_dic, lod_values,
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
                if (nb_under_LOD >= r.size - 1) or (nb_nan <= r.size - 1):
                    auto_bad_metabolites[co].append(i)

        for tab_name in frames_dic.keys():
            for co in compartments:
                tmpdf = frames_dic[tab_name][co]
                if not "isotopolog" in tab_name.lower():
                    to_drop_now = auto_bad_metabolites[co]
                    tmpdf = tmpdf.drop(columns=to_drop_now)
                    frames_dic[tab_name][co] = tmpdf
                else:
                    to_drop_now_isos = list()
                    for i in tmpdf.columns:
                        for j in auto_bad_metabolites[co]:
                            if i.startswith(j):
                                to_drop_now_isos.append(i)
                    tmpdf = tmpdf.drop(columns=to_drop_now_isos)
                frames_dic[tab_name][co] = tmpdf
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

        assert material_df.shape[1] == 1, "amountMaterial table must have 1 single column"
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


def set_samples_names(frames_dic, todrop_x_y):
    # do smartest as possible: detect if dfs' rownames match to former_name or sample or none,
    # if match sample do nothing, if match former_name set sample as rownames, finally if none, error stop
    trans_dic = dict()
    for k in frames_dic.keys():
        df = frames_dic[k]
        df = df.loc[:, ~df.columns.isin(todrop_x_y["x"])]
        df = df.loc[~df.index.isin(todrop_x_y["y"]), :]
        # df.reset_index(inplace=True)
        # df.rename(columns={ df.columns[0]: "former_name" }, inplace = True)
        # careful_samples_order = pd.merge(df.iloc[:, 0], metadata[['sample', 'former_name']],
        #                                  how="left", on="former_name")
        # df = df.assign(sample=careful_samples_order['sample'])
        # df = df.set_index('sample')
        # df = df.drop(columns=['former_name'])

        compartments = metadata['short_comp'].unique().tolist()
        trans_dic[k] = dict()
        for co in compartments:
            metada_co = metadata.loc[metadata['short_comp'] == co, :]
            # df_co = df.loc[metada_co['sample'], :]
            df_co = df.loc[metada_co['former_name'], :]
            trans_dic[k][co] = df_co

    frames_dic = trans_dic.copy()
    return


def drop_metabolites_infile(frames_dic, workingdir, exclude_list: [str, None]):
    if exclude_list is not None:
        # a file name to open
        try:
            file = workingdir + "data/" + exclude_list
            exclude_df = pd.read_csv(file, index_col=0)
        except FileNotFoundError as err_file:
            print(err_file)
        except Exception as e:
            print(e)

    return 0


def propagate_nan(confidic, frames_dic):
    # propagates nan from abundance to isotopologues and fractional contributions
    for co in confidic['compartments']:
        abu_co = frames_dic[confidic['name_abundance']][co]
        frac_co = frames_dic[confidic['name_meanE_or_fracContrib']][co]
        frac_co.mask(abu_co.isnull())
        # pending the propagation to isotopologues, think !

    return frames_dic


def prep_isotopologues_df(isos_df, isosprop_min_admitted: float, stomp_values: bool):
    # threshold for negatives
    # stomp values
    # recognize style and make names as said in somewhere function
    return 0



def do_vib_prep(workingdir, args, confidic):
    # the order of the steps is the one recommended by VIB
    frames_dic = excel2frames_dic(workingdir,  confidic)
    metadata = fg.open_metadata(workingdir, confidic)
    abundance_df = frames_dic[confidic['name_abundance']]
    lod_values, blanks_df, internal_standards_df, bad_x_y = pull_LOD_blanks_IS(abundance_df)
    frames_dic = reshape_frames_dic_elems(frames_dic, metadata, bad_x_y)
    frames_dic = abund_under_lod_set_nan(confidic, frames_dic, lod_values,
                                         args.under_detection_limit_set_nan)
    frames_dic = auto_drop_metabolites(confidic, frames_dic,  lod_values,
                                       args.auto_drop_metabolite_LOD_based)
    frames_dic = abund_subtract_blankavg(frames_dic, confidic,
                                              blanks_df, args.subtract_blankavg)
    frames_dic = abund_divideby_amount_material(frames_dic, workingdir, confidic['amountMaterial'],
                                                args.amount_material_div_alternative)
    frames_dic = abund_divideby_internalStandard(frames_dic, confidic, internal_standards_df,
                                                 args.use_internal_standard)

    # # set_sample_names(df, metadata)
    # drop_metabolites_infile(df, exclude_list)
    # propagate_nan
    # prep_isotopologues_df(isos_df, isosprop_min_admitted, stomp_values) **
    # savetables_bycompartment()
    return 0


def generictype_div_intern_stand(frames_dic, confidic, picked_internal_standard):
    if picked_internal_standard is not None:
        # make a here_IS_df (single column)
        hereISdf = 0
        # take away that column from abund and everything else in frames_dic
        # ...
        # and finally do :
        frames_dic =  abund_divideby_internalStandard(frames_dic, confidic, hereISdf,
                                                       picked_internal_standard)
    return frames_dic

def do_isocorOutput_prep(workingdir, args, confidic):
    # frames_dic =
    #frames_dic = abund_divideby_amount_material(frames_dic, workingdir, confidic['amountMaterial'],
    #                                                     args.amount_material_div_alternative)
    frames_dic = generictype_div_intern_stand(frames_dic, confidic, args.use_internal_standard)

    # # set_sample_names(df, metadata)
    # drop_metabolites_infile(df, exclude_list)
    # propagate_nan
    # prep_isotopologues_df(isos_df, isosprop_min_admitted, stomp_values) **
    # savetables_bycompartment()


    
    return 0 #

def do_generic_prep(): #  because user has a 'generic' set of csv tables, and does not have blanks nor LOD infos
    # frames_dicto = funcitonharvesting the tables into a dictionary
    # frames_dicto = functionthat makes the frames to transposed or not, looking the metadata.samples in index and cols
    # frames_dic_trans = abund_divideby_amount_material(frames_dic_trans, workingdir, confidic['amountMaterial'],
    #                                                     args.amount_material_div_alternative)
    # frames_dic = generictype_div_intern_stand(frames_dic, confidic, args.use_internal_standard)

    # prep_isotopologues_df(isos_df, isosprop_min_admitted, stomp_values)
    # frames_dicto = drop_metabolites_infile(df, exclude_list)
    return 0


def perform_type_prep(workingdir, args, confidic):
    if confidic['typeprep'] == 'vib_prep':
        do_vib_prep(workingdir, args, confidic)
    elif confidic['typeprep'] == 'isocorOutput_prep':
        do_isocorOutput_prep(workingdir, args, confidic)
    elif confidic['typeprep'] == 'generic_prep':
        do_generic_prep(workingdir, args, confidic)
    # for saving, transpose as finally desired: met x sam
    outputdir = args.config.split("/")[-2] + "/results/prepared_tables/"
    fg.detect_and_create_dir(workingdir + outputdir)
    # save all tables :

def save_advanced_options_used(args):
    return 0  # see better how to save a 'log'



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








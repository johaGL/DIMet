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

    parser.add_argument("--substract_blankavg", action=argparse.BooleanOptionalAction, default=True,
                        help="On VIB results. To samples' abundances, substract the average of the blanks")

    parser.add_argument("--amount_material_div_alternative", action=argparse.BooleanOptionalAction, default=False,
                        help="On VIB results, when dividing abundances by amount of material, \
                             stay in the abundance scale: (abundance_i/amountMaterial_i)* mean(amountMaterial)")

    parser.add_argument("--use_internal_standard", default=None, type=str,
                        help='The name of the Internal Standard for performing the division: abundances/internal_standard, \
                        example : "Myristic_acid_d27". By default is not performed')

    # for isotopologues
    parser.add_argument("--isosprop_min_admitted", default=-0.05, type=float, # TODO: change this default to -1 ! keep-0.05 for my cmd file
                        help="Negative proportions can occur, for example, after El-Maven isotopologues correction,  \
                        metabolites whose isotopologues are less or equal this cutoff are removed") # only isotopologues

    parser.add_argument("--stomp_values", action=argparse.BooleanOptionalAction, default=True, \
                        help="Stomps isotopologues' proportions to max 1.0 and min 0.0")

    # for all
    parser.add_argument("--exclude_list", default=None, type=str,
                        help="path to .csv file with columns: compartment, metabolite. Metabolites to be excluded")  # all tables affected

    return parser


def excel2individual_dic(workingdir: str, confidic: dict) -> dict:
    individual_dic = dict()
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
        tmp = tmp.loc[:,~tmp.columns.isin(badcols)]
        tmp.columns = tmp.columns.str.replace(" ", "_")
        tmp.index = tmp.index.str.replace(" ", "_")
        tmp = tmp.replace(" ", "_", regex=False)
        tmp = tmp.dropna(axis=0, how="all")
        individual_dic[i] = tmp

    return individual_dic


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


def reshape_individual_dic_elems(individual_dic: dict, metadata: pd.DataFrame, todrop_x_y: dict):
    trans_dic = dict()
    for k in individual_dic.keys():
        df = individual_dic[k]
        df = df.loc[:, ~df.columns.isin(todrop_x_y["x"])]
        df = df.loc[~df.index.isin(todrop_x_y["y"]), :]
        compartments = metadata['short_comp'].unique().tolist()
        trans_dic[k] = dict()
        for co in compartments:
            metada_co = metadata.loc[metadata['short_comp'] == co, :]
            df_co = df.loc[metada_co['former_name'], :]
            trans_dic[k][co] = df_co

    individual_dic = trans_dic.copy()
    return individual_dic


def abund_under_lod_set_nan(confidic, individual_dic,
                       lod_values, under_detection_limit_set_nan) -> dict:

    if under_detection_limit_set_nan:

        for co in confidic['compartments']:
            abund_co = individual_dic[confidic['name_abundance']][co]
            abund_coT = abund_co.T
            for i, r in abund_coT.iterrows():
                abund_coT.loc[i, :].loc[abund_coT.loc[i, :] < lod_values[i]] = np.nan
            individual_dic[confidic['name_abundance']][co] = abund_coT.T #tmp_coT.T

    return individual_dic


def auto_drop_metabolites(compartments, individual_dic, lod_values,
                          auto_drop_metabolite_LOD_based: bool) -> dict:
    # affects all the datasets in individual_dic
    auto_bad_metabolites = dict()
    for k in compartments:
        auto_bad_metabolites[k] = list()

    if auto_drop_metabolite_LOD_based:

        for co in compartments:
            abund_co = individual_dic[confidic['name_abundance']][co]
            abund_coT = abund_co.T
            for i, r in abund_coT.iterrows():
                # nbzeros = (abund_coT.loc[i, :] == 0).sum()
                nb_nan = abund_coT.loc[i, :].isna().sum()
                nb_under_LOD = (abund_coT.loc[i, :] <= lod_values[i]).sum()
                if (nb_under_LOD >= r.size - 1) or (nb_nan <= r.size -1):
                    auto_bad_metabolites[co].append(i)

        for tab_name in individual_dic.keys():
            for co in compartments:
                tmpdf = individual_dic[tab_name][co]
                if not "isotopolog" in tab_name.lower():
                    to_drop_now = auto_bad_metabolites[co]
                    tmpdf = tmpdf.drop(columns=to_drop_now)
                    individual_dic[tab_name][co] = tmpdf
                else:
                    to_drop_now_isos = list()
                    for i in tmpdf.columns:
                        for j in auto_bad_metabolites[co]:
                            if i.startswith(j):
                                to_drop_now_isos.append(i)
                    tmpdf = tmpdf.drop(columns=to_drop_now_isos)
                individual_dic[tab_name][co] = tmpdf
    return individual_dic


def abund_substract_blankavg(individual_dic: dict, confidic: dict,
                             blanks_df: pd.Series, substract_blankavg: bool):
    abund_dic = individual_dic[confidic['name_abundance']].copy()
    if substract_blankavg:
        for compartment in abund_dic.keys():
            blanks_df_s = blanks_df[list(abund_dic[compartment].columns)]
            blanksAvg_s = blanks_df_s.mean(axis=0)
            abu_df_T = abund_dic[compartment].T
            tmp = abu_df_T.subtract(blanksAvg_s, axis = 'index')
            tmp[tmp < 0] = 0
            abund_dic[compartment] = tmp.T

        individual_dic[confidic['name_abundance']] = abund_dic

    return individual_dic


def abund_divideby_amount_material(individual_dic: dict, workingdir: str,
                                   amount_material: str, alternative_method: bool):
    if amount_material is None:
        return individual_dic
    else:
        try:
            file = workingdir + "data/" + amount_material
            material_df = pd.read_csv(file, index_col=0)
        except FileNotFoundError as err_file:
            print(err_file)
        except Exception as e:
            print(e)

        assert material_df.shape[1] == 1, "amountMaterial table must have 1 single column"
        assert (material_df.iloc[:, 0] <= 0).sum() == 0, "amountMaterial table must not contain zeros nor negative numbers"
        abund_dic = individual_dic[confidic['name_abundance']].copy()
        for compartment in abund_dic.keys():
            material_df_s = material_df.loc[list(abund_dic[compartment].index), :]
            if alternative_method:
                material_avg = material_df_s.iloc[:, 0].mean()
                material_avg_ser = pd.Series([float(material_avg) for i in range(material_df_s.shape[0])],
                                             index = material_df_s.index)
                tmp = abund_dic[compartment].div(material_df_s.iloc[:, 0], axis=0)
                tmp = tmp.mul(material_avg_ser, axis = 0)
            else:
                tmp = abund_dic[compartment].div(material_df_s.iloc[:, 0], axis=0)

            individual_dic[confidic['name_abundance']][compartment] = tmp

        return individual_dic


def abund_divideby_internalStandard(individual_dic, confidic,
                                    internal_standards_df, use_internal_standard: [str, None]):
    if use_internal_standard is None:
        return individual_dic
    else:
        picked_internal_standard = use_internal_standard
        assert picked_internal_standard in internal_standards_df.columns, "\
               Error, that internal standard is not present in the data"
        abund_dic = individual_dic[confidic['name_abundance']].copy()
        for compartment in abund_dic.keys():
            inte_sta_co = internal_standards_df.loc[abund_dic[compartment].index, :]
            inte_sta_co_is = inte_sta_co[picked_internal_standard]
            # replace zeros to avoid zero division, uses min of the pd series :
            inte_sta_co_is[inte_sta_co_is == 0] = inte_sta_co_is[inte_sta_co_is > 0].min()
            tmp = abund_dic[compartment].div(inte_sta_co_is, axis = 0)
            individual_dic[confidic['name_abundance']][compartment] = tmp
        return individual_dic


def set_samples_names(individual_dic, todrop_x_y):
    trans_dic = dict()
    for k in individual_dic.keys():
        df = individual_dic[k]
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

    individual_dic = trans_dic.copy()
    return

def propagate_nan(confidic, individual_dic):
    # propagates nan from abundance to isotopologues and fractional contributions
    for co in confidic['compartments']:
        abu_co = individual_dic[confidic['name_abundance']][co]
        frac_co = individual_dic[confidic['name_meanE_or_fracContrib']][co]
        frac_co.mask(abu_co.isnull())
        # pending the propagation to isotopologues, think !

    return individual_dic
def drop_metabolites_infile( df ,exclude_list: [str, None]):
    print(exclude_list)
    return 0

def prep_isotopologues_df(isos_df, isosprop_min_admitted: float, stomp_values: bool):
    return 0

def prep_fractionalcontrib_df():
    tableFC = confidic["name_fractional_contribs"].split(".")[0]  # no extension



def do_vib_prep(workingdir, args, confidic):
    # the order of the steps is the one recommended by VIB
    individual_dic = excel2individual_dic(workingdir,  confidic)
    metadata = fg.open_metadata(workingdir, confidic)
    abundance_df = individual_dic[confidic['name_abundance']]
    lod_values, blanks_df, internal_standards_df, bad_x_y = pull_LOD_blanks_IS(abundance_df)
    individual_dic = reshape_individual_dic_elems(individual_dic, metadata, bad_x_y)

    individual_dic = abund_under_lod_set_nan(confidic, individual_dic, lod_values,
                                             args.under_detection_limit_set_nan)

    individual_dic = auto_drop_metabolites(confidic['compartments'],
                                           individual_dic,  lod_values,
                                           args.auto_drop_metabolite_LOD_based)

    individual_dic = abund_substract_blankavg(individual_dic, confidic,
                                              blanks_df, args.substract_blankavg)
    individual_dic = abund_divideby_amount_material(individual_dic, workingdir, confidic['amountMaterial'],
                                                    args.amount_material_div_alternative)
    individual_dic = abund_divideby_internalStandard(individual_dic, confidic, internal_standards_df,
                                                     args.use_internal_standard)

    # prep_isotopologues_df
    # # set_sample_names(df, metadata)
    # drop_metabolites_infile(df, exclude_list)
    # propagate_nan
    # savetables_bycompartment()
    return 0

def do_isocorOutput_prep():
    # individual_dic =

    #individual_dic = abund_divideby_amount_material(individual_dic, workingdir, confidic['amountMaterial'],
    #                                                     args.amount_material_div_alternative)

    # set_sample_names(df, metadata)
    # drop_metabolites_infile(df, exclude_list)
    return 0 #

def do_simple_prep(): #  because user did pretreatment his own

    # individual_dic =

    # individual_dic = abund_divideby_amount_material(individual_dic, workingdir, confidic['amountMaterial'],
    #                                                     args.amount_material_div_alternative)
    # set_sample_names(df, metadata)

    # drop_metabolites_infile(df, exclude_list)
    return 0


def perform_type_prep(workingdir, args, confidic):
    if confidic['typeprep'] == 'vib_prep':
        do_vib_prep(workingdir, args, confidic)
    elif confidic['typeprep'] == 'isocorOutput_prep':
        do_isocorOutput_prep(workingdir, args, confidic)
    elif confidic['typeprep'] == 'simple':
        do_simple_prep(workingdir, args, confidic)

    outputdir = args.config.split("/")[-2] + "/results/prepared_tables/"
    print(outputdir)
    fg.detect_and_create_dir(workingdir + outputdir)

def save_advanced_options_used(args):
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

    save_advanced_options_used







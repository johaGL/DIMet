#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Prepare tables

@author: johanna
"""
import os
import argparse
import pandas as pd
import functions_general as fg


def prep_args():
    show_defaults = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(prog="python -m DIMet.src.prepare",
                                     formatter_class=show_defaults)

    parser.add_argument('config', type=str,
                        help="Configuration file in absolute path")

    return parser


def df_to__dic_bycomp(df: pd.DataFrame, metadata: pd.DataFrame) -> dict:
    # splits df into dictionary of dataframes, each for one compartment:
    out_dic = dict()
    for co in metadata['short_comp'].unique():
        metada_co = metadata.loc[metadata['short_comp'] == co, :]
        df_co = df.loc[:, metada_co['original_name']]
        out_dic[co] = df_co
    return out_dic


def tabs_2_frames_dic(confidic, data_dir) -> dict:
    frames_dic = dict()
    list_config_tabs = [confidic['name_abundance'],
                        confidic['name_meanE_or_fracContrib'],
                        confidic['name_isotopologue_prop'],
                        confidic['name_isotopologue_abs']]
    list_config_tabs = [i for i in list_config_tabs if i is not None]

    for fi_name in list_config_tabs:
        try:
            measu_file = f'{data_dir}{fi_name}.csv'
            tmp = pd.read_csv(measu_file, sep='\t', header=0, index_col=0)
        except FileNotFoundError:
            measu_file = f'{data_dir}{fi_name}.tsv'
            tmp = pd.read_csv(measu_file, sep='\t', header=0, index_col=0)
        except FileNotFoundError:
            measu_file = f'{data_dir}{fi_name}.TSV'
            tmp = pd.read_csv(measu_file, sep='\t', header=0, index_col=0)
        except FileNotFoundError:
            measu_file = f'{data_dir}{fi_name}.CSV'
            tmp = pd.read_csv(measu_file, sep='\t', header=0, index_col=0)
        except Exception as e:
            print("Error in tabs_2_frames_dic : ", e)

        badcols = [i for i in list(tmp.columns) if i.startswith("Unnamed")]
        tmp = tmp.loc[:, ~tmp.columns.isin(badcols)]
        tmp.columns = tmp.columns.str.replace(" ", "_")
        tmp.index = tmp.index.str.replace(" ", "_")
        tmp = tmp.replace(" ", "_", regex=False)
        tmp = tmp.dropna(axis=0, how="all")
        frames_dic[fi_name] = tmp
    return frames_dic


def set_samples_names(frames_dic, metadata):
    compartments = metadata['short_comp'].unique().tolist()

    for tab in frames_dic.keys():
        for co in compartments:
            metada_co = metadata.loc[metadata['short_comp'] == co, :]
            df = frames_dic[tab][co]
            df = df.T
            df.reset_index(inplace=True)
            df.rename(columns={df.columns[0]: "original_name"},
                      inplace=True)
            careful_samples_order = pd.merge(df.iloc[:, 0],
                                             metada_co[['name_to_plot',
                                                        'original_name']],
                                             how="left", on="original_name")
            df = df.assign(name_to_plot=careful_samples_order['name_to_plot'])
            df = df.set_index('name_to_plot')
            df = df.drop(columns=['original_name'])
            frames_dic[tab][co] = df.T

    return frames_dic


def drop_all_nan_metabolites_on_comp_frames(frames_dic, metadata):
    """ metabolites must be in rows """
    compartments = metadata['short_comp'].unique().tolist()
    for tab in frames_dic.keys():
        for co in compartments:
            tmp = frames_dic[tab][co]
            tmp = tmp.dropna(how="all", axis=0)
            frames_dic[tab][co] = tmp
    return frames_dic


def do_prep(args, confidic, meta_path):
    metadata = fg.open_metadata(meta_path)
    fg.verify_metadata_sample_not_duplicated(metadata)
    elems_path_meta = meta_path.split("/")[:-1]
    data_dir = "/".join(elems_path_meta) + "/"

    frames_dic = tabs_2_frames_dic(confidic, data_dir)

    tabs_isotopologues = [s for s in frames_dic.keys() if
                          "isotopol" in s.lower()]
    assert len(
        tabs_isotopologues) >= 1, "\nError, bad or no isotopologues input"

    for k in frames_dic.keys():
        tmp_co_dic = df_to__dic_bycomp(frames_dic[k],
                                       metadata)  # split by compartment
        frames_dic[k] = tmp_co_dic

    frames_dic = drop_all_nan_metabolites_on_comp_frames(frames_dic, metadata)
    frames_dic = set_samples_names(frames_dic, metadata)

    return frames_dic


def perform_prep(args, confidic, meta_path, out_path) -> None:
    output_dir = out_path + "results/"
    fg.detect_and_create_dir(output_dir)

    output_tabs_dir = out_path + "results/prepared_tables/"
    fg.detect_and_create_dir(output_tabs_dir)

    frames_dic = do_prep(args, confidic, meta_path)

    suffix_str = confidic['suffix']
    for k in frames_dic.keys():
        for compartment in frames_dic[k].keys():
            tmp = frames_dic[k][compartment]
            tmp.index.name = "metabolite_or_isotopologue"
            tmp = tmp.reset_index()
            tmp = tmp.drop_duplicates()
            tmp.to_csv(
                f"{output_tabs_dir}{k}--{compartment}--{suffix_str}.tsv",
                sep='\t', header=True, index=False)


if __name__ == "__main__":
    parser = prep_args()
    args = parser.parse_args()
    configfile = os.path.expanduser(args.config)
    confidic = fg.open_config_file(configfile)
    fg.auto_check_validity_configuration_file(confidic)
    fg.verify_good_extensions_measures(confidic)
    confidic = fg.remove_extensions_names_measures(confidic)
    meta_path = os.path.expanduser(confidic['metadata_path'])
    out_path = os.path.expanduser(confidic['out_path'])
    perform_prep(args, confidic, meta_path, out_path)

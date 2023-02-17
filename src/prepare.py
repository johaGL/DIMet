import os
import sys
sys.path.append(os.path.dirname(__file__))
import argparse
import functions_general as fg
import shutil
import yaml
# from .prep_steps import *
#     from .pca_fun import massage_datadf_4pca, calcPCAand2Dplot
#
#     from .abundances_bars import *
#     from .fun_fm import countnan_samples, add_alerts, add_metabolite_column, \
#         add_isotopologue_type_column, save_heatmap_sums_isos, \
#         givelevels, table_minimalbymet, save_rawisos_plot


def prep_args():
    parser = argparse.ArgumentParser(prog="python -m DIMet.src.prepare")

    parser.add_argument('wdir', type=str,
                        help="working directory, absolute path")
    parser.add_argument('config', type=str,
                        help="configuration file, also absolute path")


    parser.add_argument('--nesaispas', type=str,
                        help="sfskfjks")
    parser.add_argument('--option2', type=str,
                        help=":::::")
    return parser

def bonitox(confidic):
    tableFC = confidic["name_fractional_contribs"].split(".")[0]  # no extension




if __name__ == "__main__":

    parser = prep_args()  # combine global and advanced prep args
    args = parser.parse_args()

    fg.wdir_config_confirmation(args.wdir, args.config)

    confifile = os.path.expanduser(args.config)
    confidic = fg.open_config_file(confifile)




    for k in confidic.keys():
        print(k, ":", confidic[k])

    print(args.nesaispas)




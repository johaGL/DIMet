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


def prep_advanced_args(global_args_parser):
    parser = global_args_parser

    parser.add_argument('--nesaispas', type=str,
                        help="sfskfjks")
    parser.add_argument('--option2', type=str,
                        help=":::::")
    return parser

def bonitox(confidic):
    tableFC = confidic["name_fractional_contribs"].split(".")[0]  # no extension




if __name__ == "__main__":

    adv_parser = prep_advanced_args(fg.global_args())  # combine global and advanced prep args
    adv_args = adv_parser.parse_args()

    confifile = os.path.expanduser(adv_args.config)
    with open(confifile, "r") as f:
        confidic = yaml.load(f, Loader=yaml.Loader)

    for k in confidic.keys():
        print(k, ":", confidic[k])

    print(adv_args.nesaispas)




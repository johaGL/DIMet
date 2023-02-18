import os
import sys
sys.path.append(os.path.dirname(__file__))
import argparse
import functions_general as fg
import shutil
import yaml

#

#     from .fun_fm import countnan_samples, add_alerts, add_metabolite_column, \
#         add_isotopologue_type_column, save_heatmap_sums_isos, \
#         givelevels, table_minimalbymet, save_rawisos_plot


def prep_args():
    parser = argparse.ArgumentParser(prog="python -m DIMet.src.prepare")

    parser.add_argument('wdir', type=str,
                        help="working directory, absolute path")
    parser.add_argument('config', type=str,
                        help="configuration file, also absolute path")

    #abundances  prep config IN ORDER ! :  # separate internal standards in another table !
    parser.add_argument("--under_LOD_setnan", action=argparse.BooleanOptionalAction, default=True,
                        help="On VIB results. Abundances under the Limit Of Detection (LOD) become NaN (np.nan). \
                           To note, when all samples are NaN for a metabolite, the metabolite is excluded")  # abus,: NO to internal standards  ?the other t: metabolites blacklisted and excluded?

    parser.add_argument("--substract_blankavg", action=argparse.BooleanOptionalAction, default=True,
                        help="On VIB results. For each metabolite abundances, substract the average of the blanks") # abus only !!!!! NO to internal standards

    parser.add_argument("--amountMaterial", default=None, type=str,
                        help="path to a .tsv file containing the amount of material (number of cells, mg tissue) by sample.\
                         The abundances will be divided by these values (sample specific)")  # NO to internal standards

    parser.add_argument("--useInternalStandard", default=None, type=str,
                        help='The name of the Internal Standard for performing abundances correction (division), \
                        example : "Myristic acid d27". By default is not performed (None)')  # TODO: transform to Myristic_acid_d27  #only abus!!!!!

    # isotopologues :

    parser.add_argument("--isosprop_min_admitted", default=-0.05, type=float, # TODO: change this default to -1 ! keep-0.05 for my cmd file
                        help="Negative proportions can occur (after ElMaven isotopologues correction, for example), \
                        metabolites whose isotopologues are less or equal this cutoff are removed ") # isos only !!!!! !!!!!

    parser.add_argument("--stomp_values", action=argparse.BooleanOptionalAction, default=True, \
                        help="Stomps isotopologues' proportions to max 1.0 and min 0.0. ")

    # both

    parser.add_argument("--exclude_list", default=None, type=str,
                        help="path to .csv file with columns: compartment, metabolite. Metabolites to be excluded")  # all tables affected !!!!!

   # TODO: think if there is a way to split into optional AND advanced arguments ?????

    # parser.add_argument("--zero_to_nan", default="F", type=int,
    #                     help="Optional, with mode prepare. Only keep zero if entire group is zeros\
    #                     i.e. given 3 bio replicates:[0 22 0] change to [NaN 22 NaN].\
    #                     Will be performed in all tables ")

    return parser

def bonitox(confidic):
    tableFC = confidic["name_fractional_contribs"].split(".")[0]  # no extension




if __name__ == "__main__":

    parser = prep_args()  # combine global and advanced prep args
    args = parser.parse_args()

    # systematically README says that metadata has 'former_name','sample' and the other columns, and it is 'sample' what will replace former_name in the tables !
    # and it can be the same if you want to keep the spirit.
    # other_prep NO, better 'universal_prep' :  and set a toyUni1/ with which I call universal data

    fg.wdir_config_confirmation(args.wdir, args.config)

    confifile = os.path.expanduser(args.config)
    confidic = fg.open_config_file(confifile)

    #parser.add_argument("--absolute-isotopologues-available", default="N",
     #                   dest="absolute_isos_avail", help="[Y/N].
     #                     Y if you have  absolute values (not in percentage) for isotopologues. Default : N")
    # this stupidity will be replaced by type_isotopologues : absolute or proportions or ...both ?
    # the advantage of absolute is that, if more than 5 samples, T-test is allowed, whereas proportions allways disfit or KW


    """
    for k in confidic.keys():
        print(k, ":", confidic[k])
    """
    print(confidic['typeprep'])




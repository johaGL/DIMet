"""
needed for DIMet/__main__.py to be able to import .py files in this location (DIMet/src/)
"""
import argparse
import os

import yaml
from .differential_univariate import *
from .abund_frompercentages import *
from .extruder import *

def dimet_message():
    return "DIMet: *D*ifferential *I*sotopically-labeled targeted *Met*abolomics\n"

parser = argparse.ArgumentParser()
parser.add_argument("--mywdir")
parser.add_argument("--config", help = "configuration in yaml file")
parser.add_argument("--mode", help = "timeseriesplots or abundplots or diffabund")
args = parser.parse_args()



confifile = os.path.expanduser(args.config)
with open(confifile, "r") as f:
    confidic = yaml.load(f, Loader=yaml.Loader)

namesuffix = confidic["namesuffix"]
datadi = confidic["datadi"]
extrulist_fi = confidic["extrulist_fi"]
names_compartments = confidic["names_compartments"]
metadata_fi = confidic["metadata_fi"]
levelstimepoints_ = confidic["levelstime"]
condilevels = confidic["conditions"] # <= locate where it is used
tableIC = confidic["name_isotopologue_contribs"]
time_sel = confidic["time_sel"]  # locate where it is used
selectedmetsD = confidic["selectedmets_forbars"] # locate where it is used
tableAbund = confidic["name_abundances"]
max_m_species = confidic["max_m_species"]

# set working directory as current
os.chdir(os.path.expanduser(args.mywdir))

metadata = pd.read_csv(datadi + metadata_fi, index_col=False)

# whatever the option is, prepare output tmp data folder for intermediary files
dirtmpdata = "tmp/"
abunda_species_4diff_dir = dirtmpdata + "abufromperc/"
allfi = os.listdir(datadi)
print()

if not os.path.exists(dirtmpdata):
    print(" 0. Preparing dataset for analysis\n")
    tsvfi = [i for i in allfi if ".tsv" in i]
    print("Your .tsv files in data folder: ", tsvfi, "\n")

    # using list of metabolites to exclude, compartment aware:
    print("using list of undesired metabolites to drop off from data")
    for filename in tsvfi:
        save_new_dfs(datadi, names_compartments,
                     filename, metadata, extrulist_fi, dirtmpdata)

    print("splited (by compartment) and clean files in tmp/ ready for analysis\n")

    # NOTE : for abundances bars and Differential,
    # compulsory step: calculate isotopologues abundances from IC percentages
    saveabundfrompercentagesIC(
        dirtmpdata,
        tableAbund,
        tableIC,
        metadata,
        names_compartments,
        namesuffix,
        abunda_species_4diff_dir,
        max_m_species,
    )
else:
    spefiles = [i for i in os.listdir(abunda_species_4diff_dir)]

if args.mode == "diffabund":
    print("\n 4. Differentially Abundant Metabolites [or Isotopologues] : DAM\n")

    whichtest = confidic["whichtest"]
    newcateg = confidic["newcateg"]  # see yml in example/configs/
    technical_toexclude = confidic["technical_toexclude"]
    contrast = confidic["contrast"]

    outdiffdir = "results/tables/"
    if not os.path.exists(outdiffdir):
        os.makedirs(outdiffdir)

    for co in names_compartments.values():
        rundiffer(
            dirtmpdata,
            tableAbund,
            namesuffix,
            metadata,
            newcateg,
            contrast,
            whichtest,
            technical_toexclude,
            co,
            outdiffdir,
            "TOTAL",
        )

        tableabuspecies_co_ = [i for i in spefiles if co in i]
        for tabusp in tableabuspecies_co_:
            rundiffer(
                abunda_species_4diff_dir,
                tabusp,
                namesuffix,
                metadata,
                newcateg,
                contrast,
                whichtest,
                technical_toexclude,
                co,
                outdiffdir,
                "species",
            )
            # end for tabusp
        # end for co

    print("\nended analysis")
    # end if args.mode == "diffabund"

if args.mode == "abundplots":
    print("!n!n!n!")


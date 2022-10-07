"""
needed for DIMet/__main__.py to be able to import .py files in this location (DIMet/src/)
"""
import argparse
import os
import shutil
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
extrudf_fi = confidic["extrulist_fi"]
names_compartments = confidic["names_compartments"]
metadata_fi = confidic["metadata_fi"]
levelstimepoints_ = confidic["levelstime"]

tableIC = confidic["name_isotopologue_contribs"].split(".")[0] # no extension


tableAbund = confidic["name_abundances"].split(".")[0] # no extension
max_m_species = confidic["max_m_species"]

# set working directory as current
os.chdir(os.path.expanduser(args.mywdir))

metadata = pd.read_csv(datadi + metadata_fi, index_col=False)

isotopolog_style = autodetect_isotop_nomenclature(datadi, tableIC, namesuffix)

print("\nPreparing dataset for analysis\n")
# whatever the option is, prepare output n data folder for intermediary files
allfi = os.listdir(datadi)
dirtmpdata = "tmp/"
abunda_species_4diff_dir = dirtmpdata + "abufromperc/"
print(" [any temporary files (tmp) are being deleted by default]")
if os.path.exists(dirtmpdata):
    shutil.rmtree(dirtmpdata) # clear at each run

os.makedirs(dirtmpdata)

tsvfi = [i for i in allfi if ".tsv" in i]
print("Your .tsv files in data folder: ", tsvfi, "\n")

# using list of metabolites to exclude, compartment aware:
print("using list of undesired metabolites to drop off from data")
extrudf = pd.read_csv(datadi + extrudf_fi, sep=",")

for filename in tsvfi:
    if extrudf.shape[0] > 1:
        save_new_dfs(datadi, names_compartments,
                     filename, metadata, extrudf, dirtmpdata, isotopolog_style)
    else: # table of mets to exclude is empty
        print(save_new_dfs_simple(datadi, names_compartments, filename,
                      metadata,  dirtmpdata))

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
spefiles = [i for i in os.listdir(abunda_species_4diff_dir)]    


if args.mode == "diffabund":
    print("\n testing for Differentially Abundant Metabolites [or Isotopologues] : DAM\n")

    whichtest = confidic["whichtest"]
    newcateg = confidic["newcateg"]  # see yml in example/configs/
    technical_toexclude = confidic["technical_toexclude"]
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
                os.makedirs(x) 

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
                technical_toexclude,
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
                    technical_toexclude,
                    co,
                    outdiffdirs,
                    outkey,
                )
                # end for tabusp
            # end for co
        # end for contrast
    print("\nended analysis")
    # end if args.mode == "diffabund"

if args.mode == "abundplots":
    time_sel = confidic["time_sel"]  # locate where it is used
    selectedmetsD = confidic["selectedmets_forbars"]  # locate where it is used
    condilevels = confidic["conditions"]  # <= locate where it is used
    print("!n!n!n!")

if args.mode == "timeseriesplots":
    tableFC = confidic["name_fractional_contribs"].split(".")[0] # no extension

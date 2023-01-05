#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 19:44:26 2022

@author: johanna
"""

import argparse

import yaml

from abund_frompercentages import *
from abundances_bars import *
from differential_univariate import *
from extruder import *
from frac_contrib_lineplot import *
from isotopologcontrib_stacked import *

parser = argparse.ArgumentParser()
parser.add_argument("--mywdir")
parser.add_argument("--config", help="configuration in yaml file")
args = parser.parse_args()

confifile = os.path.expanduser(args.config)
with open(confifile, "r") as f:
    confidic = yaml.load(f, Loader=yaml.Loader)

print(" 0. Preparing dataset for analysis\n")

os.chdir(os.path.expanduser(args.mywdir))
namesuffix = confidic["namesuffix"]
datadi = confidic["datadi"]
extrulist_fi = confidic["extrulist_fi"]
names_compartments = confidic["names_compartments"]
metadata_fi = confidic["metadata_fi"]

# by default tmp is temporary data storage
dirtmpdata = "tmp/"
metadata = pd.read_csv(datadi + metadata_fi, index_col=False)

allfi = os.listdir(datadi)
tsvfi = [i for i in allfi if ".tsv" in i]
print("Your .tsv files in data folder: ", tsvfi, "\n")


# using list of metabolites to exclude, compartment aware:
print("using list of undesired metabolites to drop off from data")
for filename in tsvfi:
    save_new_dfs(
        datadi, names_compartments, filename, metadata, extrulist_fi, dirtmpdata
    )

print("splited (by compartment) and clean files in tmp/ ready for analysis\n")


print(" 1. Isotopologue contributions : stacked bars\n")

levelstimepoints_ = confidic["levelstime"]
condilevels = confidic["conditions"]
tableIC = confidic["name_isotopologue_contribs"]

# NOTE: Leuven called "CorrectedIsotopologues" to isotopologueContributions (%) !
assert (
    "isotopol" in tableIC.lower()
), "Error!: your table here must have \
    isotopologues percentages, we try to do the stacked bars"

darkbarcolor, palsD = default_colors_stacked()

selbycompD = confidic["groups_toplot_isotopol_contribs"]
# saveisotopologcontriplot(dirtmpdata, tableIC, names_compartments,
#             levelstimepoints_, namesuffix, metadata, selbycompD,
#            darkbarcolor, palsD, condilevels )


print("\n 2. Fractional contributions : line plots\n")
tableFC = confidic["name_fractional_contribs"]
gbycompD = confidic["groups_toplot_frac_contribs"]
coloreachmetab = yieldcolorsbymet()

# savefraccontriplots(dirtmpdata, names_compartments,
#                   metadata, tableFC, namesuffix,
#           gbycompD, coloreachmetab)


# NOTE : for abundances bars and Differential,
# preliminary step: calculate isotopologues abundances from IC percentages
tableAbund = confidic["name_abundances"]
abunda_species_4diff_dir = dirtmpdata + "abufromperc/"
max_m_species = confidic["max_m_species"]
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
all_m_species_ = ["m+" + str(i) for i in range(max_m_species + 1)]
all_m_species_.append("totmk")
spefiles = [i for i in os.listdir(abunda_species_4diff_dir)]

print("\n 3. Abundances : bar plots\n")

time_sel = confidic["time_sel"]


selectedmetsD = confidic["selectedmets_forbars"]

odirbars = "results/plots/abundbars/"
if not os.path.exists(odirbars):
    os.makedirs(odirbars)

SMX = "TOTAL"  # marked and unmarked
for CO in names_compartments.values():
    file_total_co_ = [i for i in os.listdir(dirtmpdata) if tableAbund in i and CO in i]
    assert (
        len(file_total_co_) == 1
    ), "error, multiple abundance files for same comparment"
    abutab = pd.read_csv(dirtmpdata + file_total_co_[0], sep="\t", index_col=0)
    metada_sel = metadata.loc[metadata["sample"].isin(abutab.columns), :]

    # metadata and abundances time of interest
    metada_sel = metada_sel.loc[metadata["timepoint"].isin(time_sel), :]
    abu_sel = abutab[metada_sel["sample"]]

    # total piled-up data:
    piled_sel = stackallabundace(abu_sel, metada_sel)
    piled_sel["condition"] = pd.Categorical(piled_sel["condition"], condilevels)
    piled_sel["timepoint"] = pd.Categorical(piled_sel["timepoint"], time_sel)

    plotwidth = 3.5 * len(selectedmetsD[CO])
    print(f"sending to plot file  :  {selectedmetsD[CO]}")

    printabundbarswithdots(piled_sel, selectedmetsD[CO], CO, SMX, plotwidth, odirbars)

# the legend alone :
oD = mean_sd_D(abu_sel, metada_sel)
piled_alt = tmpstack(oD)
piled_alt["condition"] = pd.Categorical(piled_alt["condition"], condilevels)
piled_alt["timepoint"] = pd.Categorical(piled_alt["timepoint"], time_sel)
printtestpluslegend(piled_alt, selectedmetsD[CO][:2], CO, SMX, 10, odirbars)

# plot bars for the different isotopologues:
for myfi in spefiles:
    SMX, CO = fi2_smx_compartment(myfi, 2, 3)  # also used as titles of the plot
    print(f"running abundances in {CO} : {SMX} ")
    abutab = pd.read_csv(abunda_species_4diff_dir + myfi, sep="\t", index_col=0)
    # metadata selection
    metada_sel = metadata.loc[metadata["sample"].isin(abutab.columns), :]

    # metadata and abundances time of interest
    metada_sel = metada_sel.loc[metadata["timepoint"].isin(time_sel), :]
    abu_sel = abutab[metada_sel["sample"]]

    # total piled-up data:
    piled_sel = stackallabundace(abu_sel, metada_sel)
    piled_sel["condition"] = pd.Categorical(piled_sel["condition"], condilevels)
    piled_sel["timepoint"] = pd.Categorical(piled_sel["timepoint"], time_sel)

    plotwidth = 3.5 * len(selectedmetsD[CO])
    print(f"sending to plot file  :  {selectedmetsD[CO]}")
    print("")
    printabundbarswithdots(piled_sel, selectedmetsD[CO], CO, SMX, plotwidth, odirbars)


##############################################################################
############              DIFFERENTIAL  ANALYSIS            ##################
##############################################################################

print("\n 4. Differentially Abundant Metabolites [or Isotopologues] : DAM\n")

whichtest = confidic["whichtest"]
newcateg = confidic["newcateg"]  # example 'Control_T0'
contrast = confidic["contrast"]
technical_toexclude = confidic["technical_toexclude"]


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
print("\nended analysis")
########## END

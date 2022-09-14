#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 19:44:26 2022

@author: johanna
"""

import yaml
import matplotlib as plt
from extruder import *
from isotopologcontrib_stacked import *
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("--mywdir")
parser.add_argument("--config", help="configuration in yaml file")
args = parser.parse_args()

confifile = os.path.expanduser(args.config)
with open(confifile, 'r') as f:
    confidic = yaml.load(f, Loader=yaml.Loader)

print(" 1. Preparing dataset for analysis\n")

os.chdir(os.path.expanduser(args.mywdir))
namesuffix = confidic['namesuffix']
datadi = confidic['datadi']
extrulist_fi = confidic['extrulist_fi']
names_compartments = confidic['names_compartments']
metadata_fi = confidic['metadata_fi']

# by default tmp is temporary data storage
dirtmpdata = "tmp/"
metadata = pd.read_csv(datadi + metadata_fi, index_col=False) 

allfi = os.listdir(datadi)
tsvfi = [ i for i in allfi if ".tsv" in i ]
print("Your .tsv files in data folder: ", tsvfi, "\n")


# using list of metabolites to exclude, compartment aware:
print("using list of undesired metabolites to drop off from data")
for filename in tsvfi:
    save_new_dfs(datadi, names_compartments, filename, 
             metadata, extrulist_fi, dirtmpdata)

print("splited (by compartment) and clean files in tmp/ ready for analysis\n")


print(" 2. Isotopologue contributions : stacked bars\n")

levelstimepoints_ = confidic['levelstime']
cnds_ = confidic['conditions']
tablePicked = confidic['Config_Isotopologue_Contribs']['pickedTable']  

# NOTE: Leuven called "CorrectedIsotopologues" to isotopologueContributions (%) !
assert "isotopol" in tablePicked.lower(), "Error!: your table here must have \
    isotopologues percentages, we try to do the stacked bars"
    
darkbarcolor, palsD = default_colors_stacked()  
 
selbycompD = confidic['Config_Isotopologue_Contribs']['groups_toplot_isotopol_contribs']
saveisotopologcontriplot(dirtmpdata, tablePicked, names_compartments,
             levelstimepoints_, 
             namesuffix, metadata, selbycompD, palsD, cnds_ )

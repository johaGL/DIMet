#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 19:44:26 2022

@author: johanna
"""

import argparse

import yaml

from extruder import *

parser = argparse.ArgumentParser()
parser.add_argument("--mywdir")
parser.add_argument("--config", help="configuration in yaml file")

args = parser.parse_args()

print(args.mywdir)
print(args.config) # the yaml file path


# yf = os.path.expanduser(args.config)

# with open(yf, 'r') as f:
#     doc = yaml.load(f, Loader=yaml.Loader)


#yaml.reader()

print(" 1. Preparing dataset for analysis\n")

os.chdir(os.path.expanduser("~/DIMet/example/"))
namesuffix = "cycloser"
datadi = "data/"

extrulist_fi = "mets_toexclude.csv"

names_compartments = {"Cell_extracts" : "cell", 
                      "Medium" : "med" }
# by default tmp is temporary data storage
dirtmpdata = "tmp/"
metadata = pd.read_csv(datadi+"metadata_"+namesuffix+".csv", index_col=False) 

allfi = os.listdir(datadi)
tsvfi = [ i for i in allfi if ".tsv" in i ]
print("Your .tsv files in data folder: ", tsvfi, "\n")


# using list of metabolites to exclude, compartment aware:

print("using list of undesired metabolites to drop off from data")
for filename in tsvfi:
    save_new_dfs(datadi, names_compartments, filename, 
             metadata, extrulist_fi, dirtmpdata)

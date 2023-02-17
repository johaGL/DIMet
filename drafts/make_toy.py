#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  9 11:20:09 2022

@author: johanna


"""
import os
import pandas as pd
import numpy as np

os.chdir(os.path.expanduser("~/toy1/smpls_raw/"))

def splitsampledescription(metadata, mycolname):
    time = []
    cond = []
    hour = []
    for s in metadata[mycolname]:
        elems = s.split("_")
        condhere = elems[0]
        timehere = elems[1].split("-")[0]
        hourhere = int(timehere.replace("T", "").replace("h", ""))
        print((condhere, timehere, hourhere))
        time.append(timehere)
        cond.append(condhere)
        hour.append(hourhere)
    metadata["timepoint"] = time
    metadata["condition"] = cond
    metadata["timenum"] = hour
    return metadata

namesuffix = "toy1"


names_compartments = {"Cell_extracts": "cell",
                      "Medium": "med"}


print("formatting metadata")

metadata = pd.read_csv("premeta.csv", index_col=0)
metadata = metadata.replace(" ", "_", regex=True)  # all dataframe " " replaced
metadata.columns = metadata.columns.str.replace(" ", "_")
metadata.reset_index(inplace=True)  # former_name : set as column instead index
print(metadata.columns)
metadata["short_comp"] = metadata["compartment"]
# shorten compartment names:
for nm in names_compartments.keys():
    newnm = names_compartments[nm]
    metadata["short_comp"] = metadata["short_comp"].str.replace(nm, newnm)

metadata = splitsampledescription(metadata, "Sample_description")

# build a column by zipping values of 3 columns
metadata["workingcol"] = metadata["condition"].map(str) + "_" + \
                         metadata["short_comp"].map(str) + "_" + \
                         metadata["timepoint"].map(str)

# def setreplicates_seqs()
metadata["replicate"] = [0 for i in range(metadata.shape[0])]
metadata = metadata.sort_values(axis=0, by="workingcol")
combiuni = set(metadata.workingcol)
# repsD = dict()
for ci in combiuni:
    subtmp = metadata.loc[metadata["workingcol"] == ci]
    orsa = sorted(subtmp["former_name"])
    nreps = subtmp.shape[0]
    aseq = [x for x in range(1, nreps + 1)]
    # ot = zip(salab, aseq)
    # print([i for i in ot])
    for k in range(0, len(orsa)):
        # print(aseq[k])
        # print(metadata.loc[metadata["former_name"] == orsa[k], "replicate" ])
        metadata.loc[metadata["former_name"] == orsa[k], "replicate"] = aseq[k]
    # print(":::::::::::::::\n")

metadata = metadata.sort_values(axis=0, by="former_name")
metadata["sample"] = metadata["condition"].map(str) + "_" + \
                     metadata["short_comp"].map(str) + "_" + \
                     metadata["timepoint"].map(str) + "-" + \
                     metadata["replicate"].map(str)
metadata["sample_descrip"] = metadata["condition"].map(str) + "_" + \
                             metadata["timepoint"].map(str) + "-" + \
                             metadata["replicate"].map(str)
metadata = metadata.drop(columns=["workingcol", "compartment",
                                  "Sample_description", "sample_descrip"])
metadata = metadata[ ['sample', 'condition', 'timepoint', 'timenum', 'short_comp',
            'replicate', 'former_name', 'Sample_label' ] ] # reorder
print("ended metadata formatting, the only unique columns are former_name \
      and sample")
metadata.to_csv("metadata_toy1.csv", index=False)

print("======================End metadata generation ===========================")




finame = "MCF001089_6668_Cycloserine-metabolomics_results.xlsx"

shD = dict()

xl = pd.ExcelFile(finame)
sheetnames = xl.sheet_names

def set_r_todrop(index: list):
    todrop = list()
    for s in index:
        if not s is np.nan:
            try:
                if s.startswith("MCF"):
                    if s not in metadata['former_name'].tolist():
                        todrop.append(s)
            except Exception as e:
                print(str(s))
                print("errr in set_r_todrop   ", e)
    return todrop


for i in sheetnames:
    shD[i] = pd.read_excel(finame, sheet_name=i,
                           engine='openpyxl',
                           header=0,
                           index_col=0)

writer = pd.ExcelWriter("results_toy1.xlsx")

for i in shD.keys():
    tmp = shD[i]
    try:
        rowstodrop = set_r_todrop(list(tmp.index))
        finaltmp = tmp.drop(rowstodrop)
        tmp = finaltmp
    except Exception as e:
        print(e)
    tmp.to_excel(writer, sheet_name=i)
writer.save()



#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  9 11:20:09 2022

@author: johanna
"""
import os
import pandas as pd
import numpy as np

os.chdir(os.path.expanduser("~/example2/smpls_raw/"))
namesuffix = "cycloser"
finame = "data/modifiedfile_Cycloserine.xlsx"
datadi = "data/"
rows_todrop = ["Blank01", "Blank02", "Blank03", "LOD",
               "Cell extracts", "Medium"]

names_compartments = {"Cell_extracts" : "cell", 
                      "Medium" : "med" }


shD = dict()
sheetsnames = ['rawAbundances',
               'FracContribution_C',
               'CorrectedIsotopologues',
               'rawIntensity',
               'samples_sheet']
for i in sheetsnames:
    shD[i] = pd.read_excel(finame, sheet_name = i ,
                           engine='openpyxl',
                           header = 0,
                           index_col = 0)

print("formatting metadata")

def splitsampledescription(metadata, mycolname):
    time = []
    cond = []
    hour = []
    for s in metadata[mycolname]:
        elems = s.split("_")
        condhere = elems[0]
        timehere = elems[1].split("-")[0]
        hourhere = int(timehere.replace("T", "").replace("h",""))
        print((condhere,timehere, hourhere))
        time.append(timehere)
        cond.append(condhere)
        hour.append(hourhere)
    metadata["timepoint"] = time
    metadata["condition"] = cond
    metadata["Hours"] = hour
    return metadata


#get metadata : 
metadata =  shD["samples_sheet"]  
metadata = metadata.replace(" ", "_", regex=True) # all dataframe " " replaced
metadata.columns = metadata.columns.str.replace(" ", "_")
metadata.reset_index(inplace=True)  # former_name : set as column instead index
print(metadata.columns)
metadata["short_comp"] = metadata["compartment"]
#shorten compartment names:
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
#repsD = dict()
for ci in combiuni:
    subtmp = metadata.loc[metadata["workingcol"]==ci]
    orsa = sorted(subtmp["former_name"])
    nreps = subtmp.shape[0]
    aseq = [ x for x in range(1,nreps+1)]    
    #ot = zip(salab, aseq)
    #print([i for i in ot])
    for k in range(0,len(orsa)):
        #print(aseq[k])
        #print(metadata.loc[metadata["former_name"] == orsa[k], "replicate" ])
        metadata.loc[metadata["former_name"] == orsa[k], "replicate" ] = aseq[k]
    #print(":::::::::::::::\n")
    
metadata = metadata.sort_values(axis=0, by="former_name")
metadata["sample"] = metadata["condition"].map(str) + "_" + \
                    metadata["short_comp"].map(str) + "_" + \
                  metadata["timepoint"].map(str) + "-" + \
                      metadata["replicate"].map(str)
metadata["sample_descrip"] = metadata["condition"].map(str) + "_" + \
                  metadata["timepoint"].map(str) + "-" + \
                      metadata["replicate"].map(str)
metadata = metadata.drop(columns=["workingcol", "compartment", "Sample_description"])
print("ended metadata formatting, the only unique columns are former_name \
      and sample")
print("=================================================")
    
print("as a samples_sheet is available for all Sample_Names details (metadata)\
      clear all dataframes from Blank rows, subtitle rows , and so on")
      
def setnewcolnames(adf, metadata):
    nii = pd.DataFrame(data={ 'foname':adf.columns})
    metadahere =  pd.merge(nii, metadata, left_on="foname", right_on="former_name")
    alltrue = 0
    for h in range(len(nii.foname)):
        if metadahere.foname[h] == adf.columns[h]:
            alltrue += 1
    if alltrue == len(adf.columns):
        adf.columns = metadahere["sample"]
        return adf
    else:
        print("ERROR ! when running function setnewcolnames in formatter_be.py")
        return None
    
        
metadata.former_name
metadata.sample
tmpD = dict()
for i in sheetsnames[:-1]: # all sheets except the last need big reformatting
    print(i)    
    tmp = shD[i]
    todrop_here = set(tmp.index.to_list()).intersection(rows_todrop)
    #print(todrop_here)
    tmp = tmp.drop(todrop_here)    
    tmp = tmp.dropna(axis=0, how="all") # if all cells in row == NaN 
    tmp = tmp.loc[:, ~tmp.columns.str.contains("^Unnamed")]     
    tmp = tmp.dropna(axis=1, how="all") # if all ccell in col == NaN             
    print("transposing the dataframe, so samples are columns")
    tmp = tmp.T  
    tmp.index = tmp.index.str.replace(" ", "_")
    #tmp.columns = tmp.columns.str.replace(" ")
    columnsorder = tmp.columns
    tmp = setnewcolnames(tmp, metadata)
    tmpD[i] = tmp
    print("==================")

    """
    print(tmp.columns[:5])
    print(tmp.index[:5])  
    print(tmp.iloc[0:2, 0:2])
    print(tmp.index[:5])
    print(tmp.columns[:5])
    #print(tmp.head)
    """
 
metadata.to_csv(datadi+"metadata_"+namesuffix+".csv", index=False)   
for k in tmpD.keys():
    ofik = datadi + k + "_" +namesuffix+".tsv"
    tmpD[k].to_csv(ofik, sep="\t") 


for k in tmpD.keys():
    hasmyristic = False
    tmp = tmpD[k]
    #print("glu" in tmp.index[0].lower())
    r = 0; indexmyr = None
    while r < len(tmp.index) and hasmyristic == False:
        if "myristic" in tmp.index[r].lower():
            hasmyristic = True
            indexmyr = r
        r += 1
    if hasmyristic:
        print(f" {k} has myristic acid measurements")
        nparr = np.array(tmp )
        arrdiv = nparr / nparr[indexmyr]
        arrdiv = pd.DataFrame(arrdiv)
        arrdiv.index = tmp.index
        arrdiv.columns = tmp.columns
        ofik = datadi + k + "_" + namesuffix + "_normalized.tsv"
        arrdiv.to_csv(ofik, sep="\t") 
    else:
        print(f"{k} has not")
    

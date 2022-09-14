#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  9 11:20:09 2022

@author: johanna
"""
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np


def yieldcolorsbymet():
    # coloreachmetab dictionary: contains many metabolites, just defines colors.
    coloreachmetab = {'L-Lactic_acid': 'blue',
     'Citric_acid': "green",
     'Oxoglutaric_acid': '#903ee0',
     'Succinic_acid': '#019085',
     'Fumaric_acid': '#01668a',
     'L-Malic_acid': '#afc98f',
     'L-Alanine':  "blueviolet" ,
     'Pyruvic_acid':    '#33d3e3' ,
     'L-Glutamic_acid': '#643561',
     'L-Glutamine': '#950065' ,
     'L-Aspartic_acid': '#9d6d00',
     'L-Asparagine': '#ff7c7c' # greenfluo : '#11dc79'
     }
    return coloreachmetab

def yieldfraccountrib(dicos, tablePicked, co):
    """
    description TODO
    """    
    assert "contribution" in tablePicked.lower(), "this funcion only works \
        with fractionnal contribution table or similar"
    dfcompartment = dicos[co][tablePicked].T
    
    metabolites = dfcompartment.columns
    dfcompartment["sample"] = dfcompartment.index
    dfcompartment = pd.merge(dfcompartment, dicos[co]["metadata"], on="sample")
    
    # empty dataframe to fill
    df4plot = pd.DataFrame(columns=["Hours", "condition",  "metabolite", 
                                   "Fractional Contribution (%)"])
    
    for z in range(len(metabolites)):    
        subdf = dfcompartment.loc[:, [metabolites[z], "Hours", "condition"]]
        subdf["metabolite"] = metabolites[z]
        subdf["Fractional Contribution (%)"] = subdf[metabolites[z]] * 100
        subdf = subdf.drop(columns=[metabolites[z]])
        df4plot = df4plot.append(subdf, ignore_index=True)
    return df4plot

def nestedDi2list(anestedD):
    elems_ = []; 
    for p in anestedD.values():
        for q in p:
            elems_.append(q)
    return elems_

def df2musddf(one_m):
    """
    input: dataframe with "condition", "Hours", "mean", "sd", "metabolite"
    returns dataframe with  """
    print("-*-*-*-*-*-*-*-")
    print(set(one_m["metabolite"]))
    m_s = one_m[["condition", "Hours"]].drop_duplicates()
    print(m_s.columns)
    m_s["mean"] = 0
    m_s["sd"] = 0
    for kc in m_s["condition"]:
        for hou in one_m["Hours"]:
            ss = one_m.loc[(one_m["Hours"] == hou ) & (one_m["condition"]==kc)]
            mymean = np.mean(ss["Fractional Contribution (%)"])
            mysd = np.std(ss["Fractional Contribution (%)"])
            m_s.loc[(m_s["condition"] == kc) & (m_s["Hours"] == hou), "mean"] = mymean
            m_s.loc[(m_s["condition"] == kc) & (m_s["Hours"] == hou), "sd"] = mysd
    return m_s
    
def complextimetracer(co, df4plot, grmetsD,  mycolorsD, outfile):
    themets = nestedDi2list(grmetsD)
    somem = df4plot.loc[df4plot["metabolite"].isin(themets)]
    #print(somem.columns) # ['Hours', 'condition', 'metabolite', 'Fractional Contribution (%)']
    m_s = pd.DataFrame(columns = ["condition", "Hours", "mean", "sd", "metabolite"] )
    for k in set(somem["metabolite"]):  
        print(k)
        one_m = somem[somem["metabolite" ] == k]
        m_s1 = df2musddf(one_m)
        m_s1["metabolite"] = k
        print(m_s1.iloc[0])
        m_s = pd.concat([m_s, m_s1])  # mean and sd ( bio replicates)
    figziz = 5.2 * len(grmetsD)    
    fig, axs = plt.subplots(1, len(grmetsD), sharey=False, figsize = (figziz, 6))
    handless = []
    labels = []
    for z in range(len(grmetsD)):
        sns.lineplot( ax = axs[z] ,
            x="Hours", y="Fractional Contribution (%)", 
                     hue="metabolite",
                     style = "condition",
                     err_style=None,
                     alpha=0.9, palette = mycolorsD,
                     data= somem.loc[somem["metabolite"].isin(grmetsD[z])],
                     legend=True ## !!!!!!!!!!!!!!!!attention here   <=====
                     )
        #for k in set(somem["metabolite"]):
        m_s1 = m_s.loc[m_s["metabolite"].isin(grmetsD[z])]
        print(m_s1)
        axs[z].scatter(m_s1["Hours"],
                    m_s1["mean"], s=16,
                    facecolors='none', edgecolors='black')
        axs[z].errorbar(m_s1["Hours"], m_s1["mean"], 
                        yerr = m_s1["sd"],
                      fmt = 'none', capsize = 3, 
                      ecolor = 'black', 
                      zorder = 1)
        axs[z].set(ylabel = None)
        
        ha , la = axs[z].get_legend_handles_labels()
        handless.append(ha)
        labels.append(la)
        #plt.title(co)
        #plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        #plt.savefig(outfile, format="pdf", bbox_inches='tight')
    fig.suptitle(co)
    fig.text(0.015, 0.24, "Fractional Contribution (%)", va = "center", rotation = 'vertical')
    #loc="upper left"
    #fig.legend(handless, labels, bbox_to_anchor=(0.65,1.05))
    fig.legend(handles = handless, labels = labels)
    fig.savefig(outfile, format="pdf")
    return "hungry3"



os.chdir(os.path.expanduser("~/fluxdiff/"))
namesuffix = "cycloser"
ometadata = "data/metadata.csv"
datadi = "data/"

tablePicked = 'FracContribution_C'  # desired table

# do a dictionary with groups:  try to group them by similarity in values !!! 

gbycompD = {    
   "cell" : { 0 : ["L-Lactic_acid" , "Pyruvic_acid" ],
           1 : [ "L-Alanine"] ,
           2:  [  "Fumaric_acid", "L-Malic_acid"] ,
           3 : ["L-Glutamic_acid", "L-Glutamine"], 
           4 :  ["L-Aspartic_acid" ] },
   "med" :  {  0 : ["L-Lactic_acid" , "Pyruvic_acid" ],
           1 : [ "L-Alanine"] ,
           2: [ "L-Glutamine","L-Asparagine"]   }
    }


names_compartments = {"Cell_extracts" : "cell", 
                      "Medium" : "med" }
coloreachmetab = yieldcolorsbymet()

for co in names_compartments.values():
    print(co)
    #co = "cell"
    adf = pd.read_csv(datadi+ tablePicked + "_"+namesuffix+"_"+co+".tsv",
                          sep="\t", index_col=0 )
    metadata = pd.read_csv(datadi+"metadata_"+namesuffix+".csv", index_col=False) 
    
    # check which tables have myristic acid : 
    # DONE in fromatter_be.py !!!! 
    dicos = dict()
    dicos[co] = {}
    dicos[co]["metadata"] = metadata.loc[metadata.short_comp == co]
    dicos[co][tablePicked] = adf[dicos[co]["metadata"]["sample"]]
    
    df4plot = yieldfraccountrib(dicos, tablePicked, co)
    
    grmetsD = gbycompD[co]
    
    sns.set_style({'font.family'  :'sans-serif',
                   'font.sans-serif' : 'Liberation Sans'})
    plt.rcParams.update({'font.size' : 22})
    # https://stackoverflow.com/questions/53137983/define-custom-seaborn-color-palette    
    
    ofile = 'fc_plot_{}.pdf'.format(co)
    print(complextimetracer(co, df4plot, grmetsD, coloreachmetab, ofile))
    
    #### uninteresting metabolites : quick ugly print : 
    pickedmets_ = nestedDi2list(grmetsD)
    thiscompartment_mets = set(df4plot['metabolite'])
    
    trashmets = set(thiscompartment_mets) - set(pickedmets_)
    ppp = df4plot.loc[df4plot["metabolite"].isin(trashmets)]
    g = sns.FacetGrid(ppp, row="condition",  
                      col="metabolite" ,
                      sharey=False,
                      margin_titles=True)
    g.map(sns.lineplot, "Hours", "Fractional Contribution (%)",
                 marker = "o", err_style = "bars",
                 ci=95, color="black", alpha=0.4, size=5)    
    g.savefig("fc_trash"+co+".pdf", format = "pdf")



#######################
# interestingmets = ["L-Lactic_acid" , "Citric_acid" , "Oxoglutaric_acid",       
#                 "Succinic_acid",           
#                  "Fumaric_acid", "L-Malic_acid", 
#                 "L-Alanine",               
#                 "Pyruvic_acid"  , "L-Glutamic_acid", "L-Glutamine", 
#                 "L-Aspartic_acid", 
#                 "L-Asparagine"]

# colors3 = ["#9d6d00", "#903ee0", "#11dc79", "#f568ff", "#419500", "#013fb0", 
#           "#f2b64c", "#007ae4", "#ff905a", "#33d3e3", "#9e003a", "#019085", 
#           "#950065", "#afc98f", "#ff9bfa", "#83221d", "#01668a", "#ff7c7c", 
#           "#643561", "#75608a"]
# colors4 = ["blue", "#33d3e3", "#903ee0", "#019085", 
# "#950065", "#afc98f", "#9d6d00",  "#007ae4"  , "#01668a", "#ff7c7c", 
# "#643561", "#11dc79" ] # "#75608a"
# colors4D
# set(colors3) - set(colors4)
# colors4D = dict()
# for l in range(len(colors4)):
#     print(l)
#     colors4D[interestingmets[l]] = colors4[l]



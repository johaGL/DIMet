#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 31 12:58:52 2022

@author: johanna

Calculates absolute abundances of isotopologues
"""


from fun_fm import *


################# fun

def correction_prop_df(prop_df):
    """"
    does two corrections:
        a) all multiply 100 (percentage)
        b )if # if % < 0 , or > 100, cut 
    """
    prop_df = prop_df * 100
    prop_df[prop_df < 0 ] = 0
    prop_df[prop_df > 100] = 100
    return prop_df

def makematch_abund_rowdata(abund, rowdata):
    """
    check lines that are in abundance but not in corrected %
     NOTE: rowdata has been obtained from prop_df (%) 
    """
    set(abund.index) - set(rowdata['metabolite'] )
    rtodrop = list(set(abund.index) - set(rowdata['metabolite'] ))
    abund = abund.drop(rtodrop,axis=0)
    return abund

def getspecificmk(prop_df, rowdata,selmk):
    """
    obtain proportion of CI for specified label selmk
    example : selmk = "m+0"
    output: columns are [samples], rownames metabolites
    """
    tmp = prop_df.copy()
    tmp = tmp.assign(isotopolgFull = tmp.index)
    tmp = pd.merge( tmp, rowdata, on="isotopolgFull", how="inner")
    tmp = tmp.loc[tmp["m+x"] == selmk,:]
    tmp = tmp.set_index("metabolite") # make metabolite to be the index
    tmp = tmp.drop(["isotopolgFull", "m+x"], axis=1)
    return tmp

def yieldmarkedabu(prop_mx, abund, *totalmarked):
    """  
    * calculates abundance using proportion specific m+x dataframe *
    NOTE: if option totalmarked is 'totalmarked' (not case sentitive):
    * * calculates by substracting mzero abundance from total abundance * *
    Parameters
    ----------
    prop_mx : pandas
        example : prop_mx = getspecificmk(prop_df,rowdata, "m+0")
          so prop_mx :    sampleConditionA-1   sampleConditionA-2 ....
                   metabolites           0                      0
              
    abund : pandas

    Returns
    -------
    abu_mx : pandas
    """
    total_opt = [o for o in totalmarked]
    #make sure same order in columns
    ordercol = prop_mx.columns
    abund = abund[ordercol]
    markedout = dict()
    # calculate by rows, across mx marked proportions
    # note: if prop_zero as input, row will be zero marked proportions
    for i,row in prop_mx.iterrows():
        nameindex = row.name # the metabolite in row 
        totabu = abund.loc[nameindex,:]
        if len(total_opt) == 0 : 
            abumk = (totabu * row) / 100
        else :
            if total_opt[0].lower() == "totalmarked":
                zeroabu = (totabu * row) / 100 # zero marked abundance, 
                abumk = totabu - zeroabu  # total - zero marked abundance 
            else : 
                # wrote string not recognized (total only accepted)
                print("3dr arg not recognized, must be 'totalmarked' or empty")
                return 1
        markedout[nameindex] = abumk
        
    abu_mx = pd.DataFrame(markedout).T
    return abu_mx

def saveabundance(df, nameout):
    df.to_csv(nameout,  sep="\t") 

def callfuns_perc2abu(odir, dicos, co, tableAbund, tableIC, \
          selmk, *totalmarked):
    """
    Parameters
    ----------
    dicos : TYPE
        DESCRIPTION.
    co : TYPE
        DESCRIPTION.
    tableAbund : TYPE
        DESCRIPTION.
    tableIC : TYPE
        DESCRIPTION.
    selmk : TYPE
        DESCRIPTION.
    *totalmarked : string
        "totalmarked" (if willing totalMARKED abundance, set selmk = m+0)

    Returns
    -------
    str
        DESCRIPTION.

    """
    total_opt = [o for o in totalmarked]
    abund = dicos[co][tableAbund]
    prop_df = dicos[co][tableIC]
    prop_df = correction_prop_df(prop_df)
    rowdata = yieldrowdata(prop_df)  # this funciton is from fun_fm
  
    prop_mx = getspecificmk(prop_df, rowdata, selmk)
    if len(total_opt) == 0:
        abu_mx = yieldmarkedabu(prop_mx, abund)
        nameout =  f"{odir}abux_byProp_{selmk}_{co}.tsv"
        saveabundance(abu_mx,  nameout)
    else:
        if total_opt[0].lower() == "totalmarked":
            abu_mkall = yieldmarkedabu(prop_mx, abund, "totalmarked")
            nameout =  f"{odir}abux_byProp_totmk_{co}.tsv"
            saveabundance(abu_mkall, nameout)
    return 0
   

def saveabundfrompercentagesIC(datadi, tableAbund, tableIC, metadata,
                               names_compartments, namesuffix, odir, max_m_species):
    if not os.path.exists(odir):
        os.makedirs(odir)

    for co in names_compartments.values():
        abun = pd.read_csv(datadi+ tableAbund + "_"+  namesuffix +\
                           "_"+co+".tsv", sep="\t", index_col=0 )
        # note that pandas automatically transform any 99.9% in decimal 0.999

        propc = pd.read_csv(datadi+ tableIC + "_"+ namesuffix +\
                            "_"+co+".tsv",  sep="\t", index_col=0 )

        dicos = dict()
        dicos[co] = {}
        dicos[co]["metadata"] = metadata.loc[metadata.short_comp == co]
        selecols = dicos[co]["metadata"]["sample"]
        dicos[co][tableAbund] = abun[selecols]
        dicos[co][tableIC] = propc[selecols]

        # by default calculates total marked (taking away m+0 unmarked)
        callfuns_perc2abu(odir, dicos, co, tableAbund, tableIC, "m+0", "totalmarked")
        callfuns_perc2abu(odir, dicos, co, tableAbund, tableIC, "m+0")   # true unmarked

        for k in range(1,max_m_species+1):
            callfuns_perc2abu(odir, dicos, co, tableAbund, tableIC,  f'm+{k}')


    return 0



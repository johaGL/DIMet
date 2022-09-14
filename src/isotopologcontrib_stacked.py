#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 18 11:14:34 2022

@author: johanna
"""

import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches



def icontrib_2df4plot(dicos, tablePicked, co, levelstimepoints_):
    """
    Imitates de behaviour of a 'melt', but this function is more secure
    example:
      pd.merge output, where 1st column has the values we want:
    
        L-Phenylalanine_C13-label-2    timepoint   condition ...  sample descr...
        0.01                           T0          Control        xxxxxx
        
      is transformed into:
          timepoint    condition    isotopolgFull        Isotopologue Contrib  
          T0           control      L-Phenylala...    0.01
        
    """    
    dfcompartment = dicos[co][tablePicked].T
    
    metabolites = dfcompartment.columns
    dfcompartment["sample"] = dfcompartment.index
    dfcompartment = pd.merge(dfcompartment, dicos[co]["metadata"], on="sample")
    # empty dataframe to fill
    df4plot = pd.DataFrame(columns=["timepoint", "condition",  "isotopolgFull", 
                                   "Isotopologue Contribution (%)"])

    df4plot['timepoint'] = pd.Categorical( df4plot['timepoint'], 
                          levelstimepoints_ )
    # iteratively pile up
    for z in range(len(metabolites)):    # corrected! error I had : -30)): now ok :)     
        subdf = dfcompartment.loc[:,[metabolites[z], "timepoint", "condition"]]
        subdf["isotopolgFull"] = metabolites[z] # 1st colname as cell value, reps
        subdf["Isotopologue Contribution (%)"] = subdf[metabolites[z]] * 100
        subdf = subdf.drop(columns=[metabolites[z]]) # 1st col no longer needed 
        df4plot = pd.concat([df4plot, subdf])
        del(subdf)
    return df4plot


def massageisotopologues(df4plot):
    """
    # dealing with name style : isotopolog(variable)_C13_label-xx and .._PARENT
    C13 and PARENT are mutually exclusive
    and also correcting weird values
    """
    xu = {"name" : [ ] , "m+x" : []}
    for ch in df4plot["isotopolgFull"]:
        if "_C13-" in ch:
            elems = ch.split("_C13-")
            xu['name'].append(elems[0])
            xu['m+x'].append('m+{}'.format(elems[-1].split("-")[-1]))        
        elif "_PARENT" in ch:
            elems = ch.split("_PARENT")
            xu['name'].append(elems[0])
            xu['m+x'].append('m+0')    
    df4plot["metabolite"] = xu['name']
    df4plot["m+x"] = xu["m+x"]     
    # dealing with weird values: bigger than 100 and less than 0 :
    df4plot.loc[df4plot['Isotopologue Contribution (%)'] > 100 , \
        'Isotopologue Contribution (%)' ]  = 100
    #amin = min([ x for x in df4plot['Isotopologue Contribution (%)'].tolist() if x > 0])
    #print(amin)
    df4plot.loc[df4plot['Isotopologue Contribution (%)'] < 0 , \
        'Isotopologue Contribution (%)' ]  = 0  

    return df4plot
    

def preparemeansreplicates(df4plot, cnd, selectedmets, levelstimepoints_):
    """
    condition : 'Control', or another
    """
    secnd = df4plot.loc[df4plot["condition"] == cnd]
    secnd = secnd.groupby(["condition","metabolite", "m+x", "timepoint"]).mean()
    secnd = secnd.reset_index()

    ohmeh = dict()
    for i in selectedmets:
        tmp = secnd.loc[secnd["metabolite"] == i, ].reset_index()
        tmp['timepoint'] = pd.Categorical( tmp['timepoint'], 
                             levelstimepoints_ )
        ohmeh[i] = tmp.sort_values(by = 'm+x', axis=0, ascending=True, inplace=False)
    return ohmeh
        

 
def complexstacked(co, cnd, selectedmets, ohmeh, palsD, outfilename, figuziz):
    """plot highly custom, recommended that selectedmets <= 6 subplots"""
    ### set font style
    sns.set_style({ 'font.family': 'sans-serif', 
                     'font.sans-serif' : 'Liberation Sans'   }) 
    f, axs = plt.subplots(1, len(selectedmets), 
                          sharey=True, figsize = (figuziz, 7))
    plt.rcParams.update({'font.size': 24})
    meh = 0
    for z in range(len(selectedmets)): 
        #sns.set_style({ 'font.family': 'sans-serif', 
         #                'font.sans-serif' : 'Liberation Sans'   }) 
        axs[z].set_title(selectedmets[meh])
        sns.histplot(
            ax = axs[z],
                data=ohmeh[selectedmets[meh]],
                x='timepoint',
                # Use the value variable here to turn histogram counts into weighted
                # values.
                weights='Isotopologue Contribution (%)',
                hue='m+x',
                multiple='stack',
                palette= palsD,   #['#0000c0',  '#410257' ] , # '#440154FF'],
                # Add white borders to the bars.
                edgecolor='black',
                # Shrink the bars a bit so they don't touch.
                shrink=0.8, 
                alpha = 1,
                legend = False
            ) # end sns.histplot
        axs[z].tick_params(axis="x",labelrotation=45)
        axs[z].tick_params(axis='y', length=11, labelsize = 28)
    
        for bar in axs[z].patches:
            selcol = "black" # black, but if color is in whitelist, do different
            #print(bar.get_facecolor())
            herergba = bar.get_facecolor()
            if herergba in darkbarcolor:
                selcol = "white"
            thebarvalue = round(bar.get_height(), 1)
            if thebarvalue >= 100:
                thebarvalue = 100  # no decimals if 100
            if round(bar.get_height(), 1) > 4:
                axs[z].text(
                    # Put the text in the middle of each bar. get_x returns the start
                    # so we add half the width to get to the middle.
                    bar.get_x() + bar.get_width() / 2,
                    # Vertically, add the height of the bar to the start of the bar,
                    # along with the offset.
                    (bar.get_height() / 2) + (bar.get_y() ) + 2 ,#
                    # This is actual value we'll show.
                    thebarvalue,
                    # Center the labels and style them a bit.
                    ha='center',
                    color=selcol,
                    size= int( ( figuziz/len(selectedmets)) * 3) + 1
                    ) # end axs[z].text
            else:
                continue
            # end if round(...)
        # end for bar
        plt.ylim(0,100)
        plt.gca().invert_yaxis() # invert, step1
        plt.yticks(np.arange(0,100,20), np.arange(100,0,-20)) # invert, step2
        axs[z].set_ylabel('Isotopologue\nContribution (%)', size=26)
        axs[z].set_xlabel('', size=22)
        
        # ax.legend() no legend because goes into grid
        meh += 1
    # end for z
    f.subplots_adjust(hspace=0)
    f.suptitle(f'{co.upper()} {cnd.upper()}\n', fontsize=  18)
    
    f.savefig(outfilename, format="pdf")
    return 0


def default_colors_stacked():
    """
    Returns colors defined by default
    """
    darkbarcolor = [ (0.0, 0.0, 0.7529411764705882, 1.0),
                (0.2549019607843137, 0.00784313725490196, 
                 0.3411764705882353, 1.0)]
    palsD = {'m+0': '#410257' , 'm+1' : '#0000c0' , 'm+2': '#55a0fb' ,
         'm+3' :'#41f0ae',  'm+4' :'#addead',  'm+5' :'#eec046', 
         'm+6' : '#ffa040'  }
    return darkbarcolor, palsD


def saveisotopologcontriplot(datadi, tablePicked, names_compartments,
             levelstimepoints_, 
             namesuffix, metadata, selbycompD, palsD,  cnds_ ) :
    for co in names_compartments.values(): #
        print(co)  
        
        adf = pd.read_csv(datadi+ tablePicked + "_"+namesuffix+"_"+co+".tsv",
                              sep="\t", index_col=0 )
        # note that pandas automatically transform any 99.9% in decimal 0.999
     
        dicos = dict()
        dicos[co] = {}
        dicos[co]["metadata"] = metadata.loc[metadata.short_comp == co]
        dicos[co][tablePicked] = adf[dicos[co]["metadata"]["sample"]]
        
        # call complicated functions
        
        df4plot = icontrib_2df4plot(dicos, tablePicked, co, levelstimepoints_)
        df4plot = massageisotopologues(df4plot)
        ####
        # conditions to plot in desired order : 
        ####
       
        cnd = cnds_[0]  # test with only the first one
        
        metscustomgroups = selbycompD[co]
        print(sorted(list(set(df4plot["timepoint"]))))
        for cnd in cnds_:
            for j in range(len(metscustomgroups)):
                selectedmets = metscustomgroups[j]
                print(selectedmets)
                outfname = 'ic_{}_{}_group{}.pdf'.format(co,cnd,j)
                print(outfname)
                ohmeh = preparemeansreplicates(df4plot, cnd, selectedmets, levelstimepoints_)       
                ohmeh.keys()  # just the metabolites subframes, one co, one cnd
                figuziz = 7.5 * len(selectedmets) # note, change width  
                complexstacked(co, cnd, selectedmets, ohmeh, palsD, outfname, figuziz)
        
        #plt.show()
        # legend alone
        myhandless = []
        for c in palsD.keys():   
            paobj = mpatches.Patch(color = palsD[c], label= c)
            myhandless.append(paobj)
        plt.legend(handles=myhandless)  
        plt.savefig("ic_legendTEST.pdf", format="pdf")
        print("ended blahblahcontributio")
    return 0



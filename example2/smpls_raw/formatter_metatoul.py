#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 12:10:53 2022

Inputs : MetaToul style files corrected with isocor
Output : abundance matrix, metadata, row_data, and any other related

REcall: T0 == NORMOXIA.
       T24 and T48 Hypoxia

@author: johanna
"""

import os
import pandas as pd
import numpy as np



def giveEndo_metadata(samples_vec):
    metadata = pd.DataFrame({'name_old': samples_vec})
    metadata = metadata.assign(timepoint = metadata['name_old'].str.split("_").str[0])
    metadata = metadata.assign(genotype = metadata['name_old'].str.split("_").str[1])

    metadata = metadata.assign(replicate_bio = metadata['name_old'].str.split("_").str[2])
    metadata['genoplust'] = metadata.genotype.str.cat(metadata.timepoint, sep="_")
    metadata['compartment'] = "cell"
    metadata['presample'] = metadata.genoplust.str.cat(metadata.compartment, sep="_")
    metadata['sample'] = metadata.presample.str.cat(metadata.replicate_bio, sep="-")
    metadata = metadata.drop(columns = ['presample'])
    return metadata

def securecolnameschange(metadata, abutot_mat):
    goodorder_colnames = list()
    for i in abutot_mat.columns.tolist():
        tmpok = metadata.loc[metadata['name_old'] == i, 'sample'].tolist()[0]
        goodorder_colnames.append(tmpok)
    abutot_mat.columns = goodorder_colnames
    return abutot_mat

def giveExoPRE_metadata(samples_vec):
    metadata = pd.DataFrame({'sample': samples_vec})
    metadata = metadata.assign(replicate_bio = metadata['sample'].str.split("-").str[1])
    metadata['genoplust'] = metadata['sample'].str.split("-").str[0]
    metadata['compartment'] = "med"
    return metadata

mywd = "~/example2/smpls_raw/"
endometabo_infi1 = "dec2020_ALL_DATA_res_isocor.xlsx"
# more endometabo files: glioblastome/Metabo_glioblastome/Metatoul_analysis/
exometabo_infi = "Appendix_Data_RMN_quanti&marquageMOD.xlsx"
sheetsnames_exometabo = ['Data Quantification', 'Data Marquage Enrichment'] # 'Acceptance Criteria'
# endometabo_infi1 is a single sheet file, with corrected areas only
# note : areas == abundance
# see ownCloud/glioblastome/Metabo_glioblastome/Metatoul_analysis/README.txt
aberrant_samples_exometabo = ['CÆ 24h - 2', 'CÆ 48h - 2']

# start
os.chdir(os.path.expanduser(mywd))

# ##############################
# ENDOMETABO
# ######################


# endo open
tmp_table = pd.read_excel(endometabo_infi1,
                           sheet_name = 'Feuil1',
                           header=0, index_col=False)
endo_table = tmp_table[['sample', 'metabolite', 'isotopologue', 'corrected_area']]
trueiso_v = endo_table.metabolite.str.cat(endo_table.isotopologue.astype(str),sep = "_label")
endo_table = endo_table.assign(trueiso = trueiso_v)


# build metadata 
samples_endo = sorted(list(set(endo_table['sample'].tolist())))
metadata_endo = giveEndo_metadata(samples_endo)

metadata_endo.to_csv("metadata_endo_ldh.csv", header=True, index=False)

# build matrices : type pandas, for practicity

# fill isotopologues :
# # de-activated as it takes long (moreover, huge decimal parts in those floats!)
# endoisotopologue_outfi = "endoCorrectedIsotopologues_ldh.tsv"
# if not os.path.exists(endoisotopologue_outfi):
#     endo_trueiso_mat = pd.DataFrame(0, index=endo_table['trueiso'].tolist(),
#                                     columns=metadata_endo['name_old'].tolist())
#     print("printing isotopologues")
#     for i in endo_trueiso_mat.index:
#         for j in endo_trueiso_mat.columns:
#             vx = endo_table.loc[(endo_table['trueiso'] == i) & (endo_table['sample']== j),: ]
#             endo_trueiso_mat.loc[i,j] = float(vx['corrected_area'])
#     endo_trueiso_mat.columns = metadata_endo['sample'].tolist()    # change colnames
#     endo_trueiso_mat.to_csv(endoisotopologue_outfi, header=True, sep="\t")
#     print("ended isotopologues")


# fill total abundance
endo_metabolites_orig = list(set(endo_table['metabolite'].tolist()))
endo_abutot_mat = pd.DataFrame(0, index = endo_metabolites_orig, 
                              columns = metadata_endo['name_old'].tolist())
for i in endo_metabolites_orig:
    for j in endo_abutot_mat.columns:
        vx = endo_table.loc[(endo_table['metabolite'] == i) & (endo_table['sample']== j),: ]
        endo_abutot_mat.loc[i,j] = float(vx['corrected_area'].sum())
endo_abutot_mat = endo_abutot_mat.astype(float)


# change colnames
#endo_abutot_mat.columns = metadata_endo['sample'].tolist()
endo_abutot_mat = securecolnameschange(metadata_endo, endo_abutot_mat)
endo_abutot_mat.to_csv("endo_abundance_ldh.tsv", header=True, sep="\t")

# rowdata
# rowdata_endo = pd.DataFrame({'met_orig': endo_abutot_mat.index.tolist(),
#                              'compartment': 'cell'})
# rowdata_endo.to_csv("rowdata_endo.tsv", header=True, sep='\t', index=False)

# ##############################
# EXOMETABO
# ######################

# exo open
exometabo = dict()
for i in sheetsnames_exometabo:
    tmp = pd.read_excel(exometabo_infi, sheet_name = i ,
                       engine='openpyxl',
                        header = 0, index_col = 0)
    #print(tmp.columns)
    #tmp.set_index('Unnamed: 0', inplace=True) 
    try:
        tmp = tmp.drop(columns=['Unnamed: 13'])
    except:
        pass
    
    try: # delete aberrant_samples_exometabo 
        tmp = tmp.loc[~tmp.index.isin(aberrant_samples_exometabo),:]
    except:
        pass
    exometabo[i] = tmp
koo = exometabo['Data Marquage Enrichment']
koo.shape
print("original samples names [data received as it], note samples nomenclature must be improved")
print(koo.index)

############## improving data samples nomenclature ###########################
print(" *** improving data samples nomenclature ***")
for i in exometabo.keys():
    print(exometabo[i].shape)
    
exo_true = exometabo['Data Quantification']

# index has weird characters, correct:
indexgood = exo_true.index.str.replace(" - ", "-")
indexgood = indexgood.str.replace(" -", "-")
indexgood = indexgood.str.replace(" ", "_")
exo_true.index = indexgood


### visualize to see the weird samples
from sklearn import decomposition
import seaborn as sns
import matplotlib.pyplot as plt
forpca = exo_true*1000
forpca = forpca.replace(0,1)
pca = decomposition.PCA(n_components=3)
pc = pca.fit_transform(forpca)
pc_df = pd.DataFrame(data= pc, columns = ['PC1','PC2', 'PC3'])
fff = exo_true.index.str[:-2]
pc_df['genoplust'] = fff

dfvar = pd.DataFrame({'var':pca.explained_variance_ratio_,
             'PC':['PC1','PC2','PC3']})
# sns.barplot(x='PC',y="var", data=dfvar, color="c")
# plt.show()
# plt.close()
print()
sns.set_style("whitegrid")
sns.scatterplot( x="PC1", y="PC2",
  data=pc_df,
  hue='genoplust', style='genoplust',
  legend=True)
plt.axhline(0,ls="--", color="gray" )  
plt.axvline(0,ls="--", color="gray" ) 
plt.xlabel(f'PC1 { round( dfvar.iloc[0,:]["var"] , 2) *100 } %')
plt.ylabel(f'PC2 { round( dfvar.iloc[1,:]["var"] , 2) *100 } %')
plt.title("exometabolome")
plt.savefig("plots/exometabo_weird.png", dpi=300)
## end visualize to see the weird samples

# *** assigning CÆ_ columns to AB_T0 and CONT_T0 (decision from PCA plot) ***
# note: 'T0' alone is not possible to know if belongs to A or B genotype
# 'T0' alone: drop off
indexgoodII = exo_true.index.str.replace("CÆ_24h", "AB_T0")
indexgoodII = indexgoodII.str.replace("CÆ_48h", "CONT_T0")
exo_true.index = indexgoodII
exo_true = exo_true.drop(labels=['T0-1', 'T0-2', 'T0-3'], axis=0)
exo_true.index = exo_true.index.str.replace("CONT", "Cont")
exo_true.index = exo_true.index.str.replace("H", "")
exo_true.index = exo_true.index.str.replace("h", "")
exo_true.index = exo_true.index.str.replace("24", "T24")
exo_true.index = exo_true.index.str.replace("48", "T48")

exo_true = exo_true.T
exo_true.to_csv("exo_abundance_totalPRE.tsv", header=True, sep="\t")

# metadata exo
metadata_exo = giveExoPRE_metadata(exo_true.columns)
print("[in these particular data] exometabolome has been problematic, samples missing compared to endometabolome")
print(metadata_exo.shape)
print(exo_true.shape)


print("filling the missing samples with NAs in exometabo abundances")
missingsamples = set(metadata_endo.genoplust.str.cat(metadata_endo.replicate_bio,"-")) - set(exo_true.columns)
# reliminary steps : concerning metadata
# fill metadata with those samples, even if unexistant
# for simplicity, instead filling existing metadata, take advantage of
# already created endometabo metadata:  just change the compartment
mdexo = metadata_endo.copy()
mdexo = mdexo.assign(compartment = "med")

print()
mdexo = mdexo.assign(presample = mdexo.genoplust.str.cat(mdexo.compartment, sep= "_"))
mdexo = mdexo.assign(sample = mdexo.presample.str.cat(mdexo.replicate_bio, sep="-"))
mdexo = mdexo.drop(columns = ['presample'])
mdexo.to_csv("metadata_exo.csv", header=True)
# join in a single metadata for entire dataset:
bigmetadata = pd.concat([metadata_endo,mdexo], axis=0)

bigmetadata.to_csv("metadataLDH_endo_exo.csv", header=True, index = False)

# now really replace missing samples with NAs in exometabo abundances
print(missingsamples)
for sa in missingsamples:
    exo_true[sa] = np.nan
# reorder the columns to be identical of those in metadata, using name_old
#print(exo_true.head)

# style : AB_T24_med-1
desstyle = list()
for i in exo_true.columns:
    geno = i.split("_")[0]
    time = i.split("_")[1].split("-")[0]
    biorep = i.split("-")[-1]
    desstyle.append(geno + "_" + time + "_med" + "-" + biorep)
exo_true.columns = desstyle

exo_true.to_csv("exo_abundance_total.tsv", sep ="\t", header=True)


# row_data
rowdata_ex = pd.DataFrame({'met_orig': exo_true.index.tolist(),
                           'compartment': 'med'})
rowdata_ex.to_csv("rowdata_exo.tsv", header=True, sep='\t', index=False)
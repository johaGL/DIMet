"""
detecting aberrant isotopologues %
note: it is written as proportion, 1.0 == 100%
first multiply all by 100
 a)  print isotopologues which sum > 100
 b) print negative! isotopologues
"""
import os
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--opt", help = "minidemo | sup100 | negall  | supsum100  ", )
args= parser.parse_args()


os.chdir(os.path.expanduser("~/myelin_metabo/"))
co = "cell"

idf = pd.read_csv(f"tmp/correctedIsotopologues_myelin_{co}.tsv", sep='\t', index_col=0)
m = pd.read_csv("data/metadata_myelin.csv")
m = m.loc[m.short_comp == co,:]

idf = idf * 100
auxiliary = pd.DataFrame({'isotopologue': idf.index,
                          'metabolite' : ['' for i in range(len(idf.index))]})
auxiliary.metabolite = auxiliary.isotopologue.str.split("_m+", regex=False).str[0]
mets = auxiliary.metabolite.unique()


"""a ) isotopologues or having negative       """
if args.opt == "minidemo"  :
    print("\n    ABERRANT negative values:")
    for s in ['Vehicle_cell_T0h-1']: # idf.columns:
        print("\n    ======>  ", s)
        tmp = auxiliary.copy()
        tmp[s] = idf[s].tolist()
        for met in mets:
            minit = tmp.loc[tmp.metabolite == met,:]
            negs = minit.loc[minit[s] < 0,:]
            if negs.shape[0] > 0 :
                print(negs)

    print("endded aberrant negative values, one single sample")


if args.opt == "negall"  :
    print("\n    ABERRANT negative values:")
    for s in idf.columns: # idf.columns:
        print("\n    ======>  ", s)
        tmp = auxiliary.copy()
        tmp[s] = idf[s].tolist()
        for met in mets:
            minit = tmp.loc[tmp.metabolite == met,:]
            negs = minit.loc[minit[s] < 0,:]
            if negs.shape[0] > 0 :
                print(negs)


if args.opt == "sup100":
    print("\n    ABERRANT superior to 100 values:")
    for s in idf.columns: # idf.columns:
        print("\n    ======>  ", s)
        tmp = auxiliary.copy()
        tmp[s] = idf[s].tolist()
        for met in mets:
            minit = tmp.loc[tmp.metabolite == met,:]
            p100 = minit.loc[minit[s] > 101,:]
            if p100.shape[0] > 0 :
                print(p100)



if args.opt == "supsum100":
    """b) isotopologues summing up more than 100  %"""
    for s in idf.columns : # idf.columns:
        print("\n    ======>  ", s)
        tmp = auxiliary.copy()
        tmp[s] = idf[s].tolist()
        #print(tmp)
        sumdf = tmp.groupby('metabolite').sum()
        print("    superior to 100:")
        aberrant_over = sumdf[sumdf[s]>100]
        print("    InFerior! to 100:")
        aberrant_over = sumdf[sumdf[s] > 100]
        print(aberrant_over)



# END
# note: no need to use directly the excel cause already splitted and colnames ok by formatter_be.py

def noneed():
    os.chdir(os.path.expanduser("~/myelin_metabo/smpls_raw"))
    namesuffix = "myelin"
    finame = "MCF001289_8510_metabolomics_results.xlsx"

    rows_todrop = ["Blank01", "Blank02", "Blank03", "LOD",
                   "Cell extracts", "Medium"]

    names_compartments = {"Cell_extracts" : "cell",
                          "Medium" : "med" }

    namesheet =  'correctedIsotopologues' # this time lower case first char !
    idf = pd.read_excel(finame, sheet_name = namesheet,
                               engine='openpyxl',
                               header = 0,
                               index_col = 0)

    todrop_here = set(idf.index.to_list()).intersection(rows_todrop)
    #print(todrop_here)
    idf = idf.drop(todrop_here)
    idf = idf.dropna(axis=0, how="all") # if all cells in row == NaN
    idf = idf.loc[:, ~idf.columns.str.contains("^Unnamed")]
    idf = idf.dropna(axis=1, how="all") # if all ccell in col == NaN
    print("transposing the dataframe, so samples are columns")
    idf = idf.T
    idf = idf * 100







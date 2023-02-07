import os
import sys

import numpy as np
import pandas as pd
import numpy as np
#infile = sys.argv[1]
#outfile = sys.argv[2]


lo = os.path.expanduser("~/myelin_metabo/results/tables/extended/")



infile = "cell_TOTAL_Myelin_T0h_Pellets_T0h_Tt.tsv"

df = pd.read_csv(lo + "TOTAL/" + infile, sep='\t' )
tmp = df.copy()
tmp['abslfc'] = np.absolute(tmp['log2FC'])
print(tmp[['distance/span', 'log2FC', 'abslfc']])
dfgood = tmp.loc[(tmp['distance/span'] >= 0.2) & (tmp['abslfc'] >= 0.3),:]
try:
    dfgood = dfgood[['metabolite', 'distance/span', 'ratio', 'log2FC', 'p-value', 'padj']]
except:
    dfgood = dfgood[['metabolite', 'distance/span', 'ratio', 'log2FC', 'pvalue', 'padj']]

print(dfgood.shape)

odi = os.path.expanduser("~/myelin_metabo/results/tables/interesting/")

if not os.path.exists(odi):
    os.makedirs(odi)

outfile = "cell_TOTAL_Myelin_T0h_Pellets_T0h_Tt_good.tsv"
dfgood.to_csv(odi + outfile, sep='\t', index=False)

##### again

lo = os.path.expanduser("~/myelin_metabo/results/tables/extended/")


subdirs = ["isos", "TOTAL"]

sx = ["_disfit", "_MW", "_Tt"]
odi = os.path.expanduser("~/myelin_metabo/results/tables/interesting/")

if not os.path.exists(odi):
    os.makedirs(odi)

for SUFIX in sx:
    for subdi in subdirs:
        files = os.listdir(os.path.expanduser("~/myelin_metabo/results/tables/extended/") + subdi)
        for fi in files:
            if fi.split(".")[0].endswith(SUFIX):
                df = pd.read_csv(lo + subdi + "/" +  fi, sep='\t' )
                tmp = df.copy()
                tmp['abslfc'] = np.absolute(tmp['log2FC'])
                print(tmp[['distance/span', 'log2FC', 'abslfc']])
                dfgood = tmp.loc[(tmp['distance/span'] >= 0.2) & (tmp['abslfc'] >= 0.2),:]
                try:
                    dfgood = dfgood[['metabolite', 'distance/span', 'ratio', 'log2FC', 'p-value', 'padj']]
                except:
                    dfgood = dfgood[['metabolite', 'distance/span', 'ratio', 'log2FC', 'pvalue', 'padj']]

                print(dfgood.shape)
                if dfgood.shape[0] >= 1 :
                    outfile = fi.split(".")[0] + "_filtered.tsv"
                    dfgood.to_csv(odi + outfile, sep='\t', index=False)
                else:
                    print(dfgood.shape)

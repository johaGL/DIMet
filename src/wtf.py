import numpy

tsv = "/home/johanna/padj_addmethod.tsv"
import pandas as pd
import numpy as np
import statsmodels.stats.multitest as smm
import statsmodels
df = pd.read_csv(tsv,sep='\t', header=0, index_col=None)
print(df.head())
import matplotlib.pyplot as plt

x = [4000,4000,4000]
y = [10585,11000,10900]
from scipy import stats
yes = stats.kruskal(x, y)
print(yes)


import seaborn as sns
sns.histplot(df, x="pvalue")
plt.show()

df["pvalue_round"] = df["pvalue"].round(decimals=2)
df = df.sort_values(by="pvalue")

order = np.array(df['pvalue']).argsort()
ranks = order.argsort()

df['rank'] = ranks
for method in ['bonferroni', 'sidak',
               'holm-sidak',  'holm', 'simes-hochberg',
               'hommel',  'fdr_bh', 'fdr_by',
                'fdr_tsbh', 'fdr_tsbky']:
    print(method)
    (sgs, corrP, _, _) = smm.multipletests(df["pvalue_round"], alpha=float(0.05),
                                           method=method)
    df[method] = corrP

"""
def give_ranks(values):
    original_df = pd.DataFrame({'value': values.unique()})
    uniq_df = pd.DataFrame({'value': values.unique()})
    uniq_df = uniq_df.sort_values(by="pvalue", ascending=True)
    uniq_df['rank'] = range(uniq_df.shape[0])
    original_df = pd.merge(original_df, uniq_df, on='value',  how='left')

foo = give_ranks(np.array([9,12,4,4]))
"""
ohno = df['rank'].div(df.shape[0])
print(ohno)
df['mymethod'] = ohno * 0.95


print(statsmodels.__version__)

#df.to_csv(tsv,sep='\t' )

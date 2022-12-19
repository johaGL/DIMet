import numpy as np
import pandas as pd
import statsmodels.stats.multitest as smm

def compute_p_adjusted(df: pd.DataFrame, correction_method: str) -> pd.DataFrame:
    rej, pval_corr = smm.multipletests(df['pvalue'].values, alpha=float('0.05'), method=correction_method)[:2]
    df['padj'] = pval_corr
    return df

print("using a dummy array of values")
#kn = [0.04, 0.004, 0.1, 0.5, 0.9]
#kn = [0.24, 0.5, 0.5, 1, 0.0001]
kn = [1,  0.290034529038187, 0.320932951534461, 0.2065564914483, 0.174633900731984]


rej_bool, pval_arr = smm.multipletests(kn, alpha=0.05, method="fdr_bh")[:2]
print(kn)
print(pval_arr)
# or using the function
df = pd.DataFrame({"rownb": [i for i in range(len(kn))],
                   "pvalue" : kn})
df = compute_p_adjusted(df, "fdr_bh")
print(df)

### ####
print("\nUsing some real values")
prekn = [0.586519014724443 , 0.609126394687426 , 0.490666135722899, 0.425655265804451, \
    0.239121611248683, 0.516810415015153, 1, 0.585822774747661 ]

kn = np.array(prekn).round(2)
print(kn)
rej_bool, pval_arr = smm.multipletests(kn, alpha=0.05, method="fdr_bh")[:2]
print(pval_arr)
print("using 'is_sorted'")
# now function but ordered
def compute_p_adjustedV2(df: pd.DataFrame, correction_method: str, is_sorted: bool) -> pd.DataFrame:
    rej, pval_corr = smm.multipletests(df['pvalue'].values,
                                       alpha=float('0.05'),
                                       method=correction_method,
                                       is_sorted=is_sorted)[:2]
    df['padj'] = pval_corr
    return df
df = pd.DataFrame({"rownb": [i for i in range(len(kn))],
                   "pvalue" : kn})
df2 = df.sort_values(by='pvalue', ascending=True)
df2 = compute_p_adjustedV2(df2, "fdr_bh", is_sorted=True)
print(df)
print("useless is_sorted\n")

print("using fdrcorrection function")

rej, pval_corr = smm.fdrcorrection(df['pvalue'].values, method='p',alpha= 0.05)
print(rej)
print(pval_corr)


## ## ##
print("\nUsing several method in pure multipletests function")
print("pvalues :", df['pvalue'])
for imethod in ['bonferroni', 'sidak', 'homl-sidak', 'holm','hommel',
               'fdr_bh', 'fdr_by', 'fdr_tsbh', 'fdr_tsbky']:
    try:
        print(imethod)
        rejtmp, pvaltmp = smm.multipletests(df['pvalue'].values, alpha=float('0.05'), method=imethod)[:2]
        print(pvaltmp)
    except:
        pass

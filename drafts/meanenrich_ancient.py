


# def calculate_meanEnri(
#                 datadi,
#                 tableAbund,
#                 tableIC,
#                 metadata,
#                 names_compartments,
#                 namesuffix,
#                 odir ):
#     """
#     saves a tsv file :
#        meanEnrich_ :  weighted mean as in Castano-cerezo, Bellvert et al, 2019.
#     """
#     if not os.path.exists(odir):
#         os.makedirs(odir)
#
#     for co in names_compartments.values():
#         abun = pd.read_csv(
#             datadi + tableAbund + "_" + namesuffix + "_" + co + ".tsv",
#             sep="\t",
#             index_col=0,
#         )
#         # note that pandas automatically transform any 99.9% in decimal 0.999
#
#         propc = pd.read_csv(
#             datadi + tableIC + "_" + namesuffix + "_" + co + ".tsv",
#             sep="\t",
#             index_col=0,
#         )
#         dicos = dict()
#         dicos[co] = {}
#         dicos[co]["metadata"] = metadata.loc[metadata.short_comp == co]
#         selecols = dicos[co]["metadata"]["sample"]
#         # order  columns (samples) , same order in both df:
#         dicos[co][tableAbund] = abun[selecols]
#         dicos[co][tableIC] = propc[selecols]
#         prop_df = dicos[co][tableIC]
#         prop_df = correction_prop_df(prop_df)
#         isotosrowdata = yieldrowdataB(prop_df)  # this funciton is from fun_fm
#
#         isotosrowdata = isotosrowdata.assign(
#                             coefficient=isotosrowdata["m+x"].str.replace("m+", "", regex=False))
#         isotosrowdata.coefficient = isotosrowdata.coefficient.astype('int')
#
#         outdf = dicos[co][tableAbund].copy()
#         proxy_correctabu = prop_df.copy()
#         metabolites = outdf.index
#
#         for met_i in metabolites:
#             isos_rowdata_i = isotosrowdata.loc[isotosrowdata["metabolite"] == met_i, :]
#             prop_df_i = prop_df.loc[isos_rowdata_i.isotopolgFull.tolist(), :]
#             abu_i = dicos[co][tableAbund].loc[met_i, :]
#             abus_isotops_i = prop_df_i.multiply(abu_i.T) / 100
#             # take advantage of this corrected  abus_isotops_i, to save into proxy_...
#             proxy_correctabu.loc[abus_isotops_i.index, :] = abus_isotops_i
#
#             # prepare to apply formula (sum(coeffs*abus_isotops_i)) / n
#             # from Castano-cerezo, Bellvert et al, 2019.
#             nmax = max(isos_rowdata_i.coefficient)
#             abu_coef_prod_df = abus_isotops_i.multiply(isos_rowdata_i.coefficient.tolist(),
#                                                        axis="index")
#             abu_coef_sum_df = abu_coef_prod_df.sum(axis="index")
#             meanenrich_i = abu_coef_sum_df / nmax
#             meanenrich_i = meanenrich_i[selecols]
#             outdf.loc[met_i, ] = meanenrich_i
#
#
#
#         ofile = odir + "meanEnrich_" + namesuffix + "_" + co + ".tsv"
#         outdf.to_csv(ofile, sep='\t', header =True, index=True)
#         print("\nSaved mean enrichment to tmp/. Compartment:", co, "\n")
#     return 0

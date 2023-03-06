
# this is just a draft, to be integrated into
# the downstream analysis
def compute_uptake_secretion(individual_dic: dict, confidic: dict, metadata: pd.DataFrame):
    # following the recommendations of VIB 'Medium Analysis' : (medium + cells)/ medium
    def recognize_cell_med_words(words: list):
        cellwords = list()
        medwords = list()
        for w in words:
            if ("cell" in w.lower()) or ("endo" in w.lower()) :
                cellwords.append(w)
            elif ("med" in w.lower()) or ("exo" in w.lower()) or \
                    ("supernatant" in w.lower() or ("plasma" in w.lower())):
                medwords.append(w)
        return cellwords, medwords

    if len(confidic['compartments']) == 2:
        cellwords, medwords = recognize_cell_med_words(confidic['compartments'])
        if (len(cellwords) == 1) and (len(medwords) == 1):
            cellword = cellwords[0]
            medword = medwords[0]
            if (cellword != "") and (medword != ""):
                abu_cell_df = individual_dic[confidic['name_abundance']][cellword]
                abu_med_df = individual_dic[confidic['name_abundance']][medword]
                inmed_not_incell = set(abu_med_df.columns) - set(abu_cell_df.columns)
                abu_cell_df[list(inmed_not_incell)] = np.nan # (was deleted because zero all or all <LOD)
                # per condition and per time, create new rows with averages over replicates (use metadata)
                # indexes must match !
                # abu_cell_df_avg =  yield_avg_bycond_bytime(abu_cell_df)
                # abu_med_df_avg = yield_avg_bycond_bytime(abu_med_df)
                #uptake_secretion_df = abu_med_df_avg.sum(abu_cell_df_avg[abu_med_df_avg.columns])
                #uptake_secretion_df = uptake_secretion_df.div(abu_med_df_avg)
                return uptake_secretion_df
            else:
                print("unable to compute uptake_secretion table, check compartments")
                return None
        else:
            print("2 or more ambiguous compartments names, unable to compute uptake_secretion table")
            return None

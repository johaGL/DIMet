
# this is just a draft, to be integrated into
# the downstream analysis
def compute_uptake_secretion(individual_dic: dict, confidic: dict, metadata: pd.DataFrame):
    # following the recommendations of VIB 'Medium Analysis' : (medium + cells)/ medium
    def recognize_cell_med_words(words: list):
        cellword = ""
        medword = ""
        for w in words:
            if ("cell" in w.lower()) or ("endo" in w.lower()) :
                cellword = w
            elif ("med" in w.lower()) or ("exo" in w.lower()) or \
                    ("supernatant" in w.lower() or ("plasma" in w.lower())):
                medword = w
        return cellword, medword

    if len(confidic['compartments']) == 2:
        cellword, medword = recognize_cell_med_words(confidic['compartments'])
        if (cellword != "") and (medword != ""):
            abu_cell_df = individual_dic[confidic['name_abundance']][cellword]
            abu_med_df = individual_dic[confidic['name_abundance']][medword]
            inmed_not_incell = set(abu_med_df.columns) - set(abu_cell_df.columns)
            abu_cell_df[list(inmed_not_incell)] = 0 # equilibrate (was deleted because zero all or all <LOD)
            # per condition and per time, create new rows with averages over replicates (use metadata)
            # indexes must match !
            #uptake_secretion_df = abu_med_df_avg.sum(abu_cell_df_avg[abu_med_df_avg.columns])
            #uptake_secretion_df = uptake_secretion_df.div(abu_med_df_avg)
            return uptake_secretion_df
        else:
            print("unable to compute uptake_secretion, check compartments")
            return None
    else:
        print("too many compartments for computing uptake_secretion (max 2)")
        return None
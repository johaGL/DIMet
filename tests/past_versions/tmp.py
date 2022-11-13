# interesting functions but no longer used :

def mean_sd_D(newdf, ametadasel):
    """
    returns dico (oD) with keys: metabolite,
    mu|condition_timepoint ,
    sd|condition_timepoint

    All in same order as metabolites
    example : sd|L-Cycloserine_T0
    """
    oD = dict()
    condis = set(ametadasel["condition"])
    times = set(ametadasel["timepoint"])
    for cd in condis:
        for t in times:
            sasas = ametadasel.loc[
                (ametadasel["condition"] == cd) & (ametadasel["timepoint"] == t)
            ]

            tmp = newdf[sasas["sample"]]
            means_here = []
            sd_here = []
            for i, row in tmp.iterrows():  # iterates metabolites
                means_here.append(np.mean(row))
                sd_here.append(np.std(row))
            oD["metabolite"] = tmp.index
            oD[f"sd|{cd}_{t}"] = sd_here
            oD[f"mu|{cd}_{t}"] = means_here
    return oD


def tmpstack(oD):
    oudf = pd.DataFrame(
        data={"metabolite": [], "condition": [], "timepoint": [], "mean": [], "sd": []}
    )
    mets = oD["metabolite"]
    for k in list(oD.keys())[1:]:
        cd, t = k.split("|")[1].split("_")  # this separator "_" by default
        tmp = pd.DataFrame(
            data={
                "metabolite": mets,
                "condition": [cd for i in range(len(mets))],
                "timepoint": [t for i in range(len(mets))],
                "mean": [x for x in oD[f"mu|{cd}_{t}"]],
                "sd": [x for x in oD[f"sd|{cd}_{t}"]],
            }
        )
        oudf = pd.concat([oudf, tmp])
    oudf = oudf.drop_duplicates()
    return oudf
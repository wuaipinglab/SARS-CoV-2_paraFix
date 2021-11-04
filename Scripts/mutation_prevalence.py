import os
import json
import datetime

import pandas as pd
from matplotlib import pyplot as plt


PROTEIN_NAME = "Spike"

SITESMAPPING_FILE = "Data/sitesMapping.csv"
SURVEILLANCE_FILE = "Data/variant_surveillance.tsv"
PERCENTAGE_SUM_PLOT = "Output/percentage_sum.pdf"

print("Load data...")

df = pd.read_csv(SURVEILLANCE_FILE, sep="\t", low_memory=False)
df = df[df["Collection date"].str.len() == 10]
df["Collection date"] = pd.to_datetime(df["Collection date"])
df = df[df["Collection date"] > datetime.datetime(2019, 11, 30)]

background = df["Collection date"].value_counts()
background = background.sort_index()

sitesMapping = pd.read_csv(SITESMAPPING_FILE, index_col=0)
allSites = sitesMapping[sitesMapping["product"] == PROTEIN_NAME]
allSites = allSites["aaPos"].drop_duplicates().values

print("Summarize percentage sum...")

mutSites = []
for aaPos in allSites:
    print(aaPos)
    c_date = df.loc[df["AA Substitutions"].str.contains(
        "Spike_[A-Z]{}".format(aaPos),
        na=False,
        regex=True
    ), "Collection date"].value_counts()
    c_date = c_date.sort_index()
    
    percentage_sum = 0
    for d in background.index:
        if d in c_date.index:
            percentage_sum += c_date[d] / background[d]
    
    mutSites.append({
        "aaPos": aaPos,
        "percentage_sum": percentage_sum
    })
mutSites = pd.DataFrame.from_records(mutSites)

plt.clf()
plt.hist(mutSites["percentage_sum"])
plt.xlabel("percentage_sum")
plt.ylabel("Frequency")
plt.savefig(PERCENTAGE_SUM_PLOT, bbox_inches="tight")
plt.show()

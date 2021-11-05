import os
import json
import datetime

import pandas as pd


PROTEIN_NAME = "Spike"

SITESMAPPING_FILE = "Data/sitesMapping.csv"
SURVEILLANCE_FILE = "Data/variant_surveillance.tsv"
BACKGROUND_NUM_FILE = "Data/background_num.json"
MUTATION_NUM_FILE = "Data/mutation_num.json"

print("Load data...")

df = pd.read_csv(SURVEILLANCE_FILE, sep="\t", low_memory=False)
df = df[df["Collection date"].str.len() == 10]
df["Collection date"] = pd.to_datetime(df["Collection date"])
df = df[df["Collection date"] > datetime.datetime(2019, 11, 30)]

background = df["Collection date"].value_counts()
background = background.sort_index()


bgNum = {}
for d in background.index:
    d_str = d.strftime("%Y-%m-%d")
    bgNum[d_str] = int(background[d])
    
with open(BACKGROUND_NUM_FILE, "w") as f:
    json.dump(bgNum, f)

sitesMapping = pd.read_csv(SITESMAPPING_FILE, index_col=0)
allSites = sitesMapping[sitesMapping["product"] == PROTEIN_NAME]
allSites = allSites["aaPos"].drop_duplicates().values

print("Summarize percentage sum...")

mutation_num = {}
for aaPos in allSites:
    print(aaPos)
    c_date = df.loc[df["AA Substitutions"].str.contains(
        "Spike_[A-Z]{}".format(aaPos),
        na=False,
        regex=True
    ), "Collection date"].value_counts()
    c_date = c_date.sort_index()
    
    percentage_daily = {}
    for d in background.index:
        d_str = d.strftime("%Y-%m-%d")
        if d in c_date.index:
            percentage_daily[d_str] = int(c_date[d])
        else:
            percentage_daily[d_str] = 0
    
    mutation_num[int(aaPos)] = percentage_daily
    
with open(MUTATION_NUM_FILE, "w") as f:
    json.dump(mutation_num, f)

import os
import gc
import json
import datetime
from collections import defaultdict

import pandas as pd


PROTEIN_NAME = "Spike"
# PROTEIN_NAME = "N"

SITESMAPPING_FILE = "Data/sitesMapping.csv"
SURVEILLANCE_FILE = "Data/variant_surveillance.tsv"
BACKGROUND_NUM_FILE = "Data/background_num.json"
MUTATION_NUM_FILE = "Data/mutation_num_" + PROTEIN_NAME + ".json"

print("Load data...")

df = pd.read_csv(SURVEILLANCE_FILE, sep="\t", low_memory=False)
df = df[df["Collection date"].str.len() == 10]
df["Collection date"] = pd.to_datetime(df["Collection date"])
df = df[df["Collection date"] > datetime.datetime(2019, 11, 30)]
df["Continent"] = df["Location"].str.split(" / ").str[0]


globalDates = pd.to_datetime(df["Collection date"].unique())

bgNum = defaultdict(dict) 
for continent, group in df.groupby("Continent"):
    background = group["Collection date"].value_counts()
    background = background.sort_index()

    for d in globalDates:
        d_str = d.strftime("%Y-%m-%d")
        if d in background.index:
            bgNum[continent][d_str] = int(background[d])
        else:
            bgNum[continent][d_str] = 0
    
with open(BACKGROUND_NUM_FILE, "w") as f:
    json.dump(bgNum, f)
    print(BACKGROUND_NUM_FILE, "saved")
    
del bgNum

sitesMapping = pd.read_csv(SITESMAPPING_FILE, index_col=0)
allSites = sitesMapping[sitesMapping["product"] == PROTEIN_NAME]
allSites = allSites["aaPos"].drop_duplicates().values

print("Summarize percentage sum...")

mutation_num = defaultdict(dict)
# for aaPos in allSites[612:615]:
for aaPos in allSites:
    print(aaPos)
    for continent, group in df.groupby("Continent"):
        c_date = group.loc[group["AA Substitutions"].str.contains(
            f"{PROTEIN_NAME}_[A-Z]{aaPos}[A-Z]",
            na=False,
            regex=True
        ), "Collection date"].value_counts()
        c_date = c_date.sort_index()
        
        percentage_daily = {}
        for d in globalDates:
            d_str = d.strftime("%Y-%m-%d")
            if d in c_date.index:
                percentage_daily[d_str] = int(c_date[d])
            else:
                percentage_daily[d_str] = 0
        
        mutation_num[continent][int(aaPos)] = percentage_daily

    gc.collect()
    
with open(MUTATION_NUM_FILE, "w") as f:
    json.dump(mutation_num, f)
    print(MUTATION_NUM_FILE, "saved")

# with open("../test.json", "w") as f:
#     json.dump(mutation_num, f)

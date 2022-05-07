#!/usr/bin/env python

import gc
import datetime

import pandas as pd


PROTEIN_NAME = "Spike"
# PROTEIN_NAME = "N"

SITESMAPPING_FILE = "data/sitesMapping.csv"
SURVEILLANCE_FILE = "data/variant_surveillance.tsv"
BACKGROUND_NUM_FILE = "output/background_num.csv"
MUTATION_NUM_FILE = "output/mutation_num_" + PROTEIN_NAME + ".csv"

print("Load data...")

df: pd.DataFrame = pd.read_csv(SURVEILLANCE_FILE, sep="\t", low_memory=False)
df = df[df["Collection date"].str.len() == 10]
df["Collection date"] = pd.to_datetime(df["Collection date"])
df = df[df["Collection date"] > datetime.datetime(2019, 11, 30)]
df["Continent"] = df["Location"].str.split(" / ").str[0]


globalDates = pd.to_datetime(df["Collection date"].unique())

bgNum = []
group: pd.DataFrame
d: datetime.datetime
for continent, group in df.groupby("Continent"):
    background = group["Collection date"].value_counts()
    background = background.sort_index()

    for d in globalDates:
        d_str = d.strftime("%Y-%m-%d")
        total_num = 0
        if d in background.index:
            total_num = int(background[d])
        bgNum.append({
            "Area": continent,
            "Date": d_str,
            "Total_num": total_num
        })

pd.DataFrame.from_records(bgNum).to_csv(BACKGROUND_NUM_FILE, index=False)
print(BACKGROUND_NUM_FILE, "saved!")

del bgNum

sitesMapping = pd.read_csv(SITESMAPPING_FILE, index_col=0)
allSites = sitesMapping[sitesMapping["product"] == PROTEIN_NAME]
allSites = allSites["aaPos"].drop_duplicates().values

print("Summarize percentage sum...")

mutation_num = []
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

        for d in globalDates:
            num = 0
            d_str = d.strftime("%Y-%m-%d")
            if d in c_date.index:
                num = int(c_date[d])
            mutation_num.append({
                "Area": continent,
                "Date": d_str,
                "Mut_num": num,
                "Site": int(aaPos)
            })
    gc.collect()

pd.DataFrame.from_records(mutation_num).to_csv(MUTATION_NUM_FILE, index=False)
print(MUTATION_NUM_FILE, "saved!")

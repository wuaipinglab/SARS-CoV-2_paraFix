import os
import json

import pandas as pd


SAMPLED_DIR = "Data/sampled_trees_with_MSA/"
NEXTSTRAIN_DIR = "Data/nextstrain_trees_with_MSA/"
SAMPLED_DATE_FILE = "Data/sampled_dates.json"
NEXTSTRAIN_DATE_FILE = "Data/nextstrain_dates.json"

sampled_dates = pd.DataFrame.from_dict({"date": os.listdir(SAMPLED_DIR)})
nextstrain_dates = pd.DataFrame.from_dict({"date": os.listdir(NEXTSTRAIN_DIR)})

sampled_dates["date"] = pd.to_datetime(sampled_dates["date"])
nextstrain_dates["date"] = pd.to_datetime(nextstrain_dates["date"])

startDate = max(min(sampled_dates["date"]), min(nextstrain_dates["date"]))
endDate = min(max(sampled_dates["date"]), max(nextstrain_dates["date"]))

reference_dates = set()
adjacent_dates = set()

for d1, _ in sampled_dates.groupby("date"):
    d2 = nextstrain_dates.loc[abs(nextstrain_dates["date"] - d1).idxmin(), "date"]
    if (d1 <= endDate or d2 <= endDate) and (d1 >= startDate or d2 >= startDate):
        reference_dates.add(d1.strftime("%Y-%m-%d"))
        adjacent_dates.add(d2.strftime("%Y-%m-%d"))


with open(SAMPLED_DATE_FILE, "w") as f:
    json.dump(list(reference_dates), f)

with open(NEXTSTRAIN_DATE_FILE, "w") as f:
    json.dump(list(adjacent_dates), f)


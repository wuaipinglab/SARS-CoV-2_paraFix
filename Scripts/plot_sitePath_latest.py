import os
import json
import datetime
from collections import Counter

import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.dates import DateFormatter


PROTEIN_NAME = "Spike"
# PROTEIN_NAME = "N"

BACKGROUND_NUM_FILE = "Data/background_num.json"
MUTATION_NUM_FILE = "Data/mutation_num_" + PROTEIN_NAME + ".json"
SITES_PREVALENCE_FILE = "Data/sitesPrevalence_" + PROTEIN_NAME + ".json"
PREVALENCE_INTO_FILE = "Data/prevalenceInfo_" + PROTEIN_NAME + ".csv"

PARAFIXSITES_FILE = "Data/nextstrain_sitePath_results.csv"
PARAFIXMONTHLY_PLOT = "Output/nextstrain_sitePath_results.pdf"
DATES_FILE = "Data/nextstrain_dates.json"

# PARAFIXSITES_FILE = "Data/sampled_sitePath_results.csv"
# PARAFIXMONTHLY_PLOT = "Output/sampled_sitePath_results.pdf"
# DATES_FILE = "Data/sampled_dates.json"


sitesMapping = pd.read_csv("Data/sitesMapping.csv", index_col=0)
siteLabel = sitesMapping[["peptidePos", "gene"]].drop_duplicates()

with open(DATES_FILE) as f:
    allDates = json.load(f)
    allDates = pd.to_datetime(allDates).values

paraFixSites = pd.read_csv(PARAFIXSITES_FILE)
paraFixSites["date"] = pd.to_datetime(paraFixSites["date"])
paraFixSites = paraFixSites[paraFixSites["type"] == "paraFix"]
paraFixSites = paraFixSites[paraFixSites["date"].isin(allDates)]

plt.rcParams.update({'font.size': 20, 'font.weight': 'bold'})
fontsize = 20

typeColors = {
    "fixation": "gold",
    "paraFix": "red",
    "parallel": "blue"
}

fig, axes = plt.subplots(
    nrows=2,
    ncols=1,
    sharex=True,
    figsize=(30, 8),
    gridspec_kw={ "height_ratios": [1, 10], "hspace": 0 }
)


ax = axes[0]
for name, group in siteLabel.groupby("gene", sort=False):
    ax.fill_between(group["peptidePos"], 0, 1, label=name)
ax.axis("off")
ax.legend(bbox_to_anchor=(1.1, 1.1))


ax = axes[1]
for mutType, group in paraFixSites.groupby("type"):
    ax.scatter(
        group["site"],
        group["date"],
        s=100,
        alpha=0.4,
        c=typeColors[mutType],
        label=mutType
    )
ax.set_xlim([-100, 9900])
# ax.set_xticks(paraFixSites["site"])
ax.set_xticks([])
ax.set_xticklabels([])
ax.yaxis.set_major_formatter(DateFormatter('%b %Y'))
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)

plt.savefig(PARAFIXMONTHLY_PLOT, bbox_inches="tight")
plt.show()

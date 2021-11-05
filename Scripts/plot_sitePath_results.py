import os
import json
import datetime
from collections import Counter

import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.dates import DateFormatter


PROTEIN_NAME = "Spike"

SITESMAPPING_FILE = "Data/sitesMapping.csv"
BACKGROUND_NUM_FILE = "Data/background_num.json"
MUTATION_NUM_FILE = "Data/mutation_num.json"

# PARAFIXSITES_FILE = "Data/nextstrain_sitePath_results.csv"
# PARAFIXMONTHLY_PLOT = "Output/nextstrain_detections.pdf"
# PREVALENCE_PLOT = "Output/nextstrain_first_detection.pdf"
# DATES_FILE = "Data/nextstrain_dates.json"

PARAFIXSITES_FILE = "Data/sampled_sitePath_results.csv"
PARAFIXMONTHLY_PLOT = "Output/sampled_detections.pdf"
PREVALENCE_PLOT = "Output/sampled_first_detection.pdf"
DATES_FILE = "Data/sampled_dates.json"


sitesMapping = pd.read_csv(SITESMAPPING_FILE, index_col=0)
siteLabel = sitesMapping[["peptidePos", "gene"]].drop_duplicates()

with open(DATES_FILE) as f:
    allDates = json.load(f)
    allDates = pd.to_datetime(allDates).values

paraFixSites = pd.read_csv(PARAFIXSITES_FILE)
paraFixSites["date"] = pd.to_datetime(paraFixSites["date"])
# paraFixSites.loc[paraFixSites["type"] == "paraFix", "type"] = "parallel"
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
for name, group in siteLabel.groupby("gene"):
    ax.fill_between(group["peptidePos"], 0, 1)
ax.axis("off")

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
ax.legend(bbox_to_anchor=(1, 1))

plt.savefig(PARAFIXMONTHLY_PLOT, bbox_inches="tight")
plt.show()

# The mutation trend plot

selectedSites = paraFixSites[
    # (paraFixSites["type"] == "paraFix") &
    (paraFixSites["gene"] == PROTEIN_NAME)
]

paraFixSitesDate = {
    aaPos: sorted(group["date"].values)
    for aaPos, group in selectedSites.groupby("aaPos")
}

# sortedSites = sorted(
#     paraFixSitesDate,
#     key=lambda site: len(paraFixSitesDate[site]),
#     reverse=True
# )

sortedSites = sorted(
    paraFixSitesDate,
    key=lambda site: paraFixSitesDate[site][0],
    reverse=True
)

with open(BACKGROUND_NUM_FILE) as f:
    background = json.load(f)
    background = pd.DataFrame(background.items(), columns=["date", "num"])
    background["date"] = pd.to_datetime(background["date"])
    background = background.set_index("date")
    background = background.sort_index()

with open(MUTATION_NUM_FILE) as f:
    mutation_num = json.load(f)

mutSites = {}
for aaPos in sortedSites[:10]:
    daily_num = mutation_num[str(aaPos)]
    daily_num = pd.DataFrame(daily_num.items(), columns=["date", "num"])
    daily_num["date"] = pd.to_datetime(daily_num["date"])
    daily_num = daily_num.set_index("date")
    daily_num = daily_num.sort_index()
    mutSites[aaPos] = daily_num

x_pos = background.index.values
x_pos.sort()

nrows = 2
ncols = 5

fig, axes = plt.subplots(
    nrows=nrows,
    ncols=ncols,
    sharex=True,
    sharey=True,
    figsize = (3 * ncols * 2.4, 24)
)
axes2 = []
for ax in axes:
    axes2.extend(ax)

n = 0
bgNum = background["num"].values
for site, meta in mutSites.items():
    mutNum = [meta.loc[d, "num"] if d in meta.index else 0 for d in x_pos]

    ax = axes2[n]
    if n >= nrows * ncols:
        break
    n += 1
    
    ax.fill_between(x_pos, 0, bgNum, label = "Backgroup", facecolor='#AFDAE8')
    ax.fill_between(x_pos, 0, mutNum, label = "mutation", facecolor='#F7E15F')
    ax.vlines(
        x=paraFixSitesDate[site][0],
        colors="red",
        ymin=0,
        ymax=max(bgNum),
        linewidth=3
    )
    ax.tick_params(axis='x', labelrotation=60)
    ax.set_xlim([x_pos[1], x_pos[-1]])
    ax.xaxis.set_major_formatter(DateFormatter('%b %Y'))
    ax.set_xlabel(site, fontsize=fontsize, fontweight='bold')

plt.savefig(PREVALENCE_PLOT, bbox_inches="tight")
plt.show()

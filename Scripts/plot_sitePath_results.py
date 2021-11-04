import os
import json
import datetime
from collections import Counter

import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.dates import DateFormatter


PROTEIN_NAME = "Spike"

SITESMAPPING_FILE = "Data/sitesMapping.csv"
SURVEILLANCE_FILE = "Data/variant_surveillance.tsv"

# PARAFIXSITES_FILE = "Data/sampled_sitePath_results.csv"
# PARAFIXMONTHLY_PLOT = "Output/paraFixedMonthly_sampled.pdf"
# PARALLELSITES_PLOT = "Output/parallelSites_sampled_recent.pdf"

PARAFIXSITES_FILE = "Data/nextstrain_sitePath_results.csv"
PARAFIXMONTHLY_PLOT = "Output/paraFixedMonthly_nextstrain.pdf"
PARALLELSITES_PLOT = "Output/parallelSites_nextstrain.pdf"


sitesMapping = pd.read_csv(SITESMAPPING_FILE, index_col=0)
siteLabel = sitesMapping[["peptidePos", "gene"]].drop_duplicates()

df = pd.read_csv("Data/variant_surveillance.tsv", sep="\t", low_memory=False)
df = df[df["Collection date"].str.len() == 10]
df["Collection date"] = pd.to_datetime(df["Collection date"])
df = df[df["Collection date"] > datetime.datetime(2019, 11, 30)]

background = df["Collection date"].value_counts()
background = background.sort_index()

paraFixSites = pd.read_csv(PARAFIXSITES_FILE)
paraFixSites["date"] = pd.to_datetime(paraFixSites["date"])
# paraFixSites.loc[paraFixSites["type"] == "paraFix", "type"] = "parallel"


plt.rcParams.update({'font.size': 20, 'font.weight': 'bold'})
fontsize = 20

typeColors = {
    "fixation": "gold",
    "paraFix": "red"
}

y_pos = paraFixSites["date"].unique()
y_pos.sort()

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
    if mutType == "parallel":
        continue
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
    (paraFixSites["type"] == "paraFix") &
    (paraFixSites["gene"] == PROTEIN_NAME)
]

paraFixSitesDate = {
    aaPos: sorted(group["date"].values)
    for aaPos, group in selectedSites.groupby("aaPos")
}

sortedSites = sorted(
    paraFixSitesDate,
    key=lambda site: len(paraFixSitesDate[site]),
    reverse=True
)

# sortedSites = sorted(
#     paraFixSitesDate,
#     key=lambda site: paraFixSitesDate[site][0],
#     reverse=True
# )

mutSites = {}
for aaPos in sortedSites[:10]:
    c_date = df.loc[df["AA Substitutions"].str.contains(
        "Spike_[A-Z]{}".format(aaPos),
        na=False,
        regex=True
    ), "Collection date"]
    c_date = c_date.value_counts()
    mutSites[aaPos] = c_date.sort_index()


x_pos = background.index.unique()
x_pos = x_pos.sort_values()

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
for site, meta in mutSites.items():
    bgNum = background.values
    mutNum = [meta[d] if d in meta.index else 0 for d in x_pos]

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

plt.savefig(PARALLELSITES_PLOT, bbox_inches="tight")
plt.show()

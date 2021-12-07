import os
import json
import datetime
from collections import Counter

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.dates import DateFormatter
from matplotlib_venn import venn2, venn2_circles
import matplotlib.colors as mcolors


PROTEIN_NAME = "Spike"

BACKGROUND_NUM_FILE = "Data/background_num.json"
MUTATION_NUM_FILE = "Data/mutation_num_" + PROTEIN_NAME + ".json"
SITES_PREVALENCE_FILE = "Data/sitesPrevalence_" + PROTEIN_NAME + ".json"
PREVALENCE_INTO_FILE = "Data/prevalenceInfo_" + PROTEIN_NAME + ".csv"

PARAFIXSITES_OLD_FILE = "Data/nextstrain_sitePath_results.csv"
PARAFIXSITES_FILE = "Data/latest_sitePath_results.csv"
PARAFIXSITES_COMBINDED_FILE = "Data/nextstrain_sitePath_results_combined.csv"

TARGET_SITES_FILE = "Data/SARS-CoV-2_spike/targetSites.json"
PARAFIXMONTHLY_PLOT = "Output/latest_sitePath_results.pdf"


DAY_TOTAL_THRESHOLD = 20
PREVALENCE_PERCENTAGE_THRESHOLD = 0.5

PROTEIN_DOMAIN = {
    "RBD": [i for i in range(319, 541 + 1)],
    "NTD": [i for i in range(13, 303 + 1)],
    "S2": [i for i in range(686, 1273 + 1)],
    "FP": [i for i in range(816, 855 + 1)],
    "HR1": [i for i in range(920, 970 + 1)],
    "HR2": [i for i in range(1163, 1202 + 1)]
}

DOMAIN_COLOR = {
    "Seq": "#BCBDBD",
    "NTD": "#9BCAC8",
    "RBD": "#C6716B",
    "S2": "#9EC9A1",
    "FP": "#BAD3E1",
    "HR1": "#F0D0CE",
    "HR2": "#F0D0CE"
}

sitesMapping = pd.read_csv("Data/sitesMapping.csv", index_col=0)
siteLabel = sitesMapping[["peptidePos", "aaPos", "gene", "genomePos", "product"]].drop_duplicates()

siteLabel["domain"] = "Seq"
siteLabel = siteLabel[siteLabel["gene"] == "Spike"]
for domain_name, sites in PROTEIN_DOMAIN.items():
    siteLabel.loc[siteLabel["aaPos"].isin(sites), "domain"] = domain_name


with open(BACKGROUND_NUM_FILE) as f:
    background = json.load(f)

with open(MUTATION_NUM_FILE) as f:
    mutation_num = json.load(f)

with open(SITES_PREVALENCE_FILE) as f:
    sites_prevalence = json.load(f)


paraFixSitesOld = pd.read_csv(PARAFIXSITES_OLD_FILE)
paraFixSitesOld["date"] = pd.to_datetime(paraFixSitesOld["date"])

paraFixSites = pd.read_csv(PARAFIXSITES_FILE)
paraFixSites["date"] = pd.to_datetime(paraFixSites["date"])


combinedParaFixSites = pd.concat([paraFixSitesOld, paraFixSites])
combinedParaFixSites = combinedParaFixSites[combinedParaFixSites["gene"] == "Spike"]
combinedParaFixSites.to_csv(PARAFIXSITES_COMBINDED_FILE)

allPara = combinedParaFixSites[combinedParaFixSites["type"] != "fixation"]["aaPos"].unique()
allFix = combinedParaFixSites[combinedParaFixSites["type"] != "parallel"]["aaPos"].unique()
allParaFix = combinedParaFixSites[combinedParaFixSites["type"] == "paraFix"]["aaPos"].value_counts()

targetSites = {
    "fixation": list(int(i) for i in allFix),
    "parallel": list(int(i) for i in allPara),
    "paraFix": list(allParaFix[allParaFix > 10].index)
}

with open(TARGET_SITES_FILE, "w") as f:
    json.dump(targetSites, f)

paraOnly = set(allPara).difference(allFix)
fixOnly = set(allFix).difference(allPara)
paraFix = set(allFix).intersection(allPara)


plt.rcParams.update({'font.size': 40, 'font.weight': 'bold'})
fontsize = 40

typeColors = {
    "fixation": "teal",
    "parallel": "pink",
    "paraFix": "red",
}

plt.figure(figsize=(20, 20))

venn2(
    subsets = [set(allFix), set(allPara)],
    set_labels = ("fixation", "parallel"),
    set_colors=(typeColors["fixation"], typeColors["parallel"]),
    alpha = 0.4
)
plt.savefig("Output/nextstrain_combined_Spike.pdf", bbox_inches="tight")
plt.close()

# combinedParaFixSites = combinedParaFixSites[combinedParaFixSites["type"] != "parallel"]
combinedParaFixSites = combinedParaFixSites.sort_values("date")
combinedParaFixSites.reset_index(inplace=True)


allDates = list(combinedParaFixSites["date"].unique())

fig, axes = plt.subplots(
    nrows=1,
    ncols=2,
    sharey=True,
    figsize=(20, 60),
    gridspec_kw={ "width_ratios": [1, 20], "hspace": 0 }
)

ax = axes[0]
for name, group in siteLabel.groupby("domain", sort=False):
    if name == "Seq":
        ax.fill_between(
            [0.4, 0.6],
            min(group["aaPos"]),
            max(group["aaPos"]),
            # label=name,
            color=DOMAIN_COLOR[name]
        )
    else:
        ax.fill_between(
            [0, 1],
            min(group["aaPos"]),
            max(group["aaPos"]),
            label=name,
            color=DOMAIN_COLOR[name],
            # capstyle="round",
            # joinstyle="round"
        )
ax.axis("off")
ax.legend(bbox_to_anchor=(-0.5, 0.6))
ax.invert_yaxis()

plt.subplots_adjust(wspace = .001)

combinedParaFixSites["typeOrder"] = 0
combinedParaFixSites.loc[combinedParaFixSites["type"] == "fixation", "typeOrder"] = 1
combinedParaFixSites.loc[combinedParaFixSites["type"] == "parallel", "typeOrder"] = 2
combinedParaFixSites = combinedParaFixSites.sort_values("typeOrder")

ax = axes[1]
for mutType, group in combinedParaFixSites.groupby("type", sort=False):
    ax.scatter(
        group["date"],
        group["aaPos"],
        s=100 if mutType != "paraFix" else 150,
        alpha=0.4 if mutType != "paraFix" else 1,
        c=typeColors[mutType],
        label=mutType
    )
for site, info in combinedParaFixSites.groupby("aaPos"):
    info = info.sort_values("date")
    prev_index = -float("inf")
    prev_d = None
    for d in info["date"].values:
        curr_index = allDates.index(d)
        if curr_index - prev_index > 1:
            linestyles="dashed"
        else:
            linestyles="solid"
        if prev_d:
            ax.hlines(
                y=site,
                xmin=prev_d,
                xmax=d,
                colors="grey",
                alpha=0.4,
                linewidth=1,
                linestyles=linestyles
            )
        prev_d = d
        prev_index = curr_index
ax.legend(bbox_to_anchor=(-0.05, 0.7))
ax.yaxis.set_label_position("right")
ax.yaxis.tick_right()
ax.set_yticks([])
ax.set_yticklabels([])
ax.xaxis.set_major_formatter(DateFormatter('%b %Y'))
ax.tick_params(axis='x', labelrotation=60)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)

plt.savefig(PARAFIXMONTHLY_PLOT, bbox_inches="tight")
plt.close()

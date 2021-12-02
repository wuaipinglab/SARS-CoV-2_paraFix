import os
import json
import datetime
from collections import Counter

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.dates import DateFormatter
import matplotlib.colors as mcolors


PROTEIN_NAME = "Spike"

BACKGROUND_NUM_FILE = "Data/background_num.json"
MUTATION_NUM_FILE = "Data/mutation_num_" + PROTEIN_NAME + ".json"
SITES_PREVALENCE_FILE = "Data/sitesPrevalence_" + PROTEIN_NAME + ".json"
PREVALENCE_INTO_FILE = "Data/prevalenceInfo_" + PROTEIN_NAME + ".csv"

PARAFIXSITES_OLD_FILE = "Data/nextstrain_sitePath_results.csv"

PARAFIXSITES_FILE = "Data/latest_sitePath_results.csv"
PARAFIXMONTHLY_PLOT = "Output/latest_sitePath_results.pdf"


DAY_TOTAL_THRESHOLD = 20
PREVALENCE_PERCENTAGE_THRESHOLD = 0.5


sitesMapping = pd.read_csv("Data/sitesMapping.csv", index_col=0)
siteLabel = sitesMapping[["peptidePos", "aaPos", "gene", "genomePos", "product"]].drop_duplicates()


with open(BACKGROUND_NUM_FILE) as f:
    background = json.load(f)

with open(MUTATION_NUM_FILE) as f:
    mutation_num = json.load(f)

with open(SITES_PREVALENCE_FILE) as f:
    sites_prevalence = json.load(f)


paraFixSitesOld = pd.read_csv(PARAFIXSITES_OLD_FILE)
paraFixSitesOld["date"] = pd.to_datetime(paraFixSitesOld["date"])
paraFixSitesOld = paraFixSitesOld.sort_values("date")
paraFixSitesOld = paraFixSitesOld[paraFixSitesOld["type"] == "paraFix"]

paraFixSites = pd.read_csv(PARAFIXSITES_FILE)
paraFixSites["date"] = pd.to_datetime(paraFixSites["date"])
paraFixSites = paraFixSites.sort_values("date")
paraFixSites = paraFixSites[paraFixSites["type"] == "paraFix"]


allRecentSites = paraFixSites[paraFixSites["product"] == PROTEIN_NAME]["aaPos"].unique()
allOldSites = paraFixSitesOld[paraFixSitesOld["product"] == PROTEIN_NAME]["aaPos"].unique()

combinedParaFixSites = pd.concat([paraFixSitesOld, paraFixSites])
combinedParaFixSites = combinedParaFixSites.sort_values("date")
combinedParaFixSites.reset_index(inplace=True)

allDates = list(combinedParaFixSites["date"].unique())

plt.rcParams.update({'font.size': 40, 'font.weight': 'bold'})
fontsize = 40

typeColors = {
    "fixation": "gold",
    "paraFix": "red",
    "parallel": "blue"
}

fig, axes = plt.subplots(
    nrows=1,
    ncols=2,
    sharey=True,
    figsize=(30, 50),
    gridspec_kw={ "width_ratios": [1, 20], "hspace": 0 }
)

ax = axes[0]
for name, group in siteLabel.groupby("gene", sort=False):
    ax.fill_between([0, 1], min(group["peptidePos"]), max(group["peptidePos"]), label=name)
ax.axis("off")
ax.legend(bbox_to_anchor=(-0.5, 0.6))
ax.invert_yaxis()

plt.subplots_adjust(wspace = .001)

ax = axes[1]
for mutType, group in combinedParaFixSites.groupby("type"):
    ax.scatter(
        group["date"],
        group["site"],
        s=100,
        alpha=0.4,
        c=typeColors[mutType],
        label=mutType
    )
    for site, info in group.groupby("site"):
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
                    colors="red",
                    linewidth=2,
                    linestyles=linestyles
                )
            prev_d = d
            prev_index = curr_index
ax.set_ylim([-100, 9900])
# ax.set_xticks(paraFixSites["site"])
ax.set_yticks([])
ax.set_yticklabels([])
ax.xaxis.set_major_formatter(DateFormatter('%b %Y'))
ax.tick_params(axis='x', labelrotation=60)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.invert_yaxis()

plt.savefig(PARAFIXMONTHLY_PLOT, bbox_inches="tight")

"""
domainGroup = {}
for row in sitesMapping[["gene", "product"]].drop_duplicates().itertuples():
    _, geneName, productName = row
    domainGroup[productName] = geneName

domainNames = list(dict.fromkeys(domainGroup.values()))
cm = plt.get_cmap('rainbow')
domainColors = { name: cm(1.*i/len(domainNames)) for i, name in enumerate(domainNames) }


protein_names = []
protein_sites = {}
for name, group in siteLabel.groupby("product", sort=False):
    protein_names.append(name)
    protein_sites[name] = group['aaPos'].values

x_pos = combinedParaFixSites["date"].values
x_pos.sort()

fig, axes = plt.subplots(
    nrows=len(protein_names),
    ncols=2,
    sharex="col",
    sharey="row",
    figsize=(30, 30),
    gridspec_kw={
        "width_ratios": [1, 10],
        "height_ratios": [len(protein_sites[i]) for i in protein_names],
        "hspace": 0
    }
)


domain_labels = []
existing_labels = set()
for (domain_ax, site_ax), name in zip(axes, protein_names):
    labelName = domainGroup[name]
    domain_ax.fill_between([0, 1], 0, max(protein_sites[name]), 
                           color=domainColors[labelName],
                           label=labelName)
    domain_ax.axis("off")
    # domain_ax.spines['top'].set_visible(False)
    # domain_ax.spines['left'].set_visible(False)
    # domain_ax.spines['right'].set_visible(False)
    # domain_ax.spines['bottom'].set_visible(False)
    # domain_ax.set_ylabel(name, rotation=0, labelpad=30, loc="center")
    # domain_ax.yaxis.set_label_coords(-.2, .1)
    domain_ax.invert_yaxis()
    if labelName not in existing_labels:
        domain_labels.append(domain_ax.get_legend_handles_labels())
        existing_labels.add(labelName)
    # ax.legend(bbox_to_anchor=(1.1, 1.5))
    siteInfo = combinedParaFixSites[combinedParaFixSites["product"] == name]
    for mutType, group in siteInfo.groupby("type"):
        site_ax.scatter(
            group["date"],
            group["aaPos"],
            s=100,
            alpha=0.4,
            c=typeColors[mutType],
            label=mutType
        )
    site_ax.set_yticks([])
    site_ax.set_yticklabels([])
    site_ax.set_ylim([min(protein_sites[name]), max(protein_sites[name])])
    site_ax.set_xlim([x_pos[1] - np.timedelta64(15,'D'), x_pos[-1] + np.timedelta64(15,'D')])
    site_ax.xaxis.set_major_formatter(DateFormatter('%b %Y'))
    # site_ax.set_xticks(paraFixSites["date"].unique())
    site_ax.spines['top'].set_visible(False)
    site_ax.spines['left'].set_visible(False)
    site_ax.spines['right'].set_visible(False)
    site_ax.spines['bottom'].set_visible(False)
    site_ax.invert_yaxis()


# ax.set_xticks(paraFixSites["site"])
# ax.set_xticklabels(paraFixSites["aaPos"])
# ax.tick_params(axis='x', labelrotation=90)
handles, labels = [sum(i, []) for i in zip(*domain_labels)]
fig.legend(handles, labels, loc='center left')
plt.savefig(PARAFIXMONTHLY_PLOT, bbox_inches="tight")
# plt.show()
"""



# all_prevalent_sites = set()
# for continent, sites in sites_prevalence.items():
#     for site in sites["prevalent"]:
#         all_prevalent_sites.add(site)
# all_prevalent_sites = list(all_prevalent_sites)

sites_rep_continent = {}
for site in set([*allRecentSites, *allOldSites]):
    rep_continent = None
    earliest_dates = None
    for continent in mutation_num:
        prevalent_dates = []
        daily_num = mutation_num[continent][str(site)]
        for d in daily_num.keys():
            day_total = background[continent][d]
            if day_total > DAY_TOTAL_THRESHOLD:
                percentage = daily_num[d] / day_total
                if percentage > PREVALENCE_PERCENTAGE_THRESHOLD:
                    prevalent_dates.append(d)
        prevalent_dates = pd.to_datetime(prevalent_dates)
        if len(prevalent_dates) == 0:
            continue
        
        if rep_continent is None:
            rep_continent = continent
            earliest_dates = prevalent_dates
        else:
            is_earlier = False
            if (prevalent_dates[0] < earliest_dates[0] and 
                len(prevalent_dates) > len(earliest_dates)):
                   rep_continent = continent
                   earliest_dates = prevalent_dates
    if rep_continent is None:
            rep_continent = continent
    sites_rep_continent[site] = rep_continent

plt.rcParams.update({'font.size': 20, 'font.weight': 'bold'})
fontsize = 20

nrows = 2
ncols = 1

fig, axes = plt.subplots(
    nrows=nrows,
    ncols=ncols,
    sharex=True,
    sharey=False,
    figsize = (10 * ncols, 10 * nrows)
)

axes2 = []
for ax in axes:
    axes2.append(ax)


groupBackgroup = {}
for continent, group in background.items():
    group = pd.DataFrame(group.items(), columns=["date", "num"])
    group["date"] = pd.to_datetime(group["date"])
    group = group.set_index("date")
    group = group.sort_index()
    groupBackgroup[continent] = group

x_pos = groupBackgroup[continent].index.values
x_pos.sort()

newSites = list(set(allRecentSites).difference(set(allOldSites)))

for i in range(nrows * ncols):
    site = newSites[i]
    continent = sites_rep_continent[site]
    
    res_date = paraFixSites.loc[(paraFixSites["aaPos"] == site) &
                                (paraFixSites["gene"] == PROTEIN_NAME), "date"]
    res_date = pd.to_datetime(res_date.values[0]) if len(res_date) else None

    bgNum = groupBackgroup[continent]["num"].values
    daily_num = pd.DataFrame(
        mutation_num[continent][str(site)].items(),
        columns=["date", "num"]
    )
    daily_num["date"] = pd.to_datetime(daily_num["date"])
    daily_num = daily_num.set_index("date")
    daily_num = daily_num.sort_index()
    daily_num = daily_num["num"].values
    
    ax = axes2[i]
    ax.fill_between(x_pos, 0, bgNum, label = "Backgroup", facecolor='#AFDAE8')
    ax.fill_between(x_pos, 0, daily_num, label = "mutation", facecolor='#F7E15F')
    ax.tick_params(axis='x', labelrotation=60)
    ax.set_xlim([x_pos[1], x_pos[-1]])
    ax.xaxis.set_major_formatter(DateFormatter('%b %Y'))
    ax.set_xlabel(f"{site} ({continent})", fontsize=fontsize, fontweight='bold')
    
    if res_date is not None:
    
        ax.vlines(
            x=res_date,
            colors="red",
            ymin=0,
            ymax=max(bgNum),
            linewidth=2,
            # linestyles="dashed"
        )

plt.savefig(f"Output/latest_sitePath_new_detection.pdf", bbox_inches="tight")


plt.rcParams.update({'font.size': 30, 'font.weight': 'bold'})
fontsize = 30

nrows = 4
ncols = 2

fig, axes = plt.subplots(
    nrows=nrows,
    ncols=ncols,
    sharex=True,
    sharey=False,
    figsize = (10 * ncols, 10 * nrows)
)

axes2 = []
for ax in axes:
    axes2.extend(ax)


groupBackgroup = {}
for continent, group in background.items():
    group = pd.DataFrame(group.items(), columns=["date", "num"])
    group["date"] = pd.to_datetime(group["date"])
    group = group.set_index("date")
    group = group.sort_index()
    groupBackgroup[continent] = group

x_pos = groupBackgroup[continent].index.values
x_pos.sort()

oldSites = list(set(allOldSites).difference(set(allRecentSites)))

for i in range(nrows * ncols):
    site = oldSites[i]
    continent = sites_rep_continent[site]
    
    res_date = paraFixSitesOld.loc[(paraFixSitesOld["aaPos"] == site) &
                                   (paraFixSitesOld["gene"] == PROTEIN_NAME), "date"]
    res_date = pd.to_datetime(res_date.values[0]) if len(res_date) else None

    bgNum = groupBackgroup[continent]["num"].values
    daily_num = pd.DataFrame(
        mutation_num[continent][str(site)].items(),
        columns=["date", "num"]
    )
    daily_num["date"] = pd.to_datetime(daily_num["date"])
    daily_num = daily_num.set_index("date")
    daily_num = daily_num.sort_index()
    daily_num = daily_num["num"].values
    
    ax = axes2[i]
    ax.fill_between(x_pos, 0, bgNum, label = "Backgroup", facecolor='#AFDAE8')
    ax.fill_between(x_pos, 0, daily_num, label = "mutation", facecolor='#F7E15F')
    ax.tick_params(axis='x', labelrotation=60)
    ax.set_xlim([x_pos[1], x_pos[-1]])
    ax.xaxis.set_major_formatter(DateFormatter('%b %Y'))
    ax.set_xlabel(f"{site} ({continent})", fontsize=fontsize, fontweight='bold')
    
    if res_date is not None:
    
        ax.vlines(
            x=res_date,
            colors="red",
            ymin=0,
            ymax=max(bgNum),
            linewidth=2,
            # linestyles="dashed"
        )

plt.savefig(f"Output/latest_sitePath_first_detection.pdf", bbox_inches="tight")

import os
import json
from collections import defaultdict

import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.dates import DateFormatter


# PROTEIN_NAME = "Spike"
PROTEIN_NAME = "N"

BACKGROUND_NUM_FILE = "Data/background_num.json"
MUTATION_NUM_FILE = "Data/mutation_num_" + PROTEIN_NAME + ".json"
SITES_PREVALENCE_FILE = "Data/sitesPrevalence_" + PROTEIN_NAME + ".json"
PREVALENCE_INTO_FILE = "Data/prevalenceInfo_" + PROTEIN_NAME + ".csv"
PREVALENCE_PLOT = "Output/prevalence_" + PROTEIN_NAME + ".pdf"
PREVALENCE_LOWER_PLOT = "Output/prevalence_lower_" + PROTEIN_NAME + ".pdf"
PERCENTAGE_SUM_PLOT = "Output/percentage_sum_" + PROTEIN_NAME + ".pdf"

DAY_TOTAL_THRESHOLD = 20
PREVALENCE_PERCENTAGE_THRESHOLD = 0.5
PREVALENCE_THRESHOLD = 30

# PERCENTAGE_SUM_THRESHOLD = 70
# PREVALENCE_THRESHOLD = 50

with open(BACKGROUND_NUM_FILE) as f:
    background = json.load(f)

with open(MUTATION_NUM_FILE) as f:
    mutation_num = json.load(f)


allMutSites = []
for continent, group in mutation_num.items():
    for aaPos, daily_num in group.items():
        percentage_sum = 0
        prevalent_days = 0
        prevalent_dates = []
        for d in daily_num:
            day_total = background[continent][d]
            if day_total > DAY_TOTAL_THRESHOLD:
                percentage = daily_num[d] / day_total
                percentage_sum += percentage
                if percentage > PREVALENCE_PERCENTAGE_THRESHOLD:
                    prevalent_days += 1
                    prevalent_dates.append(d)
        
        prevalent_dates = pd.to_datetime(prevalent_dates)
        # if len(prevalent_dates):
        #     break
        
        allMutSites.append({
            "aaPos": aaPos,
            "percentage_sum": percentage_sum,
            "prevalent_days": prevalent_days,
            "continent": continent
        })
allMutSites = pd.DataFrame.from_records(allMutSites)

allMutSites.to_csv(PREVALENCE_INTO_FILE, index=False)
print(PREVALENCE_INTO_FILE, "saved")

plt.rcParams.update({'font.size': 20, 'font.weight': 'bold'})
fontsize = 20

nrows = 3
ncols = 2

fig, axes = plt.subplots(
    nrows=nrows,
    ncols=ncols,
    sharex=True,
    sharey=True,
    figsize = (5 * ncols * 2.4, 24)
)
axes2 = []
for ax in axes:
    axes2.extend(ax)

i = 0
for continent, mutSites in allMutSites.groupby("continent"):
    
    ax = axes2[i]
    if i >= nrows * ncols:
        break
    i += 1
    
    
    mutSites = mutSites[mutSites["prevalent_days"] < 300]
    (n, bins, patches) = ax.hist(mutSites["prevalent_days"], bins=100)
    # ax.vlines(
    #     x=PERCENTAGE_SUM_THRESHOLD,
    #     colors="red",
    #     ymin=0,
    #     ymax=max(n),
    #     linewidth=2,
    #     linestyles="dashed"
    #     
    # )
    ax.set_xlabel(continent, fontsize=fontsize, fontweight='bold')
    ax.set_yscale('log')
    
plt.savefig(PERCENTAGE_SUM_PLOT, bbox_inches="tight")
print(PERCENTAGE_SUM_PLOT, "saved")
# plt.show()

# The plot for the prevalent sites

sitesPrevalence = {}
groupBackgroup = {}
topMutSites = defaultdict(dict)
topMutSitesOrder = defaultdict(list)
topNum = 20

for continent, group in background.items():
    group = pd.DataFrame(group.items(), columns=["date", "num"])
    group["date"] = pd.to_datetime(group["date"])
    group = group.set_index("date")
    group = group.sort_index()
    groupBackgroup[continent] = group

    mutSites = allMutSites[allMutSites["continent"] == continent]
    mutSites = mutSites.sort_values("prevalent_days", ascending=False)
    selectedSites = mutSites.head(topNum)["aaPos"].values
    
    prevalence_sites = list()
    for aaPos in selectedSites:
        daily_num = mutation_num[continent][str(aaPos)]
        daily_num = pd.DataFrame(daily_num.items(), columns=["date", "num"])
        daily_num["date"] = pd.to_datetime(daily_num["date"])
        daily_num = daily_num.set_index("date")
        daily_num = daily_num.sort_index()
        topMutSites[continent][aaPos] = daily_num
        topMutSitesOrder[continent].append(aaPos)
        (percentage_sum, ) = mutSites.loc[mutSites["aaPos"] == aaPos, "prevalent_days"].values
        if percentage_sum > PREVALENCE_THRESHOLD:
            prevalence_sites.append(int(aaPos))

    
    allSitesPos = pd.to_numeric(mutSites["aaPos"].unique())
    non_prevalent_sites = set(allSitesPos).difference(prevalence_sites)
    
    print(continent, prevalence_sites)
    
    sitesPrevalence[continent] = {
        "prevalent": prevalence_sites,
        "non_prevalent": [int(i) for i in non_prevalent_sites]
    }
    
with open(SITES_PREVALENCE_FILE, "w") as f:
    json.dump(sitesPrevalence, f)
    print(SITES_PREVALENCE_FILE, "saved")


# Plot for prevalent sites

regions = list(background.keys())
nregion = len(regions)
nsites = int(topNum / 2)

x_pos = groupBackgroup[regions[0]].index.values
x_pos.sort()

fig, axes = plt.subplots(
    nregion,
    nsites,
    sharex=True,
    sharey='row',
    figsize = (3 * nsites * 2.4, 24)
)

for j in range(nregion):
    region = regions[j]
    sites_daily_num = topMutSites[region]
    site_names = topMutSitesOrder[region]
    for i in range(nsites):
        site = site_names[i]
        meta = sites_daily_num[site]
        mutNum = [meta.loc[d, "num"] if d in meta.index else 0 for d in x_pos]
        bgNum = groupBackgroup[region]["num"].values

        ax = axes[j][i]
        ax.fill_between(x_pos, 0, bgNum, label = "Backgroup", facecolor='#AFDAE8')
        ax.fill_between(x_pos, 0, mutNum, label = "mutation", facecolor='#F7E15F')
        ax.tick_params(axis='x', labelrotation=60)
        ax.set_xlim([x_pos[1], x_pos[-1]])
        ax.xaxis.set_major_formatter(DateFormatter('%b %Y'))
        ax.set_xlabel(site, fontsize=fontsize, fontweight='bold')
        
        if int(site) in sitesPrevalence[region]["prevalent"]:
            ax.xaxis.label.set_color('red')
        
        # if j == nregion - 1:
        #     ax.set_xlabel(site, fontsize=fontsize, fontweight='bold')
        if i == 0:
            ax.set_ylabel(region, fontsize=fontsize, fontweight='bold')
#             if i == 0:
#                 axes[i][j].legend(loc='upper left')

plt.savefig(PREVALENCE_PLOT, bbox_inches="tight")
print(PREVALENCE_PLOT, "saved")
# plt.show()


fig, axes = plt.subplots(
    nregion,
    nsites,
    sharex=True,
    sharey='row',
    figsize = (3 * nsites * 2.4, 24)
)

for j in range(nregion):
    region = regions[j]
    sites_daily_num = topMutSites[region]
    site_names = topMutSitesOrder[region]
    for i in range(nsites):
        site = site_names[i + nsites]
        meta = sites_daily_num[site]
        mutNum = [meta.loc[d, "num"] if d in meta.index else 0 for d in x_pos]
        bgNum = groupBackgroup[region]["num"].values

        ax = axes[j][i]
        ax.fill_between(x_pos, 0, bgNum, label = "Backgroup", facecolor='#AFDAE8')
        ax.fill_between(x_pos, 0, mutNum, label = "mutation", facecolor='#F7E15F')
        ax.tick_params(axis='x', labelrotation=60)
        ax.set_xlim([x_pos[1], x_pos[-1]])
        ax.xaxis.set_major_formatter(DateFormatter('%b %Y'))
        ax.set_xlabel(site, fontsize=fontsize, fontweight='bold')
        
        if int(site) in sitesPrevalence[region]["prevalent"]:
            ax.xaxis.label.set_color('red')
        
        # if j == nregion - 1:
        #     ax.set_xlabel(site, fontsize=fontsize, fontweight='bold')
        if i == 0:
            ax.set_ylabel(region, fontsize=fontsize, fontweight='bold')
#             if i == 0:
#                 axes[i][j].legend(loc='upper left')

plt.savefig(PREVALENCE_LOWER_PLOT, bbox_inches="tight")
print(PREVALENCE_LOWER_PLOT, "saved")

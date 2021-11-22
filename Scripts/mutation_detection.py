import os
import json
from collections import defaultdict

import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.dates import DateFormatter


PROTEIN_NAME = "Spike"
# PROTEIN_NAME = "N"

SAMPLING_METHOD = "nextstrain"
# SAMPLING_METHOD = "sampled"


DAY_TOTAL_THRESHOLD = 20
PREVALENCE_PERCENTAGE_THRESHOLD = 0.5

BACKGROUND_NUM_FILE = "Data/background_num.json"
MUTATION_NUM_FILE = "Data/mutation_num_" + PROTEIN_NAME + ".json"
SITES_PREVALENCE_FILE = "Data/sitesPrevalence_" + PROTEIN_NAME + ".json"

PARAFIXSITES_FILE = "Data/" + SAMPLING_METHOD + "_sitePath_results.csv"
HOMOPLASYFINDER_RES_FILE = "Data/" + SAMPLING_METHOD + "_homoplasyFinder.csv"
HYPHY_RES_FILE = "Data/" + SAMPLING_METHOD + "_hyphy_results.csv"


# sitesMapping = pd.read_csv("Data/sitesMapping.csv", index_col=0)
# siteLabel = sitesMapping[["peptidePos", "gene"]].drop_duplicates()

with open(BACKGROUND_NUM_FILE) as f:
    background = json.load(f)

with open(MUTATION_NUM_FILE) as f:
    mutation_num = json.load(f)

with open(SITES_PREVALENCE_FILE) as f:
    sites_prevalence = json.load(f)


all_prevalent_sites = set()
for continent, sites in sites_prevalence.items():
    for site in sites["prevalent"]:
        all_prevalent_sites.add(site)

all_prevalent_sites = list(all_prevalent_sites)

sites_rep_continent = {}
for site in all_prevalent_sites:
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
    sites_rep_continent[site] = rep_continent

paraFixSites = pd.read_csv(PARAFIXSITES_FILE)
paraFixSites = paraFixSites[paraFixSites["product"] == PROTEIN_NAME]

hyphySites = pd.read_csv(HYPHY_RES_FILE)
hyphySites = hyphySites[hyphySites["protein"] == PROTEIN_NAME]
hyphySites["aaPos"] = hyphySites["site"]

homoplasySites = pd.read_csv(HOMOPLASYFINDER_RES_FILE)
homoplasySites = homoplasySites[homoplasySites["product"] == PROTEIN_NAME]


plt.rcParams.update({'font.size': 20, 'font.weight': 'bold'})
fontsize = 20

method_names = ["sitePath", "hyphy", "homoplasy"]
method_res = {
    "sitePath": paraFixSites,
    "hyphy": hyphySites,
    "homoplasy": homoplasySites
}

method_color = {
    "sitePath": "purple",
    "hyphy": "green",
    "homoplasy": "orange"
}

nrows = len(method_names)
ncols = 9

for n in range(3):

    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        sharex=True,
        sharey=False,
        figsize = (5 * ncols * 2.4, 24)
    )
    
    groupBackgroup = {}
    for continent, group in background.items():
        group = pd.DataFrame(group.items(), columns=["date", "num"])
        group["date"] = pd.to_datetime(group["date"])
        group = group.set_index("date")
        group = group.sort_index()
        groupBackgroup[continent] = group
    
    x_pos = groupBackgroup[continent].index.values
    x_pos.sort()
    
    for i in range(nrows):
        method = method_names[i]
        res_info = method_res[method]
        res_info = res_info.sort_values("date")
        for j in range(ncols):
            site = all_prevalent_sites[j + n * ncols]
            continent = sites_rep_continent[site]
            
            res_date = res_info.loc[res_info["aaPos"]  == site, "date"]
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
            
            ax = axes[i][j]
            ax.fill_between(x_pos, 0, bgNum, label = "Backgroup", facecolor='#AFDAE8')
            ax.fill_between(x_pos, 0, daily_num, label = "mutation", facecolor='#F7E15F')
            ax.tick_params(axis='x', labelrotation=60)
            ax.set_xlim([x_pos[1], x_pos[-1]])
            ax.xaxis.set_major_formatter(DateFormatter('%b %Y'))
            ax.set_xlabel(f"{site}_{continent}", fontsize=fontsize, fontweight='bold')
            
            if res_date is not None:
            
                ax.vlines(
                    x=res_date,
                    colors=method_color[method],
                    ymin=0,
                    ymax=max(bgNum),
                    linewidth=2,
                    # linestyles="dashed"
            
                )
    
            if j == 0:
                ax.set_ylabel(method, fontsize=fontsize, fontweight='bold')
    
    plt.savefig(f"Output/{SAMPLING_METHOD}_detection_{PROTEIN_NAME}_{n + 1}.pdf",
                bbox_inches="tight")

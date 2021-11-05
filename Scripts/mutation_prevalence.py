import os
import json

import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.dates import DateFormatter


BACKGROUND_NUM_FILE = "Data/background_num.json"
MUTATION_NUM_FILE = "Data/mutation_num.json"
PREVALENCE_SITES_FILE = "Data/prevalence_sites.json"
PREVALENCE_PLOT = "Output/prevalence.pdf"
PERCENTAGE_SUM_PLOT = "Output/percentage_sum.pdf"

PERCENTAGE_SUM_THRESHOLD = 70

with open(BACKGROUND_NUM_FILE) as f:
    background = json.load(f)

with open(MUTATION_NUM_FILE) as f:
    mutation_num = json.load(f)

mutSites = []
for aaPos, daily_num in mutation_num.items():
    percentage_sum = 0
    
    for d in daily_num:
        percentage_sum += daily_num[d] / background[d]
    
    mutSites.append({
        "aaPos": aaPos,
        "percentage_sum": percentage_sum
    })
mutSites = pd.DataFrame.from_records(mutSites)

plt.rcParams.update({'font.size': 20, 'font.weight': 'bold'})
plt.figure(figsize=(10,10))
fontsize = 20
(n, bins, patches) = plt.hist(mutSites["percentage_sum"], bins=100)
plt.vlines(
    x=PERCENTAGE_SUM_THRESHOLD,
    colors="red",
    ymin=0,
    ymax=max(n),
    linewidth=2,
    linestyles="dashed"
    
)
plt.xlabel("percentage_sum")
plt.ylabel("Frequency (log)")
plt.yscale('log')
plt.savefig(PERCENTAGE_SUM_PLOT, bbox_inches="tight")
plt.show()


with open(BACKGROUND_NUM_FILE) as f:
    background = json.load(f)
    background = pd.DataFrame(background.items(), columns=["date", "num"])
    background["date"] = pd.to_datetime(background["date"])
    background = background.set_index("date")
    background = background.sort_index()

mutSites = mutSites.sort_values("percentage_sum", ascending=False)

sites_daily_num = {}
prevalence_sites = []
for aaPos in mutSites[mutSites["percentage_sum"] > 50]["aaPos"].values:
    daily_num = mutation_num[str(aaPos)]
    daily_num = pd.DataFrame(daily_num.items(), columns=["date", "num"])
    daily_num["date"] = pd.to_datetime(daily_num["date"])
    daily_num = daily_num.set_index("date")
    daily_num = daily_num.sort_index()
    sites_daily_num[aaPos] = daily_num
    prevalence_sites.append(int(aaPos))

x_pos = background.index.values
x_pos.sort()

nrows = 5
ncols = 8

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
for site, meta in sites_daily_num.items():
    mutNum = [meta.loc[d, "num"] if d in meta.index else 0 for d in x_pos]

    ax = axes2[n]
    if n >= nrows * ncols:
        break
    n += 1
    
    ax.fill_between(x_pos, 0, bgNum, label = "Backgroup", facecolor='#AFDAE8')
    ax.fill_between(x_pos, 0, mutNum, label = "mutation", facecolor='#F7E15F')
    ax.tick_params(axis='x', labelrotation=60)
    ax.set_xlim([x_pos[1], x_pos[-1]])
    ax.xaxis.set_major_formatter(DateFormatter('%b %Y'))
    ax.set_xlabel(site, fontsize=fontsize, fontweight='bold')

plt.savefig(PREVALENCE_PLOT, bbox_inches="tight")
plt.show()


with open(PREVALENCE_SITES_FILE, "w") as f:
    json.dump(prevalence_sites, f)




sites_daily_num = {}
for aaPos in mutSites[(mutSites["percentage_sum"] <= 50) &
                      (mutSites["percentage_sum"] > 10)]["aaPos"].values:
    daily_num = mutation_num[str(aaPos)]
    daily_num = pd.DataFrame(daily_num.items(), columns=["date", "num"])
    daily_num["date"] = pd.to_datetime(daily_num["date"])
    daily_num = daily_num.set_index("date")
    daily_num = daily_num.sort_index()
    sites_daily_num[aaPos] = daily_num

x_pos = background.index.values
x_pos.sort()

nrows = 4
ncols = 4

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
for site, meta in sites_daily_num.items():
    mutNum = [meta.loc[d, "num"] if d in meta.index else 0 for d in x_pos]

    ax = axes2[n]
    if n >= nrows * ncols:
        break
    n += 1
    
    ax.fill_between(x_pos, 0, bgNum, label = "Backgroup", facecolor='#AFDAE8')
    ax.fill_between(x_pos, 0, mutNum, label = "mutation", facecolor='#F7E15F')
    ax.tick_params(axis='x', labelrotation=60)
    ax.set_xlim([x_pos[1], x_pos[-1]])
    ax.xaxis.set_major_formatter(DateFormatter('%b %Y'))
    ax.set_xlabel(site, fontsize=fontsize, fontweight='bold')

plt.savefig("Output/prevalence_lower.pdf", bbox_inches="tight")
plt.show()


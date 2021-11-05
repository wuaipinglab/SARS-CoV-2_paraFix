import os
import json

from matplotlib import pyplot as plt
import pandas as pd


BACKGROUND_NUM_FILE = "Data/background_num.json"
MUTATION_NUM_FILE = "Data/mutation_num.json"
PERCENTAGE_SUM_PLOT = "Output/percentage_sum.pdf"

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

plt.clf()
plt.hist(mutSites["percentage_sum"])
plt.xlabel("percentage_sum")
plt.ylabel("Frequency (log)")
plt.yscale('log')
plt.savefig(PERCENTAGE_SUM_PLOT, bbox_inches="tight")
plt.show()

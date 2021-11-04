import os
import copy
import json
from collections import defaultdict

import pandas as pd
from Bio import AlignIO, Phylo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

PROTEIN = "Spike"

TREES_DIR = "Data/nextstrain_trees_with_MSA/"
HOMOPLASYFINDER_DIR = "Data/nextstrain_homoplasyFinder_results/"
DATES_FILE = "Data/nextstrain_dates.json"

# TREES_DIR = "Data/sampled_trees_with_MSA/"
# HOMOPLASYFINDER_DIR = "Data/sampled_homoplasyFinder_results/"
# DATES_FILE = "Data/sampled_dates.json"

if not os.path.exists(HOMOPLASYFINDER_DIR):
    os.mkdir(HOMOPLASYFINDER_DIR)

with open(DATES_FILE) as f:
    allDates = json.load(f)

sitesMapping = pd.read_csv("Data/sitesMapping.csv", index_col=0)

matchedTreeSeq = defaultdict(dict)
for fn in allDates:
    baseName = os.path.join(TREES_DIR, fn, fn)
    matchedTreeSeq[fn]["tree"] = Phylo.read(baseName + ".nwk", "newick")
    matchedTreeSeq[fn]["seq"] = AlignIO.read(baseName + ".fasta", "fasta")


for d in matchedTreeSeq:
    outSeqs = []
    toRemove = []
    for record in matchedTreeSeq[d]["seq"]:
        # record2 = SeqRecord(
        #     Seq(str(record.seq).replace("-", "N")),
        #     id=record.id,
        #     description=""
        # )
        outSeqs.append(record)
    
    print(toRemove)
    outDir = os.path.join(HOMOPLASYFINDER_DIR, d)
    if not os.path.exists(outDir):
        os.mkdir(outDir)
    AlignIO.write(
        MultipleSeqAlignment(outSeqs),
        os.path.join(outDir, d + "_genome.fasta"),
        "fasta"
    )
    tree = copy.deepcopy(matchedTreeSeq[d]["tree"])
    for ac in toRemove:
        ac = tree.find_any(ac)
        tree.prune(ac)
        
    for tip in tree.get_terminals():
        if tip.name is None:
            tree.prune(tip)
            
    Phylo.write(
        tree,
        os.path.join(outDir, d + "_genome.newick"),
        "newick"
    )

# # Prepare data for a single protein
# 
# proteinSites = sitesMapping[(sitesMapping["product"] == PROTEIN) &
#                             (sitesMapping["aa"] != "*")]["aaPos"].drop_duplicates().values
# 
# prevSite = 0
# for site in proteinSites:
#     if site != prevSite + 1:
#         print("Discontinued at", prevSite)
#     prevSite = site
# 
# matchedTreeSeq = defaultdict(dict)
# for fn in allDates:
#     baseName = os.path.join(TREES_DIR, fn, fn)
#     matchedTreeSeq[fn]["tree"] = Phylo.read(baseName + ".nwk", "newick")
#     matchedTreeSeq[fn]["seq"] = AlignIO.read(baseName + "_aa.fasta", "fasta")
# 
# 
# for d in matchedTreeSeq:
#     outSeqs = []
#     toRemove = []
#     for record in matchedTreeSeq[d]["seq"][:, (proteinSites[0] - 1):proteinSites[-1]]:
#         if "*" in record.seq:
#             toRemove.append(record.id)
#         else:
#             outSeqs.append(record)
#     
#     print(toRemove)
#     outDir = os.path.join(HOMOPLASYFINDER_DIR, d)
#     if not os.path.exists(outDir):
#         os.mkdir(outDir)
#     AlignIO.write(
#         MultipleSeqAlignment(outSeqs),
#         os.path.join(outDir, d + "_" + PROTEIN + ".fasta"),
#         "fasta"
#     )
#     tree = copy.deepcopy(matchedTreeSeq[d]["tree"])
#     for ac in toRemove:
#         ac = tree.find_any(ac)
#         tree.prune(ac)
#         
#     for tip in tree.get_terminals():
#         if tip.name is None:
#             tree.prune(tip)
#     Phylo.write(
#         tree,
#         os.path.join(outDir, d + "_" + PROTEIN + ".newick"),
#         "newick"
#     )


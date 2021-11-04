import os
import copy
import json
from collections import defaultdict

import pandas as pd
from Bio import AlignIO, Phylo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Nexus import Nexus, Trees

PROTEIN = "Spike"

TREES_DIR = "Data/nextstrain_trees_with_MSA/"
HYPHY_DIR = "Data/nextstrain_hyphy_results/"
DATES_FILE = "Data/nextstrain_dates.json"

# TREES_DIR = "Data/sampled_trees_with_MSA/"
# HYPHY_DIR = "Data/sampled_hyphy_results/"
# DATES_FILE = "Data/sampled_dates.json"

if not os.path.exists(HYPHY_DIR):
    os.mkdir(HYPHY_DIR)

with open(DATES_FILE) as f:
    allDates = json.load(f)
    
sitesMapping = pd.read_csv("Data/sitesMapping.csv", index_col=0)
proteinSites = sitesMapping[(sitesMapping["product"] == PROTEIN) &
                            (sitesMapping["aa"] != "*")]["genomePos"].values
prevSite = 0
for site in proteinSites:
    if site != prevSite + 1:
        print("Discontinued at", prevSite)
    prevSite = site


matchedTreeSeq = defaultdict(dict)
for fn in allDates:
    baseName = os.path.join(TREES_DIR, fn, fn)
    matchedTreeSeq[fn]["tree"] = Phylo.read(baseName + ".nwk", "newick")
    matchedTreeSeq[fn]["seq"] = AlignIO.read(baseName + ".fasta", "fasta")


for d in matchedTreeSeq:
    outNexus = Nexus.Nexus()
#     outSeqs = []
    toRemove = []
    for record in matchedTreeSeq[d]["seq"][:, (proteinSites[0] - 1):proteinSites[-1]]:
        record = SeqRecord(
            Seq(str(record.seq).replace("-", "n")),
            id=record.id,
            description=""
        )
        if "*" in str(record.translate().seq):
            toRemove.append(record.id)
        else:
#             outSeqs.append(record)
            outNexus.add_sequence(record.id, record.seq.upper())
    # AlignIO.write(
    #     MultipleSeqAlignment(outSeqs),
    #     os.path.join(HYPHY_DIR, d + "_" + PROTEIN + ".fasta"),
    #     "fasta"
    # )
    tree = copy.deepcopy(matchedTreeSeq[d]["tree"])
    for ac in toRemove:
        ac = tree.find_any(ac)
        tree.prune(ac)
    # Phylo.write(
    #     tree,
    #     os.path.join(HYPHY_DIR, d + "_" + PROTEIN + "_test.newick"),
    #     "newick"
    # )
    outNexus.trees.append(Trees.Tree(tree.format("newick")))
    
    outDir = os.path.join(HYPHY_DIR, d)
    if not os.path.exists(outDir):
        os.mkdir(outDir)
    outFileName = os.path.join(outDir, d + "_" + PROTEIN + ".nexus")
    outNexus.write_nexus_data(outFileName)
    with open(outFileName, "a") as f:
        f.write("begin trees;\n")
        f.write(str(outNexus.trees[0]))
        f.write("\nend;\n")

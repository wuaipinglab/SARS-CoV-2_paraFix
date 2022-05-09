#!/usr/bin/env python

import os
import copy
import json
from collections import defaultdict

import pandas as pd
from Bio import AlignIO, Phylo
from Bio.Align import MultipleSeqAlignment


SITESMAPPING_DIR = "data/sitesMapping.csv"
DATES_FILE = "data/nextstrain_dates.json"
TREES_DIR = "data/nextstrain_trees_with_MSA/"

HOMOPLASYFINDER_DIR = "output/nextstrain_homoplasyFinder_results"

if not os.path.exists(HOMOPLASYFINDER_DIR):
    os.mkdir(HOMOPLASYFINDER_DIR)


sitesMapping = pd.read_csv(SITESMAPPING_DIR, index_col=0)

with open(DATES_FILE) as f:
    allDates = json.load(f)


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

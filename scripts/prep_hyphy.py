#!/usr/bin/env python

import os
import copy
import json
from collections import defaultdict

import pandas as pd
from Bio import AlignIO, Phylo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Nexus import Nexus, Trees


PROTEIN_NAMES = ("Spike", "N")

SITESMAPPING_DIR = "data/sitesMapping.csv"
DATES_FILE = "data/nextstrain_dates.json"
TREES_DIR = "data/nextstrain_trees_with_MSA/"

HYPHY_DIR = "output/nextstrain_hyphy_results/"

if not os.path.exists(HYPHY_DIR):
    os.mkdir(HYPHY_DIR)


sitesMapping = pd.read_csv(SITESMAPPING_DIR, index_col=0)
    
with open(DATES_FILE) as f:
    allDates = json.load(f)


for protein in PROTEIN_NAMES:
    proteinSites = sitesMapping[(sitesMapping["product"] == protein) &
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
                outNexus.add_sequence(record.id, record.seq.upper())
        tree = copy.deepcopy(matchedTreeSeq[d]["tree"])
        for ac in toRemove:
            ac = tree.find_any(ac)
            tree.prune(ac)
            
        for tip in tree.get_terminals():
            if tip.name is None:
                tree.prune(tip)
                
        outNexus.trees.append(Trees.Tree(tree.format("newick")))
        
        outDir = os.path.join(HYPHY_DIR, d)
        if not os.path.exists(outDir):
            os.mkdir(outDir)
        outFileName = os.path.join(outDir, d + "_" + protein + ".nexus")
        outNexus.write_nexus_data(outFileName)
        with open(outFileName, "a") as f:
            f.write("begin trees;\n")
            f.write(str(outNexus.trees[0]))
            f.write("\nend;\n")

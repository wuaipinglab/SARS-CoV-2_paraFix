import os
import shutil
import json
from io import StringIO
from urllib import request
from subprocess import Popen, PIPE

import pandas as pd
from Bio import SeqIO
from Bio.PDB import PDBList

REFERENCE = "EPI_ISL_402124"
PROTEIN_NAME = "Spike"
# SAMPLING_METHOD = "nextstrain"
SAMPLING_METHOD = "sampled"

STRUCTURE_DIR = "Data/structure_" + PROTEIN_NAME
PARAFIXSITES_FILE = "Data/" + SAMPLING_METHOD + "_sitePath_results.csv"
HOMOPLASYFINDER_RES_FILE = "Data/" + SAMPLING_METHOD + "_homoplasyFinder.csv"
HYPHY_RES_FILE = "Data/" + SAMPLING_METHOD + "_hyphy_results.csv"

PROTEIN_DOMAIN = {
    "Spike": {
        "RBD": [s for s in range(319, 541 + 1)],
        "NTD": [s for s in range(13, 303 + 1)],
        "S2": [s for s in range(686, 1273 + 1)],
        "FP": [s for s in range(816, 855 + 1)],
        "HR1": [s for s in range(920, 970 + 1)],
        "HR2": [s for s in range(1163, 1202 + 1)]
    }
}

DOMAIN_COLOR = {
    "Seq": "0xBCBDBD",
    "NTD": "0x9BCAC8",
    "RBD": "0xC6716B",
    "S2": "0x9EC9A1",
    "FP": "0xBAD3E1",
    "HR1": "0xF0D0CE",
    "HR2": "0xF0D0CE"
}

SITE_TYPE_COLOR = {
    "sitePath": "gold",
    "common": "red",
    "hyphy": "blue"
}
PDB_ID = { "Spike": "7DF3" }

if not os.path.exists(STRUCTURE_DIR):
    os.mkdir(STRUCTURE_DIR)

# Download the corresponding PDB file
pdbFilePath = PDBList(verbose=False).retrieve_pdb_file(
    pdb_code=PDB_ID[PROTEIN_NAME],
    file_format="pdb",
    pdir=STRUCTURE_DIR
)

# Find the 'DBREF' entry in the PDB file for sites mapping
sitesMapping = []
with open(pdbFilePath) as f:
    monomer = True
    monomer_ids = set()
    previous_e = -1
    completeRef = False
    for line in f:
        if line.startswith("DBREF"):
            if line.startswith("DBREF1 "):
                completeRef = False
                chain_id, pdb_s, pdb_e = (
                    line[12],
                    int(line[14:18]),
                    int(line[20:24]) + 1,
                )
            elif line.startswith("DBREF2 "):
                completeRef = True
                chain_id, ref_s, ref_e, uniprot_id = (
                    line[12],
                    int(line[45:55]),
                    int(line[57:67]) + 1,
                    line[18:40].strip()
                )
            if line.startswith("DBREF "):
                completeRef = True
                chain_id, pdb_s, pdb_e, ref_s, ref_e, uniprot_id = (
                    line[12],
                    int(line[14:18]),
                    int(line[20:24]) + 1,
                    int(line[55:60]),
                    int(line[62:67]) + 1,
                    line[33:41].strip(),
                )
            if completeRef:
                if previous_e > ref_s:
                    monomer = False
                previous_e = ref_e
                if monomer:
                    monomer_ids.add(chain_id)
                    for pdb_p, ref_p in zip(range(pdb_s, pdb_e), range(ref_s, ref_e)):
                        sitesMapping.append({
                            "uniprotSite": ref_p,
                            "chain": chain_id,
                            "pdbSite": pdb_p
                        })
    sitesMappingAll = pd.DataFrame.from_records(sitesMapping)

# Download the uniprot sequence
res = request.urlopen("https://www.uniprot.org/uniprot/" + uniprot_id + ".fasta")
uniprotSeq = SeqIO.read(StringIO(res.read().decode("utf-8")), "fasta")
uniprotSeq.id = uniprotSeq.description.split('|')[1]
uniprotSeq.description = "Reference"

# Combine the sequences for multiple sequence alignment
sequences = [record for record in SeqIO.parse("Data/selected_Spike.fasta", "fasta")]
refAddedFilePath = os.path.join(STRUCTURE_DIR, "refAdded.fasta")
SeqIO.write(
    sequences=[*sequences, uniprotSeq],
    handle=refAddedFilePath,
    format="fasta"
)

# Multiple sequence alignment using muscle
process = Popen(
    ["muscle", "-in", refAddedFilePath, "-diags"],
    text=True,
    stdout=PIPE,
    stderr=PIPE
)
stdout, stderr = process.communicate()
if process.returncode:
    raise BaseException("muscle: " + stderr)
    
# Map the site numbering between PDB DBREF (uniprot) and the reference
seqs = SeqIO.to_dict(SeqIO.parse(StringIO(stdout), "fasta"))
refSite, uniprotSite = 0, 0
sitesMapping = []
for n, (ref_aa, uniprot_aa) in enumerate(zip(seqs[REFERENCE].seq, seqs[uniprot_id].seq), start=1):
    validSite = False
    if ref_aa != "-":
        refSite += 1
        validSite = True
    if uniprot_aa != "-":
        uniprotSite += 1
        validSite = True
    if validSite:
        sitesMapping.append({"refSite": refSite, "uniprotSite": uniprotSite, "alignedSite": n})
sitesMapping = pd.DataFrame.from_records(sitesMapping)

sitesMappingAll = pd.merge(sitesMappingAll, sitesMapping, on="uniprotSite", how="left")


# The sites (numbering based on EPI_ISL_402125) to be marked on the structure

hyphy_res = pd.read_csv(HYPHY_RES_FILE)
sitePath_res = pd.read_csv(PARAFIXSITES_FILE)

hyphy_sites = set(hyphy_res["site"].unique())
sitePath_sites = set(sitePath_res["aaPos"].unique())

targetSites = {
    "sitePath": list(sitePath_sites.difference(hyphy_sites)),
    "hyphy": list(hyphy_sites.difference(sitePath_sites)),
    "common": list(sitePath_sites.intersection(hyphy_sites))
}
    
chains = set()
for siteType, sites in targetSites.items():
    sites = sitesMappingAll[sitesMappingAll["refSite"].isin(sites)]
    sites = sites[["chain", "pdbSite"]].drop_duplicates()
    targetSites[siteType] = sites
    chains = chains.union(sites["chain"].unique())



pmlSetView = "set_view (\
     0.787351847,   -0.214461818,    0.577996790,\
     0.602337480,    0.067778848,   -0.795357168,\
     0.131397858,    0.974372327,    0.182546288,\
     0.000000000,    0.000000000, -539.911254883,\
   163.201431274,  163.200393677,  166.706909180,\
  -16535.880859375, 17615.703125000,  -20.000000000 )"

pmlScript = [
    "delete *",
    "fetch {}".format(PDB_ID[PROTEIN_NAME]),
    "hide",
    "bg_color white",
    "show surface",
    "set transparency, 0.9",
    "set surface_color, grey",
    "; ".join(["show cartoon, c. " + c for c in chains]),
    "color grey"
]

for domain, sites in PROTEIN_DOMAIN[PROTEIN_NAME].items():
    sites = sitesMappingAll.loc[sitesMappingAll["refSite"].isin(sites), ["chain", "pdbSite"]]
    sel = " or ".join([f"c. {c} and i. {i}" for _, c, i in sites.itertuples()])
    col = DOMAIN_COLOR[domain]
    pmlScript.extend([
        f"select {domain}, {sel}",
        f"color {col}, {domain}",
        f"set surface_color, {col}, {domain}",
        f"set transparency, 0.75, {domain}"
#         f"create Iso_{domain}, {domain}",
#         f"show surface, Iso_{domain}",
#         f"set surface_color, {col}, Iso_{domain}",
#         f"set transparency, 0.1, Iso_{domain}"
    ])

# Script for fixation sites
for category, sites in targetSites.items():
    sel = " or ".join([f"c. {c} and i. {i}" for _, c, i in sites.itertuples()])
    col = SITE_TYPE_COLOR[category]
    pmlScript.extend([
        f"select {category}, {sel}",
        f"color {col}, {category}",
        f"show sphere, {category}",
        f"disable {category}"
    ])
    pmlScript.append("reset")
    pmlScript.append(pmlSetView)
    # # For the site's label
    # label = " or ".join([f"c. {c} and i. {i} and n. CA" for _, c, i in sites.itertuples()])
    # pmlScript.extend([
    #     "set label_position,(5,7,100)",
    #     "set label_size, 10",
    #     "set label_connector, true",
    #     "set label_bg_color, white",
    #     "set label_bg_transparency, 0.2",
    #     f"label {label}, resi\n"
    # ])

with open(f"Output/{SAMPLING_METHOD}_{PROTEIN_NAME}.pml", 'w') as f:
    f.write("; ".join(pmlScript))
  

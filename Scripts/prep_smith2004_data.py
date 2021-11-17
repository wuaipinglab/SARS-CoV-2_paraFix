import os
import re
import json
from io import StringIO
from collections import defaultdict
from subprocess import Popen, PIPE

from Bio import Entrez, SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

DATA_DIR = os.path.join("Data", "H3N2_HA1_Smith2004")
PMID = "15218094"
ADDITIONAL = (
    "AF008665", "AF008697", "AF008711",
    "AF008725", "AF008755", "AF008769",
    "AF008828", "AF008867", "AF008886",
    "AF008888", "AF008903", "AF008905",
    "AF092062", "AF131997", "AF180570",
    "AF180602", "AF180643", "AF201874",
    "AF368444", "AF368446", "D21173",
    "D49961", "ISDN38157", "ISDNCDA001",
    "ISDNENG72", "ISDNHK71", "ISDNTX77",
    "ISDNVIC75",
#     "M16739",
    "U08858", "Z46405", "Z46408", "Z46413", "Z46414",
    # Manually add missing records
    "CY112289", "AF201875", "EF626609",
    "CY120992", "KM821316", "KM821317",
    "KM821324",
#     "CY113525"
)

REGION_NAMES = {
    "Akita": 'AK',
    "Amsterdam": 'AM',
    "Atlanta": 'AT',
    "Auckland": 'AU',
    "Beijing": 'BE',
    "Bilthoven": 'BI',
    "Brisbane": 'BR',
    "Canberra": 'CA',
    "Christchurch": 'CC',
    "Caen": 'CE',
    "Colorado": 'CO',
    "England": 'EN',
    "Enschede": 'ES',
    "Finland": 'FI',
    "Fujian": 'FU',
    "Guangdong": 'GD',
    "Geneva": 'GE',
    "Guildford": 'GF',
    "Guizhou": 'GU',
    "Hong Kong": 'HK',
    "Houston": 'HO',
    "Johanesburg": 'JO',
    "Johannesburg": "JO",
    "Leningrad": 'LE',
    "Lyon": 'LY',
    "Madrid": 'MA',
    "Memphis": 'ME',
    "Moscow": 'MW',
    "Nanchang": 'NA',
    "Nice": 'NE',
    "Nijmegen": 'NI',
    "Netherlands": 'NL',
    "Oslo": 'OS',
#     "Madrid": 'OV', # Not really
    "Paris": 'PA',
    "Port Chalmers": 'PC',
    "Phillipines": 'PH',
    "Philippines": 'PH',
    "Panama": 'PM',
    "Rotterdam": 'RD',
    "South Austalia": 'SA',
    "South Australia": 'SA',
    "Shandong": 'SD',
    "Sendai": 'SE',
    "Seoul": "SU",
    "Shiga": 'SG',
    "Shanghai": 'SH',
    "Sichuan": 'SI',
    "Scotland": 'SL',
    "Singapore": 'SP',
    "Stockholm": 'ST',
#     "Suita": 'SU',
    "Sydney": 'SY',
    "Texas": 'TE',
    "Tilburg": 'TI',
    "Umea": 'UM',
    "Victoria": 'VI',
    "Wellington": 'WE',
    "Wuhan": 'WU',
    "Yamagata": 'YA'
}

# Internet access to NCBI entrez service is required
Entrez.email = "youremail@domain.com"

(res,) = Entrez.read(Entrez.elink(linkname="pubmed_nuccore", id=PMID, idtype="acc"))
(res,) = res["LinkSetDb"]

IdList = [*[acc["Id"] for acc in res["Link"]], *ADDITIONAL]

# Download the sequences from Smith's paper
handle = Entrez.efetch(db="nuccore", id=IdList, rettype="gb", retmode="text")
records = [record for record in SeqIO.parse(handle, "gb")]

refSeq = SeqIO.read(os.path.join(DATA_DIR, "HA1_reference.fasta"), "fasta")

outSeqs = [refSeq]
seqname2ac = defaultdict(set)

for record in records:
    organism = None
    organism2 = None
    seq = None
    for feature in record.features:
        if feature.type == "source":
            # Get the organism name
            (organism,) = feature.qualifiers["organism"]
            m = re.search(r"/[A-Za-z ]+/[A-Za-z0-9]+/[0-9]+", organism)
            organism = m.group(0)
            _, region, name, year = organism.split("/")
            m = re.search(r"[0-9]+", name)
            name = m.group(0)
            # The renamed organism for mapping cluster name
            organism = (REGION_NAMES[region], name, year[-2:])
            organism2 = organism
            # Get the isolate/strain name of the virus
            if "strain" in feature.qualifiers:
                (strain,) = feature.qualifiers["strain"]
            elif "isolate" in feature.qualifiers:
                (strain,) = feature.qualifiers["isolate"]
            # Correct the year if possible
            year2 = strain.split("/")[-1]
            if year2:
                year2 = year2.split('(')[0].strip()
                organism2 = (REGION_NAMES[region], name, year2[-2:])
        if feature.type == "CDS":
            (seq,) = feature.qualifiers["translation"]
    # Map the renamed organism with accession id
    seqname2ac[organism].add(record.id)
    seqname2ac[organism2].add(record.id)
    outSeqs.append(SeqRecord(seq=Seq(seq), id=record.id, description=""))

sequencesFilePath = os.path.join(DATA_DIR, "sequences.fasta")
SeqIO.write(outSeqs, sequencesFilePath, "fasta")

clusters = {}

with open(os.path.join(DATA_DIR, "metadata_copied.txt")) as f:
    for row in f:
        row = row.split(" ")
        seqname = row[1].split("/")
        m = re.search(r"[0-9]+", seqname[1])
        seqname[1] = m.group(0)
        if seqname[0] == "OV":
            seqname[0] = "MA"
        seqname = tuple(seqname)
        for ac in seqname2ac[seqname]:
            clusters[ac] = row[0]

with open(os.path.join(DATA_DIR, "metadata.json"), "w") as f:
    json.dump(clusters, f)

all(record.id in clusters for record in records)

# Multiple sequence alignment using muscle
process = Popen(
    ["muscle", "-in", sequencesFilePath, "-diags"],
    stdout=PIPE,
    stderr=PIPE
    
)
stdout, stderr = process.communicate()
if process.returncode:
    raise BaseException("muscle: " + stderr.decode("UTF-8"))
    
# Read in the output MSA
seqs = AlignIO.read(StringIO(stdout.decode("UTF-8")), "fasta")
refSeq_index = None
for n, record in enumerate(seqs):
    if record.id == refSeq.id:
        refSeq_index = n
        break

assert refSeq_index is not None
prev_gap_index = -1
prev_aa_index = -1
start_index = None
end_index = None

reference_indexes = []

for index, aa in enumerate(seqs[refSeq_index]):
    if aa == "-":
        if index != prev_gap_index + 1:
            print(index, aa)
            reference_indexes.append((start_index, end_index))
        prev_gap_index = index
    else:
        if index != prev_aa_index + 1:
            start_index = index
        end_index = index + 1
        prev_aa_index = index
        

non_gap_seqs = None

for start_index, end_index in reference_indexes:
    if non_gap_seqs:
        non_gap_seqs += seqs[:, start_index:end_index]
    else:
        non_gap_seqs = seqs[:, start_index:end_index]
out_aligned_seqs = []

for n, record in enumerate(non_gap_seqs):
    if n != refSeq_index:
        out_aligned_seqs.append(record)
        
out_aligned_seqs = MultipleSeqAlignment(out_aligned_seqs)

AlignIO.write(out_aligned_seqs, os.path.join(DATA_DIR, "HA1_aligned.fasta"), "fasta")


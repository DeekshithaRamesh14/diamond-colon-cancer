"""
prepare_data.py   (step 1 of 2)

Makes the two input files that the authors' DIAMOnD.py needs:
    data/network.txt   the protein-protein interaction network: one interaction (two gene IDs) per line
    data/seeds.txt     the known colon cancer genes ("seeds"): one gene ID per line
It also saves data/gene_symbols.tsv (gene ID -> gene name), only so that people can read the results.

Run:   python prepare_data.py
Next:  python DIAMOnD.py data/network.txt data/seeds.txt 100 1 out.txt
"""

import csv                                  # read and write tables
import os                                   # folders and file paths
import urllib.request                       # download files from the internet

# ----------------------------- settings ------------------------------------
DISEASE = "DOID:219"                        # colon cancer, as written in the DisGeNET table
CURATED = {"UNIPROT", "CLINVAR"}            # a gene is a seed only if one of these curated databases reports it

# the three data files on GitHub; every link points to one fixed version (commit), so everybody gets the same data
PPI_URL = "https://raw.githubusercontent.com/dhimmel/ppi/9ae2ca12c20e0566e6644dc05b908de7f37ce15a/data/ppi-hetio-ind.tsv"
DISGENET_URL = "https://raw.githubusercontent.com/dhimmel/disgenet/a636cc0912105e2b8d8c6f2a99ea0b9444fdd412/data/consolidated.tsv"
NAMES_URL = "https://raw.githubusercontent.com/hetio/hetionet/8a6cc0c5604d0e2908786e4de38af67b8e46ee4a/hetnet/tsv/hetionet-v1.0-nodes.tsv"

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")   # the folder "data" next to this file


def download(url, path):
    """download one file, unless it is already there"""
    if not os.path.exists(path):                                              # skip if we already have it
        print("downloading", os.path.basename(path))                          # tell the user
        request = urllib.request.Request(url, headers={"User-Agent": "student-project"})   # GitHub wants a name
        with urllib.request.urlopen(request, timeout=120) as response:        # open the web address
            with open(path, "wb") as f:                                       # open the local file
                f.write(response.read())                                      # copy the content


os.makedirs(DATA_DIR, exist_ok=True)                                          # create the folder "data"
ppi_file = os.path.join(DATA_DIR, "ppi.tsv")                                  # raw interaction table
disgenet_file = os.path.join(DATA_DIR, "disgenet.tsv")                        # raw gene-disease table
names_file = os.path.join(DATA_DIR, "gene_names.tsv")                         # raw table with gene names
download(PPI_URL, ppi_file)                                                   # 1. interactions
download(DISGENET_URL, disgenet_file)                                         # 2. gene-disease associations
download(NAMES_URL, names_file)                                               # 3. gene names

# ---- network.txt: every interaction once, as "geneA,geneB" ----
links = set()                                                                 # a set keeps every pair only once
with open(ppi_file) as f:
    next(f)                                                                   # skip the header line
    for line in f:
        a, b = line.rstrip("\n").split("\t")[:2]                              # the first two columns = the two genes
        if a != b:                                                            # ignore a protein linked to itself
            links.add((min(a, b), max(a, b)))                                 # sort the pair, so A-B and B-A are the same
with open(os.path.join(DATA_DIR, "network.txt"), "w") as f:
    for a, b in sorted(links):
        f.write(a + "," + b + "\n")                                           # one interaction per line

# ---- seeds.txt: colon cancer genes reported by a curated database (UniProt or ClinVar) ----
seeds = set()                                                                 # a set, so no gene appears twice
with open(disgenet_file) as f:
    for row in csv.DictReader(f, delimiter="\t"):                             # read the table row by row
        if row["doid_code"] == DISEASE:                                       # only colon cancer rows
            if CURATED & set(row["source"].split("|")):                       # reported by UniProt or ClinVar?
                seeds.add(row["geneId"])                                      # keep this gene
with open(os.path.join(DATA_DIR, "seeds.txt"), "w") as f:
    for gene in sorted(seeds):
        f.write(gene + "\n")                                                  # one gene ID per line

# ---- gene_symbols.tsv: gene ID -> gene name (only for reading the results) ----
with open(names_file) as f, open(os.path.join(DATA_DIR, "gene_symbols.tsv"), "w") as out:
    for row in csv.DictReader(f, delimiter="\t"):
        if row["kind"] == "Gene":                                             # the table also has other node types
            out.write(row["id"].split("::")[1] + "\t" + row["name"] + "\n")   # "Gene::7157" -> "7157", then its name

print("network: %d interactions | seeds: %d" % (len(links), len(seeds)))

# DIAMOnD on colon cancer with curated seed genes and five values of `alpha`

## Question

DIAMOnD (Ghiassian, Menche, Barabási 2015) starts from known disease genes known as seeds in a protein-protein interaction network and adds proteins one by one that are connected to the seeds more than expected by chance. The `alpha` is the weight given to the seeds. The author used default 1 as alpha value, as bigger number counts the seeds more.

My questions:

1. If I run the authors' DIAMOnD on my own data, which has a different network and curated **colon cancer** genes from DisGeNET, which proteins does it predict?
2. Does changing `alpha` (1, 2, 3, 5, 10) change the predictions, and does it find hidden seed genes better or worse?
3. Is the connectivity that DIAMOnD relies on real, or would random networks like the two null models (ER and MS model) give the same result?

What I changed compared with the authors' tutorial:

* **(a) a new dataset** - a different network and a different gene list for colon cancer
* **(b) a parameter changes over 5 values of** `alpha`; I also hide some seeds to see if DIAMOnD finds them again
* **(c) a comparison with a method from the course**: the two null models of practical P3 (Erdős-Rényi and Maslov-Sneppen), applied to the seeds. I also added a fixed hash seed, so that every run gives exactly the same numbers.

## Data

All data is downloaded by a separate small program (`python prepare_data.py`) from GitHub.

| what | source | use |
|---|---|---|
| protein-protein interactions | `dhimmel/ppi`, file `ppi-hetio-ind.tsv` (a human network compiled from several resources for the Hetionet project) | the network: 15,475 proteins, 148,778 links |
| gene-disease table | `dhimmel/disgenet`, file `consolidated.tsv` | colon cancer genes reported by UniProt or ClinVar = the seeds |
| gene names | `hetio/hetionet`, file `hetionet-v1.0-nodes.tsv` | only to show gene symbols |
| the authors' code file | `dinaghiassian/DIAMOnD`, file `DIAMOnD.py` | my `DIAMOnD.py` is this file with one changed line (in `logchoose()`) and an added block; MIT licence, see `LICENSE_DIAMOnD_authors.txt` |

So, the idea is the same: a physical-interaction network plus known disease genes, but my network is a different mixture, and my seeds are a different, shorter list.

**Why curated seeds, and how I chose the disease.**

The authors only used diseases with at least 20 known genes from OMIM and from genome-wide significant GWAS hits. My first attempt was periodontitis, but only about 5 of its genes come from a curated source in this DisGeNET, too few to hide some of them for the robustness test. So, I counted, for every disease in the file, the genes reported by UniProt or ClinVar that are in my network.

**Seed rule:** a gene is a seed if UniProt or ClinVar report it for colon cancer. This gives 22 genes, of which **21 are in the network**.

## How to run

There are two programs. The first prepares the two input files; the second is the authors' `DIAMOnD.py` with my changes and is given those files.

In Windows Command Prompt, run:

```
python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
python prepare_data.py
python DIAMOnD.py data/network.txt data/seeds.txt 100 1 diamond_output.txt
python DIAMOnD.py data/network.txt data/seeds.txt --analysis
```

* `python prepare_data.py` downloads the data into `data/` and writes `data/network.txt` (the interactions) and `data/seeds.txt` (the colon cancer genes).
* `python DIAMOnD.py data/network.txt data/seeds.txt 100 1 diamond_output.txt`: network file, seed file, number of proteins to add (100), `alpha` (1), output file. It is one plain DIAMOnD run on my data.
* `python DIAMOnD.py data/network.txt data/seeds.txt --analysis`: the second `DIAMOnD.py` command (`--analysis`) runs my three changes: the `alpha` sweep with the hide-seeds test, the null models, and the figures. Tables and figures are written into `results/`. It takes about 30 minutes.
* The program sets `PYTHONHASHSEED=0` by restarting itself once.
* A harmless `SyntaxWarning` may appear; it comes from a comment in the authors' code.

## What is the parameter change, and why does the hide-seeds test belong to it?

The parameter I change is `alpha`, the seed weight (1, 2, 3, 5, 10). But `alpha` only changes which proteins are predicted; to say whether a different list is better, I need a test with a known answer. So I give DIAMOnD only 75% of the seeds (the "seeds" input is reduced) and check how many of the 25% hidden ones it finds again.

**Why the fixed hash seed:** When several proteins have exactly the same p-value, the authors' code keeps one of them, and which one depends on Python's internal ordering of sets. That ordering changes between runs unless `PYTHONHASHSEED` is fixed. I checked this on an earlier version of this program (with the earlier periodontitis data): two runs, one of them started with a different `PYTHONHASHSEED` in the environment, gave identical result files. In this two-file version, the hold-out and `alpha` results of the full run were identical to the full run of that earlier version.

**Comparison method inside the hide-seeds test.** I also rank proteins by "number of links to the visible seeds" (a very simple method), and I show the chance level (picking proteins at random).

**Course method (`run_null_model()`): the null models of practical P3.** A null model is a random version of the real network used as a control. Practical P3 has two:

* **Erdős–Rényi (ER):** the same number of proteins and links, everything else random. It asks: is the statistic explained by size and density alone?
* **Maslov–Sneppen (M–S):** the real network with its links swapped, so that every protein keeps exactly its number of links (degree). One swap takes two links A–B and C–D and rewires them to A–D and C–B. It asks: is the statistic explained by the degrees alone?

## Results

The tables, the two figures and the log of my full run are in the `results/` folder. The interpretation is in a separate Word document.

## Credit

DIAMOnD: Ghiassian SD, Menche J, Barabási A-L (2015), *A DIseAse MOdule Detection (DIAMOnD) algorithm derived from a systematic analysis of connectivity patterns of disease proteins in the human interactome*, PLoS Comput Biol 11(4): e1004120. Code: github.com/dinaghiassian/DIAMOnD (MIT licence; my `DIAMOnD.py` contains their code, see `LICENSE_DIAMOnD_authors.txt`). Network and DisGeNET files: the `dhimmel/ppi` and `dhimmel/disgenet` repositories; gene names: Hetionet.

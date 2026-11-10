# DIAMOnD on colon cancer

## Question

DIAMOnD (Ghiassian, Menche, Barabási 2015) starts from known disease genes (the seeds) in a protein-protein interaction network and adds proteins one at a time, always the one that is connected to the module more than expected by chance. The parameter `alpha` is the weight given to the seeds (the authors' default is 1).

My questions:

1. If I run the authors' DIAMOnD on a new dataset (colon cancer), which proteins does it predict?
2. Does changing `alpha` (1, 2, 3, 5, 10) change the predictions, and does it find hidden seed genes better or worse?
3. How is DIAMOnD different from a module method taught in the course - Louvain community detection?

## Data source

All data is downloaded automatically by `prepare_data.py` from GitHub:

* **Network:** protein-protein interactions from `dhimmel/ppi` (file `ppi-hetio-ind.tsv`): 15,475 proteins and 148,778 links.
* **Seed genes:** colon cancer genes from `dhimmel/disgenet` (file `consolidated.tsv`), taken only if UniProt or ClinVar reports them: 22 genes, 21 of them in the network.
* **Gene names:** `hetio/hetionet` (file `hetionet-v1.0-nodes.tsv`), only to read the results.
* **Authors' code:** `dinaghiassian/DIAMOnD` (file `DIAMOnD.py`). It is kept unchanged in `authors_original/`.

## How to run

In Windows Command Prompt:

```
python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
python prepare_data.py
python run_authors_original.py data/network.txt data/seeds.txt 100 1 results/original_output.txt
python DIAMOnD.py data/network.txt data/seeds.txt 100 1,2,3,5,10 --hide-seeds
python compare_with_coursework.py
```

1. `prepare_data.py` downloads the data and writes `data/network.txt` and `data/seeds.txt`.
2. `run_authors_original.py` runs the authors' original code with no change (alpha = 1, 100 proteins). It only replaces one line that fails on new NumPy, in memory.
3. `DIAMOnD.py` is the authors' code with the `alpha` list added (lines marked `CHANGED`). It runs `alpha` = 1, 2, 3, 5, 10 and the hide-seeds test.
4. `compare_with_coursework.py` compares DIAMOnD with Louvain communities (a few seconds).

All results are written to the `results/` folder.

## Credit

* DIAMOnD: Ghiassian SD, Menche J, Barabási A-L (2015), PLoS Comput Biol 11(4): e1004120. Code: github.com/dinaghiassian/DIAMOnD (MIT licence, see `LICENSE_DIAMOnD_authors.txt`).
* Data: `dhimmel/ppi`, `dhimmel/disgenet` and `hetio/hetionet` on GitHub.
* Tools: Python, NetworkX, NumPy, SciPy and python-igraph.
* Course: CB 5113 Systems Biology & Network Modeling (Louvain communities and the hypergeometric test).
* AI assistance: Claude (Anthropic) helped write and test the code and the documents.

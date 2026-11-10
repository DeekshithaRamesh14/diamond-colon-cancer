"""
compare_with_coursework.py   (step 3 of 3)

Compares DIAMOnD with ONE method from the course (CB 5113): Louvain community detection (week 6, modularity and communities).
Why this one: DIAMOnD is a disease-module method, and the course method that also looks for modules in a network is community detection.
The DIAMOnD paper itself uses Louvain as its comparison (Fig. 1F), so this is the fair course counterpart.

It does NOT change DIAMOnD. It reads the list of 100 proteins that the authors' unchanged DIAMOnD.py made (Run 1) and asks:
  1. Louvain is run on the whole network WITHOUT looking at the seeds. Where do the colon cancer seeds fall?
     (hypergeometric test + Benjamini-Hochberg correction, as in practical P5, Part 4)
  2. Do DIAMOnD's 100 picks fall in the community that holds most seeds?
  3. Hide-seeds test (the same hidden seeds as in DIAMOnD.py): can the best Louvain community find the hidden seeds,
     compared with DIAMOnD?

Run:
    python compare_with_coursework.py                      (uses results/original_output.txt)
    python compare_with_coursework.py path/to/authors_output.txt
Time: a few seconds.   Needs: python-igraph, numpy, scipy (all in requirements.txt)
"""

import os                                   # file paths
import sys                                  # command-line arguments
import csv                                  # read and write tables
import random                               # random numbers (igraph uses Python's random)
import numpy as np                          # numbers and arrays
import igraph as ig                         # igraph, the network tool of the course
from scipy.stats import hypergeom, rankdata   # hypergeometric test; ranks

# ----------------------------- settings ------------------------------------
NETWORK_FILE = "data/network.txt"           # the network (made by prepare_data.py)
SEED_FILE = "data/seeds.txt"                # the known colon cancer genes
SYMBOL_FILE = "data/gene_symbols.tsv"       # gene ID -> gene name (only for reading)
RESULTS_DIR = "results"                     # output folder
N_PICKS = 100                               # how many DIAMOnD picks we look at
HIDE_FRACTION = 0.25                        # same hide-seeds test as DIAMOnD.py: hide 25% of the seeds ...
N_REPEATS = 30                              # ... 30 times
RNG_SEED = 2026                             # same random seed as DIAMOnD.py, so the SAME seeds are hidden in each round

log_lines = []                              # everything printed is also saved in results/course_comparison_log.txt


def say(text=""):
    """print a line and remember it for the log file"""
    print(text)
    log_lines.append(text)


def write_csv(name, header, rows):
    """save a table in the results folder"""
    with open(os.path.join(RESULTS_DIR, name), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)


def mean_and_se(values):
    """mean and standard error of a list of numbers"""
    v = np.array(values, dtype=float)
    return v.mean(), v.std(ddof=1) / np.sqrt(len(v))


def benjamini_hochberg(p):
    """Benjamini-Hochberg correction for many tests (as p.adjust(method='BH') in R)"""
    p = np.array(p, dtype=float)
    n = len(p)
    order = np.argsort(p)                                       # smallest p first
    adjusted = p[order] * n / (np.arange(n) + 1)                # p x n / rank
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]      # keep the order of the adjusted values
    out = np.empty(n)
    out[order] = np.minimum(adjusted, 1.0)
    return out


# ============================================================================
# 0. Load everything
# ============================================================================
os.makedirs(RESULTS_DIR, exist_ok=True)
diamond_file = sys.argv[1] if len(sys.argv) > 1 else os.path.join(RESULTS_DIR, "original_output.txt")

edges = [tuple(line.strip().split(",")) for line in open(NETWORK_FILE) if line.strip()]   # one interaction per line
g = ig.Graph.TupleList(edges, directed=False)                   # the network as an igraph graph
names = g.vs["name"]                                            # igraph numbers proteins 0, 1, 2 ...; this keeps their gene IDs
index = {n: i for i, n in enumerate(names)}                     # gene ID -> igraph number
N = g.vcount()                                                  # number of proteins

symbols = {}                                                    # gene ID -> gene name
for line in open(SYMBOL_FILE):
    parts = line.rstrip("\n").split("\t")
    if len(parts) > 1:
        symbols[parts[0]] = parts[1]

seed_genes = sorted(s.strip() for s in open(SEED_FILE) if s.strip() in index)   # seeds that are in the network (sorted = same order every run)
seeds = [index[s] for s in seed_genes]                          # igraph numbers of the seeds
seed_set = set(seeds)

picks = []                                                      # DIAMOnD's picks, in the order they were added
for line in open(diamond_file):
    if not line.startswith("#") and line.strip():
        picks.append(index[line.split()[1]])
picks = picks[:N_PICKS]

say("Network: %d proteins, %d links | seeds in the network: %d | DIAMOnD picks read from %s: %d" %
    (N, g.ecount(), len(seeds), diamond_file, len(picks)))

# the simple baseline of DIAMOnD.py: the 100 non-seed proteins with the most links to seeds (ties by gene ID)
links_to_seeds = {}
for s in seeds:
    for nb in g.neighbors(s):
        if nb not in seed_set:
            links_to_seeds[nb] = links_to_seeds.get(nb, 0) + 1
baseline = sorted(links_to_seeds, key=lambda v: (-links_to_seeds[v], names[v]))[:N_PICKS]

# ============================================================================
# Louvain community detection (course week 6)
# ============================================================================
say("\nLOUVAIN (course method): what does it say about the colon cancer genes?")
random.seed(RNG_SEED)                                           # Louvain has a random part; fixed so everybody gets the same result
louvain = g.community_multilevel()                              # Louvain communities of the whole network (seeds are NOT used)
member = louvain.membership                                     # community number of each protein
sizes = [len(c) for c in louvain]                               # size of each community
say("  Louvain finds %d communities, modularity %.3f, sizes from %d to %d proteins" % (len(louvain), louvain.modularity, min(sizes), max(sizes)))

n_seeds = len(seeds)
seeds_in = np.array([sum(member[s] == c for s in seeds) for c in range(len(louvain))])   # seeds per community
p_raw = [hypergeom.sf(k - 1, N, n_seeds, sizes[c]) if k > 0 else 1.0 for c, k in enumerate(seeds_in)]   # P(at least k seeds by chance)
p_adj = benjamini_hochberg(p_raw)                               # correct for testing every community (as in practical P5)
rows_b = [[c, sizes[c], int(seeds_in[c]), float("%.3g" % p_raw[c]), float("%.3g" % p_adj[c])] for c in range(len(louvain)) if seeds_in[c] > 0]
rows_b.sort(key=lambda r: -r[2])
write_csv("course_communities.csv", ["community", "size", "seeds_inside", "p_hypergeometric", "p_BH_adjusted"], rows_b)
say("  communities that contain seeds (of %d seeds):" % n_seeds)
for r in rows_b:
    say("    community %-3d size %5d  seeds %2d  p = %.2g  BH-adjusted p = %.2g" % tuple(r))
say("  seeds are spread over %d of %d communities" % (len(rows_b), len(louvain)))

in_seed_comm = {c for c in range(len(louvain)) if seeds_in[c] > 0}   # communities that hold at least one seed
share_chance = sum(1 for v in range(N) if member[v] in in_seed_comm and v not in seed_set) / (N - n_seeds)
share_diamond = sum(member[v] in in_seed_comm for v in picks) / len(picks)
say("  share of DIAMOnD's 100 picks that lie in a community with a seed: %.2f (for a random protein: %.2f)" % (share_diamond, share_chance))
top_comm = int(np.argmax(seeds_in))                              # the community with the most seeds
say("  DIAMOnD picks inside the community with most seeds (community %d, %d proteins = %.0f%% of the network): %d of %d" %
    (top_comm, sizes[top_comm], 100 * sizes[top_comm] / N, sum(member[v] == top_comm for v in picks), len(picks)))

# hide-seeds test for the best community (same hidden seeds as in DIAMOnD.py: same random seed and same sorted seed list)
random.seed(RNG_SEED)
n_hide = round(HIDE_FRACTION * len(seeds))
hold = {"all": [], "top100": [], "size": []}
for rep in range(N_REPEATS):
    hidden = random.sample(seed_genes, n_hide)                  # the same call as in DIAMOnD.py -> the same hidden seeds
    hidden_ids = {index[h] for h in hidden}
    visible_ids = [s for s in seeds if s not in hidden_ids]
    count = [0] * len(louvain)                                  # visible seeds per community
    for s in visible_ids:
        count[member[s]] += 1
    best = max(range(len(louvain)), key=lambda c: (count[c], -sizes[c]))   # community with most visible seeds (smaller one if equal)
    cand = [v for v in range(N) if member[v] == best and v not in visible_ids]   # all its members except the visible seeds
    hold["all"].append(len(hidden_ids & set(cand)) / n_hide)    # recall if we predict the WHOLE community
    links = {}                                                  # links of each candidate to visible seeds
    for s in visible_ids:
        for nb in g.neighbors(s):
            links[nb] = links.get(nb, 0) + 1
    cand_sorted = sorted(cand, key=lambda v: (-links.get(v, 0), names[v]))[:N_PICKS]   # best 100 of the community by links to seeds
    hold["top100"].append(len(hidden_ids & set(cand_sorted)) / n_hide)
    hold["size"].append(len(cand))

diamond_recall = None                                           # DIAMOnD alpha=1 recall at k=100, from the earlier hold-out run
summary_file = os.path.join(RESULTS_DIR, "holdout_summary.csv")
if os.path.exists(summary_file):
    for row in csv.DictReader(open(summary_file)):
        if row["method"] == "alpha=1" and row["k"] == "100":
            diamond_recall = (float(row["mean_recall"]), float(row["standard_error"]))
m_all, s_all = mean_and_se(hold["all"])
m_top, s_top = mean_and_se(hold["top100"])
say("\n  Hide-seeds test (hide %d of %d seeds, %d rounds, same hidden seeds as DIAMOnD.py):" % (n_hide, len(seeds), N_REPEATS))
say("    whole best Louvain community as the prediction: %.3f +- %.3f of hidden seeds found, but it holds about %d proteins" %
    (m_all, s_all, round(np.mean(hold["size"]))))
say("    best 100 of that community (most links to seeds): %.3f +- %.3f" % (m_top, s_top))
if diamond_recall:
    say("    DIAMOnD alpha=1, 100 picks (from holdout_summary.csv): %.3f +- %.3f" % diamond_recall)
write_csv("course_holdout_communities.csv", ["method", "recall_mean", "recall_se", "proteins_predicted"],
          [["whole best Louvain community", round(m_all, 4), round(s_all, 4), round(np.mean(hold["size"]))],
           ["best 100 of the community (links to seeds)", round(m_top, 4), round(s_top, 4), N_PICKS],
           ["DIAMOnD alpha=1 (100 picks)", diamond_recall[0] if diamond_recall else "", diamond_recall[1] if diamond_recall else "", N_PICKS]])

with open(os.path.join(RESULTS_DIR, "course_comparison_log.txt"), "w") as f:
    f.write("\n".join(log_lines) + "\n")
say("\nSaved: course_communities.csv, course_holdout_communities.csv, course_comparison_log.txt")

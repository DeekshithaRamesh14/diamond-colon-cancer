#! /usr/bin/env python


"""
# -----------------------------------------------------------------------
# encoding: utf-8
# DIAMOnD.py
# Joerg Menche, Susan D. Ghiassian
# Last Modified: 2020-22-09
# This code runs the DIAMOnD algorithm as described in
#
# A DIseAse MOdule Detection (DIAMOnD) Algorithm derived from a
# systematic analysis of connectivity patterns of disease proteins in
# the Human Interactome
#
# by Susan Dina Ghiassian, Joerg Menche & Albert-Laszlo Barabasi
#
#
# -----------------------------------------------------------------------
"""

import time
import networkx as nx
import numpy as np
import copy
import scipy.stats
from collections import defaultdict
import csv
import sys


# =============================================================================
def print_usage():

    print(' ')
    print('        usage: python3 DIAMOnD.py network_file seed_file n alpha(optional) outfile_name (optional)')
    print('        -----------------------------------------------------------------')
    print('        network_file : The edgelist must be provided as any delimiter-separated')
    print('                       table. Make sure the delimiter does not exit in gene IDs')
    print('                       and is consistent across the file.')
    print('                       The first two columns of the table will be')
    print('                       interpreted as an interaction gene1 <==> gene2')
    print('        seed_file    : table containing the seed genes (if table contains')
    print('                       more than one column they must be tab-separated;')
    print('                       the first column will be used only)')
    print('        n            : desired number of DIAMOnD genes, 200 is a reasonable')
    print('                       starting point.')
    print('        alpha        : an integer representing weight of the seeds,default')
    print('                       value is set to 1')
    print('        outfile_name : results will be saved under this file name')
    print('                       by default the outfile_name is set to "first_n_added_nodes_weight_alpha.txt"')
    print(' ')


# =============================================================================
def check_input_style(input_list):
    try:
        network_edgelist_file = input_list[1]
        seeds_file = input_list[2]
        max_number_of_added_nodes = int(input_list[3])
    # if no input is given, print out a usage message and exit
    except:
        print_usage()
        sys.exit(0)
        return

    alpha = 1
    outfile_name = 'first_%d_added_nodes_weight_%d.txt' % (max_number_of_added_nodes, alpha)

    if len(input_list) == 5:
        try:
            alpha = int(input_list[4])
            outfile_name = 'first_%d_added_weight_%d.txt' % (max_number_of_added_nodes, alpha)
        except:
            outfile_name = input_list[4]

    if len(input_list) == 6:
        try:
            alpha = int(input_list[4])
            outfile_name = input_list[5]
        except:
            print_usage()
            sys.exit(0)
            return
    return network_edgelist_file, seeds_file, max_number_of_added_nodes, alpha, outfile_name


# =============================================================================
def read_input(network_file, seed_file):
    """
    Reads the network and the list of seed genes from external files.
    * The edgelist must be provided as a tab-separated table. The
    first two columns of the table will be interpreted as an
    interaction gene1 <==> gene2
    * The seed genes mus be provided as a table. If the table has more
    than one column, they must be tab-separated. The first column will
    be used only.
    * Lines that start with '#' will be ignored in both cases
    """

    sniffer = csv.Sniffer()
    line_delimiter = None
    for line in open(network_file, 'r'):
        if line[0] == '#':
            continue
        else:
            dialect = sniffer.sniff(line)
            line_delimiter = dialect.delimiter
            break
    if line_delimiter == None:
        print
        'network_file format not correct'
        sys.exit(0)

    # read the network:
    G = nx.Graph()
    for line in open(network_file, 'r'):
        # lines starting with '#' will be ignored
        if line[0] == '#':
            continue
        # The first two columns in the line will be interpreted as an
        # interaction gene1 <=> gene2
        # line_data   = line.strip().split('\t')
        line_data = line.strip().split(line_delimiter)
        node1 = line_data[0]
        node2 = line_data[1]
        G.add_edge(node1, node2)

    # read the seed genes:
    seed_genes = set()
    for line in open(seed_file, 'r'):
        # lines starting with '#' will be ignored
        if line[0] == '#':
            continue
        # the first column in the line will be interpreted as a seed
        # gene:
        line_data = line.strip().split('\t')
        seed_gene = line_data[0]
        seed_genes.add(seed_gene)

    return G, seed_genes


# ================================================================================
def compute_all_gamma_ln(N):
    """
    precomputes all logarithmic gammas
    """
    gamma_ln = {}
    for i in range(1, N + 1):
        gamma_ln[i] = scipy.special.gammaln(i)

    return gamma_ln


# =============================================================================
def logchoose(n, k, gamma_ln):
    if n - k + 1 <= 0:
        return scipy.infty
    lgn1 = gamma_ln[n + 1]
    lgk1 = gamma_ln[k + 1]
    lgnk1 = gamma_ln[n - k + 1]
    return lgn1 - (lgnk1 + lgk1)   # CHANGED BY STUDENT: round brackets. Authors' line: return lgn1 - [lgnk1 + lgk1]  (crashes with NumPy 2)


# =============================================================================
def gauss_hypergeom(x, r, b, n, gamma_ln):
    return np.exp(logchoose(r, x, gamma_ln) +
                  logchoose(b, n - x, gamma_ln) -
                  logchoose(r + b, n, gamma_ln))


# =============================================================================
def pvalue(kb, k, N, s, gamma_ln):
    """
    -------------------------------------------------------------------
    Computes the p-value for a node that has kb out of k links to
    seeds, given that there's a total of s sees in a network of N nodes.

    p-val = \sum_{n=kb}^{k} HypergemetricPDF(n,k,N,s)
    -------------------------------------------------------------------
    """
    p = 0.0
    for n in range(kb, k + 1):
        if n > s:
            break
        prob = gauss_hypergeom(n, s, N - s, k, gamma_ln)
        # print prob
        p += prob

    if p > 1:
        return 1
    else:
        return p

    # =============================================================================


def get_neighbors_and_degrees(G):
    neighbors, all_degrees = {}, {}
    for node in G.nodes():
        nn = set(G.neighbors(node))
        neighbors[node] = nn
        all_degrees[node] = G.degree(node)

    return neighbors, all_degrees


# =============================================================================
# Reduce number of calculations
# =============================================================================
def reduce_not_in_cluster_nodes(all_degrees, neighbors, G, not_in_cluster, cluster_nodes, alpha):
    reduced_not_in_cluster = {}
    kb2k = defaultdict(dict)
    for node in not_in_cluster:

        k = all_degrees[node]
        kb = 0
        # Going through all neighbors and counting the number of module neighbors
        for neighbor in neighbors[node]:
            if neighbor in cluster_nodes:
                kb += 1

        # adding wights to the the edges connected to seeds
        k += (alpha - 1) * kb
        kb += (alpha - 1) * kb
        kb2k[kb][k] = node

    # Going to choose the node with largest kb, given k
    k2kb = defaultdict(dict)
    for kb, k2node in kb2k.items():
        min_k = min(k2node.keys())
        node = k2node[min_k]
        k2kb[min_k][kb] = node

    for k, kb2node in k2kb.items():
        max_kb = max(kb2node.keys())
        node = kb2node[max_kb]
        reduced_not_in_cluster[node] = (max_kb, k)

    return reduced_not_in_cluster


# ======================================================================================
#   C O R E    A L G O R I T H M
# ======================================================================================
def diamond_iteration_of_first_X_nodes(G, S, X, alpha):
    """
    Parameters:
    ----------
    - G:     graph
    - S:     seeds
    - X:     the number of iterations, i.e only the first X gened will be
             pulled in
    - alpha: seeds weight
    Returns:
    --------

    - added_nodes: ordered list of nodes in the order by which they
      are agglomerated. Each entry has 4 info:
      * name : dito
      * k    : degree of the node
      * kb   : number of +1 neighbors
      * p    : p-value at agglomeration
    """

    N = G.number_of_nodes()

    added_nodes = []

    # ------------------------------------------------------------------
    # Setting up dictionaries with all neighbor lists
    # and all degrees
    # ------------------------------------------------------------------
    neighbors, all_degrees = get_neighbors_and_degrees(G)

    # ------------------------------------------------------------------
    # Setting up initial set of nodes in cluster
    # ------------------------------------------------------------------

    cluster_nodes = set(S)
    not_in_cluster = set()
    s0 = len(cluster_nodes)

    s0 += (alpha - 1) * s0
    N += (alpha - 1) * s0

    # ------------------------------------------------------------------
    # precompute the logarithmic gamma functions
    # ------------------------------------------------------------------
    gamma_ln = compute_all_gamma_ln(N + 1)

    # ------------------------------------------------------------------
    # Setting initial set of nodes not in cluster
    # ------------------------------------------------------------------
    for node in cluster_nodes:
        not_in_cluster |= neighbors[node]
    not_in_cluster -= cluster_nodes

    # ------------------------------------------------------------------
    #
    # M A I N     L O O P
    #
    # ------------------------------------------------------------------

    all_p = {}

    while len(added_nodes) < X:

        # ------------------------------------------------------------------
        #
        # Going through all nodes that are not in the cluster yet and
        # record k, kb and p
        #
        # ------------------------------------------------------------------

        info = {}

        pmin = 10
        next_node = 'nix'
        reduced_not_in_cluster = reduce_not_in_cluster_nodes(all_degrees,
                                                             neighbors, G,
                                                             not_in_cluster,
                                                             cluster_nodes, alpha)

        for node, kbk in reduced_not_in_cluster.items():
            # Getting the p-value of this kb,k
            # combination and save it in all_p, so computing it only once!
            kb, k = kbk
            try:
                p = all_p[(k, kb, s0)]
            except KeyError:
                p = pvalue(kb, k, N, s0, gamma_ln)
                all_p[(k, kb, s0)] = p

            # recording the node with smallest p-value
            if p < pmin:
                pmin = p
                next_node = node

            info[node] = (k, kb, p)

        # ---------------------------------------------------------------------
        # Adding node with smallest p-value to the list of aaglomerated nodes
        # ---------------------------------------------------------------------
        added_nodes.append((next_node,
                            info[next_node][0],
                            info[next_node][1],
                            info[next_node][2]))

        # Updating the list of cluster nodes and s0
        cluster_nodes.add(next_node)
        s0 = len(cluster_nodes)
        not_in_cluster |= (neighbors[next_node] - cluster_nodes)
        not_in_cluster.remove(next_node)

    return added_nodes


# ===========================================================================
#
#   M A I N    D I A M O n D    A L G O R I T H M
#
# ===========================================================================
def DIAMOnD(G_original, seed_genes, max_number_of_added_nodes, alpha, outfile=None):
    """
    Runs the DIAMOnD algorithm
    Input:
    ------
     - G_original :
             The network
     - seed_genes :
             a set of seed genes
     - max_number_of_added_nodes:
             after how many added nodes should the algorithm stop
     - alpha:
             given weight to the sees
     - outfile:
             filename for the output generates by the algorithm,
             if not given the program will name it 'first_x_added_nodes.txt'
     Returns:
     --------
      - added_nodes: A list with 4 entries at each element:
            * name : name of the node
            * k    : degree of the node
            * kb   : number of neighbors that are part of the module (at agglomeration)
            * p    : connectivity p-value at agglomeration
      -
    """

    # 1. throwing away the seed genes that are not in the network
    all_genes_in_network = set(G_original.nodes())
    seed_genes = set(seed_genes)
    disease_genes = seed_genes & all_genes_in_network

    if len(disease_genes) != len(seed_genes):
        print("DIAMOnD(): ignoring %s of %s seed genes that are not in the network" % (
            len(seed_genes - all_genes_in_network), len(seed_genes)))

    # 2. agglomeration algorithm.
    added_nodes = diamond_iteration_of_first_X_nodes(G_original,
                                                     disease_genes,
                                                     max_number_of_added_nodes, alpha)
    # 3. saving the results
    with open(outfile, 'w') as fout:

        fout.write('\t'.join(['#rank', 'DIAMOnD_node', 'p_hyper']) + '\n')
        rank = 0
        for DIAMOnD_node_info in added_nodes:
            rank += 1
            DIAMOnD_node = DIAMOnD_node_info[0]
            p = float(DIAMOnD_node_info[3])

            fout.write('\t'.join(map(str, ([rank, DIAMOnD_node, p]))) + '\n')

    return added_nodes


# ===========================================================================
#
# "Hey Ho, Let's go!" -- The Ramones (1976)
#
# ===========================================================================


# ===========================================================================
#
#   A D D E D   B Y   T H E   S T U D E N T
#
#   Everything between this banner and the line   if __name__ == '__main__':
#   is new, and so are the first lines inside that block.
#   The authors' code above is not changed, except ONE marked line in logchoose() (a NumPy fix).
#
#   The data files are made by a separate program (prepare_data.py), then fed to the authors' code:
#
#     python prepare_data.py                                       makes data/network.txt and data/seeds.txt
#     python DIAMOnD.py data/network.txt data/seeds.txt 100 1 out.txt   the authors' own command (unchanged)
#     python DIAMOnD.py data/network.txt data/seeds.txt --analysis      the three changes below; results go to results/
#
#   The three changes (assignment):
#     1. new dataset      colon cancer with curated seed genes instead of the authors' own gene lists (made by prepare_data.py)
#     2. parameter sweep  DIAMOnD with five values of alpha, judged by a hide-seeds test
#     3. course method    null models of practical P3 (Erdos-Renyi and Maslov-Sneppen)
#
# ===========================================================================
import os                                   # folders and file paths
import random                               # random choices (seeds to hide, shuffling of the network)
import subprocess                           # to restart the program with a fixed hash seed
import igraph as ig                         # igraph: builds and shuffles networks, as in practical P3

if not hasattr(scipy, "infty"):             # new SciPy removed scipy.infty, which the authors' logchoose() uses
    scipy.infty = np.inf                    # give the old name back (same value)

# ----------------------------- settings ------------------------------------
ALPHAS = [1, 2, 3, 5, 10]                   # the 5 seed weights we compare (1 = the authors' default)
N_ADD = 100                                 # how many proteins DIAMOnD adds in every run
KS = [25, 50, 100]                          # we look at the first 25, 50 and 100 proteins added
HIDE_FRACTION = 0.25                        # share of the seeds that is hidden in the hide-seeds test
N_REPEATS = 30                              # number of random hide-seeds rounds (use 5 for a quick test)
N_NULL = 1000                               # random networks per null model, as B = 1000 in practical P3 (use 20 for a quick test)
SWAPS_PER_LINK = 10                         # swap attempts per link when a network is shuffled, as 10 x m in practical P3
RNG_SEED = 2026                             # fixed random seed, so every run gives the same numbers
HERE = os.path.dirname(os.path.abspath(__file__))     # the folder this file is in
RESULTS_DIR = os.path.join(HERE, "results")           # tables and figures are saved here


# =============================================================================
def ensure_fixed_hash_seed():
    """
    When several proteins have exactly the same p-value, the authors' code keeps one of them.
    Which one depends on Python's ordering of sets, and that ordering changes between runs unless
    PYTHONHASHSEED is fixed BEFORE python starts. So: if it is not fixed, run the same command again with it fixed.
    """
    if os.environ.get("PYTHONHASHSEED") != "0":                              # is the hash seed already fixed?
        env = dict(os.environ, PYTHONHASHSEED="0")                           # copy the environment and fix the seed
        code = subprocess.call([sys.executable, "-W", "ignore::SyntaxWarning"] + sys.argv, env=env)   # run the same command again
        sys.exit(code)                                                       # this first copy stops here


def load_symbols(path):
    """gene ID -> gene name, read from the table that prepare_data.py saved (empty if the file is missing)"""
    symbols = {}                                                             # the lookup table
    if os.path.exists(path):                                                 # only if the file exists
        for line in open(path):                                              # one gene per line
            gene_id, name = line.rstrip("\n").split("\t")                    # "7157<TAB>TP53"
            symbols[gene_id] = name                                          # remember the name
    return symbols


def mean_and_se(values):
    """average and standard error (= standard deviation / sqrt(number of values))"""
    v = np.array(values, dtype=float)                                     # list -> numpy array
    return v.mean(), v.std(ddof=1) / np.sqrt(len(v))                      # mean and standard error


def write_csv(name, header, rows):
    """save a table in results/"""
    os.makedirs(RESULTS_DIR, exist_ok=True)                               # create results/ if needed
    with open(os.path.join(RESULTS_DIR, name), "w", newline="") as f:     # open the file
        writer = csv.writer(f)                                            # csv writer
        writer.writerow(header)                                           # first line = column names
        writer.writerows(rows)                                            # the rows


# =============================================================================
def run_alpha_sweep(G, seeds, symbols):
    """parameter sweep: DIAMOnD() for alpha = 1, 2, 3, 5, 10, first on all seeds, then with hidden seeds"""
    N = G.number_of_nodes()                                               # number of proteins in the network

    # ---- Part A: all seeds, how much does the list change with alpha? ----
    print("\nPART A: all %d seeds, alpha = %s" % (len(seeds), ALPHAS))     # progress message
    runs = {a: DIAMOnD(G, set(seeds), N_ADD, a, outfile=os.devnull) for a in ALPHAS}   # the authors' DIAMOnD(), once per alpha
    first = [x[0] for x in runs[1]]                                       # genes predicted with alpha = 1 (the default)
    overlap_rows = []                                                     # rows of the overlap table
    for a in ALPHAS:                                                      # for each alpha ...
        genes = [x[0] for x in runs[a]]                                   # ... the predicted genes
        shared = len(set(genes) & set(first))                             # how many are also in the alpha = 1 list
        zeros = sum(float(x[3]) == 0 for x in runs[a])                    # how many picks have a p-value of exactly 0
        overlap_rows.append([a, shared, zeros])                           # store the numbers
        print("  alpha=%-2d  genes shared with alpha=1: %3d of 100 | picks with p-value exactly 0: %3d" % (a, shared, zeros))
    write_csv("overlap_with_alpha1.csv", ["alpha", "genes_shared_with_alpha1_of_100", "picks_with_p_exactly_0"], overlap_rows)
    write_csv("top100_per_alpha.csv", ["rank"] + ["alpha_%d_gene" % a for a in ALPHAS] + ["alpha_1_p_value"],
              [[i + 1] + [symbols.get(runs[a][i][0], runs[a][i][0]) for a in ALPHAS] + [float(runs[1][i][3])] for i in range(N_ADD)])
    print("  top 10 for alpha=1:", ", ".join(symbols.get(g, g) for g in first[:10]))   # a quick look

    # ---- Part B: hide seeds, how many come back? ----
    # Why this belongs to the parameter change: alpha only matters through the seeds, so to judge which alpha is better we need a test with a
    # known answer. We give DIAMOnD only 75% of the seeds (the input 'seeds' is changed) and count how many of the hidden 25% it finds again.
    # The hiding is NOT a parameter of DIAMOnD; it is the yardstick (the authors did the same with 25%, 50%, 75% hidden, Fig. 3E-F).
    random.seed(RNG_SEED)                                                 # fixed random seed
    n_hide = round(HIDE_FRACTION * len(seeds))                            # how many seeds to hide (for example 5 of 21)
    print("\nPART B: hide %d seeds, %d rounds, %d DIAMOnD runs..." % (n_hide, N_REPEATS, N_REPEATS * len(ALPHAS)))
    methods = ["alpha=%d" % a for a in ALPHAS] + ["most seed neighbours"]  # DIAMOnD for each alpha + one simple comparison
    recall = {(m, k): [] for m in methods for k in KS}                    # recall for each method and k, one value per round
    curves = {m: [] for m in methods}                                     # recall after 1, 2, ..., N_ADD proteins, one curve per round
    ceiling = []                                                          # share of hidden seeds with a direct link to a visible seed
    for rep in range(N_REPEATS):                                          # repeat with different hidden seeds
        hidden = random.sample(seeds, n_hide)                             # choose the seeds to hide
        visible = [s for s in seeds if s not in hidden]                   # the seeds DIAMOnD may see
        ceiling.append(sum(bool(set(G.neighbors(s)) & set(visible)) for s in hidden) / n_hide)   # ceiling check
        for a in ALPHAS:                                                  # for each alpha ...
            genes = [x[0] for x in DIAMOnD(G, set(visible), N_ADD, a, outfile=os.devnull)]   # ... run DIAMOnD on the visible seeds
            curves["alpha=%d" % a].append(np.cumsum([g in hidden for g in genes]) / n_hide)   # recall curve of this round
            for k in KS:                                                  # for each list length ...
                found = len(set(genes[:k]) & set(hidden))                 # ... how many hidden seeds are in the first k
                recall["alpha=%d" % a, k].append(found / n_hide)          # recall = share of hidden seeds found
        count = {}                                                        # simple comparison: links to the visible seeds
        for s in visible:                                                 # for each visible seed ...
            for nb in G.neighbors(s):                                     # ... and each of its neighbours
                if nb not in visible:                                     # (not the visible seeds themselves)
                    count[nb] = count.get(nb, 0) + 1                      # one more link from a seed
        ranked = sorted(count, key=lambda g: (-count[g], g))              # most links first, ties by gene id
        curves["most seed neighbours"].append(np.cumsum([g in hidden for g in ranked[:N_ADD]]) / n_hide)   # recall curve
        for k in KS:                                                      # for each list length ...
            recall["most seed neighbours", k].append(len(set(ranked[:k]) & set(hidden)) / n_hide)   # ... same recall measure
        print("  round %d of %d done" % (rep + 1, N_REPEATS))             # progress message
    random_recall = {k: k / (N - (len(seeds) - n_hide)) for k in KS}      # chance level: k random proteins

    print("\n  recall = share of the hidden seeds found in the first k predictions (mean +- standard error)")
    print("  %-22s" % "method" + "".join("   k=%-3d        " % k for k in KS) + "  diff vs alpha=1 (k=50)")
    summary = []                                                          # rows for the csv file
    for m in methods:                                                     # one line per method
        line = "  %-22s" % m                                              # start of the printed line
        for k in KS:                                                      # one column per k
            mean, se = mean_and_se(recall[m, k])                          # average over the rounds
            line += "   %.3f +- %.3f" % (mean, se)                        # add it to the line
            summary.append([m, k, round(mean, 4), round(se, 4)])          # and to the csv rows
        d_mean, d_se = mean_and_se(np.array(recall[m, 50]) - np.array(recall["alpha=1", 50]))   # difference to alpha = 1, same rounds
        print(line + ("   %+.3f +- %.3f" % (d_mean, d_se) if m != "alpha=1" else "   -"))
    print("  %-22s" % "random guessing" + "".join("   %.3f         " % random_recall[k] for k in KS))
    c_mean, c_se = mean_and_se(ceiling)                                   # average over the rounds
    print("  ceiling: hidden seeds with a direct link to a visible seed: %.3f +- %.3f" % (c_mean, c_se))
    summary += [["random guessing", k, round(random_recall[k], 4), ""] for k in KS]   # chance level rows
    summary.append(["ceiling: hidden seeds with a direct link to a visible seed", "", round(c_mean, 4), round(c_se, 4)])
    write_csv("holdout_summary.csv", ["method", "k", "mean_recall", "standard_error"], summary)
    return curves, N - (len(seeds) - n_hide)                              # curves and number of proteins a random guess can pick from


# =============================================================================
# Course method: the null models of practical P3 (ER and Maslov-Sneppen), applied to the seeds
def z_score(observed, null):
    """z = (observed - mean of the null) / standard deviation of the null  (same formula as practical P3)"""
    null = np.asarray(null, dtype=float)                                  # list -> numpy array
    return (observed - null.mean()) / null.std(ddof=1)                    # ddof=1: sample standard deviation


def empirical_p(observed, null):
    """p = (1 + number of random networks with a value >= the observed one) / (B + 1)  (same formula as practical P3)"""
    null = np.asarray(null, dtype=float)                                  # list -> numpy array
    return (1 + np.sum(null >= observed)) / (len(null) + 1)               # one-sided: "as high or higher than observed"


def run_null_model(G, seeds):
    """
    Practical P3 asked: is a statistic of the real network surprising, compared with MANY random networks?
    Here the statistics describe the disease genes, and the same two null models are used:
      ER  (Erdos-Renyi, same number of proteins and links, everything else random)  -> is it explained by size and density alone?
      M-S (Maslov-Sneppen, links swapped, every protein keeps its number of links) -> is it explained by the degrees alone?
    """
    g = ig.Graph.TupleList(G.edges(), directed=False)                     # the real network as an igraph graph
    n, m = g.vcount(), g.ecount()                                         # number of proteins and of links
    names = g.vs["name"]                                                  # igraph numbers the proteins 0, 1, 2 ...; this list keeps their gene IDs
    seed_ids = [names.index(s) for s in seeds]                            # the igraph numbers of the seeds
    print("\nNULL MODELS: %d random networks per null model (ER and Maslov-Sneppen)" % N_NULL)

    def statistics(h):
        """the two numbers measured on one network h"""
        links = h.induced_subgraph(seed_ids).ecount()                     # 1. links between two seeds: do the disease genes touch each other?
        net = nx.Graph([(names[a], names[b]) for a, b in h.get_edgelist()])   # same network for the authors' DIAMOnD() (needs networkx)
        picks = DIAMOnD(net, set(seeds), 10, 1, outfile=os.devnull)       # the authors' DIAMOnD(): first 10 proteins, alpha = 1
        strength = np.median([-np.log10(max(float(p[3]), 1e-300)) for p in picks])   # 2. typical -log10(p-value) of those 10 (bigger = more significant)
        return [links, strength]

    observed = statistics(g)                                              # the real network (like "obs" in practical P3)
    null = {}                                                             # the random networks' numbers, for "ER" and "M-S"
    for kind in ("ER", "M-S"):                                            # the two null models
        random.seed(RNG_SEED)                                             # fixed random seed (igraph uses Python's random)
        rows = []                                                         # one line per random network
        for b in range(N_NULL):                                           # B random networks (an ensemble, not just one)
            if kind == "ER":
                h = ig.Graph.Erdos_Renyi(n=n, m=m)                        # ER: same n and m, links placed at random
            else:
                h = g.copy()                                              # M-S: start from the real network ...
                h.rewire(n=SWAPS_PER_LINK * m)                            # ... and swap links (degrees stay the same), about 10 x m attempts
            rows.append(statistics(h))                                    # measure the random network
            if (b + 1) % 100 == 0:                                        # a progress message every 100 networks
                print("  %-3s %d of %d done" % (kind, b + 1, N_NULL))
        null[kind] = np.array(rows)                                       # columns: links between seeds, strength of the first 10 picks

    labels = ["links between seeds", "strength of first 10 picks"]         # names of the two statistics
    print("\n  %-28s %9s | %-22s %8s %7s | %-22s %8s %7s" % ("statistic", "observed", "ER mean +- sd", "z", "p", "M-S mean +- sd", "z", "p"))
    table = []                                                            # rows for the csv file
    for j, label in enumerate(labels):                                    # one line per statistic
        line = [label, round(observed[j], 3)]                             # start of the table row
        text = "  %-28s %9.2f" % (label, observed[j])                     # start of the printed line
        for kind in ("ER", "M-S"):                                        # the two null models
            values = null[kind][:, j]                                     # this statistic in all random networks
            z, p = z_score(observed[j], values), empirical_p(observed[j], values)
            line += [round(values.mean(), 3), round(values.std(ddof=1), 3), round(z, 2), round(p, 4)]
            text += " | %8.2f +- %-9.2f %8.1f %7.4f" % (values.mean(), values.std(ddof=1), z, p)
        table.append(line)
        print(text)
    write_csv("null_model_comparison.csv",
              ["statistic", "observed", "ER_mean", "ER_sd", "ER_z", "ER_p", "MS_mean", "MS_sd", "MS_z", "MS_p"], table)
    return observed, null, labels


# =============================================================================
def make_figures(curves, pool, observed, null, labels):
    """two figures in results/: recall against iteration (like the authors' recall plots) and the null-model histograms (practical P3, part 5)"""
    import matplotlib                                                     # plotting library (only needed here)
    matplotlib.use("Agg")                                                 # draw into files, no window
    import matplotlib.pyplot as plt                                       # plotting functions
    x = np.arange(1, N_ADD + 1)                                           # iteration = number of proteins added so far
    ink = "#0b0b0b"                                                       # colour of the text and axes

    # ---- Figure 1: recall against iteration for the five alpha values (parameter sweep) ----
    fig, ax = plt.subplots(figsize=(6.5, 4.5))                            # one panel
    colours = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]     # one colour per alpha
    for a, colour, marker in zip(ALPHAS, colours, "osD^v"):               # one line per alpha
        mean = np.mean(curves["alpha=%d" % a], axis=0)                    # average over the hide-seeds rounds
        ax.plot(x, mean, color=colour, marker=marker, markevery=10, ms=5, lw=1.6, label="DIAMOnD alpha = %d" % a)
    ax.plot(x, np.mean(curves["most seed neighbours"], axis=0), color="gray", ls="--", lw=1.4, label="most seed neighbours")   # simple comparison
    ax.plot(x, x / pool, color="gray", ls=":", lw=1.4, label="random guessing")   # chance level
    ax.set_xlabel("iteration (proteins added)", color=ink)                # x-axis title
    ax.set_ylabel("recall of the hidden seeds", color=ink)                # y-axis title
    ax.set_title("Colon cancer: %d%% of the seeds hidden, mean of %d rounds" % (HIDE_FRACTION * 100, N_REPEATS), fontsize=10)
    ax.set_ylim(bottom=0)                                                 # start at zero
    ax.spines[["top", "right"]].set_visible(False)                        # clean frame
    ax.set_xlim(1, N_ADD)                                                 # no empty margin
    ax.legend(frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(1.02, 1))   # legend with all lines, outside the plot
    fig.savefig(os.path.join(RESULTS_DIR, "fig_alpha_sweep_recall.png"), dpi=200, bbox_inches="tight")   # save
    plt.close(fig)                                                        # free the memory

    # ---- Figure 2: observed value (red line) over the histogram of the random networks, as in practical P3 part 5 ----
    fig, axes = plt.subplots(2, 2, figsize=(9, 6))                        # rows = the two statistics, columns = the two null models
    for j, label in enumerate(labels):                                    # one row per statistic
        for k, kind in enumerate(("ER", "M-S")):                          # one column per null model
            ax, values = axes[j, k], null[kind][:, j]                     # this panel and its random-network values
            if j == 0:                                                    # links between seeds are whole numbers: one bar per number
                bins = np.arange(-0.5, max(values.max(), observed[j]) + 1.5)
            else:                                                         # strength is continuous: 30 bars
                bins = 30
            ax.hist(values, bins=bins, color="#b5b4ad")                   # histogram of the random networks
            ax.axvline(observed[j], color="#e34948", lw=2)                # the real network = red line
            ax.set_title("%s: %s null" % (label, kind), fontsize=10)      # panel title
            ax.set_xlabel(label, color=ink)                               # x-axis title
            ax.set_ylabel("random networks", color=ink)                   # y-axis title
            ax.text(0.55, 0.92, "red line = real network (%.2f), z = %.1f" % (observed[j], z_score(observed[j], values)),
                    transform=ax.transAxes, ha="center", va="top", fontsize=8, color="#b02a29")   # label of the red line, away from the bars
            ax.spines[["top", "right"]].set_visible(False)                # clean frame
    fig.tight_layout()                                                    # avoid overlaps
    fig.savefig(os.path.join(RESULTS_DIR, "fig_null_model.png"), dpi=200, bbox_inches="tight")   # save
    plt.close(fig)                                                        # free the memory


# =============================================================================
def run_analysis(network_file, seed_file):
    """the three changes in one go: alpha sweep + hide-seeds test, null models, figures"""
    G, seed_set = read_input(network_file, seed_file)                     # the authors' function reads the two files
    seeds = sorted(s for s in seed_set if s in G)                         # seeds that are in the network (sorted = same order every run)
    print("network: %d proteins, %d links | seeds: %d" % (G.number_of_nodes(), G.number_of_edges(), len(seeds)))
    symbols = load_symbols(os.path.join(os.path.dirname(network_file), "gene_symbols.tsv"))   # gene names, if available
    curves, pool = run_alpha_sweep(G, seeds, symbols)                     # change 2: alpha sweep and hide-seeds test
    observed, null, labels = run_null_model(G, seeds)                     # change 3: null models of practical P3
    make_figures(curves, pool, observed, null, labels)                    # figures
    print("\nresults saved in", RESULTS_DIR)                              # done


if __name__ == '__main__':
    # ---- added by the student: fixed hash seed and the new mode (the authors' code continues below) ----
    ensure_fixed_hash_seed()                                       # every run gives the same results
    if len(sys.argv) > 3 and sys.argv[3] == '--analysis':          # python DIAMOnD.py network.txt seeds.txt --analysis
        run_analysis(sys.argv[1], sys.argv[2])                     # alpha sweep, hide-seeds test, null models
        sys.exit(0)                                                # then stop
    # ---- end of the added lines ----
    # -----------------------------------------------------
    # Checking for input from the command line:
    # -----------------------------------------------------
    #
    # [1] file providing the network in the form of an edgelist
    #     (tab-separated table, columns 1 & 2 will be used)
    #
    # [2] file with the seed genes (if table contains more than one
    #     column they must be tab-separated; the first column will be
    #     used only)
    #
    # [3] number of desired iterations
    #
    # [4] (optional) seeds weight (integer), default value is 1
    # [5] (optional) name for the results file

    # check if input style is correct
    input_list = sys.argv
    network_edgelist_file, seeds_file, max_number_of_added_nodes, alpha, outfile_name = check_input_style(input_list)

    # read the network and the seed genes:
    G_original, seed_genes = read_input(network_edgelist_file, seeds_file)

    # run DIAMOnD
    added_nodes = DIAMOnD(G_original,
                          seed_genes,
                          max_number_of_added_nodes, alpha,
                          outfile=outfile_name)

    print("\n results have been saved to '%s' \n" % outfile_name)

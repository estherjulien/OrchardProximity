import networkx as nx
import numpy as np
import argparse
import random
import pickle
import time
import os
'''
Code for generating birth-hybridization network
- Apdated from Remie Janssen's code for the paper (with help of Celine Scornavacca for distances on arcs):
Janssen, R., & Liu, P. (2021). Comparing the topology of phylogenetic network generators. 
Journal of Bioinformatics and Computational Biology, 19(06), 2140012.

- Originally, the method is from the paper:
Zhang, C., Ogilvie, H.A., Drummond, A.J., Stadler, T.: Bayesian inference of species networks from multilocus
sequence data. Molecular biology and evolution 35(2), 504–517 (2018)
'''


def get_reticulations(net):
    return {v for v in net.nodes() if net.in_degree(v) == 2}


def get_leaves(net):
    return {u for u in net.nodes() if net.out_degree(u) == 0}


def birth_hyb(time_limit, speciation_rate, hybridization_rate, rng, taxa_goal=None, max_retics=None, distances=True):
    nw = nx.DiGraph()
    nw.add_node(0)
    leaves = {0}
    current_node = 1
    no_of_leaves = 1
    retics = 0
    if max_retics == -1:
        max_retics = None
    extra_time = rng.exponential(1/float(speciation_rate))
    current_time = extra_time
    current_speciation_rate = float(speciation_rate)
    current_hybridization_rate = hybridization_rate
    rate = current_speciation_rate + current_hybridization_rate

    while current_time < time_limit and (taxa_goal is not None and taxa_goal != no_of_leaves):
        # if rng.random() < stopping_rate and :
        #     break
        if (0 in leaves) or (rng.random() < current_speciation_rate / rate):
            # Speciate
            splitting_leaf = random.choice(list(leaves))
            if distances:
                nw.add_weighted_edges_from([(splitting_leaf, current_node, 0), (splitting_leaf, current_node + 1, 0)],
                                           weight='length')
            else:
                nw.add_edges_from([(splitting_leaf, current_node), (splitting_leaf, current_node + 1)])
            leaves.remove(splitting_leaf)
            leaves.add(current_node)
            leaves.add(current_node+1)
            current_node += 2
            no_of_leaves += 1
        elif len(leaves) >= 2:
            # Hybridize
            merging = rng.choice(tuple(leaves), 2, replace=False)
            l0 = merging[0]
            l1 = merging[1]
            pl0 = -1
            for p in nw.predecessors(l0):
                pl0 = p
            pl1 = -1
            for p in nw.predecessors(l1):
                pl1 = p
            # If pl0==pl1, the new hybridization results in parallel edges.
            if pl0 != pl1:
                nw.remove_node(l1)
                if distances:
                    nw.add_weighted_edges_from([(pl1, l0, 0), (l0, current_node, 0)], weight='length')
                else:
                    nw.add_edges_from([(pl1, l0), (l0, current_node)])
                leaves.remove(l0)
                leaves.remove(l1)
                leaves.add(current_node)
                current_node += 1
                no_of_leaves -= 1
                retics += 1
        else:
            break

        # extend all pendant edges
        if distances:
            for l in leaves:
                pl = -1
                for p in nw.predecessors(l):
                    pl = p
                nw[pl][l]['length'] += extra_time
        current_speciation_rate = float(speciation_rate*no_of_leaves)
        current_hybridization_rate = float(hybridization_rate * (no_of_leaves * (no_of_leaves - 1))/2)
        rate = current_speciation_rate + current_hybridization_rate
        extra_time = rng.exponential(1/rate)  # REMIE'S CODE THIS WAS 1/rate
        current_time += extra_time
        if max_retics is not None and retics > max_retics:
            return None, retics, no_of_leaves

    # nothing has happened yet, and there is only one node
    if len(nw) == 1:
        nw.add_edges_from([(0, 1)])
    return nw, retics, no_of_leaves


def make_zods(l, retic, seed):
    net_info = f"L{l}_R{retic}_ZODS"

    st = time.time()
    # MAKE NETWORK
    s_rate = 1.0
    rng = np.random.RandomState(seed)

    print(f"JOB {seed}: Start creating NETWORK (ZODS, L = {l}, R = {retic})")
    net = False
    while True:
        h_rate = rng.uniform(0.0001, 0.4)
        net, ret_num, num_leaves = birth_hyb(500, s_rate, h_rate, rng, taxa_goal=l, max_retics=retic)
        if net is not None and len(get_leaves(net)) == l and ret_num == retic:
            break
    if not net:
        print(f"\nJOB {seed}: FAILED {net_info}\n")
        return False

    print(f"JOB {seed}: finished creating network in {np.round(time.time() - st, 2)}s (ZODS, L = {l}, R = {retic})")
    # SAVE INSTANCE AS .el FILE
    os.makedirs("simulatedNetworks", exist_ok=True)
    file = f"simulatedNetworks/n{seed}_{net_info}.el"
    f = open(file, "w")
    for x, y in net.edges:
        f.write(f"{x} {y}\n")
    f.close()
    return True


def make_network_from_el_simulated(inst_num, l, retic):
    net_info = f"L{l}_R{retic}_ZODS"
    file = f"simulatedNetworks/n{inst_num}_{net_info}.el"
    # open edges file
    with open(file, "r") as handle:
        edges_file = handle.readlines()
    edges = np.loadtxt(edges_file, dtype=str)

    # add edges
    nw = nx.DiGraph()
    nw.add_edges_from(edges)

    # save nw
    os.makedirs("data/simulated_network", exist_ok=True)
    with open(f"data/simulated_network/network_{net_info}_{inst_num}.pkl", "wb") as handle:
        pickle.dump(nw, handle)

    return True


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-l', '--leaves',
        help='Number of leaves.',
        default=10,
        type=int
    )
    parser.add_argument(
        '-r', '--retic',
        help='Number of reticulations.',
        type=int,
        default=0,
    )
    parser.add_argument(
        '-s', '--seed',
        help='Seed number.',
        type=int,
        default=0,
    )
    args = parser.parse_args()
    # args = argparse.Namespace(leaves=10, retic=5, seed=0)

    make_zods(args.leaves, args.retic, args.seed)

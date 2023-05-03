import networkx as nx
import numpy as np
import argparse
import pickle
import os


def make_network_from_el(inst_num):
    # open edges file
    file = f"realNetworks/n{inst_num}.el"
    with open(file, "r") as handle:
        edges_file = handle.readlines()
    edges = np.loadtxt(edges_file, dtype=str)

    # add edges
    nw = nx.DiGraph()
    nw.add_edges_from(edges)

    # save nw
    os.makedirs(f"data/real_network", exist_ok=True)
    save_file = f"data/real_network/n{inst_num}_nw.pkl"
    with open(save_file, "wb") as handle:
        pickle.dump(nw, handle)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'instance',
        help='Instance number.',
        type=int,
        default=1,
    )
    args = parser.parse_args()

    make_network_from_el(args.instance)
from network_gen import make_network_from_el_simulated, get_leaves, get_reticulations
from real_network_reader import make_network_from_el
from pyscipopt import Model, SCIP_PARAMSETTING

import numpy as np
import argparse
import pickle
import os
'''
Code for writing and solving the L-distance problem.
'''


def write_lp_get_l_or(network_file, problem_file):
    """
    Code for writing an lp file using the PySCIPOpt interface for a given network.
    :param network_file: path to a networkx pkl file with the network
    :param problem_file: path to the lp problem file
    """

    # GET INSTANCE
    # get information on simulated_network
    with open(network_file, "rb") as handle:
        network = pickle.load(handle)

    # PREPROCESSING
    # set of vertices = simulated_network.nodes
    num_vers = len(network.nodes)
    # set of leaves
    lvs = get_leaves(network)
    # set of reticulations
    retics = get_reticulations(network)
    retic_arcs = []
    for v in retics:
        retic_arcs += [(p, v) for p in network.predecessors(v)]
    retic_arcs = set(retic_arcs)

    # write LP model
    m = Model("write")

    # create variables
    x = dict()
    for u, v in network.edges:
        if (u, v) not in retic_arcs:
            x[(u, v)] = 1
            continue
        x[(u, v)] = m.addVar(f"x_{u}_{v}", vtype="INTEGER", lb=0, ub=1)

    h = dict()
    for v in retics:
        h[v] = m.addVar(f"h_{v}", vtype="INTEGER", lb=0, ub=1)

    l = dict()
    for v in network.nodes:
        l[v] = m.addVar(f"l_{v}", vtype="CONTINUOUS", lb=0)

    # set objective
    m.setObjective(sum(h[v] for v in retics), sense="minimize")

    # CONSTRAINTS FOR VERTICAL ARCS
    # link h to x
    for v in retics:
        m.addCons(sum(x[u, v] for u in network.predecessors(v)) - 1 <= h[v])

    for u in network.nodes:
        if u in lvs:
            continue
        try:
            m.addCons(sum(x[u, v] for v in network.successors(u)) >= 1)
        except AssertionError:
            if all([x[u, v] == 1 for v in network.successors(u)]):
                continue
            else:
                assert True

    for v in retics:
        m.addCons(sum(x[u, v] for u in network.predecessors(v)) >= 1)

    # CONSTRAINTS FOR VERTEX LABELS
    m.addConss(l[u] <= l[v] for (u, v) in network.edges)

    # force difference of label if a tree arc
    for v in network.nodes:
        if v in retics:
            continue
        m.addConss(l[u] <= l[v] - 1 for u in network.predecessors(v))

    # CONSTRAINTS TO LINK x AND l
    for v in retics:
        m.addConss(l[u] <= l[v] - 1 + num_vers * (1 - x[u, v]) for u in network.predecessors(v))
        m.addConss(l[u] >= l[v] - num_vers*x[u, v] for u in network.predecessors(v))

    # write to LP file
    m.writeProblem(problem_file)


def solve_problem(problem, sol_name, time_limit=-1):
    '''
    Solve the given MILP problem using SCIP and output the solution and runtime.
    :param problem: path to lp file of problem
    :param sol_name: path to solution file
    :param time_limit: time limit for solving the MILP in seconds.
    '''
    m = Model("model")
    m.readProblem(problem)
    m.hideOutput()
    m.setPresolve(SCIP_PARAMSETTING.OFF)
    if time_limit > 0:
        m.setParam('limits/time', time_limit)

    m.optimize()
    if m.getStatus() == "timelimit":
        print(f"NOT OPTIMAL")
        f = open(sol_name, "w")
        f.write("TIME OUT\n")
        f.write("gap:\n")
        f.write(str(m.getGap()))
        f.close()
        return None
    else:
        sol_time = m.getSolvingTime()
        print(f"solved in {np.round(sol_time, 3)} sec")
    # get solution
    m.writeBestSol(sol_name)

    # print solve time in file
    f = open(sol_name, "a")
    f.write("runtime (in sec): \n")
    f.write(str(sol_time))
    f.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-real',
        help='Choose instance type',
        type=int,
        choices=[0, 1],
        default=1
    )
    parser.add_argument(
        '-s', '--seed',
        help='Random generator seed.',
        type=int,
        default=0,
    )
    parser.add_argument(
        '-w', '--write_lp',
        help='Write LP file',
        type=int,
        choices=[0, 1],
        default=1
    )
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
        default=-1
    )
    args = parser.parse_args()

    if args.real:
        network_file = f"data/real_network/n{args.seed}_nw.pkl"
        problem_info = f"n{args.seed}_real"
        if not os.path.exists(network_file):
            print("GET REAL NETWORK")
            make_network_from_el(args.seed)
    else:
        network_file = f"data/simulated_network/network_L{args.leaves}_R{args.retic}_ZODS_{args.seed}.pkl"
        problem_info = f"L{args.leaves}_R{args.retic}_ZODS_{args.seed}"
        if not os.path.exists(network_file):
            print("GET SIMULATED NETWORK")
            make_network_from_el_simulated(args.seed, args.leaves, args.retic)

    network_dir = "data/get_L_or"
    problem_file = network_dir + f"/lp_file/problem_{problem_info}.lp"
    if args.write_lp:
        os.makedirs(network_dir + "/lp_file", exist_ok=True)
        problem_file = f"data/get_L_or/lp_file/problem_{problem_info}.lp"
        write_lp_get_l_or(network_file, problem_file)

    os.makedirs(network_dir + "/solution", exist_ok=True)
    sol_name = network_dir + f"/solution/solution_{problem_info}.txt"

    solve_problem(problem_file, sol_name)

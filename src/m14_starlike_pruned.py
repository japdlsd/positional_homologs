import itertools
import os.path

import sys

import argh
import networkx as nx
from gurobipy import GRB
import gurobipy as gp

# import helpers
from grouping import *
import model_utils
from m13_starlike import solve


def main(input_filename, output_directory, time_limit=300, threads=1, **kwargs):
    colors, cost = model_utils.load_task(input_filename)
    pruned_cost = prune_edges_single_cost(cost)

    grouping_result = solve(
        colors,
        pruned_cost,
        time_limit=time_limit,
        logfile=os.path.join(output_directory, "gurobi.log"),
        threads=threads,
    )
    model_utils.save_results(output_directory, grouping_result)


if __name__ == "__main__":
    argh.dispatch_command(main)

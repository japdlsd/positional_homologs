import itertools
import os.path

import sys

import argh
import pulp
import networkx as nx

# import helpers
from grouping import *
import model_utils
from m03b_gurobi_native_components import solve


def main(input_filename, output_directory, time_limit=300, threads=1, **kwargs):
    data = model_utils.load_task_with_metadata(input_filename)
    colors = data["colors"]
    cost = data["cost"]
    del data["colors"]
    del data["cost"]

    pruned_cost = prune_edges_single_cost(cost)

    grouping_result = solve(
        colors,
        pruned_cost,
        time_limit=time_limit,
        logfile=os.path.join(output_directory, "gurobi.log"),
        threads=threads,
    )
    model_utils.save_results(output_directory, grouping_result, **data)


if __name__ == "__main__":
    argh.dispatch_command(main)

import itertools
import os.path
import random

import sys
from collections import defaultdict
from enum import Enum

import argh
import gurobipy as gp
from gurobipy import GRB
import networkx as nx

# import helpers
from grouping import *
import model_utils


def main(input_filename, output_directory, time_limit=300, threads=1, **kwargs):
    colors, cost = model_utils.load_task(input_filename)
    grouping_result = solve(
        colors,
        cost,
        time_limit=time_limit,
        logfile=os.path.join(output_directory, "gurobi.log"),
        threads=threads,
    )
    model_utils.save_results(output_directory, grouping_result)


def solve(colors, cost, time_limit=300, logfile="gurobi.log", threads=1):
    if len(colors) == 0:
        return GroupingResult([], 0, None)
    if len(colors) == 1:
        return GroupingResult([1], 0, None)

    N = len(colors)
    # MAX_C is the largest number of vertices with the same color
    color_sizes = defaultdict(int)
    for color in colors:
        color_sizes[color] += 1
    MAX_C = max(color_sizes.values())

    model = gp.Model("grouping")
    model.Params.LogFile = logfile
    model.Params.Threads = threads
    model.Params.TimeLimit = time_limit

    # create the component integer variables with bounds 1..MAX_C (inclusive)
    components = {}
    for i in range(N):
        components[i] = model.addVar(vtype=GRB.INTEGER, lb=1, ub=MAX_C, name=f"components_{i}")

    # add the constraints that vertices with the same color must have different components
    for i in range(N):
        for j in range(i+1, N):
            if colors[i] == colors[j]:
                model.addQConstr(components[i]*components[i] - 2 * components[i] * components[j] + components[j] * components[j] >= 1)

    # set the values of the components for the largest color to be 1..MAX_C
    largest_color = max(color_sizes, key=color_sizes.get)
    c = 1
    for i in range(N):
        if colors[i] == largest_color:
            components[i].setAttr(GRB.Attr.LB, c)
            components[i].setAttr(GRB.Attr.UB, c)
            c += 1

    # create are_connected indicator variables
    are_connected = {}
    for i in range(N):
        for j in range(i + 1, N):
            are_connected[i, j] = model.addVar(vtype=GRB.BINARY, name=f"are_connected_{i}_{j}")

    # forbid to connect vertices with different components
    for i in range(N):
        for j in range(i + 1, N):
            model.addQConstr(MAX_C**2 * are_connected[i, j] + components[i]*components[i] - 2 * components[i] * components[j] + components[j] * components[j] <= MAX_C**2)

    # set the objective: maximize the total cost of the connections
    model.setObjective(
        gp.quicksum(cost[i, j] * are_connected[i, j] for i in range(N) for j in range(i + 1, N) if not np.isnan(cost[i, j])),
        GRB.MAXIMIZE,
    )

    # solve the model
    model.optimize()

    # extract the solution
    raw_grouping = [int(components[i].X) for i in range(N)]
    grouping = model_utils.relabel_grouping(raw_grouping)
    total_value = model.objVal

    return GroupingResult(grouping, total_value, None)



if __name__ == "__main__":
    argh.dispatch_command(main)

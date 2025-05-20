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


def main(input_filename, output_directory, time_limit=300, threads=1, **kwargs):
    data = model_utils.load_task_with_metadata(input_filename)
    colors = data["colors"]
    cost = data["cost"]
    del data["colors"]
    del data["cost"]

    grouping_result = solve(
        colors,
        cost,
        time_limit=time_limit,
        logfile=os.path.join(output_directory, "gurobi.log"),
        threads=threads,
    )
    model_utils.save_results(output_directory, grouping_result, **data)


def solve(colors, cost, time_limit=300, logfile="gurobi.log", threads=1):
    if len(colors) == 0:
        return GroupingResult([], 0, None)
    if len(colors) == 1:
        return GroupingResult([1], 0, None)

    N = len(colors)
    size_by_color = {
        color: sum(1 for c in colors if c == color) for color in set(colors)
    }
    largest_color = max(size_by_color, key=size_by_color.get)
    MAX_C = size_by_color[largest_color]  # largest color size

    model = gp.Model("grouping")
    model.Params.LogFile = logfile
    model.Params.Threads = threads
    model.Params.TimeLimit = time_limit

    # creating component variables
    components = {}
    for i in range(N):
        components[i] = model.addVar(vtype=GRB.INTEGER, lb=1, ub=MAX_C, name=f"components_{i}")

    # assigning values to the nodes of the largest component directly
    for c, i in enumerate(
        [i for i, color in enumerate(colors) if color == largest_color]
    ):
        model.addConstr(components[i] == c + 1)
    # creating variable to have the value of the absolute difference between the labels (D)

    diff_label = {}
    for i in range(N):
        for j in range(i + 1, N):
            diff_label[i, j] = model.addVar(
                vtype=GRB.INTEGER, lb=-MAX_C, ub=MAX_C, name=f"diff_label_{i}_{j}"
            )
            model.addConstr(diff_label[i, j] == components[i] - components[j])

    abs_diff_label = {}
    for i in range(N):
        for j in range(i + 1, N):
            abs_diff_label[i, j] = model.addVar(
                vtype=GRB.INTEGER, lb=0, ub=MAX_C, name=f"abs_diff_label_{i}_{j}"
            )
    for i in range(N):
        for j in range(i + 1, N):
            model.addGenConstrAbs(abs_diff_label[i, j], diff_label[i, j])

    # forbidding nodes with same color to have absolute difference 0
    for i in range(N):
        for j in range(i + 1, N):
            if colors[i] == colors[j]:
                model.addConstr(abs_diff_label[i, j] >= 1)

    # creating indicators of nodes being in the same component
    are_connected = {}
    for i in range(N):
        for j in range(i + 1, N):
            are_connected[i, j] = model.addVar(vtype=GRB.BINARY, name=f"are_connected_{i}_{j}")
    # since we are maximizing those variables, it is sufficient to forbid to have value 1 if difference is positive
    for i in range(N):
        for j in range(i + 1, N):
            # prob += are_connected[i, j] * MAX_C <= MAX_C - abs_diff_label[i, j]
            model.addConstr(are_connected[i, j] * MAX_C <= MAX_C - abs_diff_label[i, j])

    # creating the target function
    model.setObjective(
        gp.quicksum(cost[i, j] * are_connected[i, j] for i in range(N) for j in range(i + 1, N) if not np.isnan(cost[i, j])),
        GRB.MAXIMIZE,
    )

    # solving the problem
    model.optimize()

    upper_bound = model.ObjBound

    # extract the solution
    raw_grouping = [int(components[i].X) for i in range(N)]
    grouping = model_utils.relabel_grouping(raw_grouping)
    total_value = model.objVal

    return GroupingResult(grouping, total_value, None, upper_bound=upper_bound)


if __name__ == "__main__":
    argh.dispatch_command(main)

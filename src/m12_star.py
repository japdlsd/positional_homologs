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

    model = gp.Model("grouping")
    model.Params.LogFile = logfile
    model.Params.Threads = threads
    model.Params.TimeLimit = time_limit

    # each present edge has an indicator of whether it is chosen from a particular end.
    # constraint: an edge can be activated only from one end
    is_edge_activated = {}
    for i in range(N):
        for j in range(i + 1, N):
            if not np.isnan(cost[i, j]) and colors[i] != colors[j]:
                is_edge_activated[i, j] = model.addVar(
                    vtype=GRB.BINARY, name=f"is_edge_activated_{i}_{j}"
                )
                is_edge_activated[j, i] = model.addVar(
                    vtype=GRB.BINARY, name=f"is_edge_activated_{j}_{i}"
                )
                model.addConstr(is_edge_activated[i, j] + is_edge_activated[j, i] <= 1)
    # each vertex can have at most one activated edge from its end
    # also, it's not possible to have an activated edge if you have incoming activated edges
    for i in range(N):
        model.addConstr(
            gp.quicksum(
                is_edge_activated[i, j] for j in range(N) if (i, j) in is_edge_activated # edge activated from i
            )
            * N
            + gp.quicksum(
                is_edge_activated[j, i] for j in range(N) if (j, i) in is_edge_activated # edge activated to i
            )
            <= N
        )

    # each vertex can have at most 1 incoming activated edge of each color
    for i in range(N):
        for color in set(colors):
            if sum(1 for j in range(N) if colors[j] == color and (j, i) in is_edge_activated) == 0:
                continue
            model.addConstr(
                gp.quicksum(
                    is_edge_activated[j, i] for j in range(N) if colors[j] == color and (j, i) in is_edge_activated
                )
                <= 1
            )

    # setting the target function
    model.setObjective(
        gp.quicksum(
            cost[i, j] * is_edge_activated[i, j]
            for i in range(N)
            for j in range(N)
            if (i, j) in is_edge_activated
        ),
        GRB.MAXIMIZE,
    )

    # solving the problem
    model.optimize()

    # extract the lower and upper bound
    upper_bound = model.ObjBound

    # extract the solution
    # we need to find the connected components
    G = nx.Graph()
    G.add_nodes_from(range(N))
    for i in range(N):
        for j in range(N):
            if not np.isnan(cost[i, j]) and (i, j) in is_edge_activated and is_edge_activated[i, j].X > 0.5:
                G.add_edge(i, j)
    components = list(nx.connected_components(G))
    # creating the grouping
    grouping = [0] * N
    for i, component in enumerate(components):
        for j in component:
            grouping[j] = i + 1

    # debug print
    # print all `is_edge_activated with value > 0.5`
    for i, j in is_edge_activated:
        if is_edge_activated[i, j].X > 0.5:
            print(f"Edge {i} -- {j} is activated (value: {is_edge_activated[i, j].X})")
    # print all edges of graph G
    print("Edges:")
    for i, j in G.edges:
        print(f"  edge {i} -- {j}")
    for i, component in enumerate(components):
        print(f"Component {i}: {component}")
        for j in component:
            print(f"  node {j} of color {colors[j]}")
            print(f"Active outgoing edges: {[(j, k) for k in range(N) if (j, k) in is_edge_activated and is_edge_activated[j, k].X > 0.5]}")
            print(f"Active incoming edges: {[(k, j) for k in range(N) if (k, j) in is_edge_activated and is_edge_activated[k, j].X > 0.5]}")

    return GroupingResult(grouping, model.objVal, None, upper_bound=upper_bound)


if __name__ == "__main__":
    argh.dispatch_command(main)

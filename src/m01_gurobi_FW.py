import os.path

import argh
import pulp
import networkx as nx

import helpers
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

    prob = pulp.LpProblem("Grouping_by_Floyd-Warshall", pulp.LpMaximize)
    all_possible_edges = [
        (i, j) for i in range(N) for j in range(i + 1, N) if colors[i] != colors[j]
    ]

    is_edge_chosen = pulp.LpVariable.dicts(
        "is_edge_chosen", all_possible_edges, 0, 1, pulp.LpInteger
    )
    # setting the target function
    prob += pulp.lpSum(
        [
            cost[i][j] * is_edge_chosen[(i, j)]
            for i, j in all_possible_edges
            if not np.isnan(cost[i][j])
        ]
    )

    # setting the FLoys-Warshall constraints
    add_floyd_warshall_constraints(prob, colors, is_edge_chosen, all_possible_edges)

    # solving the problem
    prob.solve(
        solver=pulp.GUROBI_CMD(
            timeLimit=int(time_limit), logPath=logfile, threads=threads
        )
    )

    # getting the result
    chosen_edges = [
        edge for edge in all_possible_edges if is_edge_chosen[edge].value() == 1
    ]
    # creatng a graph using networkx and enumerating all components
    G = nx.Graph()
    G.add_nodes_from(range(N))
    G.add_edges_from(chosen_edges)
    components = list(nx.connected_components(G))
    # creating the grouping
    grouping = [0] * N
    for i, component in enumerate(components):
        for j in component:
            grouping[j] = i + 1

    return GroupingResult(grouping, pulp.value(prob.objective), prob)


if __name__ == "__main__":
    argh.dispatch_command(main)

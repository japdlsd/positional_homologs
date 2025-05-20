import itertools
import os.path
import random

import sys
from collections import defaultdict
from enum import Enum

import argh
import pulp
import networkx as nx

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


class NodeType(Enum):
    NODE = 0
    COMPONENT = 1


def solve(colors, cost, time_limit=300, logfile="gurobi.log", threads=1):
    if len(colors) == 0:
        return GroupingResult([], 0, None)
    if len(colors) == 1:
        return GroupingResult([1], 0, None)

    N = len(colors)

    grouping = [i + 1 for i in range(N)]
    grouping, total_value = polish_grouping_iterative_color(
        colors, cost, grouping, logfile
    )

    return GroupingResult(grouping, total_value, None)


def polish_grouping_iterative_color(colors, cost, grouping, logfile):
    total_value = 0
    N = len(colors)
    for i in range(N):
        for j in range(i + 1, N):
            if grouping[i] == grouping[j] and (not np.isnan(cost[i, j])):
                total_value += cost[i, j]

    with open(logfile, "w") as f:
        has_improved = True
        while has_improved:
            print(f"Current solution: {grouping}", file=f)
            has_improved = False
            untried_colors = sorted(set(colors))
            # pick a random color
            while (not has_improved) and len(untried_colors) > 0:
                chosen_color = random.choice(untried_colors)
                # remove that color from the untried_colors
                untried_colors.remove(chosen_color)
                print(
                    f"Trying color {chosen_color}, {len(untried_colors)} untried remaining",
                    file=f,
                )

                # let's create a matching task between nodes of chosen color and components
                G = nx.Graph()
                # adding nodes of the chosen color
                for i in range(N):
                    if colors[i] == chosen_color:
                        G.add_node((i, NodeType.NODE))
                for j in set(
                    component
                    for i, component in enumerate(grouping)
                    if colors[i] != chosen_color
                ):
                    G.add_node((j, NodeType.COMPONENT))
                # adding weight to the edges
                edge_weights = defaultdict(float)
                for i in range(N):
                    for x in range(N):
                        if (
                            colors[i] == chosen_color
                            and colors[x] != chosen_color
                            and (not np.isnan(cost[i, x]))
                        ):
                            component = grouping[x]
                            edge_weights[(i, component)] += cost[i, x]
                # adding edges to the graph
                for (node, component), weight in edge_weights.items():
                    G.add_edge(
                        (node, NodeType.NODE),
                        (component, NodeType.COMPONENT),
                        weight=weight,
                    )
                # computing matching
                matching = nx.max_weight_matching(G)
                print(f"Matching: {matching}", file=f)
                # adding components to the nodes
                candidate_grouping = grouping[:]
                # remove existing component from nodes of the chosen color
                for i in range(N):
                    if colors[i] == chosen_color:
                        candidate_grouping[i] = None
                for x, y in matching:
                    if x[1] == NodeType.COMPONENT:
                        x, y = y, x
                    candidate_grouping[x[0]] = y[0]
                # add free components for elements with None in the candidate grouping
                for i in range(N):
                    if candidate_grouping[i] is None:
                        candidate_grouping[i] = smallest_free_number(candidate_grouping)
                print(f"Candidate grouping: {candidate_grouping}", file=f)
                # computing new cost
                candidate_total_value = 0
                for i in range(N):
                    for j in range(i + 1, N):
                        if candidate_grouping[i] == candidate_grouping[j] and (
                            not np.isnan(cost[i, j])
                        ):
                            candidate_total_value += cost[i, j]
                print(f"Candidate total value: {candidate_total_value}", file=f)

                if candidate_total_value > total_value:
                    has_improved = True
                    grouping = candidate_grouping
                    total_value = candidate_total_value
                    print(f"Improved solution to {total_value}", file=f)
        grouping = model_utils.relabel_grouping(grouping)
    return grouping, total_value


def smallest_free_number(xs):
    occupied = set(xs)
    for x in range(1, len(xs) + 1):
        if x not in occupied:
            return x


if __name__ == "__main__":
    argh.dispatch_command(main)

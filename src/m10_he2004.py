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

    color_set = sorted(set(colors))

    # first, create the groupings of the elements
    # in the beginnig, each element is in a group of its own
    clusters_by_color = {
        color: [{i} for i in range(N) if colors[i] == color] for color in color_set
    }
    total_value = 0
    order_of_color_merge = [
        (color_set[0], color_set[i]) for i in range(len(set(color_set))-1, 0, -1)
    ]
    print("Order of colors to merge:", order_of_color_merge, file=sys.stderr)

    for color1, color2 in order_of_color_merge:
        print(f"Processing colors {color1} and {color2}", file=sys.stderr)
        print(f"{clusters_by_color=}", file=sys.stderr)
        print(f"Clusters by color {color1}: {clusters_by_color[color1]}", file=sys.stderr)
        print(f"Clusters by color {color2}: {clusters_by_color[color2]}", file=sys.stderr)
        # we create a graph where the nodes are the clusters of color1 and color2
        G = nx.Graph()
        for i, cluster1 in enumerate(clusters_by_color[color1]):
            G.add_node((color1, i))
        for j, cluster2 in enumerate(clusters_by_color[color2]):
            G.add_node((color2, j))
        # now we need to compute the cost of edges between the clusters
        for i, cluster1 in enumerate(clusters_by_color[color1]):
            for j, cluster2 in enumerate(clusters_by_color[color2]):
                edge_cost = 0
                for node1 in cluster1:
                    for node2 in cluster2:
                        if not np.isnan(cost[node1, node2]):
                            edge_cost += cost[node1, node2]
                G.add_edge((color1, i), (color2, j), weight=edge_cost)
        # now we need to find the maximum matching in this graph
        matching = nx.max_weight_matching(G)
        print(f"Matching: {matching}", file=sys.stderr)
        # now we need to merge the clusters according to the matching
        # and add the cost of connection to the total cost
        merged_clusters_of_color2 = set()
        for node1, node2 in matching:
            if color1 == node2[0]:
                node1, node2 = node2, node1
            merged_clusters_of_color2.add(node2[1])
            total_value += G[node1][node2]["weight"]
            clusters_by_color[color1][node1[1]].update(clusters_by_color[color2][node2[1]])
        # disband the clusters of color2: remove those that were merged
        # and add the remaining to the clusters of color1
        for j, cluster in enumerate(clusters_by_color[color2]):
            if j not in merged_clusters_of_color2:
                clusters_by_color[color1].append(cluster)
        del clusters_by_color[color2]

    assert len(clusters_by_color) == 1, "There should be only one color left"
    final_clusters = list(clusters_by_color.values())[0]
    print(f"Final clusters: {final_clusters}", file=sys.stderr)
    # now we need to parse the final grouping
    grouping = [None] * N
    for cluster_num, cluster in enumerate(final_clusters):
        for x in cluster:
            grouping[x] = cluster_num + 1

    return GroupingResult(grouping, total_value, None)



if __name__ == "__main__":
    argh.dispatch_command(main)

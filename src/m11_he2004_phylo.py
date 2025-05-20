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

import helpers
# import helpers
from grouping import *
import model_utils

from m10_he2004 import solve as solve_without_phylogeny

def main(input_filename, phylogeny_filename, output_directory, time_limit=300, threads=1, **kwargs):
    data = model_utils.load_task_with_metadata(input_filename)
    colors = data["colors"]
    cost = data["cost"]
    del data["colors"]
    del data["cost"]

    phylogeny = helpers.load_phylogeny_tree_from_file(phylogeny_filename)
    logfile = os.path.join(output_directory, "gurobi.log")

    grouping_result = solve_hybrid(colors, cost, logfile, phylogeny, threads, time_limit)

    model_utils.save_results(output_directory, grouping_result, **data)


def solve_hybrid(colors, cost, logfile, phylogeny, threads, time_limit):
    phylogeny_leaves = [node for node in phylogeny.nodes() if phylogeny.out_degree(node) == 0]
    if set(colors).issubset(set(phylogeny_leaves)):
        grouping_result = solve(
            colors,
            cost,
            phylogeny,
            time_limit=time_limit,
            logfile=logfile,
            threads=threads,
        )
    else:
        print("The set of colors does not match the set of leaves of the phylogeny, ingoring the phylogeny",
              file=sys.stderr)
        grouping_result = solve_without_phylogeny(
            colors,
            cost,
            time_limit=time_limit,
            logfile=logfile,
            threads=threads,
        )
    return grouping_result


def solve(colors, cost, phylogeny, time_limit=300, logfile="gurobi.log", threads=1):
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

    order_of_color_merge = []
    remaining_colors = set(color_set)
    phylogeny = phylogeny.to_undirected()
    while len(remaining_colors) > 1:
        # find two closest colors in the phylogeny, remove one of them
        # and add this pair to the order of color merge
        min_distance = float("inf")
        closest_color_pair = (None, None)
        for color1, color2 in itertools.combinations(remaining_colors, 2):
            distance = nx.shortest_path_length(phylogeny, color1, color2)
            if distance < min_distance:
                min_distance = distance
                closest_color_pair = (color1, color2)
        order_of_color_merge.append(closest_color_pair)
        remaining_colors.remove(closest_color_pair[1])

    print("Order of colors to merge:", order_of_color_merge, file=sys.stderr)

    for color1, color2 in order_of_color_merge:
        print(f"Processing colors {color1} and {color2}", file=sys.stderr)
        if color1 not in clusters_by_color and color2 not in clusters_by_color:
            clusters_by_color[color1] = []
        if color1 not in clusters_by_color:
            clusters_by_color[color1] = clusters_by_color[color2]
            del clusters_by_color[color2]
            continue
        if color2 not in clusters_by_color:
            continue

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

    assert len(clusters_by_color) == 1, f"There should be only one color left\nFinal clustering: {clusters_by_color}"
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

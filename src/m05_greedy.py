import itertools
import os.path

import sys

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


def solve(colors, cost, time_limit=300, logfile="gurobi.log", threads=1):
    if len(colors) == 0:
        return GroupingResult([], 0, None)
    if len(colors) == 1:
        return GroupingResult([1], 0, None)

    N = len(colors)

    # order edges from the biggest to the smallest, add them as long as there is no conflict
    uf = UnionFindColored(colors)
    edges = sorted(
        (
            (cost[i, j], i, j)
            for i in range(N)
            for j in range(i + 1, N)
            if not np.isnan(cost[i, j])
        ),
        reverse=True,
    )
    total_value = 0
    for c, i, j in edges:
        try:
            uf.join_sets(i, j)
            total_value += c
        except ValueError:
            pass  # intentionally
    raw_grouping = [uf.find_root(i) for i in range(N)]
    grouping = model_utils.relabel_grouping(raw_grouping)

    return GroupingResult(grouping, total_value, None)


class UnionFindColored:
    def __init__(self, colors):
        self.n = len(colors)
        self.parent = [i for i in range(self.n)]
        self.size = [1 for _ in range(self.n)]
        self.colors = [{colors[i]} for i in range(self.n)]
        self.number_of_components = self.n

    def find_root(self, x):
        # looking for the root
        root = x
        while self.parent[root] != root:
            root = self.parent[root]
        # updating every node on the path
        q = x
        while q != root:
            next_q = self.parent[q]
            self.parent[q] = root
            q = next_q
        return root

    def are_connected(self, x, y):
        return self.find_root(x) == self.find_root(y)

    def join_sets(self, x, y):
        # check if already connected
        if self.are_connected(x, y):
            return

        # find the roots of the sets
        root_x = self.find_root(x)
        root_y = self.find_root(y)
        size_x = self.size[root_x]
        size_y = self.size[root_y]

        # check if the color sets are disjoint
        if self.colors[root_x] & self.colors[root_y]:
            raise ValueError(
                f"Error while adding edge {x} - {y}! Sets with roots {root_x} and {root_y} are not color-disjoint! {self.colors[root_x]=} {self.colors[root_y]=}"
            )

        # join the sets
        if size_x < size_y:
            self.parent[root_x] = root_y

            self.size[root_y] += size_x
            self.colors[root_y] |= self.colors[root_x]
        else:
            self.parent[root_y] = root_x

            self.size[root_x] += size_y
            self.colors[root_x] |= self.colors[root_y]
        self.number_of_components -= 1


if __name__ == "__main__":
    argh.dispatch_command(main)

from collections import namedtuple

import numpy as np
import pulp
import networkx as nx

GroupingResult = namedtuple("GroupingResult", ["grouping", "cost", "solver_object", "upper_bound"], defaults=[None])


def solve_grouping_floyd_warshall_single_cost(colors, cost):
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
    prob.solve(solver=pulp.GUROBI_CMD())

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


def add_floyd_warshall_constraints(prob, colors, is_edge_chosen, all_possible_edges):
    N = len(colors)
    # setting the constraints
    dp_coordinates = [
        (i, j, k) for i in range(N) for j in range(N) for k in range(N + 1)
    ]
    dp = pulp.LpVariable.dicts("dp", dp_coordinates, 0, 1, pulp.LpInteger)
    m = pulp.LpVariable.dicts("m", dp_coordinates, 0, 1, pulp.LpInteger)
    # DP[i, j, 0] = 1 if i == j
    for i in range(N):
        prob += dp[i, i, 0] == 1
    # DP[i, j, 0] = is_edge_chosen[i, j] if i <= j
    # DP[i, j, 0] = is_edge_chosen[j, i] if i > j
    for i, j in all_possible_edges:
        prob += dp[i, j, 0] == is_edge_chosen[(i, j)]
        prob += dp[j, i, 0] == is_edge_chosen[(i, j)]
    # the rest on the zero-th layer have to be zero
    for i in range(N):
        for j in range(N):
            if (
                i != j
                and (i, j) not in all_possible_edges
                and (j, i) not in all_possible_edges
            ):
                prob += dp[i, j, 0] == 0

    # the Floyd-Warshall constraints
    for k in range(N):
        for i in range(N):
            for j in range(N):
                if i == j:
                    prob += dp[i, j, k + 1] == 1
                    continue

                # m = AND(dp[i, k, k], dp[k, j, k])
                prob += m[i, j, k + 1] <= dp[i, k, k]
                prob += m[i, j, k + 1] <= dp[k, j, k]
                prob += m[i, j, k + 1] >= dp[i, k, k] + dp[k, j, k] - 1

                # dp[i, j, k+1] = OR(dp[i, j, k], m[i, j, k+1])
                prob += dp[i, j, k + 1] >= dp[i, j, k]
                prob += dp[i, j, k + 1] >= m[i, j, k + 1]
                prob += dp[i, j, k + 1] <= dp[i, j, k] + m[i, j, k + 1]
    # color constraints on the final layer of the DP
    for i in range(N):
        for j in range(N):
            if colors[i] == colors[j] and i != j:
                for k in range(N + 1):
                    prob += dp[i, j, k] == 0


def solve_grouping_one_pass(colors, cost_primary, cost_secondary):
    if len(colors) == 0:
        return GroupingResult([], 0, None)
    if len(colors) == 1:
        return GroupingResult([1], 0, None)

    cost_primary = np.array(cost_primary)
    cost_secondary = np.array(cost_secondary)

    coefficient = np.nansum(cost_secondary) * 10
    cost_together = cost_primary * coefficient + cost_secondary

    result = solve_grouping_floyd_warshall_single_cost(colors, cost_together)
    return result


class UnionFind:
    def __init__(self, n):
        self.parent = [i for i in range(n)]
        self.size = [1 for _ in range(n)]
        self.number_of_components = n

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
        # join the sets
        if size_x < size_y:
            self.parent[root_x] = root_y
            self.size[root_y] += size_x
        else:
            self.parent[root_y] = root_x
            self.size[root_x] += size_y
        self.number_of_components -= 1


def test_union_find():
    uf = UnionFind(5)
    assert uf.number_of_components == 5
    assert uf.are_connected(0, 1) == False, "0 and 1 are not connected"
    uf.join_sets(0, 1)
    assert uf.are_connected(0, 1) == True, "0 and 1 are connected"
    assert uf.number_of_components == 4, "4 components left"
    uf.join_sets(2, 3)
    assert uf.number_of_components == 3, "3 components left"
    uf.join_sets(0, 3)
    assert uf.number_of_components == 2, "2 components left"
    uf.join_sets(0, 3)
    assert uf.number_of_components == 2, "2 components left"
    assert uf.are_connected(0, 3) == True, "0 and 3 are connected"
    assert uf.are_connected(0, 2) == True, "0 and 2 are connected"
    assert uf.are_connected(0, 4) == False, "0 and 4 are not connected"


def prune_edges_single_cost(cost):
    # The goal is to remove as many cheapest edges as possible while keeping the graph connected
    # the optimal way to do this is to reverse the process and add the most expensive edges
    # until the graph becomes connected. The connectedness is checked by union-find data structure.
    N = cost.shape[0]
    sorted_edges = sorted(
        (
            (cost[i, j], i, j)
            for i in range(N)
            for j in range(N)
            if i < j and not np.isnan(cost[i, j])
        ),
        reverse=True,
    )
    edges_bucketed_by_cost = {}
    for c, i, j in sorted_edges:
        if c not in edges_bucketed_by_cost:
            edges_bucketed_by_cost[c] = []
        edges_bucketed_by_cost[c].append((i, j))

    new_cost = np.full((N, N), np.nan)

    uf = UnionFind(N)
    for c, edges in sorted(edges_bucketed_by_cost.items(), key=lambda x: x[0], reverse=True):
        if uf.number_of_components == 1:
            break
        for i, j in edges:
            uf.join_sets(i, j)
            new_cost[i, j] = c
            new_cost[j, i] = c
    assert uf.number_of_components == 1, "The graph is not connected"
    return new_cost


def test_prune_edges_single_cost():
    cases = [
        (np.array([[0]]), np.array([[np.nan]])),
        (np.array([[0, 1], [1, 0]]), np.array([[np.nan, 1], [1, np.nan]])),
        (
            np.array([[0, 6, 5], [6, 0, 2], [5, 2, 0]]),
            np.array([[np.nan, 6, 5], [6, np.nan, np.nan], [5, np.nan, np.nan]]),
        ),
    ]

    for cost, expected in cases:
        new_cost = prune_edges_single_cost(cost)
        print(f"case: {cost}")
        print(f"expected: {expected}")
        print(f"result: {new_cost}")
        for i in range(cost.shape[0]):
            for j in range(cost.shape[1]):
                assert (
                    np.isnan(new_cost[i, j])
                    and np.isnan(expected[i, j])
                    or new_cost[i, j] == expected[i, j]
                ), f"new_cost[{i}, {j}] != expected[{i}, {j}]"


def test_grouping_case1():
    colors = [1]
    cost_primary = [[0]]
    cost_secondary = [[0]]

    expected_grouping = [1]
    result = solve_grouping_floyd_warshall_single_cost(colors, np.array(cost_primary))
    assert (
        result.grouping == expected_grouping
    ), f"{result.grouping} != {expected_grouping}"


def test_grouping_case2():
    colors = [1, 2]
    cost_primary = [[0, 1], [1, 0]]
    cost_secondary = [[0, 1], [1, 0]]
    expected_grouping = [1, 1]
    result = solve_grouping_floyd_warshall_single_cost(colors, np.array(cost_primary))
    assert (
        result.grouping == expected_grouping
    ), f"{result.grouping} != {expected_grouping}"


def test_grouping_case3():
    colors = [1, 1]
    cost_primary = [[0, 1], [1, 0]]
    cost_secondary = [[0, 1], [1, 0]]
    expected_grouping = [1, 2]
    result = solve_grouping_floyd_warshall_single_cost(colors, np.array(cost_primary))
    assert (
        result.grouping == expected_grouping
    ), f"{result.grouping} != {expected_grouping}"


def test_grouping_case4():
    colors = [1, 2, 3, 1, 2]
    cost_primary = [
        [0, np.nan, 1, 1, 1],
        [1, 0, 1, 1, 1],
        [1, 1, 0, 1, 1],
        [1, 1, 1, 0, 1],
        [1, 1, 1, 1, 0],
    ]
    cost_secondary = np.zeros((5, 5))
    expected_groupings = (
        [1, 1, 1, 2, 2],
        [1, 2, 1, 2, 1],
        [1, 1, 2, 2, 2],
        [1, 2, 2, 2, 1],
    )
    result = solve_grouping_floyd_warshall_single_cost(colors, np.array(cost_primary))

    print(result.solver_object)

    assert (
        result.grouping in expected_groupings
    ), f"{result.grouping} not in {expected_groupings}"

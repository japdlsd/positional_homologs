import itertools
import os.path

import sys

import argh
import pulp
import networkx as nx

# import helpers
from grouping import *
import model_utils
from m07_iterative_color import solve as iterative_solve


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

    greedy_solution = iterative_solve(colors, cost)
    greedy_grouping = greedy_solution.grouping

    N = len(colors)
    size_by_color = {
        color: sum(1 for c in colors if c == color) for color in set(colors)
    }
    largest_color = max(size_by_color, key=size_by_color.get)
    MAX_C = size_by_color[largest_color]  # largest color size
    # we need to adjust MAX_C to incorporate the greedy solution
    MAX_C = max(MAX_C, len(set(greedy_grouping)))

    prob = pulp.LpProblem("Grouping_by_components", pulp.LpMaximize)
    # creating component variables
    components = pulp.LpVariable.dicts("components", range(N), 1, MAX_C, pulp.LpInteger)
    # setting the initial values to the one from the greedy solution
    for i in range(N):
        components[i].setInitialValue(greedy_grouping[i])

    # adding the constraints on the components: if the colors are the same, then the components have to be different
    for i in range(N):
        for j in range(i + 1, N):
            if colors[i] == colors[j]:
                prob += components[i] != components[j]

    # creating variable to have the value of the absolute difference between the labels (D)
    # zs is an auxiliary indicator which has value 1 if c_i >= c_j and 0 if c_i <= c_j (sic! it is not mutually exclusive!)
    abs_diff_label = pulp.LpVariable.dicts(
        "abs_diff_label",
        [(i, j) for i in range(N) for j in range(i + 1, N)],
        0,
        MAX_C,
        pulp.LpInteger,
    )
    zs = pulp.LpVariable.dicts(
        "zs", [(i, j) for i in range(N) for j in range(i + 1, N)], 0, 1, pulp.LpInteger
    )
    for i in range(N):
        for j in range(i + 1, N):
            add_abs_value_constraint(
                components[i],
                components[j],
                abs_diff_label[i, j],
                zs[i, j],
                MAX_C,
                prob,
            )
    # forbidding nodes with same color to have absolute difference 0
    for i in range(N):
        for j in range(i + 1, N):
            if colors[i] == colors[j]:
                prob += abs_diff_label[i, j] >= 1

    # creating indicators of nodes being in the same component
    are_connected = pulp.LpVariable.dicts(
        "are_connected",
        [(i, j) for i in range(N) for j in range(i + 1, N)],
        0,
        1,
        pulp.LpInteger,
    )
    # since we are maximizing those variables, it is sufficient to forbid to have value 1 if difference is positive
    for i in range(N):
        for j in range(i + 1, N):
            prob += are_connected[i, j] * MAX_C <= MAX_C - abs_diff_label[i, j]

    # creating the target function
    prob += pulp.lpSum(
        [
            cost[i, j] * are_connected[i, j]
            for i in range(N)
            for j in range(i + 1, N)
            if not np.isnan(cost[i, j])
        ]
    )

    # solving the problem
    prob.solve(
        solver=pulp.GUROBI_CMD(
            timeLimit=int(time_limit), logPath=logfile, threads=threads, warmStart=True
        )
    )

    # getting the result
    grouping = [int(round(components[i].value())) for i in range(N)]

    return GroupingResult(grouping, pulp.value(prob.objective), prob)


def add_abs_value_constraint(x, y, d, z, M, prob):
    """Add a rule that d = abs(x - y), where z is an auxiliary indicator variable,
    range of x and y is 1..M, and prob is the PuLP task object."""
    prob += x - y <= d
    prob += y - x <= d
    prob += x - y <= M * z
    prob += y - x <= M * (1 - z)
    prob += d <= x - y + 2 * M * (1 - z)
    prob += d <= y - x + 2 * M * z


def test_absolute_value_construction():
    M = 15
    for xval in range(1, M + 1):
        for yval in range(1, M + 1):
            prob = pulp.LpProblem("Test_absolute_value_construction", pulp.LpMaximize)
            x = xval
            y = yval
            d = pulp.LpVariable("d", 0, M + 1, pulp.LpInteger)
            z = pulp.LpVariable("z", 0, 1, pulp.LpInteger)
            add_abs_value_constraint(x, y, d, z, M, prob)

            prob += d - abs(x - y)

            prob.solve()

            assert d.value() == abs(
                x - y
            ), f"Expected {d.value()} == {abs(x - y)}, got x={x}, y={y} d={d.value()} z={z.value()}"

            prob = pulp.LpProblem("Test_absolute_value_construction", pulp.LpMinimize)
            x = xval
            y = yval
            d = pulp.LpVariable("d", 0, M + 1, pulp.LpInteger)
            z = pulp.LpVariable("z", 0, 1, pulp.LpInteger)
            add_abs_value_constraint(x, y, d, z, M, prob)

            prob += d - abs(x - y)

            prob.solve()

            assert d.value() == abs(
                x - y
            ), f"Expected {d.value()} == {abs(x - y)}, got x={x}, y={y} d={d.value()} z={z.value()}"


if __name__ == "__main__":
    argh.dispatch_command(main)

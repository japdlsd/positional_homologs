import json

import dill
import numpy as np


def load_task(input_filename):
    with open(input_filename, "r") as f:
        data = json.load(f)

    colors = data["colors"]
    cost_raw = data["cost"]
    for i, row in enumerate(cost_raw):
        for j, raw_value in enumerate(row):
            try:
                value = float(raw_value)
            except ValueError:
                value = np.nan
            cost_raw[i][j] = value
    cost = np.array(cost_raw, dtype=float)
    return colors, cost


def load_task_with_metadata(input_filename):
    with open(input_filename, "r") as f:
        data = json.load(f)

    cost_raw = data["cost"]
    for i, row in enumerate(cost_raw):
        for j, raw_value in enumerate(row):
            try:
                value = float(raw_value)
            except ValueError:
                value = np.nan
            cost_raw[i][j] = value
    cost = np.array(cost_raw, dtype=float)
    data["cost"] = cost

    return data


def save_results(output_directory, grouping_result, **kwargs):
    with open(f"{output_directory}/grouping.json", "w") as f:
        json.dump(grouping_result.grouping, f, indent=4)
    additional_data = {
        "upper_bound": grouping_result.upper_bound,
        "internal_score": grouping_result.cost,
        "grouping": grouping_result.grouping,
        **kwargs,
    }
    with open(f"{output_directory}/additional_data.json", "w") as f:
        json.dump(additional_data, f, indent=4)

    # with open(f"{output_directory}/full_result.dill", 'wb') as f:
    #     dill.dump(grouping_result, f)


def relabel_grouping(raw_grouping):
    N = len(raw_grouping)
    label_set = set(raw_grouping)
    remapping = {label: i + 1 for i, label in enumerate(label_set)}
    grouping = [remapping[raw_grouping[i]] for i in range(N)]
    return grouping

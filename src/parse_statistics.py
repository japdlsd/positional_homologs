import json
import sys
import tempfile

import argh

import model_utils
from m07_iterative_color import polish_grouping_iterative_color


def main(
    task_filename, grouping_filename, measurements_filename, additional_data_filename, task_label, model_label
):
    result = {
        "task": task_label,
        "model": model_label,
    }

    colors, cost_matrix = model_utils.load_task(task_filename)
    with open(grouping_filename, "r") as f:
        grouping = json.load(f)
    # compute the actual score from the grouping
    score = 0
    is_valid_grouping = True
    for i in range(len(colors)):
        for j in range(i + 1, len(colors)):
            if grouping[i] == grouping[j]:
                if colors[i] == colors[j]:
                    is_valid_grouping = False
                if cost_matrix[i, j] is not None:
                    score += cost_matrix[i, j]
                else:
                    score += 0
    result["number_of_vertices"] = len(colors)
    result["score"] = score
    result["valid"] = is_valid_grouping
    result["number_of_components"] = len(set(grouping))

    with tempfile.NamedTemporaryFile() as f:
        polished_grouping, polished_score = polish_grouping_iterative_color(
            colors, cost_matrix, grouping, logfile=f.name
        )

    result["polished_score"] = polished_score
    result["polished_valid"] = is_valid_grouping
    result["polished_number_of_components"] = len(set(polished_grouping))

    # parse the additional data file
    with open(additional_data_filename, "r") as f:
        additional_data = json.load(f)
        result["upper_bound"] = additional_data["upper_bound"]
        result["internal_score"] = additional_data["internal_score"]

    # parse the measurements file
    with open(measurements_filename, "r") as f:
        header = f.readline().strip().split("\t")
        colname_running_time = "s"
        colnum_running_time = header.index(colname_running_time)
        assert colnum_running_time >= 0
        colname_max_memory = "max_vms"
        colnum_max_memory = header.index(colname_max_memory)
        assert colnum_max_memory >= 0

        running_times = []
        max_memories = []

        for line in f:
            line = line.strip().split("\t")
            running_time = float(line[colnum_running_time]) if line[colnum_running_time] != "NA" else 0.0
            max_memory = float(line[colnum_max_memory]) if line[colnum_max_memory] != "NA" else 0.0
            running_times.append(running_time)
            max_memories.append(max_memory)

        avg_running_time = sum(running_times) / len(running_times)
        avg_max_memory = sum(max_memories) / len(max_memories)
        result["avg_running_time_sec"] = avg_running_time
        result["avg_max_memory_mb"] = avg_max_memory

    # save the result
    json.dump(result, sys.stdout, indent=4)


if __name__ == "__main__":
    argh.dispatch_command(main)

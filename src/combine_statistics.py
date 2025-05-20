import json
import sys

import argh


@argh.arg("input_filenames", nargs="+")
def main(input_filenames):
    data = []
    for filename in input_filenames:
        print(f"Loading {filename}", flush=True, file=sys.stderr)
        with open(filename, "r") as f:
            dato = json.load(f)
            data.append(dato)

    header = [
        "model",
        "task",
        "number_of_vertices",
        "polished_score",
        "upper_bound",
        "avg_running_time_sec",
        "avg_max_memory_mb",
        "raw_score",
        "internal_score",
        "valid",
        "number_of_components",
        "polished_valid",
        "polished_number_of_components",
    ]
    rows = []
    for dato in data:
        row = [
            dato["model"],
            dato["task"],
            dato["number_of_vertices"],
            dato["polished_score"],
            dato["upper_bound"],
            dato["avg_running_time_sec"],
            dato["avg_max_memory_mb"],
            dato["score"],
            dato["internal_score"],
            dato["valid"],
            dato["number_of_components"],
            dato["polished_valid"],
            dato["polished_number_of_components"],
        ]
        rows.append(row)

    rows = sorted(rows, key=lambda x: (x[0], x[2], x[1]))

    print("\t".join(header))
    for row in rows:
        print("\t".join(str(x) for x in row))


if __name__ == "__main__":
    argh.dispatch_command(main)

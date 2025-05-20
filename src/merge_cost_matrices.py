import json
import os

import argh
import numpy as np


def main(similiarity_filename, blast_filename, output_filename):
    with open(similiarity_filename) as f:
        similarity_data = json.load(f)
    with open(blast_filename) as f:
        blast_data = json.load(f)
    similarity_matrix = np.array(similarity_data["cost"])
    blast_matrix = np.array(blast_data["cost"])

    coefficient = np.nansum(blast_matrix) * 10
    cost_matrix = np.array(similarity_matrix) * coefficient + np.array(blast_matrix)

    output = similarity_data.copy()
    output["cost"] = cost_matrix.tolist()

    output_dirname = os.path.dirname(output_filename)
    os.makedirs(output_dirname, exist_ok=True)

    with open(output_filename, "w") as f:
        json.dump(output, f)


if __name__ == "__main__":
    argh.dispatch_command(main)

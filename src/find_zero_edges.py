import json
import os

import argh


def main(target_directory):
    if not target_directory.endswith("/"):
        target_directory += "/"
    if not os.path.exists(target_directory):
        raise IOError("Target directory does not exist")
    things_inside = os.listdir(target_directory)
    directories_inside = [
        d for d in things_inside if os.path.isdir(os.path.join(target_directory, d))
    ]
    # those are actually names of the gene groups
    gene_group_names = directories_inside

    total_zero_count = 0

    for gene_group_name in gene_group_names:
        filename = os.path.join(target_directory, gene_group_name, "single_cost.json")
        with open(filename) as f:
            data = json.load(f)
        blocks = data['blocks']
        costs = data['cost']
        for x in range(len(costs)):
            for y in range(len(costs[x])):
                if costs[x][y] is not None and  costs[x][y] < 0.05:
                    total_zero_count += 1
                    if len(blocks) > 15:
                        print(f"Small cost edge found in {filename}: {x}->{y} ({costs[x][y]:.4f})")


    print(f"Total Small cost edges found: {total_zero_count}")


if __name__ == "__main__":
    argh.dispatch_command(main)

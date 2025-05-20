import json
import os
import sys
from collections import defaultdict

import argh


def main(input_dirname, output_filename):
    # example of output:
    # {"comment": "Map family to positional ortholog lists. Each list has pairs (species, geneid).",
    # "data": {
    #  "GeneFamilyID1": [[["S12", "13"], ["S11", "12"], ["S6", "7"], ["S8", "9"], ["S10", "11"]]],
    #  "GeneFamilyID2": [[["S8", "19"], ["S6", "15"]], [["S11", "22"], ["S8", "17"], ["S12", "23"], ["S10", "21"], ["S6", "13"]]],
    # ...
    # }}
    output_obj = {"comment": "Map family to positional ortholog lists. Each list has pairs (species, geneid).",
                  "data": {}}

    for group_dirname in os.listdir(input_dirname):
        print(f"Processing {group_dirname}", file=sys.stderr, flush=True)
        if os.path.isdir(os.path.join(input_dirname, group_dirname)):
            filename = os.path.join(input_dirname, group_dirname, "additional_data.json")
            if not os.path.exists(filename):
                print(f"WARNING: File {filename} does not exist", file=sys.stderr, flush=True)
                continue
            with open(filename) as f:
                data = json.load(f)
            # example of data:
            # {
            #     "upper_bound": null,
            #     "internal_score": 0,
            #     "grouping": [
            #         1
            #     ],
            #     "taskname": "1",
            #     "blocks": [
            #         {
            #             "chromosome": "2L",
            #             "start": 896,
            #             "end": 3861,
            #             "homolog_group": "1",
            #             "strand": "+",
            #             "gene_name": "AALB014720",
            #             "organism": "Aalb"
            #         }
            #     ]
            # }
            output_obj["data"][group_dirname] = []
            genes_by_grouping = defaultdict(list)
            for grouping, block in zip(data["grouping"], data["blocks"]):
                genes_by_grouping[grouping].append([block["organism"], block["gene_name"]])
            print(f"Found {len(genes_by_grouping)} groupings in {filename}", file=sys.stderr, flush=True)
            for genes in genes_by_grouping.values():
                output_obj["data"][group_dirname].append(genes)

    with open(output_filename, "w") as f:
        json.dump(output_obj, f, indent=2)


if __name__ == "__main__":
    argh.dispatch_command(main)

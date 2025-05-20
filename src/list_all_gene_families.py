import json
import os
import sys
from collections import defaultdict
from pathlib import Path

import argh
import yaml

import helpers


def main(input_config_filename):
    with open(input_config_filename, "r") as f:
        config = yaml.safe_load(f)

    input_dir = os.path.dirname(input_config_filename)
    organism_identifiers = list(config["homolog_families"].keys())

    blocks_by_id = {}
    for organism_id in organism_identifiers:
        blocks_filename = os.path.join(
            input_dir, config["homolog_families"][organism_id]
        )

        blocks = helpers.load_bed_homolog(Path(blocks_filename), organism=organism_id)
        blocks_by_id[organism_id] = sorted(blocks)

    # split the input files into the groups
    homolog_groups_by_id = defaultdict(list)

    for organism_id, blocks in blocks_by_id.items():
        for block in blocks:
            block_with_organism = block._replace(organism=organism_id)
            homolog_groups_by_id[block.homolog_group].append(block_with_organism)

    json.dump(list(homolog_groups_by_id.keys()), sys.stdout)

if __name__ == "__main__":
    argh.dispatch_command(main)

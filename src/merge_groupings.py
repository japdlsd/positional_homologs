import json
import os
from collections import defaultdict

import argh

from helpers import BedHomologRecord


@argh.arg("groping_filenames", nargs="+")
def main(groping_filenames, output_directory="merged_groupings"):
    blocks_by_organism = defaultdict(list)
    for filename in groping_filenames:
        with open(filename, "r") as f:
            data = json.load(f)
        blocks = [BedHomologRecord(**record) for record in data["blocks"]]
        grouping = data["grouping"]

        groups = set(grouping)
        if len(groups) != 1:
            for i in range(len(blocks)):
                blocks[i] = blocks[i]._replace(homolog_group=f"{blocks[i].homolog_group}.{grouping[i]}")
        for block in blocks:
            blocks_by_organism[block.organism].append(block)

    os.makedirs(output_directory, exist_ok=True)

    for organism, blocks in blocks_by_organism.items():
        output_filename = f"{output_directory}/{organism}.bed"
        with open(output_filename, "w") as f:
            for block in sorted(blocks_by_organism[organism]):
                chromosome = f"{block.organism}.{block.chromosome}"
                start = block.start
                end = block.end
                name = block.homolog_group
                score = 0
                strand = block.strand
                gene_name = block.gene_name
                line = f"{chromosome}\t{start}\t{end}\t{name}\t{score}\t{strand}\t{gene_name}"
                print(line, file=f)


if __name__ == "__main__":
    argh.dispatch_command(main)

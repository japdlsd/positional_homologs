import json
import os.path
import sys
from collections import defaultdict

import argh
import numpy as np
import yaml
from pathlib import Path
from tqdm import tqdm
import dill

import helpers


def main(
    config_filename,
    group_id,
    output_directory,
    tempdir=None,
    window_radius=5
):
    with open(config_filename) as f:
        config = yaml.safe_load(f)
    if tempdir is None:
        tempdir = config.get("tempdir", "tmp")
        # create the temporary directory
        os.makedirs(tempdir, exist_ok=True)

    # load the input files
    print(f"Loading input files from {os.path.dirname(config_filename)}", file=sys.stderr)

    input_dir = os.path.dirname(config_filename)
    organism_identifiers = list(config["homolog_families"].keys())
    assert set(organism_identifiers) == set(
        config["genomes"].keys()
    ), f"The set of organism identifiers must match the set of genomes: {organism_identifiers} != {set(config['genomes'].keys())}"
    genomes_by_id = {}
    blocks_by_id = {}

    for organism_id in organism_identifiers:
        genome_filename = os.path.join(input_dir, config["genomes"][organism_id])
        blocks_filename = os.path.join(
            input_dir, config["homolog_families"][organism_id]
        )

        genome = helpers.load_fasta(Path(genome_filename))
        blocks = helpers.load_bed_homolog(Path(blocks_filename), organism=organism_id)
        genomes_by_id[organism_id] = genome
        blocks_by_id[organism_id] = sorted(blocks)

    # split the input files into the groups
    homolog_groups_by_id = defaultdict(list)

    for organism_id, blocks in blocks_by_id.items():
        for block in blocks:
            block_with_organism = block._replace(organism=organism_id)
            homolog_groups_by_id[block.homolog_group].append(block_with_organism)

    neighbor_calculators_by_organism = {
        organism_id: helpers.GeneOrderNeighbourCalculator(blocks)
        for organism_id, blocks in blocks_by_id.items()
    }

    chosen_group = homolog_groups_by_id[group_id]

    print(f"Computing the context similarity matrix for group {group_id}", file=sys.stderr)

    if window_radius > 0:
        context_similarity_matrix = np.array(
            helpers.compute_context_similarity_matrix_distinct_set(
                chosen_group, neighbor_calculators_by_organism, radius=window_radius
            )
        )
    elif window_radius == 0:
        context_similarity_matrix = np.ones((len(chosen_group), len(chosen_group)))
    else:
        raise ValueError("window_radius must be >= 0")

    print(f"Outputting the results to {output_directory}", file=sys.stderr)
    os.makedirs(output_directory, exist_ok=True)
    outputs = {
        'similarity_matrix': {
            'cost': context_similarity_matrix.tolist(),
        },
    }
    # adding colors to every output
    colors = [block.organism for block in chosen_group]
    for output in outputs.values():
        output['colors'] = colors
    # adding tasknames to every output
    for output in outputs.values():
        output['taskname'] = group_id
    # adding full info about the blocks to every output
    for output in outputs.values():
        output['blocks'] = [block._asdict() for block in chosen_group]

    # dump the outputs
    for output_name, output in outputs.items():
        with open(os.path.join(output_directory, f"{output_name}.json"), "w") as f:
            print(f"Dumping {output_name} to {f.name}", file=sys.stderr)
            json.dump(output, f)


if __name__ == "__main__":
    argh.dispatch_command(main)

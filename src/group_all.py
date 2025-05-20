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
    output_directory,
    method=None,
    threads=None,
    time_limit=None,
    threads_per_instance=None,
    tempdir=None,
):
    with open(config_filename) as f:
        config = yaml.safe_load(f)
    if method is None:
        method = config.get("method", "he2004-phylo")
    if threads is None:
        threads = config.get("threads", 1)
    if time_limit is None:
        time_limit = config.get("time_limit", 300)
    if threads_per_instance is None:
        threads_per_instance = config.get("threads_per_instance", 1)
    if tempdir is None:
        tempdir = config.get("tempdir", "tmp")
        # create the temporary directory
        os.makedirs(tempdir, exist_ok=True)
    storedir = Path(output_directory) / config.get("storedir", "store")
    storedir.mkdir(parents=True, exist_ok=True)

    if method == "he2004-phylo":
        from m11_he2004_phylo import solve_hybrid as solve
    else:
        raise NotImplementedError(f"method {method} is not implemented")

    # load the input files
    input_dir = os.path.dirname(config_filename)
    organism_identifiers = list(config["homolog_families"].keys())
    assert set(organism_identifiers) == set(
        config["genomes"].keys()
    ), f"The set of organism identifiers must match the set of genomes: {organism_identifiers} != {set(config['genomes'].keys())}"
    assert set(organism_identifiers) == set(
        config["exons"].keys()
    ), f"The set of organism identifiers must match the set of exons: {organism_identifiers} != {set(config['exons'].keys())}"
    genomes_by_id = {}
    blocks_by_id = {}
    exons_by_id = {}

    for organism_id in organism_identifiers:
        genome_filename = os.path.join(input_dir, config["genomes"][organism_id])
        blocks_filename = os.path.join(
            input_dir, config["homolog_families"][organism_id]
        )
        exons_filename = os.path.join(input_dir, config["exons"][organism_id])

        genome = helpers.load_fasta(Path(genome_filename))
        blocks = helpers.load_bed_homolog(Path(blocks_filename), organism=organism_id)
        genomes_by_id[organism_id] = genome
        blocks_by_id[organism_id] = sorted(blocks)

        exons_conflicting = helpers.load_exons(
            Path(exons_filename), organism=organism_id
        )
        exons_by_id[organism_id] = {}
        for gene_id, exons in sorted(exons_conflicting.items()):
            exons_by_id[organism_id][gene_id] = helpers.remove_overlapping_exons(exons)

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

    cost_matrices_by_group_id = get_cost_matrices_sequential(
        exons_by_id,
        genomes_by_id,
        homolog_groups_by_id,
        neighbor_calculators_by_organism,
        storedir=storedir,
    )

    # disambiguate all groups

    # dump the results into the output directory

    raise NotImplementedError()


# def get_cost_matrices_parallel(homolog_groups_by_id):
#
#     group_ids = list(homolog_groups_by_id.keys())
#
#     def compute_cost_matrix_for_group(group_id):
#         group = homolog_groups_by_id[group_id]
#
#         context_similarity_matrix = helpers.compute_context_similarity_matrix(
#             group, neighbor_calculators_by_organism
#         )
#         exon_score_matrix = helpers.compute_exon_score_matrix(
#             group, exons_by_id, genomes_by_id, mode="blastn"
#         )
#         coefficient = np.nansum(exon_score_matrix) * 10
#         cost = np.array(context_similarity_matrix) * coefficient + np.array(
#             exon_score_matrix
#         )
#         return cost


def get_cost_matrices_sequential(
    exons_by_id, genomes_by_id, homolog_groups_by_id, neighbor_calculators_by_organism, storedir
):
    cost_matrices_by_group_id = {}
    for gnum, (group_id, group) in tqdm(enumerate(list(homolog_groups_by_id.items())), desc="Cost matrices"):
        print(f"Group {group_id}: ({gnum}/{len(homolog_groups_by_id)})")
        store_filename = storedir / f"cost_matrix_{group_id}.dill"
        if store_filename.exists():
            with open(store_filename, "rb") as f:
                cost = dill.load(f)
                cost_matrices_by_group_id[group_id] = cost
                continue

        for block in group:
            print(f"  {block}")
        context_similarity_matrix = helpers.compute_context_similarity_matrix(
            group, neighbor_calculators_by_organism
        )
        exon_score_matrix = helpers.compute_exon_score_matrix(
            group, exons_by_id, genomes_by_id, mode="blastn"
        )
        print(f"Exon score matrix: ")
        print(np.array(exon_score_matrix))

        coefficient = np.nansum(exon_score_matrix) * 10
        cost = np.array(context_similarity_matrix) * coefficient + np.array(
            exon_score_matrix
        )

        print(f"Total cost matrix: ")
        print(cost)
        sys.stdout.flush()
        cost_matrices_by_group_id[group_id] = cost
        with open(store_filename, "wb") as f:
            dill.dump(cost, f)
    return cost_matrices_by_group_id


if __name__ == "__main__":
    argh.dispatch_command(main)

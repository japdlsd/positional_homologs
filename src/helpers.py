import os
import subprocess
import tempfile
from collections import namedtuple, defaultdict
from io import StringIO
from pathlib import Path
from typing import Iterator

import networkx as nx
import numpy as np
import scipy.optimize
from Bio import SeqIO, Blast, Phylo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import gffutils

BedHomologRecord = namedtuple(
    "BedHomologRecord",
    ["chromosome", "start", "end", "homolog_group", "strand", "gene_name", "organism"],
    defaults=[None],
)
ExonRecord = namedtuple(
    "ExonRecord",
    ["chromosome", "start", "end", "strand", "gene_id", "exon_id", "organism"],
    defaults=[None],
)


def load_bed_homolog(filename: Path, organism=None) -> Iterator[BedHomologRecord]:
    """
    The columns of BED files are as follows (all are standard except the last one):
    - (0) chromosome
    - (1) start
    - (2) end
    - (3) homolog group (ortholog group)
    - (4) score (_not used_)
    - (5) strand
    - (6) actual gene name (**NOT STANDARD**)

    :param filename:
    :return: Generator of tuples (chromosome, start, end, homologous_group, strand, gene_name)
    """
    with open(filename, "r") as f:
        for line in f:
            if line.strip() == "":
                continue
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            chromosome = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            homolog_group = fields[3]
            strand = fields[5]
            gene_name = fields[6]
            yield BedHomologRecord(
                chromosome, start, end, homolog_group, strand, gene_name, organism
            )


def load_fasta(filename: Path):
    return SeqIO.to_dict(SeqIO.parse(filename, "fasta"))


def remove_overlapping_exons(exons):
    """We just remove shorter exons for any conflicting pairs."""
    exons_by_length = sorted(exons, key=lambda x: x.end - x.start, reverse=True)
    result = [exons_by_length[0]]

    for exon in exons_by_length[1:]:
        if any(exon.start < x.end and exon.end > x.start for x in result):
            continue
        result.append(exon)

    result = sorted(result, key=lambda x: (x.chromosome, x.start, x.end))
    return result


def load_exons(filename: Path | str, organism=None):
    """Returns a dictionary with a list of exons for each gene id."""
    if isinstance(filename, str):
        filename = Path(filename)

    if filename.suffix == ".bed":
        return load_exons_bed(filename, organism=organism)
    if filename.suffix == ".gff":
        return load_exons_gff(filename, organism=organism)
    raise ValueError(f"Unknown file format: {filename}")


def load_exons_bed(filename: Path, organism=None):
    exons_by_gene_id = defaultdict(list)
    with open(filename) as f:
        for line in f:
            if line.strip() == "":
                continue
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            chromosome = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            strand = fields[5]
            exon_id = fields[3]
            assert len(exon_id.split("-")) == 3, f"Invalid exon id: {exon_id}"
            gene_id = exon_id.split("-")[0]

            exon_record = ExonRecord(
                chromosome, start, end, strand, gene_id, exon_id, organism
            )
            exons_by_gene_id[gene_id].append(exon_record)
    return dict(exons_by_gene_id)


def load_exons_gff(filename: Path, organism=None):
    exons_by_gene_id = defaultdict(list)

    sequences_names_to_our_names = {
        # for the AmerMAF2021.gff file
        "CM029764": "2R",
        "CM029763": "2L",
        "CM029766": "3R",
        "CM029765": "3L",
        "CM029762": "X",
    }

    db_filename = str(filename).replace(".gff", ".db")
    if not os.path.exists(db_filename):
        db = gffutils.create_db(
            str(filename),
            dbfn=str(db_filename),
            force=True,
            keep_order=True,
            merge_strategy="merge",
            sort_attribute_values=True,
        )

    db = gffutils.FeatureDB(str(db_filename))
    for exon in db.features_of_type("exon"):
        if exon.seqid not in sequences_names_to_our_names:
            continue
        chromosome = sequences_names_to_our_names[exon.seqid]
        gene_name = exon.attributes["gene_id"][0]
        exon_full_name = exon.attributes["ID"][0]

        start = exon.start
        end = exon.end
        strand = exon.strand
        exon_record = ExonRecord(
            chromosome, start, end, strand, gene_name, exon_full_name, organism
        )
        exons_by_gene_id[gene_name].append(exon_record)
    return dict(exons_by_gene_id)


def compute_coverage_per_position(bed_records, genome):
    coverage = 0
    for record in bed_records:
        coverage += len(genome[record.chromosome][record.start : record.end])
    return coverage


class GeneOrderNeighbourCalculator:
    def __init__(self, bed_records: Iterator[BedHomologRecord]):
        self.bed_records = sorted(
            bed_records, key=lambda x: (x.chromosome, x.start, x.end)
        )
        self.inverse_index = {record: i for i, record in enumerate(self.bed_records)}

    def __call__(self, record: BedHomologRecord, radius: int = 5):
        index = self.inverse_index[record]
        chromosome = record.chromosome
        start = max(0, index - radius)
        end = min(len(self.bed_records), index + radius + 1)
        result = [
            self.bed_records[i]
            for i in range(start, end)
            if self.bed_records[i].chromosome == chromosome and i != index
        ]
        return result

    def distinct_set(self, record: BedHomologRecord, radius: int = 5):
        index = self.inverse_index[record]
        chromosome = record.chromosome
        left_part = set()
        right_part = set()

        # finding distinct homolog groups to the left
        q = index - 1
        while q >= 0 and self.bed_records[q].chromosome == chromosome and len(left_part) < radius:
            current_record = self.bed_records[q]
            if current_record.homolog_group != record.homolog_group:
                left_part.add(current_record.homolog_group)
            q -= 1
        # finding distinct homolog groups to the right
        q = index + 1
        while q < len(self.bed_records) and self.bed_records[q].chromosome == chromosome and len(right_part) < radius:
            current_record = self.bed_records[q]
            if current_record.homolog_group != record.homolog_group:
                right_part.add(current_record.homolog_group)
            q += 1

        result = left_part | right_part
        return result


def blast_two_sequences(
    seq1: str | Seq | SeqRecord,
    seq2: str | Seq | SeqRecord,
    tempdir: str = "tmp",
    blastn_executable="blastn",
    blastx_executable="blastx",
    mode="blastn",
):
    def convert_to_seqrecord(s):
        if isinstance(s, str):
            return SeqRecord(Seq(s))
        if isinstance(s, Seq):
            return SeqRecord(s)
        if isinstance(s, SeqRecord):
            return s
        raise ValueError(f"Cannot convert {s} to SeqRecord")

    seq1 = convert_to_seqrecord(seq1)
    seq2 = convert_to_seqrecord(seq2)

    with (
        tempfile.NamedTemporaryFile(dir=tempdir, mode="w", suffix=".fa") as f1,
        tempfile.NamedTemporaryFile(dir=tempdir, mode="w", suffix=".fa") as f2,
        tempfile.NamedTemporaryFile(dir=tempdir, mode="w+b", suffix=".xml") as f_out,
    ):
        SeqIO.write(seq1, f1, "fasta")
        SeqIO.write(seq2, f2, "fasta")
        f1.flush()
        f2.flush()
        f_out.flush()

        if mode == "blastn" and len(seq1) > 1000 and len(seq2) > 1000:
            cmd = f"{blastn_executable} -query {f1.name} -subject {f2.name} -outfmt 5 -out {f_out.name} -task blastn -dust no"
        elif mode == "blastn":
            cmd = f"{blastn_executable} -query {f1.name} -subject {f2.name} -outfmt 5 -out {f_out.name} -task blastn-short -dust no"
        elif mode == "blastx":
            cmd = f"{blastx_executable} -query {f1.name} -subject {f2.name} -outfmt 5 -out {f_out.name} -task blastx"
        process_results = subprocess.run(cmd, shell=True, check=True)
        assert (
            process_results.returncode == 0
        ), f"Error running blast: {process_results.returncode}"
        result = Blast.read(f_out)
        return result


def context_similarity(record1, record2, neighbor_calculators_by_organism, radius=5):
    neighbors1 = neighbor_calculators_by_organism[record1.organism](
        record1, radius=radius
    )
    neighbors2 = neighbor_calculators_by_organism[record2.organism](
        record2, radius=radius
    )
    context1 = set(neighbor.homolog_group for neighbor in neighbors1)
    context2 = set(neighbor.homolog_group for neighbor in neighbors2)
    assert (
        len(context1 | context2) > 0
    ), f"Empty context for both {record1} and {record2}!"
    jaccard = len(context1 & context2) / len(context1 | context2)
    return jaccard


def context_similarity_distinct_set(record1, record2, neighbor_calculators_by_organism, radius=5):
    context1 = neighbor_calculators_by_organism[record1.organism].distinct_set(
        record1, radius=radius
    )
    context2 = neighbor_calculators_by_organism[record2.organism].distinct_set(
        record2, radius=radius
    )
    assert (
        len(context1 | context2) > 0
    ), f"Empty context for both {record1} and {record2}!"
    jaccard = len(context1 & context2) / len(context1 | context2)
    return jaccard


def compute_context_similarity_matrix(
    group, neighbor_calculators_by_organism, mask_same_organism=True
):
    n = len(group)
    matrix = [[0.0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            if mask_same_organism and group[i].organism == group[j].organism:
                matrix[i][j] = matrix[j][i] = np.nan
                continue

            matrix[i][j] = context_similarity(
                group[i], group[j], neighbor_calculators_by_organism
            )
            matrix[j][i] = matrix[i][j]
    return matrix


def compute_context_similarity_matrix_distinct_set(
    group, neighbor_calculators_by_organism, mask_same_organism=True, radius=5
):
    n = len(group)
    matrix = [[0.0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            if mask_same_organism and group[i].organism == group[j].organism:
                matrix[i][j] = matrix[j][i] = np.nan
                continue

            matrix[i][j] = context_similarity_distinct_set(
                group[i], group[j], neighbor_calculators_by_organism, radius=radius
            )
            matrix[j][i] = matrix[i][j]
    return matrix


def compute_blast_similarity_matrix(group, genomes_by_id, mask_same_organism=True):
    n = len(group)
    matrix = [[0.0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            if mask_same_organism and group[i].organism == group[j].organism:
                matrix[i][j] = matrix[j][i] = np.nan
                continue

            seq1 = genomes_by_id[group[i].organism][group[i].chromosome][
                group[i].start : group[i].end
            ]
            seq2 = genomes_by_id[group[j].organism][group[j].chromosome][
                group[j].start : group[j].end
            ]
            blast_results = blast_two_sequences(seq1, seq2)

            if len(blast_results) == 0:
                continue

            best_score = max(hsp.score for hsp in blast_results[0])
            similarity = best_score
            matrix[j][i] = matrix[i][j] = similarity
    return matrix


def compute_exon_score_two_genes(exon_set1, exon_set2, genomes_by_id, mode="blastn"):
    G = nx.Graph()
    # nodes are exons from both groups
    # edges are between exons from first and second set
    # the edge weight is the bit score of alignment between the two exons
    for exon1 in exon_set1:
        for exon2 in exon_set2:
            seq1 = genomes_by_id[exon1.organism][exon1.chromosome][
                exon1.start : exon1.end
            ]
            seq2 = genomes_by_id[exon2.organism][exon2.chromosome][
                exon2.start : exon2.end
            ]
            blast_results = blast_two_sequences(seq1, seq2, mode=mode)
            if len(blast_results) == 0:
                continue
            best_score = max(hsp.score for hsp in blast_results[0])
            G.add_edge(exon1, exon2, weight=best_score)
    # let's compute the maximum weight matching
    matching = nx.max_weight_matching(G)
    # the score is the sum of the weights of the matching edges
    score = sum(G[exon1][exon2]["weight"] for exon1, exon2 in matching)
    return score


def compute_exon_score_two_genes_assignment(
    exon_set1, exon_set2, genomes_by_id, mode="blastn"
):
    cost_matrix = np.zeros((len(exon_set1), len(exon_set2)))
    # nodes are exons from both groups
    # edges are between exons from first and second set
    # the edge weight is the bit score of alignment between the two exons
    for i1, exon1 in enumerate(exon_set1):
        for i2, exon2 in enumerate(exon_set2):
            seq1 = genomes_by_id[exon1.organism][exon1.chromosome][
                exon1.start : exon1.end
            ]
            seq2 = genomes_by_id[exon2.organism][exon2.chromosome][
                exon2.start : exon2.end
            ]
            blast_results = blast_two_sequences(seq1, seq2, mode=mode)
            if len(blast_results) == 0:
                continue
            best_score = max(hsp.score for hsp in blast_results[0])
            cost_matrix[i1, i2] = best_score

    # let's compute the maximum assignment
    row_indices, col_indices = scipy.optimize.linear_sum_assignment(
        cost_matrix, maximize=True
    )
    score = cost_matrix[row_indices, col_indices].sum()
    return score


def compute_exon_score_two_genes_assignment_normalized(
    exon_set1, exon_set2, genomes_by_id, mode="blastn", tempdir="tmp"
):
    cost_matrix = np.zeros((len(exon_set1), len(exon_set2)))
    # nodes are exons from both groups
    # edges are between exons from first and second set
    # the edge weight is the bit score of alignment between the two exons
    for i1, exon1 in enumerate(exon_set1):
        for i2, exon2 in enumerate(exon_set2):
            seq1 = genomes_by_id[exon1.organism][exon1.chromosome][
                exon1.start : exon1.end
            ]
            seq2 = genomes_by_id[exon2.organism][exon2.chromosome][
                exon2.start : exon2.end
            ]
            bss_xy = blast_two_sequences(seq1, seq2, mode=mode, tempdir=tempdir)
            bss_yx = blast_two_sequences(seq2, seq1, mode=mode, tempdir=tempdir)
            bss_xx = blast_two_sequences(seq1, seq1, mode=mode, tempdir=tempdir)
            bss_yy = blast_two_sequences(seq2, seq2, mode=mode, tempdir=tempdir)

            bs_xy = max(hsp.score for hsp in bss_xy[0]) if len(bss_xy) > 0 else 0
            bs_yx = max(hsp.score for hsp in bss_yx[0]) if len(bss_yx) > 0 else 0
            bs_xx = max(hsp.score for hsp in bss_xx[0]) if len(bss_xx) > 0 else 0
            bs_yy = max(hsp.score for hsp in bss_yy[0]) if len(bss_yy) > 0 else 0

            similarity = (bs_xy + bs_yx) / (bs_xx + bs_yy) if bs_xx + bs_yy > 0 else 0
            cost_matrix[i1, i2] = similarity

    # let's compute the maximum assignment
    row_indices, col_indices = scipy.optimize.linear_sum_assignment(
        cost_matrix, maximize=True
    )
    score = cost_matrix[row_indices, col_indices].sum()
    return score


def compute_exon_score_matrix_matching(
    group,
    exons_by_id,
    genomes_by_id,
    mask_same_organism=True,
    tqdm_bar=None,
    mode="blastn",
):
    n = len(group)
    matrix = [[0.0 for _ in range(n)] for _ in range(n)]
    if tqdm_bar is not None:
        tqdm_bar.total = n * (n + 1) // 2
    for i in range(n):
        for j in range(i, n):
            if tqdm_bar is not None:
                tqdm_bar.update(1)
            if mask_same_organism and group[i].organism == group[j].organism:
                matrix[i][j] = matrix[j][i] = np.nan
                continue

            exons1 = exons_by_id[group[i].organism][group[i].gene_name]
            exons2 = exons_by_id[group[j].organism][group[j].gene_name]
            score = compute_exon_score_two_genes_assignment_normalized(
                exons1, exons2, genomes_by_id, mode=mode
            )
            matrix[j][i] = matrix[i][j] = score
    return matrix


def compute_exon_score_matrix(
    group,
    exons_by_id,
    genomes_by_id,
    mask_same_organism=True,
    tqdm_bar=None,
    mode="blastn",
    tempdir="tmp",
):
    n = len(group)
    matrix = [[0.0 for _ in range(n)] for _ in range(n)]
    if tqdm_bar is not None:
        tqdm_bar.total = n * (n + 1) // 2
    for i in range(n):
        for j in range(i, n):
            if tqdm_bar is not None:
                tqdm_bar.update(1)
            if mask_same_organism and group[i].organism == group[j].organism:
                matrix[i][j] = matrix[j][i] = np.nan
                continue

            exons1 = exons_by_id[group[i].organism][group[i].gene_name]
            exons2 = exons_by_id[group[j].organism][group[j].gene_name]

            assert group[i].organism in genomes_by_id, f"Organism {group[i].organism} not found in genomes_by_id! Available organisms are: {genomes_by_id.keys()}"
            assert group[i].chromosome in genomes_by_id[group[i].organism], f"Chromosome {group[i].chromosome} not found in genomes_by_id[{group[i].organism}]! Available chromosomes are: {list(genomes_by_id[group[i].organism].keys())}"
            assert group[j].organism in genomes_by_id, f"Organism {group[j].organism} not found in genomes_by_id! Available organisms are: {genomes_by_id.keys()}"
            assert group[j].chromosome in genomes_by_id[group[j].organism], f"Chromosome {group[j].chromosome} not found in genomes_by_id[{group[j].organism}]! Available chromosomes are: {list(genomes_by_id[group[j].organism].keys())}"



            seq1 = "".join(
                str(genomes_by_id[group[i].organism][group[i].chromosome][
                    exon.start : exon.end
                ].seq)
                for exon in sorted(exons1)
            )
            seq2 = "".join(
                str(genomes_by_id[group[j].organism][group[j].chromosome][
                    exon.start : exon.end
                ].seq)
                for exon in sorted(exons2)
            )

            bss_xy = blast_two_sequences(seq1, seq2, mode=mode, tempdir=tempdir)
            bss_yx = blast_two_sequences(seq2, seq1, mode=mode, tempdir=tempdir)
            bss_xx = blast_two_sequences(seq1, seq1, mode=mode, tempdir=tempdir)
            bss_yy = blast_two_sequences(seq2, seq2, mode=mode, tempdir=tempdir)

            bs_xy = max(hsp.score for hsp in bss_xy[0]) if len(bss_xy) > 0 else 0
            bs_yx = max(hsp.score for hsp in bss_yx[0]) if len(bss_yx) > 0 else 0
            bs_xx = max(hsp.score for hsp in bss_xx[0]) if len(bss_xx) > 0 else 0
            bs_yy = max(hsp.score for hsp in bss_yy[0]) if len(bss_yy) > 0 else 0

            similarity = (bs_xy + bs_yx) / (bs_xx + bs_yy) if bs_xx + bs_yy > 0 else 0

            matrix[j][i] = matrix[i][j] = similarity
    return matrix


def load_phylogeny_tree_from_file(filename: Path):
    with open(filename, 'r') as f:
        s = f.read().strip()
        return load_phylogeny_tree_from_string(s)


def load_phylogeny_tree_from_string(s: str):
    f = StringIO(s)
    G_raw = Phylo.read(f, 'newick')
    next_unused_node = 0
    for node in G_raw.get_nonterminals():
        if node.name == "" or node.name is None:
            node.name = f"internal_{next_unused_node}"
            next_unused_node += 1
    G = Phylo.to_networkx(G_raw)

    root = list(G.nodes)[0]
    G2 = nx.dfs_tree(G, root)
    # add the lengths to the edges
    for u, v in G2.edges:
        G2.edges[(u, v)]['length'] = v.branch_length
        assert v.branch_length is not None, f"Branch length is None for edge {u} -> {v}!!"
    # add the true names to the edges
    for u, v in G2.edges:
        G2.edges[(u, v)]['true_name'] = (u.name, v.name)

    # now, let's remove those idiotic Clade() objects and replace them with strings
    G3 = nx.DiGraph()
    for u, v in G2.edges:
        G3.add_edge(u.name, v.name, length=G2.edges[(u, v)]['length'], true_name=G2.edges[(u, v)]['true_name'])
    return G3



def test_parse_newick_with_lengths():
    s = "(A:1,B:2,(C:3,D:4):5)"
    G = load_phylogeny_tree_from_string(s)
    print()
    print(G.nodes)
    print(G.edges)
    print(f"leaves: {[node for node in G.nodes if G.out_degree(node) == 0]}")
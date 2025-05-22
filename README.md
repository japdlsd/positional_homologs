# 🧬 Positional Ortholog Identifier

**Positional Ortholog Identifier** is a computational tool designed to identify **positional orthologs** within homologous gene families. Unlike traditional sequence similarity approaches, this method integrates **gene order (synteny)** information and family-specific context to determine which genes are true positional counterparts across genomes.

This tool is particularly useful for:
- Resolving ambiguous one-to-many or many-to-many homology relationships.
- Refining ortholog predictions in large gene families.
- Comparative genomics, evolutionary analysis, and genome annotation pipelines.

### 🔍 Key Features
- **Input**: Predefined homologous gene families across multiple species or genomes.
- **Output**: A set of predicted positional orthologs for each gene family.
- **Algorithms**:
  - Exact integer programming formulation (ILP) for optimal solutions.
  - Efficient heuristic methods for large families.
  - Adaptation of the method by *He et al. (2004)*.

### ⚙️ Use Cases
- Genome evolution studies.
- Functional annotation transfer.
- Phylogenomic reconstructions.

## Overview of the method

**TL;DR** Each gene family is processed individually.
For each gene pair inside the gene family, we compute two similarity scores: 1) **local context similarity**, and 2) **sequence similarity**.
Those are then combined into **a single composite score**. 
We then find a way to split the gene family into clusters of positional orthologs such that the sum of the composite scores inside those clusters is maximized, and no clusters contain two genes from the same organism,
This problem is formulated as a **Maximum Colorful Graph Partition** problem, which is solved using either an exact integer programming formulation or a heuristic method.

#### Local context similarity

Our method is based on the clustering of genes based on their genomic context; two genes that have very similar gene sets in their neighborhood are more likely to be positional orthologs.
Let $A = (a_1, \ldots, a_n)$ denote the gene order of a circular chromosome, where $a_i$ is the identifier of the homolog group of the $i$-th gene in the chromosome.
Fix a position $i$ in the gene order.
Consider the size $k+1$ set of homology group identifiers to the left of and including $a_i$.
We refer to this set, without $a_i$, as the *left $k$-local context* of the gene $i$.
Analogously, we define the right $k$-local context.
The *$k$-local context* of gene $i$ is defined as the union of its left and right $k$-local contexts.
For linear chromosomes, the left and right contexts may contain fewer than $k$ elements if the position is near an extremity.
In this context, we can also refer to $k$ as the *radius* (or *window radius*) of the local context.
We define the *$k$-local context similarity* between two genes as the Jaccard similarity between their $k$-local contexts.

As an example, consider the gene order $(2, 1, 3, 3, 4, 5, 3, 5, 7, 2, 8, 2, 4, 10, 1)$.
The left $3$-local context of the 6th gene (with homolog group 5) is $\{4, 3, 1\}$, and its right $3$-local context is $\{3, 7, 2\}$. 
Note that 5 is skipped in the right context, as it matches the homolog group of the gene under consideration.
The combined $3$-local context is therefore $\{1, 2, 3, 4, 7\}$.
The $3$-local context of the 8th gene (also from homolog group 5) is $\{1, 2, 3, 4, 7, 8\}$.
Their $3$-local context similarity is thus $5/6$.

To break ties between gene pairs with identical local context similarity, we use a secondary similarity measure based on sequence similarity.

#### Sequence similiarity

Let $a$ and $b$ denote the concatenated exon sequences of the two genes.
Let $BS(a, b)$ be the bit-score of the best BLAST alignment between $a$ and $b$.
Define the reciprocal score between $a$ and $b$ as
$$d(a, b) = \dfrac{BS(a, b) + BS(b, a)}{BS(a, a) + BS(b, b)}$$
if the denominator is positive, and zero otherwise.
We use the reciprocal score as the *sequence similarity* of the two genes.

#### Composing scores into a single composite score

The *composite similarity score* of two genes within a homolog group is 
the ordered pair consisting of their $k$-local context similarity and their sequence similarity.
Most of the methods we apply later operate on scalar similarity scores.
To this end, we define the *linearized similarity score* between two genes as $10 S \cdot x + y$, where $(x, y)$ is the composite similarity score and $S$ is the sum of all sequence similarity scores between pairs of genes in that homolog group.

### Positional ortholog identification as an optimization problem

We formalize the problem of identifying positional orthologs as a maximization problem we 
call the **Maximum Colorful Graph Partition** problem.

Consider a graph $G = (V, E)$ with a non-negative edge weight function $w: E \to \mathbf{R}_{+}$
and a vertex labeling function $\lambda: V \to \mathbf{N}$.
Let $X \subseteq V$ be a subset of vertices, and let $G_X$ be the corresponding induced subgraph 
(i.e., $G_X$ contains all vertices from $X$ and all edges of $G$ that connect vertices in $X$).
We denote the sum of the weights of edges in $G_X$ as
$$
w(X) := \sum_{e \in E(G_X)} w(e).
$$
We call a subset $X$ *colorful* if all labels of its vertices are distinct.
We call a partition
$\Xi$ of the vertices $V$ *colorful* if all its subsets are colorful.
We denote the total weight of a partition $\Xi$ as
$$
w(\Xi) := \sum_{X \in \Xi} w(X),
$$
and refer to this quantity as the *weight of the partition* $\Xi$.
The individual sets in the partition are called *clusters*.

#### Maximum Colorful Graph Partition problem:
_Given a graph $G = (V, E)$ with a non-negative edge weight function 
$w: E \to \mathbf{R}_{+}$ and a vertex labeling function $\lambda: V \to \mathbf{N}$, 
find a colorful partition $\Xi$ of $V$ that maximizes the total weight $w(\Xi)$._

In our setting, the vertices represent individual genes from a given homolog group, 
the labeling function maps each gene to its organism, and the edge weight function assigns 
the linearized similarity score between the corresponding genes.
Each set in the resulting colorful partition then corresponds to a cluster of positional orthologs.

## INSTALL

```bash
mamba env update -f requirements.yml
conda activate positional_homologs

# install the BLAST tool
sudo apt install ncbi-blast+
```


### first-time install 

```bash
mamba env update -f requirements_crude.yml
conda activate positional_homologs
conda env export | egrep -v '^prefix: ' > requirements.yml

# install the BLAST tool
sudo apt install ncbi-blast+
```

## Getting the data

```bash
cd raw_data
# getting the Anopheles data
cd anopheles2
bash download_genomes.sh
```

## Input specification

### snakemake config file

The config file is a YAML file that specifies the input data and parameters for the workflow. 
Usually located in the `disambiguation` directory.

The snakemake config file should contain the following keys:
- `output_dir`: the relative address of the output directory (to the working directory)
- `input_data_config`: the relative address of the input data config file (to the working directory)
- `all_gene_families`: the relative address of a JSON file that contains a single list of all gene families' names (to the working family).

Additionally, it may contain the following parameters:
- `window_radii`: the list of window radii for the context similarity (default=[5]).
- `time_limit`: the time limit for the disambiguation process (in seconds) **for a single gene family** (default=180).
- `selected_models`: a list of models to use to compute the stuff. Available models:
  - 00-bad-model
  - 00b-single-blob
  - 01-gurobi-FW (requires Gurobi licence!)
  - 03b-gurobi-native-components (requires Gurobi licence!)
  - 04b-gurobi-native-components-pruning (requires Gurobi licence!)
  - 05-greedy
  - 07-iterative-color
  - 11-he2004-phylo
  - 11a-he2004-phylo-filter
  - 12-star (requires Gurobi licence!)

Example: `disambiguation/anopheles2.yml`.

### input data config file

The input data config file is a YAML file that specifies the input data for the workflow.
Usually located in the `raw_data/<experiment>` directory.

The input data config file should contain the following keys:
- `homolog_families`: a dictionary that maps organism names to span BED files of those organisms.
- `exons`: a dictionary that maps organism names to exon BED files of those organisms.
- `genomes`: a dictionary that maps organism names to genome FASTA files of those organisms.
- `phylogeny`: name of the phylogeny file in the Newick format.

All file paths in this file are relative to the **input data config file**.

Example: `raw_data/anopheles2/input_config.yml`.

### All gene families file

A JSON file that contains a single list of all gene families' names.
Usually located in the `raw_data/<experiment>` directory.

Example: `raw_data/anopheles2/all_gene_families.json`.

### span BED file

A BED file where each line corresponds to a gene in a particular organism. The columns are:

- `chrom`: the chromosome name
- `start`: the start position of the gene
- `end`: the end position of the gene
- `name`: the name of the **gene family**
- `score`: ignored
- `strand`: the strand of the gene
- `gene_name`: the name of the gene inside that organism

Examples: `raw_data/anopheles2/*.span.bed`.

### exon BED file

A BED file where each line corresponds to a gene exon in a particular organism. The columns are:

- `chrom`: the chromosome name
- `start`: the start position of the exon
- `end`: the end position of the exon
- `name`: the name of the exon in the format `<gene_name>-R<transcript_name>-E<exon_number_inside_transcript>`
- `score`: ignored
- `strand`: the strand of the gene

Examples: `raw_data/anopheles2/*.exon.bed`.

## Output specification

The output of the workflow is a directory that contains the following files:

- `exons_longest_transcript/<organism_name>.bed` - a BED file that contains the exons of the longest transcript (by total basepair length) for each gene in the organism
- `blast_matrix/<gene_family_id>/blast_matrix.json` - a JSON file that contains the BLAST matrix for the gene family
- `similarity_matrix/<gene_family_id>/wr_<window_radius>/similarity_matrix.json` - a JSON file that contains the similarity matrix for the gene family
- `cost_matrix/<gene_family_id>/wr_<window_radius>/single_cost.json` - a JSON file that contains the composite cost matrix for the gene family
- `gene_families_disambiguated/<model>/wr_<window_radius>/<gene_family_id>/additional_data.json` - a JSON file that contains the data for the gene family disambiguation
- `results/<model>/wr_<window_radius>/<organism_name>.bed` - a BED file that contains the disambiguated gene spans for the organism
- `result_mapping_krister_format/<model>/wr_<window_radius>/mapping.json` - a JSON file that contains the mapping of the gene families to the disambiguated gene families


## Run

See `disambiguation/README.md` file for instructions on how to run the disambiguation workflow.
 
## Project Structure

- `raw_data/` - raw data
  - `anopheles2/` - Anopheles data
- `disambiguation/` - full disambiguation of the homologs
- `src/` - source code
- `workflows/` - snakemake workflows

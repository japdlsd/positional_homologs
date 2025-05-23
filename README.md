# 🧬 Positional Ortholog Identifier

Have you ever tried to analyze gene order evolution, but found that an ortholog gene family contains far too many genes, and now you wonder, which one of three genes from organism X matches that one out of four genes from organism Y?
Well, you are not alone. And you are in the right place!

What you need is to disambiguate the gene family into **positional orthologs**.
Positional orthologs are genes in different species that not only _evolved from a common ancestral gene_ 
(making them orthologs), but _also have retained the same position in their respective genomes_ as that ancestral gene. 
In other words, positional orthologs are orthologous genes that still occupy the same genomic neighborhood or 
context as their ancestor did, despite the evolutionary changes that may have occurred elsewhere in the genome

**Positional Ortholog Identifier** is a computational tool designed 
to identify positional orthologs within homologous gene families. 
Unlike traditional sequence similarity approaches, this method integrates **gene order (synteny)** 
information and family-specific context to determine which genes are true positional counterparts across genomes.


This tool is particularly useful for:
- Analyzing gene order evolution.
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


## Overview of the method


Each gene family is processed individually.
For each gene pair inside the gene family, we compute two similarity scores: 1) **local context similarity**, 
and 2) **sequence similarity**.
Those are then combined into **a single composite score**. 
We then find a way to split the gene family into clusters of positional orthologs such that the sum of the composite scores inside those clusters is maximized, and no cluster contains two genes from the same organism.
This problem is formulated as a **Maximum Colorful Graph Partition** problem (see the paper), which is solved using either an exact integer programming formulation or a heuristic method.


#### Local context similarity

The local context of a gene is a set of gene families identifers that are located in the vicinity of the gene.
The **window radius** of the local context is the number of distinct identifiers to the left and right of the gene that are included in the context.
The local context similarity between two genes is then defined as the **Jaccard similarity** of their local contexts.


To break ties between gene pairs with identical local context similarity, we use a secondary similarity measure based on sequence similarity.

#### Sequence similiarity

Let $a$ and $b$ denote the concatenated exon sequences of the two genes.
In case of multiple transcripts, the ones with the longest total basepair lengths are used.
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


## Overview of the workflow

The workflow is implemented in Snakemake, a system for creating reproducible and scalable data analyses.
The workflow consists of individual _rules_, which specify how to create a specific output file from a set of input files.
So, the only thing you need to do is to specify the input data and the parameters for the workflow, and then run it using Snakemake.

In terms of the input files, you will need the following things:

- FASTA files of the genomes for each organism
- [span BED files](#span-bed-file) of the genes with their gene families for each organism
- [exon BED files](#exon-bed-file) of the genes with their exons for each organism
- a phylogeny file of the organisms in the Newick format
- [a list of gene families](#all-gene-families-file) you want to analyze in JSON format
- [an input data config file](#input-data-config-file) in YAML format, which matches all aforementioned files together
- [a Snakemake config file](#snakemake-config-file) in YAML format, which specifies the input data config file location and the parameters for the workflow

The workflow will then create the BED [files with the disambiguated gene families](#output-specification) for each organism, as well as some intermediate files that are used during the workflow.

The rest of this README file will describe how to [set up the environment](#install), [prepare](#input-specification) the input files, [run](#run) the workflow, and what the [output](#output-specification) files look like.


### Install

```bash
mamba env update -f requirements.yml
conda activate positional_homologs

# install the BLAST tool
sudo apt install ncbi-blast+
```

### Getting the data (for the Anopheles demo)

```bash
cd raw_data
# getting the Anopheles data
cd anopheles2
bash download_genomes.sh
```

### Input specification

#### snakemake config file

The config file is a YAML file that specifies the input data and parameters for the workflow. 
Usually located in the `disambiguation` directory.

The snakemake config file should contain the following keys:
- `output_dir`: the relative address of the output directory (to the working directory)
- `input_data_config`: the relative address of the input data config file (to the working directory)
- `all_gene_families`: the relative address of a JSON file that contains a single list of all gene families' names (to the working family).
  - this file is used to limit the number of gene families to process in case you want to make a test run.
  - there is a script `src/list_all_gene_families.py` that can be used to create this file from the input data config file.

Additionally, it may contain the following parameters:
- `window_radii`: the list of window radii for the context similarity (default=`[5]`).
- `time_limit`: the time limit for the disambiguation process (in seconds) **for a single gene family** (default=`180`).
- `selected_models`: a list of models to use to compute the stuff. Available models:
  - `00-bad-model`
  - `00b-single-blob`
  - `01-gurobi-FW` (requires Gurobi licence!)
  - `03b-gurobi-native-components` (requires Gurobi licence!)
  - `04b-gurobi-native-components-pruning` (requires Gurobi licence!)
  - `05-greedy`
  - `07-iterative-color`
  - `11-he2004-phylo`
  - `11a-he2004-phylo-filter`
  - `12-star` (requires Gurobi licence!)

Example: `disambiguation/anopheles2.yml`.

#### input data config file

The input data config file is a YAML file that specifies the input data for the workflow.
Usually located in the `raw_data/<experiment>` directory.

The input data config file should contain the following keys:
- `homolog_families`: a dictionary that maps organism names to span BED files of those organisms.
- `exons`: a dictionary that maps organism names to exon BED files of those organisms.
- `genomes`: a dictionary that maps organism names to genome FASTA files of those organisms.
- `phylogeny`: name of the phylogeny file in the Newick format.

All file paths in this file are relative to the **input data config file**.

Example: `raw_data/anopheles2/input_config.yml`.

#### All gene families file

A JSON file that contains a single list of all gene families' names.
Usually located in the `raw_data/<experiment>` directory.

Example: `raw_data/anopheles2/all_gene_families.json`.

#### span BED file

A BED file where each line corresponds to a gene in a particular organism. The columns are:

- `chrom`: the chromosome name
- `start`: the start position of the gene
- `end`: the end position of the gene
- `name`: the name of the **gene family**
- `score`: ignored
- `strand`: the strand of the gene
- `gene_name`: the name of the gene inside that organism

Examples: `raw_data/anopheles2/*.span.bed`.

#### exon BED file

A BED file where each line corresponds to a gene exon in a particular organism. The columns are:

- `chrom`: the chromosome name
- `start`: the start position of the exon
- `end`: the end position of the exon
- `name`: the name of the exon in the format `<gene_name>-R<transcript_name>-E<exon_number_inside_transcript>`
- `score`: ignored
- `strand`: the strand of the gene

Examples: `raw_data/anopheles2/*.exon.bed`.

### Run

```bash
cd disambiguation
conda activate positional_homologs
snakemake -j<number_of_threads> --snakefile ../workflows/disambiguation.smk --configfile <your_config>.yml --printshellcmds --scheduler=greedy
```

For anopheles data, you can use the following command to run the workflow with 16 threads:

```bash
cd disambiguation
conda activate positional_homologs
snakemake -j16 --snakefile ../workflows/disambiguation.smk --configfile anopheles2.yml --printshellcmds --scheduler=greedy
```


### Output specification

These are probably the only files you are interested in:
- `results/<model>/wr_<window_radius>/<organism_name>.bed` - a BED file that contains the **disambiguated gene spans for the organism**
  - the only difference from the input span BED file is that the gene family names are replaced with the disambiguated gene family names.
  Specifically, a gene family name is extended by a suffix `.<cluster_id>`, where `<cluster_id>` is a number.
  Note that if there is only one cluster, the suffix is omitted.

Additionally, there are some other files that are generated during the workflow:

- `exons_longest_transcript/<organism_name>.bed` - a BED file that contains the exons of the longest transcript (by total basepair length) for each gene in the organism
- `blast_matrix/<gene_family_id>/blast_matrix.json` - a JSON file that contains the BLAST matrix for the gene family
- `similarity_matrix/<gene_family_id>/wr_<window_radius>/similarity_matrix.json` - a JSON file that contains the similarity matrix for the gene family
- `cost_matrix/<gene_family_id>/wr_<window_radius>/single_cost.json` - a JSON file that contains the composite cost matrix for the gene family
- `gene_families_disambiguated/<model>/wr_<window_radius>/<gene_family_id>/additional_data.json` - a JSON file that contains the data for the gene family disambiguation
- `result_mapping_krister_format/<model>/wr_<window_radius>/mapping.json` - a JSON file that contains the mapping of the gene families to the disambiguated gene families




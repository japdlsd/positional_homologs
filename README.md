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

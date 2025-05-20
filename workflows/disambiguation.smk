import json
import os.path

import yaml

output_dir = config["output_dir"]
all_gene_families_filename = config["all_gene_families"]

with open(all_gene_families_filename) as f:
    all_gene_families = json.load(f)

input_config_filename = config["input_data_config"]
DISAMBIGUATE_GROUP_TIME_LIMIT = config.get("time_limit", 180)  # Default to 3 minutes if not specified

window_radii = config.get("window_radii", [5])


def parse_input_config(input_config_filename):
    with open(input_config_filename) as f:
        input_config = yaml.safe_load(f)

    genomes = input_config["genomes"]
    blocks = input_config["homolog_families"]
    exons_raw = input_config["exons"]
    phylogeny = input_config["phylogeny"]

    def add_dir_prefix(d, prefix):
        return {k: f"{prefix}/{v}" for k, v in d.items()}

    dir_prefix = os.path.dirname(input_config_filename)
    genomes = add_dir_prefix(genomes,dir_prefix)
    blocks = add_dir_prefix(blocks,dir_prefix)
    exons_raw = add_dir_prefix(exons_raw,dir_prefix)
    phylogeny = os.path.join(dir_prefix,phylogeny)

    return genomes, blocks, exons_raw, phylogeny


genomes, blocks, exons_raw, phylogeny = parse_input_config(input_config_filename)
exons = {k: f"{output_dir}/exons_longest_transcript/{k}.bed" for k, _ in exons_raw.items()}

models = {
    "11-he2004-phylo": {
        "source": "../src/m11_he2004_phylo.py",
        "threads_limit": 1,
        "phylo": True,
    },
    "11a-he2004-phylo-filter": {
        "source": "../src/m11a_he2004_phylo_filter.py",
        "threads_limit": 1,
        "phylo": True,
    },
    "01-gurobi-FW": {
        "source": "../src/m01_gurobi_FW.py",
        "threads_limit": 8,
        "phylo": False,
    },
    "03b-gurobi-native-components": {
        "source": "../src/m03b_gurobi_native_components.py",
        "threads_limit": 8,
        "phylo": False,
    },
    "04b-gurobi-native-components-pruning": {
        "source": "../src/m04b_gurobi_native_components_pruning.py",
        "threads_limit": 8,
        "phylo": False,
    },
    "05-greedy": {
        "source": "../src/m05_greedy.py",
        "threads_limit": 1,
        "phylo": False,
    },
    "07-iterative-color": {
        "source": "../src/m07_iterative_color.py",
        "threads_limit": 1,
        "phylo": False,
    },
    "12-star": {
        "source": "../src/m12_star.py",
        "threads_limit": 8,
        "phylo": False,
    },
    "00-bad-model": {
        "source": "../src/m00_bad_model.py",
        "threads_limit": 1,
        "phylo": False,
    },
    "00b-single-blob": {
        "source": "../src/m00b_single_blob.py",
        "threads_limit": 1,
        "phylo": False,
    }
}

if "selected_models" in config:
    selected_models = config["selected_models"]
    models = {k: v for k, v in models.items() if k in selected_models}
else:
    selected_models = models.keys()



rule exons_longest_transcript:
    input: lambda w: exons_raw[w.organism]
    log: f"{output_dir}/exons_longest_transcript/{{organism}}.log"
    output: f"{output_dir}/exons_longest_transcript/{{organism}}.bed"
    shell: """python ../src/exons_longest_transcript.py {input} {wildcards.organism} {output} > {log} 2>&1"""

rule similarity_matrix:
    input:
        genomes=genomes.values(),
        blocks=blocks.values()
    params:
        target_directory=lambda w: f"{output_dir}/similarity_matrix/wr_{w.window_radius}/{w.group_id}",
        window_radius=lambda w: w.window_radius,
    log: f"{output_dir}/similarity_matrix/wr_{{window_radius}}/{{group_id}}/log"
    output: f"{output_dir}/similarity_matrix/wr_{{window_radius}}/{{group_id}}/similarity_matrix.json"
    shell: """python ../src/make_similarity_matrix.py {input_config_filename} {wildcards.group_id} {params.target_directory} --window-radius {params.window_radius} > {log} 2>&1"""

rule blast_matrix:
    input:
        genomes=genomes.values(),
        blocks=blocks.values(),
        exons=exons.values()
    params:
        target_directory=lambda w: f"{output_dir}/blast_matrix/{w.group_id}"
    log: f"{output_dir}/blast_matrix/{{group_id}}/log"
    output: f"{output_dir}/blast_matrix/{{group_id}}/blast_matrix.json"
    shell: """python ../src/make_blast_matrix.py {input_config_filename} {wildcards.group_id} {params.target_directory} > {log} 2>&1"""

rule cost_matrix:
    input:
        similarity_matrix=f"{output_dir}/similarity_matrix/wr_{{window_radius}}/{{group_id}}/similarity_matrix.json",
        blast_matrix=f"{output_dir}/blast_matrix/{{group_id}}/blast_matrix.json"
    log: f"{output_dir}/cost_matrix/wr_{{window_radius}}/{{group_id}}/log"
    output:
        single_cost=f"{output_dir}/cost_matrix/wr_{{window_radius}}/{{group_id}}/single_cost.json"
    shell: """python ../src/merge_cost_matrices.py {input.similarity_matrix} {input.blast_matrix} {output.single_cost} > {log} 2>&1"""


rule all_cost_matrix:
    input: expand(f"{output_dir}/cost_matrix/wr_{{window_radius}}/{{group_id}}/single_cost.json",group_id=all_gene_families, window_radius=window_radii)


rule disambiguate_group:
    input:
        single_cost=f"{output_dir}/cost_matrix/wr_{{window_radius}}/{{group_id}}/single_cost.json",
        phylogeny=phylogeny,
        model_script=lambda w: models[w.model]["source"]
    params:
        phylo_input=lambda w: (phylogeny if models[w.model].get("phylo", False) else ""),
        target_directory=f"{output_dir}/gene_families_disambiguated/{{model}}/wr_{{window_radius}}/{{group_id}}",
        time_limit=DISAMBIGUATE_GROUP_TIME_LIMIT
    threads: lambda w: models[w.model].get("threads_limit",1)
    output:
        f"{output_dir}/gene_families_disambiguated/{{model}}/wr_{{window_radius}}/{{group_id}}/grouping.json",
        f"{output_dir}/gene_families_disambiguated/{{model}}/wr_{{window_radius}}/{{group_id}}/additional_data.json"
    log:
        error_log=f"{output_dir}/gene_families_disambiguated/{{model}}/wr_{{window_radius}}/{{group_id}}/error_log",
        gurobi_log=f"{output_dir}/gene_families_disambiguated/{{model}}/wr_{{window_radius}}/{{group_id}}/gurobi.log",
    benchmark:
        f"{output_dir}/gene_families_disambiguated/{{model}}/wr_{{window_radius}}/{{group_id}}/measurements.txt"
    shell: """python {input.model_script} {input.single_cost} {params.phylo_input} {params.target_directory} --threads {threads} --time-limit {params.time_limit} > {log.error_log} 2>&1"""


rule parse_statistics:
    input:
        single_cost=f"{output_dir}/cost_matrix/wr_{{window_radius}}/{{group_id}}/single_cost.json",
        grouping=f"{output_dir}/gene_families_disambiguated/{{model}}/wr_{{window_radius}}/{{group_id}}/grouping.json",
        measurements=f"{output_dir}/gene_families_disambiguated/{{model}}/wr_{{window_radius}}/{{group_id}}/measurements.txt",
        additional_data=f"{output_dir}/gene_families_disambiguated/{{model}}/wr_{{window_radius}}/{{group_id}}/additional_data.json",
        code="../src/parse_statistics.py"
    params:
        task_label=lambda w: w.group_id,
        model_label=lambda w: w.model
    output: f"{output_dir}/gene_families_disambiguated/{{model}}/wr_{{window_radius}}/{{group_id}}/statistics.json"
    shell: "python ../src/parse_statistics.py {input.single_cost} {input.grouping} {input.measurements} {input.additional_data} {params.task_label} {params.model_label} > {output}"

rule merge_groupings:
    input: lambda w: expand(f"{output_dir}/gene_families_disambiguated/{w.model}/wr_{w.window_radius}/{{group_id}}/additional_data.json", group_id=all_gene_families)
    params:
        target_dir=lambda w: f"{output_dir}/results/{w.model}/wr_{w.window_radius}/"
    output:
        phony=f"{output_dir}/results/computed_{{model}}_wr_{{window_radius}}"
        #files=expand(f"{output_dir}/results/{{model}}/{{label}}.bed", label=genomes.keys())
    shell: "python ../src/merge_groupings.py {input} -o {params.target_dir} && touch {output.phony}"


rule mapping_krister_format:
    input: lambda w: expand(f"{output_dir}/gene_families_disambiguated/{w.model}/wr_{w.window_radius}/{{group_id}}/additional_data.json", group_id=all_gene_families)
    params:
        input_dir=lambda w: f"{output_dir}/gene_families_disambiguated/{w.model}/wr_{w.window_radius}/",
    log: f"{output_dir}/result_mapping_krister_format/{{model}}/wr_{{window_radius}}/log"
    output: f"{output_dir}/result_mapping_krister_format/{{model}}/wr_{{window_radius}}/mapping.json"
    shell: """python ../src/mapping_krister_format.py {params.input_dir} {output} > {log} 2>&1"""


rule all_merge_groupings:
    default_target: True
    #input: expand(f"{output_dir}/results/{{model}}/{{label}}.bed", label=genomes.keys(), model=models.keys())
    input: expand(f"{output_dir}/results/computed_{{model}}_wr_{{window_radius}}", model=models.keys(), window_radius=window_radii),
           expand(f"{output_dir}/result_mapping_krister_format/{{model}}/wr_{{window_radius}}/mapping.json", model=models.keys(), window_radius=window_radii)
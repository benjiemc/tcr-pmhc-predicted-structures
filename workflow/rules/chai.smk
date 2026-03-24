import csv
import os

include: "utils.smk"

rule create_chai_sequences:
    input: "data/interim/sequences/{name}.json",
    output: "data/interim/sequences_chai/{name}.fasta",
    shell:
        """
        python -m tcr_pmhc_predicted_structures.apps.convert_to_input \
            --log-level info \
            --format chai \
            -o $(dirname {output[0]}) \
            {input}
        """

rule run_chai_1:
    input:
        sequences="data/interim/sequences_chai/{name}.fasta",
        msas=lambda wildcards: get_sequence_ids(f"data/interim/sequences/{wildcards.name}.json", f'data/interim/msas_{wildcards.msa_method}_{wildcards.data}/{{name}}.aligned.pqt'),
        templates_summaries=lambda wildcards: get_sequence_ids(f"data/interim/sequences/{wildcards.name}.json", f'data/interim/templates_{wildcards.msa_method}_{wildcards.data}_converted_indices/{{name}}.m8'),
        templates=lambda wildcards: get_sequence_ids(f"data/interim/sequences/{wildcards.name}.json", f'data/interim/templates_{wildcards.msa_method}_{wildcards.data}_converted_indices/{{name}}')
    output:
        structures=protected(expand("data/interim/structures_chai_{{msa_method}}_{{data}}/{{name}}_seed-{{seed}}_sample-{sample}_model.cif", sample=range(5))),
        summary_confidences=protected(expand("data/interim/structures_chai_{{msa_method}}_{{data}}/{{name}}_seed-{{seed}}_sample-{sample}_summary_confidences.json", sample=range(5))),
        confidences=protected(expand("data/interim/structures_chai_{{msa_method}}_{{data}}/{{name}}_seed-{{seed}}_sample-{sample}_confidences.json", sample=range(5)))
    resources:
        runtime="10m",
        mem="40GB",
        gpu=1
    script: "../scripts/predict_structure_with_chai.py"

rule create_structures_chai_af3_tcrmodel2_sample_size_50:
    input: get_predicted_structure_names("data/interim/sequences_sample_50.csv", "data/interim/structures_chai_af3_tcrmodel2/{name}_seed-{seed}_sample-{sample}_model.cif", config['seed'])

rule create_structures_chai_af3_tcrmodel2_sample_size_1500:
    input: get_predicted_structure_names("data/interim/sequences_sample_1500.csv", "data/interim/structures_chai_af3_tcrmodel2/{name}_seed-{seed}_sample-{sample}_model.cif", config['seed'])

rule create_structures_chai_af3_tcrmodel2_benchmark:
    input: get_predicted_structure_names("data/interim/benchmark_structures.csv", "data/interim/structures_chai_af3_tcrmodel2/{name}_seed-{seed}_sample-{sample}_model.cif", config['seed'], column_name='id')

rule process_chai_structures_with_tcrmodel2_data_sample_size_1500:
    input: get_predicted_structure_names("data/interim/sequences_sample_1500.csv", "data/processed/structures_chai_af3_tcrmodel2/{name}_{seed}_{sample}.pdb", config['seed'])

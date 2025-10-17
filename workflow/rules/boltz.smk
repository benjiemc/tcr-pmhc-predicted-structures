include: "utils.smk"


rule create_boltz_yaml:
    input:
        sequences="data/interim/sequences/{name}.json",
        msas=lambda wildcards: get_sequence_ids(f"data/interim/sequences/{wildcards.name}.json", f'data/interim/msas_{wildcards.msa_method}_{wildcards.data}/{{name}}.csv'),
        templates=lambda wildcards: get_sequence_ids(f"data/interim/sequences/{wildcards.name}.json", f'data/interim/templates_{wildcards.msa_method}_{wildcards.data}/{{name}}')
    output: "data/interim/sequences_boltz_{msa_method}_{data}/{name}.yaml"
    resources:
        runtime="5m",
        mem="1GB"
    shell:
        """
        out_dir_name=$(dirname {output[0]})

        python -m tcr_pmhc_predicted_structures.apps.convert_to_input \
            --log-level {config[log_level]} \
            --format boltz \
            --add-paths \
            --relative-paths \
            --unpaired-msa-paths {input.msas} \
            --template-paths {input.templates} \
            -o $out_dir_name \
            {input.sequences}
        """

rule run_boltz_predictions:
    input:
        sequence="data/interim/sequences_boltz_{msa_method}_{data}/{name}.yaml",
        msas=lambda wildcards: get_sequence_ids(f"data/interim/sequences/{wildcards.name}.json", f'data/interim/msas_{wildcards.msa_method}_{wildcards.data}/{{name}}.csv'),
        templates=lambda wildcards: get_sequence_ids(f"data/interim/sequences/{wildcards.name}.json", f'data/interim/templates_{wildcards.msa_method}_{wildcards.data}/{{name}}')
    output:
        structures=protected(expand("data/interim/structures_boltz_{{msa_method}}_{{data}}/{{name}}_seed-{{seed}}_sample-{sample}_model.cif", sample=range(5))),
        summary_confidences=protected(expand("data/interim/structures_boltz_{{msa_method}}_{{data}}/{{name}}_seed-{{seed}}_sample-{sample}_summary_confidences.json", sample=range(5))),
        confidences=protected(expand("data/interim/structures_boltz_{{msa_method}}_{{data}}/{{name}}_seed-{{seed}}_sample-{sample}_confidences.json", sample=range(5)))
    resources:
        runtime="10m",
        mem="80GB",
        gpu=1
    script: "../scripts/predict_structure_with_boltz.py"

rule create_structures_boltz_af3_tcrmodel2_sample_size_50:
    input: get_predicted_structure_names("data/interim/sequences_sample_50.csv", "data/interim/structures_boltz_af3_tcrmodel2/{name}_seed-{seed}_sample-{sample}_model.cif", config['seed'])

rule create_structures_boltz_af3_tcrmodel2_sample_size_1500:
    input: get_predicted_structure_names("data/interim/sequences_sample_1500.csv", "data/interim/structures_boltz_af3_tcrmodel2/{name}_seed-{seed}_sample-{sample}_model.cif", config['seed'])

rule create_structures_boltz_af3_tcrmodel2_benchmark:
    input: get_predicted_structure_names("data/interim/benchmark_structures.csv", "data/interim/structures_boltz_af3_tcrmodel2/{name}_seed-{seed}_sample-{sample}_model.cif", config['seed'], column_name='id')

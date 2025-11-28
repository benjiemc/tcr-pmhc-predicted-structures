import os

include: "utils.smk"

envvars:
    "TCRMODEL2_ALPHAFOLD_DATA_DIR"

rule run_alphafold2:
    input:
        sequences="data/interim/sequences_chai/{name}.fasta",
        msas=lambda wildcards: get_sequence_ids(f"data/interim/sequences/{wildcards.name}.json", f'data/interim/msas_{wildcards.msa_method}_{wildcards.data}/{{name}}_unpaired.a3m'),
        templates=lambda wildcards: get_sequence_ids(f"data/interim/sequences/{wildcards.name}.json", f'data/interim/templates_{wildcards.msa_method}_{wildcards.data}_converted_indices/{{name}}')
    output:
        structures_unrelaxed=protected(expand("data/interim/structures_af2_{{msa_method}}_{{data}}/{{name}}_seed-{{seed}}_sample-{sample}_unrelaxed.pdb", sample=range(5))),
        summary_confidences=protected(expand("data/interim/structures_af2_{{msa_method}}_{{data}}/{{name}}_seed-{{seed}}_sample-{sample}_summary_confidences.json", sample=range(5))),
        confidences=protected(expand("data/interim/structures_af2_{{msa_method}}_{{data}}/{{name}}_seed-{{seed}}_sample-{sample}_confidences.json", sample=range(5)))
    params:
        tcrmodel2_alphafold_data_dir=os.environ['TCRMODEL2_ALPHAFOLD_DATA_DIR']
    resources:
        runtime="30m",
        mem="40GB",
        gpu=1
    script: "../scripts/predict_structure_with_alphafold_from_tcrmodel2.py"

rule relax_af2_structure:
    input: "data/interim/structures_af2_{msa_method}_{data}/{name}_seed-{seed}_sample-{sample}_unrelaxed.pdb"
    output: "data/interim/structures_af2_{msa_method}_{data}/{name}_seed-{seed}_sample-{sample}_relaxed.pdb"
    resources:
        runtime="20m",
        mem="10GB",
        gpu=1
    script: "../scripts/relax_structure.py"

rule create_structures_af2_af3_tcrmodel2_sample_size_50:
    input: get_predicted_structure_names("data/interim/sequences_sample_50.csv", "data/interim/structures_af2_af3_tcrmodel2/{name}_seed-{seed}_sample-{sample}_unrelaxed.pdb", config['seed'])

rule create_structures_af2_af3_tcrmodel2_sample_size_1500:
    input: get_predicted_structure_names("data/interim/sequences_sample_1500.csv", "data/interim/structures_af2_af3_tcrmodel2/{name}_seed-{seed}_sample-{sample}_unrelaxed.pdb", config['seed'])

rule create_structures_af2_af3_tcrmodel2_benchmark:
    input: get_predicted_structure_names("data/interim/benchmark_structures.csv", "data/interim/structures_af2_af3_tcrmodel2/{name}_seed-{seed}_sample-{sample}_unrelaxed.pdb", config['seed'], column_name='id')

rule process_af2_structure:
    input: "data/interim/structures_{model}_{msa_method}_{data_source}/{name}_seed-{seed}_sample-{sample}_relaxed.pdb"
    output: "data/processed/structures_{model}_{msa_method}_{data_source}/{name}_{seed}_{sample}.pdb"
    resources:
        runtime="10m",
        mem="1GB"
    run:
        import warnings

        import stcrpy

        with warnings.catch_warnings():
            warnings.simplefilter('ignore', UserWarning)
            tcr = stcrpy.load_TCR(input[0])

        tcr.standardise_chain_names()
        tcr.crop()

        output_dir = os.path.dirname(output[0])
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        tcr.save(output[0])

rule calculate_af2_dockq_score:
    input:
        prediction="data/interim/structures_af2_{msa_method}_{data_source}/{name}_seed-{seed}_sample-{sample}_relaxed.pdb",
        structure="data/raw/structures_no_crop/{name}.pdb"
    output: temp("data/interim/af2_{msa_method}_{data_source}_{name}_{seed}_{sample}_dockq.json")
    resources:
        runtime="5m",
        mem="1GB"
    script: "../scripts/calculate_dockq_score.py"

rule process_af2_structures_with_tcrmodel2_data_sample_size_1500:
    input: get_predicted_structure_names("data/interim/sequences_sample_1500.csv", "data/processed/structures_af2_af3_tcrmodel2/{name}_{seed}_{sample}.pdb", config['seed'])

import os

include: "utils.smk"

envvars:
    "TCRMODEL2_ALPHAFOLD_DATA_DIR"

wildcard_constraints:
    sample_size="\d+",
    seed="\d+",
    sample="\d+",
    method="[A-Za-z0-9]+",
    data="[A-Za-z0-9]+",
    name="[A-Za-z0-9_]+"

rule run_alphafold2:
    input:
        sequences="data/interim/sequences_chai/{name}.fasta",
        msas=lambda wildcards: get_sequence_ids(f"data/interim/sequences/{wildcards.name}.json", f'data/interim/msas_{wildcards.msa_method}_{wildcards.data}/{{name}}_unpaired.a3m'),
        templates=lambda wildcards: get_sequence_ids(f"data/interim/sequences/{wildcards.name}.json", f'data/interim/templates_{wildcards.msa_method}_{wildcards.data}_converted_indices/{{name}}')
    output:
        structures_unrelaxed=protected(expand("data/interim/structures_af2_{{msa_method}}_{{data}}/{{name}}_seed-{{seed}}_sample-{sample}_unrelaxed.pdb", sample=range(5))),
        structures_relaxed=protected(expand("data/interim/structures_af2_{{msa_method}}_{{data}}/{{name}}_seed-{{seed}}_sample-{sample}_relaxed.pdb", sample=range(5))),
        summary_confidences=protected(expand("data/interim/structures_af2_{{msa_method}}_{{data}}/{{name}}_seed-{{seed}}_sample-{sample}_summary_confidences.json", sample=range(5))),
        confidences=protected(expand("data/interim/structures_af2_{{msa_method}}_{{data}}/{{name}}_seed-{{seed}}_sample-{sample}_confidences.json", sample=range(5)))
    params:
        tcrmodel2_alphafold_data_dir=os.environ['TCRMODEL2_ALPHAFOLD_DATA_DIR']
    resources:
        runtime="30m",
        mem="40GB",
        gpu=1
    script: "../scripts/predict_structure_with_alphafold_from_tcrmodel2.py"

rule create_structures_af2_af3_tcrmodel2_sample_size_50:
    input: get_predicted_structure_names("data/interim/sequences_sample_50.csv", "data/interim/structures_af2_af3_tcrmodel2/{name}_seed-{seed}_sample-{sample}_unrelaxed.pdb", config['seed'])

rule create_structures_af2_af3_tcrmodel2_sample_size_1500:
    input: get_predicted_structure_names("data/interim/sequences_sample_1500.csv", "data/interim/structures_af2_af3_tcrmodel2/{name}_seed-{seed}_sample-{sample}_unrelaxed.pdb", config['seed'])

rule create_structures_af2_af3_tcrmodel2_benchmark:
    input: get_predicted_structure_names("data/interim/benchmark_structures.csv", "data/interim/structures_af2_af3_tcrmodel2/{name}_seed-{seed}_sample-{sample}_unrelaxed.pdb", config['seed'], column_name='id')

rule extract_af2_confidence_scores:
    input:
        summary_confidences="data/interim/structures_af2_{msa_method}_{data}/{name}_seed-{seed}_sample-{sample}_summary_confidences.json",
        confidences="data/interim/structures_af2_{msa_method}_{data}/{name}_seed-{seed}_sample-{sample}_confidences.json",
        sequences="data/interim/sequences/{name}.json"
    output: temp("data/interim/af2_{msa_method}_{data}_{name}_{seed}_{sample}_confidence_scores.csv")
    resources:
        runtime="20m",
        mem="1GB"
    script: "../scripts/extract_af2_confidence_metrics.py"

rule aggregate_af2_confidence_scores:
    input: lambda wildcards: get_predicted_structure_names(f"data/interim/sequences_sample_{wildcards.sample_size}.csv", f"data/interim/af2_{wildcards.method}_{wildcards.data}_{{name}}_{{seed}}_{{sample}}_confidence_scores.csv", config['seed'])
    output: "data/interim/confidence_scores_af2_{method}_{data}_{sample_size}.csv"
    resources:
        runtime="10m",
        mem="1GB"
    run:
        import pandas as pd

        entries = []
        for path in input:
            confidence_scores = pd.read_csv(path)

            _, _, name, seed, sample = os.path.basename(path).replace('_confidence_scores.csv', '').rsplit('_', 4)
            confidence_scores['name'] = name
            confidence_scores['seed'] = seed
            confidence_scores['sample'] = sample

            entries.append(confidence_scores)

        df = pd.concat(entries).reset_index(drop=True)
        df.to_csv(output[0], index=False)
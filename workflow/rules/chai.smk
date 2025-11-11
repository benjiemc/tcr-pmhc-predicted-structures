import csv
import os

include: "utils.smk"

wildcard_constraints:
    sample_size="\d+",
    seed="\d+",
    sample="\d+",
    method="[A-Za-z0-9]+",
    data="[A-Za-z0-9]+",
    name="[A-Za-z0-9_]+"

def _get_chai_names(path: str) -> list[str]:
    if not os.path.exists(path):
        print(path, 'does not exist!')
        return []

    with open(path, 'r') as fh:
        reader = csv.DictReader(fh)

        return [
            f"data/interim/sequences-and-restraints-chai/{entry['index']}{extension}"
            for entry in reader for extension in ('.fasta', '_restraints.csv')
        ]

rule create_chai_sequence_and_restraints:
    input:
        sequences="data/interim/sequences/{name}.json",
        restraints="data/processed/tcr_pmhc_distance_restaints.csv"
    output:
        sequence="data/interim/sequences_chai/{name}.fasta",
        restraint="data/interim/sequences_chai/{name}_restraints.csv",
    shell:
        """
        python -m tcr_pmhc_predicted_structures.apps.convert_to_input \
            --log-level info \
            --format chai \
            --chai-restraints {input.restraints} \
            -o $(dirname {output[0]}) \
            {input.sequences}
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

rule extract_chai_confidence_scores:
    input:
        summary_confidences="data/interim/structures_chai_{msa_method}_{data}/{name}_seed-{seed}_sample-{sample}_summary_confidences.json",
        confidences="data/interim/structures_chai_{msa_method}_{data}/{name}_seed-{seed}_sample-{sample}_confidences.json",
        sequences="data/interim/sequences/{name}.json"
    output: temp("data/interim/chai_{msa_method}_{data}_{name}_{seed}_{sample}_confidence_scores.csv")
    resources:
        runtime="20m",
        mem="1GB"
    script: "../scripts/extract_chai_confidence_metrics.py"

rule aggregate_chai_confidence_scores:
    input: lambda wildcards: get_predicted_structure_names(f"data/interim/sequences_sample_{wildcards.sample_size}.csv", f"data/interim/chai_{wildcards.method}_{wildcards.data}_{{name}}_{{seed}}_{{sample}}_confidence_scores.csv", config['seed'])
    output: "data/interim/confidence_scores_chai_{method}_{data}_{sample_size}.csv"
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

rule run_chai_1_with_restraints:
    input:
        sequence="data/interim/sequences-and-restraints-chai/{name}.fasta",
        restraint="data/interim/sequences-and-restraints-chai/{name}_restraints.csv"
    output: directory("data/interim/structures-chai-restraints/{name}")
    run:
        import logging
        import shutil
        from pathlib import Path

        import numpy as np

        from chai_lab.chai1 import run_inference

        logging.basicConfig(level=logging.INFO)

        if not os.path.exists(os.path.dirname(output[0])):
            os.mkdir(os.path.dirname(output[0]))

        output_dir = Path(output[0])
        if output_dir.exists():
            logging.warning('Removing old output directory: %s', output[0])
            shutil.rmtree(output_dir)

        output_dir.mkdir()

        candidates = run_inference(
            fasta_file=Path(input.sequence),
            output_dir=output_dir,
            num_trunk_recycles=3,
            num_diffn_timesteps=200,
            seed=42,
            constraint_path=Path(input.restraint),
            use_esm_embeddings=True,
            use_msa_server=True,
            use_templates_server=True,
        )

        for i in range(len(candidates.cif_paths)):
            np.savez(
                output_dir.joinpath(f'confidences.model_idx_{i}.npz'),
                pae=candidates.pae[i].numpy(),
                pde=candidates.pde[i].numpy(),
                plddt=candidates.plddt[i].numpy(),
            )

def _get_chai_names(path: str) -> list[str]:
    if not os.path.exists(path):
        print(path, 'does not exist!')
        return []

    with open(path, 'r') as fh:
        return [f"data/interim/structures-chai/{entry['index']}" for entry in csv.DictReader(fh)]

rule create_chai_structures:
    input: _get_chai_names("data/interim/sequences_sample_50.csv")

def _get_chai_restraints_names(path: str) -> list[str]:
    if not os.path.exists(path):
        print(path, 'does not exist!')
        return []

    with open(path, 'r') as fh:
        return [f"data/interim/structures-chai-restraints/{entry['index']}" for entry in csv.DictReader(fh)]

rule create_chai_structures_with_restraints:
    input: _get_chai_restraints_names("data/interim/sequences_sample_50.csv")

import csv
import os

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

rule create_chai_sequences_and_restraints:
    input:
        sequences="data/interim/sequences_sample.csv",
        restraints="data/processed/tcr_pmhc_distance_restaints.csv",
        mhc_sequences=expand("data/external/hla_sequences/{mhc}.json", mhc=[
            'h2_d',
            'h2_k',
            'h2_l',
            'hla_a',
            'hla_b',
            'hla_c',
            'hla_e',
            'hla_f',
            'hla_g',
        ])
    output: _get_chai_names('data/interim/sequences_sample.csv')
    shell:
        """
        python -m tcr_pmhc_predicted_structures.apps.convert_to_input \
            --log-level info \
            --format chai \
            --mhc-sequences {input.mhc_sequences} \
            --chai-restraints {input.restraints} \
            -o $(dirname {output[0]}) \
            {input.sequences}
        """

rule run_chai_1:
    input: "data/interim/sequences-and-restraints-chai/{name}.fasta"
    output: directory("data/interim/structures-chai/{name}")
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
            fasta_file=Path(input[0]),
            output_dir=output_dir,
            num_trunk_recycles=3,
            num_diffn_timesteps=200,
            seed=42,
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
    input: _get_chai_names("data/interim/sequences_sample.csv")

def _get_chai_restraints_names(path: str) -> list[str]:
    if not os.path.exists(path):
        print(path, 'does not exist!')
        return []

    with open(path, 'r') as fh:
        return [f"data/interim/structures-chai-restraints/{entry['index']}" for entry in csv.DictReader(fh)]

rule create_chai_structures_with_restraints:
    input: _get_chai_restraints_names("data/interim/sequences_sample.csv")

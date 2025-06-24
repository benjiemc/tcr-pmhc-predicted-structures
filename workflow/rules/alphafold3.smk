import csv
import os

envvars:
    "ALPHAFOLD_DATA_DIR",
    "ALPHAFOLD3_WEIGHTS_DIR"

def _get_json_names(path: str) -> list[str]:
    if not os.path.exists(path):
        print(path, 'does not exist!')
        return []

    with open(path, 'r') as fh:
        reader = csv.DictReader(fh)

        return [f"data/interim/sequences-af3/{entry['index']}.json" for entry in reader]

rule create_af3_jsons:
    input:
        sequences="data/interim/sequences_sample.csv",
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
    output: _get_json_names("data/interim/sequences_sample.csv")
    shell:
        """
        out_dir_name=$(dirname {output[0]})

        python -m tcr_pmhc_predicted_structures.apps.convert_to_input \
            --log-level info \
            --format alphafold3 \
            --mhc-sequences {input.mhc_sequences} \
            -o $out_dir_name \
            {input.sequences}
        """

rule run_alphafold_msas:
    input: "data/interim/sequences-af3/{name}.json"
    output: directory("data/interim/msas-af3/{name}")
    params:
        alphafold_data_dir=os.environ['ALPHAFOLD_DATA_DIR'],
    shell:
        """
        input_file_name="$(basename {input})"
        input_dir_abs="$(pwd)/$(dirname {input})"
        output_dir_abs="$(pwd)/$(dirname {output})"

        mkdir -p $output_dir_abs

        docker run -it \
            --volume $output_dir_abs:/tmp/af_output \
            --entrypoint bash \
            alphafold3 \
            -c "chown -R $(id -u):$(id -g) /tmp/af_output"


        docker run -it \
            --user $(id -u):$(id -g) \
            --volume $input_dir_abs:/tmp/af_input \
            --volume $output_dir_abs:/tmp/af_output \
            --volume {params.alphafold_data_dir}:/public_databases \
            alphafold3 \
            python run_alphafold.py \
                --norun_inference \
                --json_path=/tmp/af_input/$input_file_name \
                --output_dir=/tmp/af_output
        """

def _get_af3_msa_names(path: str) -> list[str]:
    if not os.path.exists(path):
        print(path, 'does not exist!')
        return []

    with open(path, 'r') as fh:
        return [f"data/interim/msas-af3/{entry['index']}" for entry in csv.DictReader(fh)]

rule create_af3_msas:
    input: _get_af3_msa_names("data/interim/sequences_sample.csv")

rule run_alphafold_predictions:
    input: "data/interim/msas-af3/{name}/{name}_data.json"
    output: directory("data/interim/structures-af3/{name}")
    params:
        alphafold3_weights_dir=os.environ['ALPHAFOLD3_WEIGHTS_DIR']
    shell:
        """
        input_file_name="$(basename {input})"
        input_dir_abs="$(pwd)/$(dirname {input})"
        output_dir_abs="$(pwd)/$(dirname {output})"

        mkdir -p $output_dir_abs

        docker run -it \
            --volume $output_dir_abs:/tmp/af_output \
            --entrypoint bash \
            alphafold3 \
            -c "chown -R $(id -u):$(id -g) /tmp/af_output"


        docker run -it \
            --user $(id -u):$(id -g) \
            --volume $input_dir_abs:/tmp/af_input \
            --volume $output_dir_abs:/tmp/af_output \
            --volume {params.alphafold3_weights_dir}:/models \
            --gpus all \
            alphafold3 \
            python run_alphafold.py \
                --norun_data_pipeline \
                --json_path=/tmp/af_input/$input_file_name \
                --model_dir=/models \
                --output_dir=/tmp/af_output
        """

def _get_af3_names(path: str) -> list[str]:
    if not os.path.exists(path):
        print(path, 'does not exist!')
        return []

    with open(path, 'r') as fh:
        return [f"data/interim/structures-af3/{entry['index']}" for entry in csv.DictReader(fh)]

rule create_af3_structures:
    input: _get_af3_names("data/interim/sequences_sample.csv")

rule aggregate_af3_confidence_scores:
    input: _get_af3_names("data/interim/sequences_sample.csv")
    output: "data/interim/af3_confidence_scores.csv"
    shell:
        """
        python -m tcr_pmhc_predicted_structures.apps.extract_confidence_from_af3 \
            --log-level info \
            -o {output} \
            {input}
        """

rule compute_ipsae_on_af3_structure:
    input:
        confidences="data/interim/structures-af3/{name}/seed-{seed}_sample-{sample}/{name}_seed-{seed}_sample-{sample}_confidences.json",
        summary_confidences="data/interim/structures-af3/{name}/seed-{seed}_sample-{sample}/{name}_seed-{seed}_sample-{sample}_summary_confidences.json",
        structure="data/interim/structures-af3/{name}/seed-{seed}_sample-{sample}/{name}_seed-{seed}_sample-{sample}_model.cif"
    output:
        "data/interim/ipsaes-af3/{name}_seed-{seed}_sample-{sample}_model_{distance_cutoff}_{pae_cutoff}_byres.txt",
        "data/interim/ipsaes-af3/{name}_seed-{seed}_sample-{sample}_model_{distance_cutoff}_{pae_cutoff}.pml",
        "data/interim/ipsaes-af3/{name}_seed-{seed}_sample-{sample}_model_{distance_cutoff}_{pae_cutoff}.txt"
    shell:
        """
        output_dir=$(dirname {output[0]})
        mkdir -p $output_dir

        cp {input.confidences} $TMPDIR/
        cp {input.summary_confidences} $TMPDIR/
        cp {input.structure} $TMPDIR/

        python third_party/IPSAE/ipsae.py \
            $TMPDIR/$(basename {input.confidences}) \
            $TMPDIR/$(basename {input.structure}) \
            {wildcards.distance_cutoff} \
            {wildcards.pae_cutoff}

        mv $TMPDIR/$(basename {output[0]}) $output_dir/
        mv $TMPDIR/$(basename {output[1]}) $output_dir/
        mv $TMPDIR/$(basename {output[2]}) $output_dir/
        """

def _get_ipsaes_names(path: str, seeds: list[int], distance_cutoff: float, pae_cutoff: float) -> list[str]:
    if not os.path.exists(path):
        print(path, 'does not exist!')
        return []

    template = "data/interim/ipsaes-af3/{{name}}_seed-{{seed}}_sample-{{sample}}_model_{distance_cutoff:02}_{pae_cutoff:02}.txt"

    with open(path, 'r') as fh:
        indices = [entry['index'] for entry in csv.DictReader(fh)]

    return expand(template.format(distance_cutoff=distance_cutoff, pae_cutoff=pae_cutoff), name=indices, seed=seeds, sample=range(5))

rule add_af3_ipsaes_to_confidence_scores:
    input:
        confidence_scores="data/interim/af3_confidence_scores.csv",
        ipaes=_get_ipsaes_names("data/interim/sequences_sample.csv", seeds=[1], distance_cutoff=10, pae_cutoff=10)
    output: "data/processed/af3_confidence_scores.csv"
    run:
        import pandas as pd

        chain_map = {'A': 'mhc_chain1', 'B': 'mhc_chain2', 'C': 'antigen_chain', 'D': 'alpha_chain', 'E': 'beta_chain'}

        ipsaes = []
        for path in input.ipaes:
            name, seed, sample, _, _, _ = os.path.basename(path).replace('.txt', '').rsplit('_', 5)
            ipsae = pd.read_csv(path, sep='\s+')

            ipsae = ipsae[ipsae['Type'] != 'max'].copy()

            ipsae['chain_type1'] = ipsae['Chn1'].map(chain_map)
            ipsae['chain_type2'] = ipsae['Chn2'].map(chain_map)
            index = ipsae['chain_type1'] + '_' + ipsae['chain_type2'] + '_ipsae'
            ipsae = ipsae.set_index(index)

            ipsae = ipsae[['ipSAE']].T

            ipsae['structure_name'] = name
            ipsae['seed'] = int(seed.replace('seed-', ''))
            ipsae['sample'] = int(sample.replace('sample-', ''))

            ipsaes.append(ipsae)

        ipsaes = pd.concat(ipsaes).reset_index(drop=True)

        ipsaes['combined_ipsaes'] = ipsaes.filter(
            regex=(
                '((antigen_chain|mhc_chain_1)_(alpha_chain|beta_chain))'
                '|((alpha_chain|beta_chain)_(antigen_chain|mhc_chain_1))_ipsae'
            )
        ).apply('mean', axis=1)

        confidence_scores = pd.read_csv(input.confidence_scores)
        confidence_scores['structure_name'] = confidence_scores['structure_name'].astype(str)
        confidence_scores = confidence_scores.merge(ipsaes, how='left', on=['structure_name', 'seed', 'sample'])
        confidence_scores.to_csv(output[0], index=False)

rule summarise_af3_structural_features:
    input: "data/interim/structures-af3/{name}/seed-{seed}_sample-{sample}/{name}_seed-{seed}_sample-{sample}_model.cif"
    output: temp("data/interim/{name}_{seed}_{sample}_structural_summary.json")
    script: "../scripts/summarise_structural_features.py"

def _get_af3_structural_features_names(path: str, seeds: list) -> list[str]:
    if not os.path.exists(path):
        print(path, 'does not exist!')
        return []

    with open(path, 'r') as fh:
        indices = [entry['index'] for entry in csv.DictReader(fh)]

    return [f"data/interim/{name}_{seed}_{sample}_structural_summary.json" for name in indices for seed in seeds for sample in range(5)]

rule aggregate_af3_structural_features:
    input: _get_af3_structural_features_names("data/interim/sequences_sample.csv", [1])
    output: "data/processed/af3_structural_features.csv"
    run:
        import json
        import pandas as pd

        entries = []
        for path in input:
            with open(path, 'r') as fh:
                structural_features = json.load(fh)

            name, seed, sample = os.path.basename(path).replace('_structural_summary.json', '').rsplit('_', 2)
            structural_features['name'] = name
            structural_features['seed'] = seed
            structural_features['sample'] = sample

            entries.append(structural_features)

        df = pd.DataFrame(entries)
        df.to_csv(output[0], index=False)

rule process_af3_structure:
    input: "data/interim/structures-af3/{name}/seed-{seed}_sample-{sample}/{name}_seed-{seed}_sample-{sample}_model.cif"
    output: "data/processed/structures-af3/{name}_{seed}_{sample}.pdb"
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

def _get_af3_processed_structure_names(path: str, seeds: list) -> list[str]:
    if not os.path.exists(path):
        print(path, 'does not exist!')
        return []

    with open(path, 'r') as fh:
        indices = [entry['index'] for entry in csv.DictReader(fh)]

    return expand("data/processed/structures-af3/{name}_{seed}_{sample}.pdb", name=indices, seed=seeds, sample=list(range(5)))

rule process_af3_structures:
    input: _get_af3_processed_structure_names("data/interim/sequences_sample.csv", [1])

rule create_contact_map:
    input: "data/processed/structures-af3/{name}_{seed}_{sample}.pdb"
    output: temp("data/interim/{name}_{seed}_{sample}_contacts.csv")
    params:
        contact_distance=5.0
    shell:
        """
        python -m tcr_pmhc_predicted_structures.apps.create_contact_map \
            --log-level info \
            --cutoff-distance {params.contact_distance} \
            -o {output} \
            {input}
        """

def _get_af3_contact_maps(path: str, seeds: list) -> list[str]:
    if not os.path.exists(path):
        print(path, 'does not exist!')
        return []

    with open(path, 'r') as fh:
        indices = [entry['index'] for entry in csv.DictReader(fh)]

    return [f"data/interim/{name}_{seed}_{sample}_contacts.csv" for name in indices for seed in seeds for sample in range(5)]

rule aggregate_contact_maps:
    input: _get_af3_contact_maps("data/interim/sequences_sample.csv", [1])
    output: "data/processed/af3_contact_maps.csv"
    run:
        import pandas as pd

        contact_maps = []
        for path in input:
            contact_map = pd.read_csv(path)

            name, seed, sample = os.path.basename(path).replace('_contacts.csv', '').rsplit('_', 2)
            contact_map['name'] = name
            contact_map['seed'] = seed
            contact_map['sample'] = sample

            contact_maps.append(contact_map)

        contact_maps = pd.concat(contact_maps)
        contact_maps.to_csv(output[0], index=False)

rule profile_af3_interactions:
    input: "data/processed/structures-af3/{name}_{seed}_{sample}.pdb"
    output: temp("data/interim/{name}_{seed}_{sample}_interactions.csv")
    script: "../scripts/profile_interactions.py"

def _get_af3_interactions(path: str, seeds: list) -> list[str]:
    if not os.path.exists(path):
        print(path, 'does not exist!')
        return []

    with open(path, 'r') as fh:
        indices = [entry['index'] for entry in csv.DictReader(fh)]

    return [f"data/interim/{name}_{seed}_{sample}_interactions.csv" for name in indices for seed in seeds for sample in range(5)]

rule aggregate_af3_interactions:
    input: _get_af3_interactions("data/interim/sequences_sample.csv", [1])
    output: "data/processed/af3_interactions.csv"
    run:
        import pandas as pd

        interactions = []
        for path in input:
            interaction = pd.read_csv(path)

            name, seed, sample = os.path.basename(path).replace('_interactions.csv', '').rsplit('_', 2)
            interaction['name'] = name
            interaction['seed'] = seed
            interaction['sample'] = sample

            interactions.append(interaction)

        interactions = pd.concat(interactions)
        interactions.to_csv(output[0], index=False)

import csv
import os

seed = config.get('seed', [1])

envvars:
    "ALPHAFOLD3_WEIGHTS_DIR"

wildcard_constraints:
    sample_size="\d+"

def _get_json_names(template: str, path: str) -> list[str]:
    if not os.path.exists(path):
        print(path, 'does not exist!')
        return []

    with open(path, 'r') as fh:
        reader = csv.DictReader(fh)
        indices = [entry['index'] for entry in reader]

    return expand(template, name=indices)

def lookup_sequence_ids(sequence_path: str, sequence_key_path: str, template: str) -> list[str]:
    if not os.path.exists(sequence_path):
        print(sequence_path, 'does not exist!')
        return []

    with open(sequence_path, 'r') as fh:
        sequences = json.load(fh)

    if not os.path.exists(sequence_key_path):
        print(sequence_key_path, 'does not exist!')
        return []

    with open(sequence_key_path, 'r') as fh:
        sequence_key = {sequence: key for key, sequence in (line.split() for line in fh if line)}

    ids = [sequence_key[sequence] for sequence in sequences.values() if sequence in sequence_key]

    return expand(template, name=ids)

rule run_create_af3_jsons:
    input:
        sequences="data/interim/sequences/{name}.json",
        sequence_key="data/interim/sequence_key.txt",
        unpaired_msas=lambda wildcards: lookup_sequence_ids(f"data/interim/sequences/{wildcards.name}.json", "data/interim/sequence_key.txt", 'data/interim/msas_af3_af3/{name}_paired.a3m'),
        paired_msas=lambda wildcards: lookup_sequence_ids(f"data/interim/sequences/{wildcards.name}.json", "data/interim/sequence_key.txt", 'data/interim/msas_af3_af3/{name}_unpaired.a3m'),
        templates=lambda wildcards: lookup_sequence_ids(f"data/interim/sequences/{wildcards.name}.json", "data/interim/sequence_key.txt", 'data/interim/templates_af3_af3/{name}')
    output: "data/interim/sequences_af3_af3/{name}_{seed}.json"
    resources:
        runtime="5m",
        mem="1GB"
    shell:
        """
        out_dir_name=$(dirname {output[0]})

        python -m tcr_pmhc_predicted_structures.apps.convert_to_input \
            --log-level {config[log_level]} \
            --format alphafold3 \
             --alphafold3-model-seeds {wildcards.seed} \
            --sequence-key {input.sequence_key} \
            --relative-paths \
            --unpaired-msa-paths {input.unpaired_msas} \
            --paired-msa-paths {input.paired_msas} \
            --template-paths {input.templates} \
            -o $out_dir_name \
            {input.sequences}
        """

def _get_required_msas(wildcards, path, sequence_key_path, template) -> list[str]:
    if not os.path.exists(path):
        print(path, 'does not exist!')
        return []

    if not os.path.exists(sequence_key_path):
        print(sequence_key_path, 'does not exist!')
        return []

    with open(sequence_key_path, 'r') as fh:
        sequence_key = {sequence: key for key, sequence in (line.split() for line in fh if line)}

    with open(path, 'r') as fh:
        reader = csv.DictReader(fh)
        relevant_entries = [entry for entry in reader if entry['id'] == wildcards.name]

    sequences = [
        value for entry in relevant_entries for key, value in entry.items()
        if key in ('tcr_alpha_chain', 'tcr_beta_chain', 'mhc_chain_1', 'mhc_chain_2')
    ]

    keys = [sequence_key[sequence] for sequence in sequences]

    return expand(template, key=keys)


rule run_alphafold_predictions:
    input:
        sequence="data/interim/sequences_af3_af3/{name}_{seed}.json",
        unpaired_msas=lambda wildcards: _get_required_msas(wildcards, "data/interim/sequences_full.csv", "data/interim/sequence_key.txt", "data/interim/msas_af3_af3/{key}_unpaired.a3m"),
        paired_msas=lambda wildcards: _get_required_msas(wildcards, "data/interim/sequences_full.csv", "data/interim/sequence_key.txt", "data/interim/msas_af3_af3/{key}_paired.a3m"),
        templates=lambda wildcards: _get_required_msas(wildcards, "data/interim/sequences_full.csv", "data/interim/sequence_key.txt", "data/interim/templates_af3_af3/{key}")
    output:
        structures=expand("data/interim/structures_af3_af3/{{name}}_seed-{{seed}}_sample-{sample}_model.cif", sample=range(5)),
        summary_confidences=expand("data/interim/structures_af3_af3/{{name}}_seed-{{seed}}_sample-{sample}_summary_confidences.json", sample=range(5)),
        confidences=expand("data/interim/structures_af3_af3/{{name}}_seed-{{seed}}_sample-{sample}_confidences.json", sample=range(5))
    params:
        alphafold3_weights_dir=os.environ['ALPHAFOLD3_WEIGHTS_DIR']
    resources:
        runtime="10m",
        mem="64GB",
        gpu=1
    shell:
        """
        input_dir_abs=$TMPDIR/{wildcards.name}/input
        temp_output_dir_abs=$TMPDIR/{wildcards.name}/output
        if [ -d $input_dir_abs ]; then
            rm -r $input_dir_abs
        fi
        if [ -d $temp_output_dir_abs ]; then
            rm -r $temp_output_dir_abs
        fi
        mkdir -p $input_dir_abs $temp_output_dir_abs

        cp {input.sequence} $input_dir_abs
        cp {input.unpaired_msas} $input_dir_abs
        cp {input.paired_msas} $input_dir_abs
        cp -r {input.templates} $input_dir_abs

        mkdir -p $(dirname {output.structures[0]})

        input_file_name="$(basename {input.sequence})"

        apptainer exec \
            --nv \
            --bind $input_dir_abs:/root/input \
            --bind $temp_output_dir_abs:/root/output \
            --bind {params.alphafold3_weights_dir}:/root/models \
            containers/alphafold3.sif \
            python third_party/alphafold3/run_alphafold.py \
                --norun_data_pipeline \
                --json_path=/root/input/$input_file_name \
                --model_dir=/root/models \
                --output_dir=/root/output

        output_dir=$(dirname {output.structures[0]})
        mkdir -p $output_dir

        for sample_num in {{0..4}}; do
            prediction_dir="${{temp_output_dir_abs}}/{wildcards.name}/seed-{wildcards.seed}_sample-${{sample_num}}"

            mv "$prediction_dir/{wildcards.name}_seed-{wildcards.seed}_sample-${{sample_num}}_model.cif" $output_dir
            mv "$prediction_dir/{wildcards.name}_seed-{wildcards.seed}_sample-${{sample_num}}_summary_confidences.json" $output_dir
            mv "$prediction_dir/{wildcards.name}_seed-{wildcards.seed}_sample-${{sample_num}}_confidences.json" $output_dir
        done
        """

def _get_af3_names(path: str, template: str, seed: int) -> list[str]:
    if not os.path.exists(path):
        print(path, 'does not exist!')
        return []

    with open(path, 'r') as fh:
        indices = [entry['index'] for entry in csv.DictReader(fh)]

    return expand(template, name=indices, seed=seed, sample=range(5))


rule create_af3_structures_sample_size_50:
    input: _get_af3_names("data/interim/sequences_sample_50.csv", "data/interim/structures_af3_af3/{name}_seed-{seed}_sample-{sample}_model.cif", seed)

rule run_create_af3_jsons_tcrmodel2_data:
    input:
        sequences="data/interim/sequences/{name}.json",
        sequence_key="data/interim/sequence_key.txt",
        unpaired_msas=lambda wildcards: lookup_sequence_ids(f"data/interim/sequences/{wildcards.name}.json", "data/interim/sequence_key.txt", 'data/interim/msas_af3_tcrmodel2/{name}_paired.a3m'),
        paired_msas=lambda wildcards: lookup_sequence_ids(f"data/interim/sequences/{wildcards.name}.json", "data/interim/sequence_key.txt", 'data/interim/msas_af3_tcrmodel2/{name}_unpaired.a3m'),
        templates=lambda wildcards: lookup_sequence_ids(f"data/interim/sequences/{wildcards.name}.json", "data/interim/sequence_key.txt", 'data/interim/templates_af3_tcrmodel2/{name}')
    output: "data/interim/sequences_af3_tcrmodel2/{name}_{seed}.json"
    resources:
        runtime="5m",
        mem="1GB"
    shell:
        """
        out_dir_name=$(dirname {output[0]})

        python -m tcr_pmhc_predicted_structures.apps.convert_to_input \
            --log-level {config[log_level]} \
            --format alphafold3 \
             --alphafold3-model-seeds {wildcards.seed} \
            --sequence-key {input.sequence_key} \
            --relative-paths \
            --unpaired-msa-paths {input.unpaired_msas} \
            --paired-msa-paths {input.paired_msas} \
            --template-paths {input.templates} \
            -o $out_dir_name \
            {input.sequences}
        """

rule run_alphafold_predictions_with_tcrmodel2_data:
    input:
        sequence="data/interim/sequences_af3_tcrmodel2/{name}_{seed}.json",
        unpaired_msas=lambda wildcards: _get_required_msas(wildcards, "data/interim/sequences_full.csv", "data/interim/sequence_key.txt", "data/interim/msas_af3_tcrmodel2/{key}_unpaired.a3m"),
        paired_msas=lambda wildcards: _get_required_msas(wildcards, "data/interim/sequences_full.csv", "data/interim/sequence_key.txt", "data/interim/msas_af3_tcrmodel2/{key}_paired.a3m"),
        templates=lambda wildcards: _get_required_msas(wildcards, "data/interim/sequences_full.csv", "data/interim/sequence_key.txt", "data/interim/templates_af3_tcrmodel2/{key}")
    output:
        structures=expand("data/interim/structures_af3_tcrmodel2/{{name}}_seed-{{seed}}_sample-{sample}_model.cif", sample=range(5)),
        summary_confidences=expand("data/interim/structures_af3_tcrmodel2/{{name}}_seed-{{seed}}_sample-{sample}_summary_confidences.json", sample=range(5)),
        confidences=expand("data/interim/structures_af3_tcrmodel2/{{name}}_seed-{{seed}}_sample-{sample}_confidences.json", sample=range(5))
    params:
        alphafold3_weights_dir=os.environ['ALPHAFOLD3_WEIGHTS_DIR']
    resources:
        runtime="10m",
        mem="64GB",
        gpu=1
    shell:
        """
        input_dir_abs=$TMPDIR/{wildcards.name}/input
        temp_output_dir_abs=$TMPDIR/{wildcards.name}/output
        if [ -d $input_dir_abs ]; then
            rm -r $input_dir_abs
        fi
        if [ -d $temp_output_dir_abs ]; then
            rm -r $temp_output_dir_abs
        fi
        mkdir -p $input_dir_abs $temp_output_dir_abs

        cp {input.sequence} $input_dir_abs
        cp {input.unpaired_msas} $input_dir_abs
        cp {input.paired_msas} $input_dir_abs
        cp -r {input.templates} $input_dir_abs

        mkdir -p $(dirname {output.structures[0]})

        input_file_name="$(basename {input.sequence})"

        apptainer exec \
            --nv \
            --bind $input_dir_abs:/root/input \
            --bind $temp_output_dir_abs:/root/output \
            --bind {params.alphafold3_weights_dir}:/root/models \
            containers/alphafold3.sif \
            python third_party/alphafold3/run_alphafold.py \
                --norun_data_pipeline \
                --json_path=/root/input/$input_file_name \
                --model_dir=/root/models \
                --output_dir=/root/output

        output_dir=$(dirname {output.structures[0]})
        mkdir -p $output_dir

        for sample_num in {{0..4}}; do
            prediction_dir="${{temp_output_dir_abs}}/{wildcards.name}/seed-{wildcards.seed}_sample-${{sample_num}}"

            mv "$prediction_dir/{wildcards.name}_seed-{wildcards.seed}_sample-${{sample_num}}_model.cif" $output_dir
            mv "$prediction_dir/{wildcards.name}_seed-{wildcards.seed}_sample-${{sample_num}}_summary_confidences.json" $output_dir
            mv "$prediction_dir/{wildcards.name}_seed-{wildcards.seed}_sample-${{sample_num}}_confidences.json" $output_dir
        done
        """


rule create_af3_structures_with_tcrmodel2_data_sample_size_50:
    input: _get_af3_names("data/interim/sequences_sample_50.csv", "data/interim/structures_af3_tcrmodel2/{name}_seed-{seed}_sample-{sample}_model.cif", seed)

rule create_af3_structures_with_tcrmodel2_data_sample_size_1500:
    input: _get_af3_names("data/interim/sequences_sample_1500.csv", "data/interim/structures_af3_tcrmodel2/{name}_seed-{seed}_sample-{sample}_model.cif", seed)

rule aggregate_af3_confidence_scores:
    input:
        summary_confidences=lambda wildcards: _get_af3_names(f"data/interim/sequences_sample_{wildcards.sample_size}.csv", "data/interim/structures_{{method}}_{{data}}/{name}_seed-{seed}_sample-{sample}_summary_confidences.json", seed),
        confidences=lambda wildcards: _get_af3_names(f"data/interim/sequences_sample_{wildcards.sample_size}.csv", "data/interim/structures_{{method}}_{{data}}/{name}_seed-{seed}_sample-{sample}_confidences.json", seed),
        sequences=lambda wildcards: expand(_get_json_names("data/interim/sequences_{{{{method}}}}_{{{{data}}}}/{name}_{{seed}}.json", "data/interim/sequences_sample_{wildcarsd.sample_size}.csv"), seed=seed)
    output: "data/interim/confidence_scores_{method}_{data}_{sample_size}.csv"
    resources:
        runtime=lambda wildcards : f"{max(-(wildcards.sample_size // -5), 300)}m",
        mem=lambda wildcards : f"{max(-(wildcards.sample_size // -25), 100)}GB"
    shell:
        """
        python -m tcr_pmhc_predicted_structures.apps.extract_confidence_from_af3 \
            --log-level {config[log_level]} \
            -o {output} \
            --summary-confidence-paths {input.summary_confidences} \
            --confidence-paths {input.confidences} \
            --sequence-paths {input.sequences}
        """

rule compute_ipsae_on_af3_structure:
    input:
        confidences="data/interim/structures_{method}_{data}/{name}_seed-{seed}_sample-{sample}_confidences.json",
        summary_confidences="data/interim/structures_{method}_{data}/{name}_seed-{seed}_sample-{sample}_summary_confidences.json",
        structure="data/interim/structures_{method}_{data}/{name}_seed-{seed}_sample-{sample}_model.cif"
    output:
        "data/interim/ipsaes_{method}_{data}/{name}_seed-{seed}_sample-{sample}_model_{distance_cutoff}_{pae_cutoff}_byres.txt",
        "data/interim/ipsaes_{method}_{data}/{name}_seed-{seed}_sample-{sample}_model_{distance_cutoff}_{pae_cutoff}.pml",
        "data/interim/ipsaes_{method}_{data}/{name}_seed-{seed}_sample-{sample}_model_{distance_cutoff}_{pae_cutoff}.txt"
    resources:
        runtime="10m",
        mem="2GB"
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

    template = "data/interim/ipsaes_{{{{method}}}}_{{{{data}}}}/{{name}}_seed-{{seed}}_sample-{{sample}}_model_{distance_cutoff:02}_{pae_cutoff:02}.txt"

    with open(path, 'r') as fh:
        indices = [entry['index'] for entry in csv.DictReader(fh)]

    return expand(template.format(distance_cutoff=distance_cutoff, pae_cutoff=pae_cutoff), name=indices, seed=seeds, sample=range(5))

rule add_ipsaes_to_confidence_scores:
    input:
        confidence_scores="data/interim/confidence_scores_{method}_{data}_{sample_size}.csv",
        ipaes=lambda wildcards: _get_ipsaes_names(f"data/interim/sequences_sample_{wildcards.sample_size}.csv", seeds=[1], distance_cutoff=10, pae_cutoff=10)
    output: "data/processed/confidence_scores_{method}_{data}_{sample_size}.csv"
    resources:
        runtime="10m",
        mem="2GB"
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
    input: "data/interim/structures_af3_af3/{name}_seed-{seed}_sample-{sample}_model.cif"
    output: temp("data/interim/{name}_{seed}_{sample}_structural_summary.json")
    resources:
        runtime="10m",
        mem="2GB"
    script: "../scripts/summarise_structural_features.py"

def _get_af3_structural_features_names(path: str, seeds: list) -> list[str]:
    if not os.path.exists(path):
        print(path, 'does not exist!')
        return []

    with open(path, 'r') as fh:
        indices = [entry['index'] for entry in csv.DictReader(fh)]

    return [f"data/interim/{name}_{seed}_{sample}_structural_summary.json" for name in indices for seed in seeds for sample in range(5)]

rule aggregate_af3_structural_features:
    input: lambda wildcards: _get_af3_structural_features_names(f"data/interim/sequences_sample_{wildcards.sample_size}.csv", [1])
    output: "data/processed/af3_structural_features_{sample_size}.csv"
    resources:
        runtime="10m",
        mem="1GB"
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
    input: "data/interim/structures_{method}_{data}/{name}_seed-{seed}_sample-{sample}_model.cif"
    output: "data/processed/structures_{method}_{data}/{name}_{seed}_{sample}.pdb"
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

def _get_processed_structure_names(path: str, template: str, seeds: list) -> list[str]:
    if not os.path.exists(path):
        print(path, 'does not exist!')
        return []

    with open(path, 'r') as fh:
        indices = [entry['index'] for entry in csv.DictReader(fh)]

    return expand(template, name=indices, seed=seeds, sample=list(range(5)))


rule process_af3_structures_sample_size_50:
    input: _get_processed_structure_names("data/interim/sequences_sample_50.csv", "data/processed/structures_af3_af3/{name}_{seed}_{sample}.pdb", [1])

rule process_af3_structures_with_tcrmodel2_data_sample_size_50:
    input: _get_processed_structure_names("data/interim/sequences_sample_50.csv", "data/processed/structures_af3_tcrmodel2/{name}_{seed}_{sample}.pdb", [1])

rule process_af3_structures_with_tcrmodel2_data_sample_size_1500:
    input: _get_processed_structure_names("data/interim/sequences_sample_1500.csv", "data/processed/structures_af3_tcrmodel2/{name}_{seed}_{sample}.pdb", [1])

rule create_contact_map:
    input: "data/processed/structures_af3_af3/{name}_{seed}_{sample}.pdb"
    output: temp("data/interim/{name}_{seed}_{sample}_contacts.csv")
    params:
        contact_distance=5.0
    resources:
        runtime="10m",
        mem="1GB"
    shell:
        """
        python -m tcr_pmhc_predicted_structures.apps.create_contact_map \
            --log-level {config[log_level]} \
            --cutoff-distance {params.contact_distance} \
            -o {output} \
            {input}
        """

def _get_af3_contact_maps(path: str, template: str, seeds: list) -> list[str]:
    if not os.path.exists(path):
        print(path, 'does not exist!')
        return []

    with open(path, 'r') as fh:
        indices = [entry['index'] for entry in csv.DictReader(fh)]

    return [template.format(name=name, seed=seed, sample=sample) for name in indices for seed in seeds for sample in range(5)]


rule aggregate_contact_maps:
    input: lambda wildcards: _get_af3_contact_maps(f"data/interim/sequences_sample_{wildcards.sample_size}.csv", "data/interim/{name}_{seed}_{sample}_contacts.csv", [1])
    output: "data/processed/af3_contact_maps_{sample_size}.csv"
    resources:
        runtime="10m",
        mem="1GB"
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
    input: "data/processed/structures_af3_af3/{name}_{seed}_{sample}.pdb"
    output: temp("data/interim/{name}_{seed}_{sample}_interactions.csv")
    resources:
        runtime="10m",
        mem="2GB"
    script: "../scripts/profile_interactions.py"

def _get_af3_interactions(path: str, seeds: list) -> list[str]:
    if not os.path.exists(path):
        print(path, 'does not exist!')
        return []

    with open(path, 'r') as fh:
        indices = [entry['index'] for entry in csv.DictReader(fh)]

    return [f"data/interim/{name}_{seed}_{sample}_interactions.csv" for name in indices for seed in seeds for sample in range(5)]

rule aggregate_af3_interactions:
    input: lambda wildcards: _get_af3_interactions(f"data/interim/sequences_sample_{wildcards.sample_size}.csv", [1])
    output: "data/processed/af3_interactions_{sample_size}.csv"
    resources:
        runtime="10m",
        mem="1GB"
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

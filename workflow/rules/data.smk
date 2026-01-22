import csv
import json
import os

wildcard_constraints:
    unbalance="(_unbalanced)?",
    sample_size="\d+"

rule create_full_length_sequences:
    input:
        sequences="data/raw/sequences_pMHC_split.csv",
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
    output: "data/interim/sequences_full.csv"
    resources:
        runtime="20m",
        mem="1GB"
    shell:
        """
        python -m tcr_pmhc_predicted_structures.apps.create_full_length_sequences \
            --log-level {config[log_level]} \
            --mhc-sequences {input.mhc_sequences} \
            -o $TMPDIR/sequences_full.csv \
            {input.sequences}

        grep -v '*' $TMPDIR/sequences_full.csv > {output}
        """

rule get_bad_sequence_ids:
    input:
        raw_sequences="data/raw/sequences_pMHC_split.csv",
        full_sequences="data/interim/sequences_full.csv"
    output: "data/interim/bad_sequences.txt"
    resources:
        runtime="1m",
        mem="100MB"
    shell:
        """
        comm -3 \
            <(sed 1d {input.raw_sequences} | cut -f1 -d, | sort) \
            <(sed 1d {input.full_sequences} | cut -f1 -d, | sort) > {output}
        """

rule sample_sequences:
    input:
        sequences="data/raw/sequences_pMHC_split.csv",
        bad_sequences="data/interim/bad_sequences.txt"
    output: "data/interim/sequences_sample_{sample_size}{unbalance}.csv"
    params:
        seed=123
    resources:
        runtime="10m",
        mem="1GB"
    script: "../scripts/sample_sequences.py"

rule create_sequence_key:
    input: "data/interim/sequences_full.csv"
    output: "data/interim/sequence_key.txt"
    resources:
        runtime="1m",
        mem="1GB"
    shell:
        """
        cat \
            <(sed 1d {input} | cut -d, -f2 | sort | uniq) \
            <(sed 1d {input} | cut -d, -f3 | sort | uniq) \
            <(sed 1d {input} | cut -d, -f5 | sort | uniq) \
            <(sed 1d {input} | cut -d, -f6 | sort | uniq) | sort | uniq > $TEMPDIR/sequences.txt

        paste -d ' ' <(cat $TEMPDIR/sequences.txt | xargs -I {{}} sh -c 'echo -n "{{}}" | sha256sum | cut -c1-8') $TEMPDIR/sequences.txt > {output}
        """

def _get_complex_names(path: str, bad_sequences_path: str, template: str) -> list[str]:
    if not os.path.exists(bad_sequences_path):
        print(bad_sequences_path, 'does not exist!')
        return []

    with open(bad_sequences_path, 'r') as fh:
        bad_sequences = [line.strip() for line in fh if line]

    if not os.path.exists(path):
        print(path, 'does not exist!')
        return []

    with open(path, 'r') as fh:
        indices = [entry['index'] for entry in csv.DictReader(fh) if entry['index'] not in bad_sequences]

    return expand(template, name=indices)

rule convert_sequences_to_jsons:
    input: "data/interim/sequences_full.csv"
    output: _get_complex_names("data/raw/sequences_pMHC_split.csv", "data/interim/bad_sequences.txt", "data/interim/sequences/{name}.json")
    resources:
        runtime="10m",
        mem="1GB"
    run:
        output_dir = os.path.dirname(output[0])
        os.makedirs(output_dir, exist_ok=True)

        with open(input[0], 'r') as fh:
            reader = csv.DictReader(fh)
            for entry in reader:
                with open(os.path.join(output_dir, f"{entry['id']}.json"), 'w') as output_fh:
                    json.dump(entry, output_fh)

rule select_benchmark_structures:
    input:
        shi_parks_smith_benchmark="data/external/benchmark_set_from_shi_parks_smith.csv",
        structures="data/raw/structures_no_crop",
        structures_summary="data/raw/structures_summary_no_crop.csv",
        stcrdab="data/external/stcrdab"
    output: "data/interim/benchmark_structures.csv"
    resources:
        runtime="10m",
        mem="1GB"
    script: "../scripts/select_benchmark_structures.py"

def get_sequence_names(path: str, template: str) -> list[str]:
    if not os.path.exists(path):
        print(path, 'does not exist!')
        return []

    with open(path, 'r') as fh:
        ids = [entry['id'] for entry in csv.DictReader(fh)]

    return expand(template, name=ids)

rule extract_benchmark_sequences:
    input: "data/interim/benchmark_structures.csv"
    output: get_sequence_names("data/interim/benchmark_structures.csv", "data/interim/sequences/{name}.json")
    resources:
        runtime="1m",
        mem="10MB"
    script: "../scripts/extract_sequences_from_csv.py"

rule select_af_training_structures:
    input:
        structures="data/raw/structures_no_crop",
        structures_summary="data/raw/structures_summary_no_crop.csv",
        stcrdab="data/external/stcrdab"
    output: "data/interim/af_training_structures.csv"
    params:
        seed=123
    resources:
        runtime="10m",
        mem="1GB"
    script: "../scripts/select_af_training_structures.py"

rule extract_af_training_sequences:
    input: "data/interim/af_training_structures.csv"
    output: get_sequence_names("data/interim/af_training_structures.csv", "data/interim/sequences/{name}.json")
    resources:
        runtime="5m",
        mem="100MB"
    script: "../scripts/extract_sequences_from_csv.py"

rule combine_af_training_sequences_with_sample:
    input:
        raw_sequences="data/raw/sequences_pMHC_split.csv",
        sequences_sample="data/interim/sequences_sample_{sample_size}.csv",
        af_training_sequences="data/interim/af_training_structures.csv"
    output: "data/interim/sequences_sample_{sample_size}_af_training.csv"
    resources:
        runtime="5m",
        mem="1GB"
    script: "../scripts/combine_sequences_with_af_training_structures.py"

def create_paths(summary_file: str, path_name: str) -> list[str]:
    if not os.path.exists(summary_file):
        print(summary_file, 'does not exist!')
        return []

    with open(summary_file, 'r') as fh:
        names = [entry['path'].replace('.pdb', '') for entry in csv.DictReader(fh)]

    return expand(path_name, name=names)

rule summarise_structural_features:
    input: "data/raw/structures/{name}.pdb"
    output: temp("data/interim/{name}_structural_summary.json")
    resources:
        runtime="5m",
        mem="1GB"
    script: "../scripts/summarise_structural_features.py"

rule aggregate_structural_features:
    input: create_paths("data/raw/structures_summary.csv", "data/interim/{name}_structural_summary.json")
    output: "data/processed/tcr_pmhc_structural_features.csv"
    resources:
        runtime="5m",
        mem="1GB"
    run:
        import json
        import pandas as pd

        entries = []
        for path in input:
            with open(path, 'r') as fh:
                structural_features = json.load(fh)

            name = os.path.basename(path).replace('_structural_summary.json', '')
            structural_features['name'] = name

            entries.append(structural_features)

        df = pd.DataFrame(entries)
        df.to_csv(output[0], index=False)

rule profile_interactions:
    input: "data/raw/structures/{name}.pdb"
    output: temp("data/interim/{name}_interactions.csv")
    resources:
        runtime="5m",
        mem="1GB"
    script: "../scripts/profile_interactions.py"

rule aggregate_interactions:
    input: create_paths("data/raw/structures_summary.csv", "data/interim/{name}_interactions.csv")
    output: "data/processed/tcr_pmhc_interactions.csv"
    resources:
        runtime="5m",
        mem="1GB"
    run:
        import pandas as pd

        interactions = []
        for path in input:
            interaction = pd.read_csv(path)

            name = os.path.basename(path).replace('_interactions.csv', '')
            interaction['name'] = name

            interactions.append(interaction)

        interactions = pd.concat(interactions)
        interactions.to_csv(output[0], index=False)

rule create_distance_restraints:
    input:
        "data/raw/structures_summary.csv",
        "data/raw/structures"
    output: "data/processed/tcr_pmhc_distance_restaints.csv"
    resources:
        runtime="15m",
        mem="2GB"
    shell:
        """
        jupyter nbconvert --to notebook --execute --inplace notebooks/measure_constraint_distances.ipynb
        """

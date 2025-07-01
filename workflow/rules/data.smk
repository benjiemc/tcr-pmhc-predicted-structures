import csv
import os

rule sample_sequences:
    input: "data/raw/sequences_pMHC_split.csv"
    output: "data/interim/sequences_sample.csv"
    params:
        seed=123,
        sample_size=50
    run:
        import pandas as pd

        sequences = pd.read_csv(input[0], index_col='index')

        num_folds = len(sequences['fold'].unique())
        folds_sample_size = params.sample_size // num_folds
        group_sample_size = folds_sample_size // 2

        sequences_sample = sequences.groupby(['fold', 'label']).sample(group_sample_size, random_state=params.seed)
        sequences_sample.to_csv(output[0])

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
    script: "../scripts/summarise_structural_features.py"

rule aggregate_structural_features:
    input: create_paths("data/raw/structures_summary.csv", "data/interim/{name}_structural_summary.json")
    output: "data/processed/tcr_pmhc_structural_features.csv"
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
    script: "../scripts/profile_interactions.py"

rule aggregate_interactions:
    input: create_paths("data/raw/structures_summary.csv", "data/interim/{name}_interactions.csv")
    output: "data/processed/tcr_pmhc_interactions.csv"
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
    shell:
        """
        jupyter nbconvert --to notebook --execute --inplace notebooks/measure_constraint_distances.ipynb
        """

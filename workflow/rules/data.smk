import csv
import json
import os

envvars:
    "ALPHAFOLD_DATA_DIR",
    "TCRMODEL2_ALPHAFOLD_DATA_DIR",
    "COLABFOLD_DATA_DIR"

wildcard_constraints:
    sample_size="\d+"

rule sample_sequences:
    input: "data/raw/sequences_pMHC_split.csv"
    output: "data/interim/sequences_sample_{sample_size}.csv"
    params:
        seed=123
    run:
        import pandas as pd

        sequences = pd.read_csv(input[0], index_col='index')

        num_folds = len(sequences['fold'].unique())
        folds_sample_size = int(wildcards.sample_size) // num_folds
        group_sample_size = folds_sample_size // 2

        sequences_sample = sequences.groupby(['fold', 'label']).sample(group_sample_size, random_state=params.seed)
        sequences_sample.to_csv(output[0])

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
    shell:
        """
        python -m tcr_pmhc_predicted_structures.apps.create_full_length_sequences \
            --log-level {config[log_level]} \
            --mhc-sequences {input.mhc_sequences} \
            -o {output} \
            {input.sequences}
        """

rule create_sequence_key:
    input: "data/interim/sequences_full.csv"
    output: "data/interim/sequence_key.txt"
    shell:
        """
        cat \
            <(sed 1d {input} | cut -d, -f2 | sort | uniq) \
            <(sed 1d {input} | cut -d, -f3 | sort | uniq) \
            <(sed 1d {input} | cut -d, -f5 | sort | uniq) \
            <(sed 1d {input} | cut -d, -f6 | sort | uniq) | sort | uniq > $TEMPDIR/sequences.txt

        paste -d ' ' <(cat $TEMPDIR/sequences.txt | xargs -I {{}} sh -c 'echo -n "{{}}" | sha256sum | cut -c1-8') $TEMPDIR/sequences.txt > {output}
        """

def _get_complex_names(path: str, template: str) -> list[str]:
    if not os.path.exists(path):
        print(path, 'does not exist!')
        return []

    with open(path, 'r') as fh:
        indices = [entry['index'] for entry in csv.DictReader(fh)]

    return expand(template, name=indices)

rule convert_sequences_to_jsons:
    input: "data/interim/sequences_full.csv"
    output: _get_complex_names(f"data/raw/sequences_pMHC_split.csv", "data/interim/sequences/{name}.json")
    run:
        output_dir = os.path.dirname(output[0])
        os.makedirs(output_dir, exist_ok=True)

        with open(input[0], 'r') as fh:
            reader = csv.DictReader(fh)
            for entry in reader:
                with open(os.path.join(output_dir, f"{entry['id']}.json"), 'w') as output_fh:
                    json.dump(entry, output_fh)

def _create_msa_names(path: str, template: list[str]) -> list[str]:
    if not os.path.exists(path):
        print(path, 'does not exist!')
        return []

    with open(path, 'r') as fh:
        ids = [line.split()[0] for line in fh.readlines() if line]

    return expand(template, name=ids)

rule create_af3_msa_inputs:
    input: "data/interim/sequence_key.txt"
    output: _create_msa_names(f"data/interim/sequence_key.txt", 'data/interim/msa-inputs/{name}.json')
    run:
        output_dir = os.path.dirname(output[0])
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        with open(input[0], 'r') as fh:
            for id_, sequence in (line.split() for line in fh.readlines() if line):
                with open(os.path.join(output_dir, id_ + '.json'), 'w') as output_fh:
                    json.dump(
                    {
                        'name': id_,
                        'modelSeeds': [1],
                        'sequences': [
                            {'protein': {'id': 'A', 'sequence': sequence}},
                        ],
                        'dialect': 'alphafold3',
                        'version': 1,
                    },
                    output_fh,
                )

rule run_create_af3_msas_and_find_templates:
    input: "data/interim/msa-inputs/{name}.json"
    output:
        unpaired_msa='data/interim/msas/{name}_unpaired.a3m',
        paired_msa='data/interim/msas/{name}_paired.a3m',
        templates=directory('data/interim/templates/{name}')
    params:
        alphafold_data_dir=os.environ['ALPHAFOLD_DATA_DIR'],
        uniref90_database_path=lambda wildcards: "${DB_DIR}/uniref90_2022_05.fa",
        mgnify_database_path=lambda wildcards: "${DB_DIR}/mgy_clusters_2022_05.fa",
        small_bfd_database_path=lambda wildcards: "${DB_DIR}/bfd-first_non_consensus_sequences.fasta",
        uniprot_cluster_annot_database_path=lambda wildcards: "${DB_DIR}/uniprot_all_2021_04.fa",
        seqres_database_path=lambda wildcards: "${DB_DIR}/pdb_seqres_2022_09_28.fasta",
        pdb_database_path=lambda wildcards: "${DB_DIR}/mmcif_files"
    benchmark:
        "data/benchmarks/msa-generation_alphafold3_{name}.tsv"
    threads: 32
    script: "../scripts/create_msas_and_find_templates.py"

def _get_msas_for_sample(sample_path: str, full_sequence_path: str, sequence_key_path: str, template: str) -> list[str]:
    flag = False
    for path in sample_path, full_sequence_path, sequence_key_path:
        if not os.path.exists(path):
            print(sample_path, 'does not exist!')
            flag = True

    if flag:
        return []

    with open(sample_path, 'r') as fh:
        reader = csv.DictReader(fh)
        ids = [entry['index'] for entry in reader]

    with open(full_sequence_path, 'r') as fh:
        reader = csv.DictReader(fh)
        sequences = [
            entry[chain_name]
            for entry in reader
            if entry['id'] in ids
            for chain_name in ('tcr_alpha_chain', 'tcr_beta_chain', 'mhc_chain_1', 'mhc_chain_2')
        ]

    with open(sequence_key_path, 'r') as fh:
        sequence_key = {sequence: key for key, sequence in (line.split() for line in fh if line)}

    names = [sequence_key[sequence] for sequence in sequences]

    return expand(template, name=names)

rule create_af3_msas_and_find_templates_sample_size_50:
    input: _get_msas_for_sample("data/interim/sequences_sample_50.csv", "data/interim/sequences_full.csv", "data/interim/sequence_key.txt", 'data/interim/msas/{name}_unpaired.a3m')

rule create_colabfold_msa_inputs:
    input: "data/interim/sequence_key.txt"
    output: _create_msa_names("data/interim/sequence_key.txt", "data/interim/colabfold-msa-inputs/{name}.fasta")
    run:
        output_dir = os.path.dirname(output[0])
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        with open(input[0], 'r') as fh:
            for id_, sequence in (line.split() for line in fh.readlines() if line):
                with open(os.path.join(output_dir, id_ + '.fasta'), 'w') as output_fh:
                    output_fh.write(f'>{id_}\n{sequence}\n')

rule run_create_colabfold_msas_and_find_templates:
    input: "data/interim/colabfold-msa-inputs/{name}.fasta"
    output:
        msa='data/interim/colabfold-msas/{name}.a3m',
        templates='data/interim/colabfold-templates/{name}.m8'
    threads: 32
    params:
        colabfold_data_dir=os.environ['COLABFOLD_DATA_DIR']
    benchmark:
        "data/benchmarks/msa-generation_mmseq_colabfold_{name}.tsv"
    shell:
        """
        working_dir=$TMPDIR/{wildcards.name}
        if [ -d $working_dir ]; then
            rm -r $working_dir
        fi
        mkdir $working_dir

        colabfold_search \
            --mmseqs $(which mmseqs) \
            --use-templates 1 \
            --threads {threads} \
            --db2 pdb100_230517 \
            {input} \
            {params.colabfold_data_dir} \
            $working_dir

        mv $working_dir/0.a3m {output.msa}
        mv $working_dir/pdb100_230517.m8 {output.templates}
        """

rule create_colabfold_msas_and_find_templates_sample_size_50:
    input: _get_msas_for_sample("data/interim/sequences_sample_50.csv", "data/interim/sequences_full.csv", "data/interim/sequence_key.txt", 'data/interim/colabfold-msas/{name}.a3m')

rule run_create_tcrmodel2_msas_and_find_templates:
    input: "data/interim/msa-inputs/{name}.json"
    output:
        unpaired_msa='data/interim/tcrmodel2-msas/{name}_unpaired.a3m',
        paired_msa='data/interim/tcrmodel2-msas/{name}_paired.a3m',
        templates=directory('data/interim/tcrmodel2-templates/{name}')
    params:
        alphafold_data_dir=os.environ['TCRMODEL2_ALPHAFOLD_DATA_DIR'],
        uniref90_database_path=lambda wildcards: "${DB_DIR}/uniref90.tcrmhc.fasta",
        mgnify_database_path=lambda wildcards: "${DB_DIR}/mgnify.fasta",
        small_bfd_database_path=lambda wildcards: "${DB_DIR}/small_bfd.tcrmhc.fasta",
        uniprot_cluster_annot_database_path=lambda wildcards: "${DB_DIR}/uniprot.tcrmhc.fasta",
        seqres_database_path=lambda wildcards: "${DB_DIR}/pdb_seqres.txt",
        pdb_database_path=lambda wildcards: "${DB_DIR}/pdb_mmcif/mmcif_files"
    threads: 32
    benchmark:
        "data/benchmarks/msa-generation_alphafold3_tcrmodel2_{name}.tsv"
    script: "../scripts/create_msas_and_find_templates.py"

rule create_tcrmodel2_msas_and_find_templates_sample_size_50:
    input: _get_msas_for_sample("data/interim/sequences_sample_50.csv", "data/interim/sequences_full.csv", "data/interim/sequence_key.txt", 'data/interim/tcrmodel2-msas/{name}_unpaired.a3m')

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

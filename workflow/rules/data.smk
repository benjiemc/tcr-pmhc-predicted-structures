import csv
import json
import os

envvars:
    "ALPHAFOLD_DATA_DIR",
    "TCRMODEL2_ALPHAFOLD_DATA_DIR",
    "COLABFOLD_DATA_DIR"

wildcard_constraints:
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
    output: "data/interim/sequences_sample_{sample_size}.csv"
    params:
        seed=123
    resources:
        runtime="10m",
        mem="1GB"
    run:
        import pandas as pd

        sequences = pd.read_csv(input.sequences, index_col='index')

        with open(input.bad_sequences, 'r') as fh:
            bad_sequences = [int(line.strip()) for line in fh if line]

        sequences = sequences.loc[~sequences.index.isin(bad_sequences)]

        num_folds = len(sequences['fold'].unique())
        folds_sample_size = int(wildcards.sample_size) // num_folds
        group_sample_size = folds_sample_size // 2

        sequences_sample = sequences.groupby(['fold', 'label']).sample(group_sample_size, random_state=params.seed)
        sequences_sample.to_csv(output[0])

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

def _get_benchmark_sequence_names(path: str, template: str) -> list[str]:
    if not os.path.exists(path):
        print(path, 'does not exist!')
        return []

    with open(path, 'r') as fh:
        ids = [entry['id'] for entry in csv.DictReader(fh)]

    return expand(template, name=ids)

rule extract_benchmark_sequences:
    input: "data/interim/benchmark_structures.csv"
    output: _get_benchmark_sequence_names("data/interim/benchmark_structures.csv", "data/interim/sequences/{name}.json")
    resources:
        runtime="1m",
        mem="10MB"
    run:
        import pandas as pd

        df = pd.read_csv(input[0])
        output_df = df[
            [
                column_name
                for column_name in df.columns
                if column_name.endswith('_sequence') or column_name in ('id', 'mhc_type')
            ]
        ]

        output_df.columns = [{
            'Achain_sequence': 'tcr_alpha_chain',
            'Bchain_sequence': 'tcr_beta_chain',
            'antigen_chain_sequence': 'peptide',
            'mhc_chain1_sequence': 'mhc_chain_1',
            'mhc_chain2_sequence': 'mhc_chain_2',
        }.get(column_name, column_name) for column_name in output_df.columns]

        output_dir = os.path.dirname(output[0])
        os.makedirs(output_dir, exist_ok=True)
        for entry in [{k: v for k, v in row.items() if pd.notna(v)} for row in output_df.to_dict(orient='records')]:
            with open(os.path.join(output_dir, f"{entry['id']}.json"), 'w') as fh:
                json.dump(entry, fh)

def _create_msa_names(path: str, template: list[str]) -> list[str]:
    if not os.path.exists(path):
        print(path, 'does not exist!')
        return []

    with open(path, 'r') as fh:
        ids = [line.split()[0] for line in fh.readlines() if line]

    return expand(template, name=ids)

rule create_af3_msa_inputs:
    input: "data/interim/sequence_key.txt"
    output: _create_msa_names(f"data/interim/sequence_key.txt", 'data/interim/msa-inputs_af3/{name}.json')
    resources:
        runtime="10m",
        mem="1GB"
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
    input: "data/interim/msa-inputs_af3/{name}.json"
    output:
        unpaired_msa=protected('data/interim/msas_af3_af3/{name}_unpaired.a3m'),
        paired_msa=protected('data/interim/msas_af3_af3/{name}_paired.a3m'),
        templates=protected(directory('data/interim/templates_af3_af3/{name}'))
    params:
        alphafold_data_dir=os.environ['ALPHAFOLD_DATA_DIR'],
        uniref90_database_path=lambda wildcards: "${DB_DIR}/uniref90_2022_05.fa",
        mgnify_database_path=lambda wildcards: "${DB_DIR}/mgy_clusters_2022_05.fa",
        small_bfd_database_path=lambda wildcards: "${DB_DIR}/bfd-first_non_consensus_sequences.fasta",
        uniprot_cluster_annot_database_path=lambda wildcards: "${DB_DIR}/uniprot_all_2021_04.fa",
        seqres_database_path=lambda wildcards: "${DB_DIR}/pdb_seqres_2022_09_28.fasta",
        pdb_database_path=lambda wildcards: "${DB_DIR}/mmcif_files"
    benchmark:
        "data/benchmarks/msa-generation_af3_af3_{name}.tsv"
    threads: 32
    resources:
        runtime="3h",
        mem="64GB"
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
    input:
        unpaired_msas=_get_msas_for_sample("data/interim/sequences_sample_50.csv", "data/interim/sequences_full.csv", "data/interim/sequence_key.txt", 'data/interim/msas_af3_af3/{name}_unpaired.a3m'),
        paired_msas=_get_msas_for_sample("data/interim/sequences_sample_50.csv", "data/interim/sequences_full.csv", "data/interim/sequence_key.txt", 'data/interim/msas_af3_af3/{name}_paired.a3m'),
        templates=_get_msas_for_sample("data/interim/sequences_sample_50.csv", "data/interim/sequences_full.csv", "data/interim/sequence_key.txt", 'data/interim/templates_af3_af3/{name}')

rule create_colabfold_msa_inputs:
    input: "data/interim/sequence_key.txt"
    output: _create_msa_names("data/interim/sequence_key.txt", "data/interim/msa-inputs_mmseqs2/{name}.fasta")
    resources:
        runtime="10m",
        mem="1GB"
    run:
        output_dir = os.path.dirname(output[0])
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        with open(input[0], 'r') as fh:
            for id_, sequence in (line.split() for line in fh.readlines() if line):
                with open(os.path.join(output_dir, id_ + '.fasta'), 'w') as output_fh:
                    output_fh.write(f'>{id_}\n{sequence}\n')

rule run_create_colabfold_msas_and_find_templates:
    input: "data/interim/msa-inputs_mmseqs2/{name}.fasta"
    output:
        msa=protected('data/interim/msas_mmseqs2_colabfold/{name}_unpaired.a3m'),
        templates=protected('data/interim/templates_mmseqs2_colabfold/{name}.m8')
    threads: 32
    resources:
        runtime="5h",
        mem="64GB"
    params:
        colabfold_data_dir=os.environ['COLABFOLD_DATA_DIR']
    benchmark:
        "data/benchmarks/msa-generation_mmseqs2_colabfold_{name}.tsv"
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

rule extract_colabfold_templates:
    input: "data/interim/templates_mmseqs2_colabfold/{name}.m8"
    output: directory("data/interim/templates_mmseqs2_colabfold/{name}")
    resources:
        runtime="3m",
        mem="1GB"
    params:
        colabfold_data_dir=os.environ['COLABFOLD_DATA_DIR'],
        number_of_templates=20
    script: '../scripts/extract_colabfold_templates.py'

rule spoof_colabfold_paired_msas:
    input: "data/interim/msa-inputs_mmseqs2/{name}.fasta"
    output: "data/interim/msas_mmseqs2_colabfold/{name}_paired.a3m"
    resources:
        runtime="1m",
        mem="1MB",
    run:
        with open(input[0], 'r') as fh:
            sequence = next(line.strip() for line in fh if line and not line.startswith('>'))

        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        with open(output[0], 'w') as fh:
            fh.write(f">query\n{sequence}\n")

rule create_colabfold_msas_and_find_templates_sample_size_50:
    input:
        unpaired_msas=_get_msas_for_sample("data/interim/sequences_sample_50.csv", "data/interim/sequences_full.csv", "data/interim/sequence_key.txt", 'data/interim/msas_mmseqs2_colabfold/{name}_unpaired.a3m'),
        paired_msas=_get_msas_for_sample("data/interim/sequences_sample_50.csv", "data/interim/sequences_full.csv", "data/interim/sequence_key.txt", 'data/interim/msas_mmseqs2_colabfold/{name}_paired.a3m'),
        templates=_get_msas_for_sample("data/interim/sequences_sample_50.csv", "data/interim/sequences_full.csv", "data/interim/sequence_key.txt", 'data/interim/templates_mmseqs2_colabfold/{name}')

rule run_create_tcrmodel2_msas_and_find_templates:
    input: "data/interim/msa-inputs_af3/{name}.json"
    output:
        unpaired_msa=protected('data/interim/msas_af3_tcrmodel2/{name}_unpaired.a3m'),
        paired_msa=protected('data/interim/msas_af3_tcrmodel2/{name}_paired.a3m'),
        templates=protected(directory('data/interim/templates_af3_tcrmodel2/{name}'))
    params:
        alphafold_data_dir=os.environ['TCRMODEL2_ALPHAFOLD_DATA_DIR'],
        uniref90_database_path=lambda wildcards: "${DB_DIR}/uniref90.tcrmhc.fasta",
        mgnify_database_path=lambda wildcards: "${DB_DIR}/mgnify.fasta",
        small_bfd_database_path=lambda wildcards: "${DB_DIR}/small_bfd.tcrmhc.fasta",
        uniprot_cluster_annot_database_path=lambda wildcards: "${DB_DIR}/uniprot.tcrmhc.fasta",
        seqres_database_path=lambda wildcards: "${DB_DIR}/pdb_seqres.txt",
        pdb_database_path=lambda wildcards: "${DB_DIR}/pdb_mmcif/mmcif_files"
    threads: 32
    resources:
        runtime="1h",
        mem="64GB"
    benchmark:
        "data/benchmarks/msa-generation_af3_tcrmodel2_{name}.tsv"
    script: "../scripts/create_msas_and_find_templates.py"

rule create_tcrmodel2_msas_and_find_templates_sample_size_50:
    input:
        unpaired_msas=_get_msas_for_sample("data/interim/sequences_sample_50.csv", "data/interim/sequences_full.csv", "data/interim/sequence_key.txt", 'data/interim/msas_af3_tcrmodel2/{name}_unpaired.a3m'),
        paired_msas=_get_msas_for_sample("data/interim/sequences_sample_50.csv", "data/interim/sequences_full.csv", "data/interim/sequence_key.txt", 'data/interim/msas_af3_tcrmodel2/{name}_paired.a3m'),
        templates=_get_msas_for_sample("data/interim/sequences_sample_50.csv", "data/interim/sequences_full.csv", "data/interim/sequence_key.txt", 'data/interim/templates_af3_tcrmodel2/{name}')

rule create_tcrmodel2_msas_and_find_templates_sample_size_1500:
    input:
        unpaired_msas=_get_msas_for_sample("data/interim/sequences_sample_1500.csv", "data/interim/sequences_full.csv", "data/interim/sequence_key.txt", 'data/interim/msas_af3_tcrmodel2/{name}_unpaired.a3m'),
        paired_msas=_get_msas_for_sample("data/interim/sequences_sample_1500.csv", "data/interim/sequences_full.csv", "data/interim/sequence_key.txt", 'data/interim/msas_af3_tcrmodel2/{name}_paired.a3m'),
        templates=_get_msas_for_sample("data/interim/sequences_sample_1500.csv", "data/interim/sequences_full.csv", "data/interim/sequence_key.txt", 'data/interim/templates_af3_tcrmodel2/{name}')

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

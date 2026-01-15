import csv
import hashlib
import json
import os

envvars:
    "ALPHAFOLD_DATA_DIR",
    "TCRMODEL2_ALPHAFOLD_DATA_DIR",
    "COLABFOLD_DATA_DIR"

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

def get_sequence_ids(
    path: str, template: str, sequence_key_path: str | None = None, additional_structures_path: str | None = None
) -> list[str]:
    past_ids = set()

    if sequence_key_path:
        if not os.path.exists(sequence_key_path):
            print(sequence_key_path, 'does not exist!')
            return []

        with open(sequence_key_path, 'r') as fh:
            past_ids = past_ids | {line.strip().split()[0] for line in fh if line}

    if additional_structures_path:
        if not os.path.exists(additional_structures_path):
            print(additional_structures_path, 'does not exist!')
            return []

        with open(additional_structures_path, 'r') as fh:
            for entry in csv.DictReader(fh):
                for sequence_name in (
                    'Achain_sequence', 'Bchain_sequence', 'mhc_chain1_sequence', 'mhc_chain2_sequence'
                ):
                    if sequence_name not in entry:
                        continue

                    sequence = entry[sequence_name]
                    if not sequence:
                        continue

                    past_ids.add(hashlib.sha256(sequence.encode()).hexdigest()[:8])

    if not os.path.exists(path):
        print(path, 'does not exist!')
        return []

    with open(path, 'r') as fh:
        ids = {
            sequence_id
            for entry in csv.DictReader(fh)
            for sequence_name in ('Achain_sequence', 'Bchain_sequence', 'mhc_chain1_sequence', 'mhc_chain2_sequence')
            if (sequence := entry[sequence_name])
            and (sequence_id := hashlib.sha256(sequence.encode()).hexdigest()[:8]) not in past_ids
        }

    return expand(template, name=ids)

rule create_af3_msa_inputs_from_benchmark:
    input:
        structures="data/interim/benchmark_structures.csv",
        sequence_key="data/interim/sequence_key.txt"
    output: get_sequence_ids('data/interim/benchmark_structures.csv', 'data/interim/msa-inputs_af3/{name}.json', sequence_key_path='data/interim/sequence_key.txt')
    resources:
        runtime="10m",
        mem="10MB"
    script: "../scripts/create_af3_msa_inputs_from_structures.py"

rule create_af3_msa_inputs_from_af_training_structures:
    input:
        structures="data/interim/af_training_structures.csv",
        additional_structures="data/interim/benchmark_structures.csv",
        sequence_key="data/interim/sequence_key.txt"
    output: get_sequence_ids('data/interim/af_training_structures.csv', 'data/interim/msa-inputs_af3/{name}.json', sequence_key_path='data/interim/sequence_key.txt', additional_structures_path='data/interim/benchmark_structures.csv')
    resources:
        runtime="10m",
        mem="10MB"
    script: "../scripts/create_af3_msa_inputs_from_structures.py"

def get_dash_sequence_ids(path: str, sequence_key_path: str, template: str) -> list[str]:
    if not os.path.exists(sequence_key_path):
        print(sequence_key_path, 'does not exist!')
        return []

    with open(sequence_key_path, 'r') as fh:
        sequence_key = dict(line.strip().split() for line in fh if line)

    if not os.path.exists(path):
        print(path, 'does not exist!')
        return []

    with open(path, 'r') as fh:
        reader = csv.DictReader(fh)

        ids = {
            sequence_id
            for entry in reader
            for sequence_name in ('tcr_alpha_chain', 'tcr_beta_chain', 'mhc_chain_1', 'mhc_chain_2')
            if (sequence := entry[sequence_name])
            and (sequence_id := hashlib.sha256(sequence.encode()).hexdigest()[:8]) not in sequence_key
        }

    return expand(template, name=ids)

rule create_af3_msa_inputs_from_dash_dataset:
    input:
        sequences="data/interim/sequences_dash_with_negatives.csv",
        sequence_key="data/interim/sequence_key.txt"
    output: get_dash_sequence_ids("data/interim/sequences_dash_with_negatives.csv", "data/interim/sequence_key.txt", "data/interim/msa-inputs_af3/{name}.json")
    resources:
        runtime="10m",
        mem="500MB"
    run:
        output_dir = os.path.dirname(output[0])
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        with open(input.sequence_key, 'r') as fh:
            sequence_key = dict(line.strip().split() for line in fh if line)

        past_ids = set()
        with open(input.sequences, 'r') as fh:
            reader = csv.DictReader(fh)
            for entry in reader:
                for sequence_name in (
                    'tcr_alpha_chain', 'tcr_beta_chain', 'mhc_chain_1', 'mhc_chain_2'
                ):
                    if sequence_name not in entry:
                        continue

                    sequence = entry[sequence_name]
                    if not sequence:
                        continue

                    sequence_id = hashlib.sha256(sequence.encode()).hexdigest()[:8]
                    if sequence_id in sequence_key or sequence_id in past_ids:
                        continue

                    with open(os.path.join(output_dir, sequence_id + '.json'), 'w') as output_fh:
                        json.dump(
                            {
                                'name': sequence_id,
                                'modelSeeds': [1],
                                'sequences': [
                                    {'protein': {'id': 'A', 'sequence': sequence}},
                                ],
                                'dialect': 'alphafold3',
                                'version': 1,
                            },
                            output_fh,
                        )

                    past_ids.add(sequence_id)

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

def _get_msas_for_benchmark(path: str, template: str) -> list[str]:
    if not os.path.exists(path):
        print(path, 'does not exist!')
        return []

    with open(path, 'r') as fh:
        reader = csv.DictReader(fh)
        ids = {
            hashlib.sha256(sequence.encode()).hexdigest()[:8]
            for entry in reader
            for sequence_name in ('Achain_sequence', 'Bchain_sequence', 'mhc_chain1_sequence', 'mhc_chain2_sequence')
            if (sequence := entry[sequence_name])
        }

    return expand(template, name=ids)

rule create_af3_msas_and_find_templates_benchmark_structures:
    input:
        unpaired_msas=_get_msas_for_benchmark('data/interim/benchmark_structures.csv', 'data/interim/msas_af3_af3/{name}_unpaired.a3m'),
        paired_msas=_get_msas_for_benchmark('data/interim/benchmark_structures.csv', 'data/interim/msas_af3_af3/{name}_paired.a3m'),
        templates=_get_msas_for_benchmark('data/interim/benchmark_structures.csv', 'data/interim/templates_af3_af3/{name}')

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

rule create_colabfold_msa_inputs_from_benchmark:
    input:
        benchmark_structures="data/interim/benchmark_structures.csv",
        sequence_key="data/interim/sequence_key.txt"
    output: get_sequence_ids('data/interim/benchmark_structures.csv', 'data/interim/msa-inputs_mmseqs2/{name}.fasta', sequence_key_path='data/interim/sequence_key.txt')
    resources:
        runtime="10m",
        mem="10MB"
    run:
        output_dir = os.path.dirname(output[0])
        if not os.path.exists(output_dir):
            print('Creating output dir')
            os.mkdir(output_dir)

        print('Loading sequence key')
        with open(input.sequence_key, 'r') as fh:
            sequence_key = dict(line.strip().split() for line in fh if line)

        past_ids = set()
        with open(input.benchmark_structures, 'r') as fh:
            reader = csv.DictReader(fh)
            for entry in reader:
                for sequence_name in (
                    'Achain_sequence', 'Bchain_sequence', 'mhc_chain1_sequence', 'mhc_chain2_sequence'
                ):
                    if sequence_name not in entry:
                        continue

                    sequence = entry[sequence_name]
                    if not sequence:
                        continue

                    sequence_id = hashlib.sha256(sequence.encode()).hexdigest()[:8]
                    if sequence_id in sequence_key or sequence_id in past_ids:
                        continue

                    with open(os.path.join(output_dir, sequence_id + '.fasta'), 'w') as output_fh:
                        output_fh.write(f'>{sequence_id}\n{sequence}\n')

                    past_ids.add(sequence_id)

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

rule create_colabfold_msas_and_find_templates_benchmark_structures:
    input:
        unpaired_msas=_get_msas_for_benchmark('data/interim/benchmark_structures.csv', 'data/interim/msas_mmseqs2_colabfold/{name}_unpaired.a3m'),
        paired_msas=_get_msas_for_benchmark('data/interim/benchmark_structures.csv', 'data/interim/msas_mmseqs2_colabfold/{name}_paired.a3m'),
        templates=_get_msas_for_benchmark('data/interim/benchmark_structures.csv', 'data/interim/templates_mmseqs2_colabfold/{name}')

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

rule create_tcrmodel2_msas_and_find_templates_benchmark_structures:
    input:
        unpaired_msas=_get_msas_for_benchmark('data/interim/benchmark_structures.csv', 'data/interim/msas_af3_tcrmodel2/{name}_unpaired.a3m'),
        paired_msas=_get_msas_for_benchmark('data/interim/benchmark_structures.csv', 'data/interim/msas_af3_tcrmodel2/{name}_paired.a3m'),
        templates=_get_msas_for_benchmark('data/interim/benchmark_structures.csv', 'data/interim/templates_af3_tcrmodel2/{name}')

rule create_tcrmodel2_msas_and_find_templates_af_training_structures:
    input:
        unpaired_msas=_get_msas_for_benchmark('data/interim/af_training_structures.csv', 'data/interim/msas_af3_tcrmodel2/{name}_unpaired.a3m'),
        paired_msas=_get_msas_for_benchmark('data/interim/af_training_structures.csv', 'data/interim/msas_af3_tcrmodel2/{name}_paired.a3m'),
        templates=_get_msas_for_benchmark('data/interim/af_training_structures.csv', 'data/interim/templates_af3_tcrmodel2/{name}')

def _get_msas_for_dash(path: str, template: str, peptide: str) -> list[str]:
    if not os.path.exists(path):
        print(path, 'does not exist!')
        return []

    with open(path, 'r') as fh:
        reader = csv.DictReader(fh)
        ids = {
            hashlib.sha256(sequence.encode()).hexdigest()[:8]
            for entry in reader
            for sequence_name in ('tcr_alpha_chain', 'tcr_beta_chain', 'mhc_chain_1', 'mhc_chain_2')
            if (sequence := entry[sequence_name]) and entry['peptide'] == peptide
        }

    return expand(template, name=ids)

rule create_tcrmodel2_msas_and_find_templates_dash_bmlf:
    input:
        unpaired_msas=_get_msas_for_dash("data/interim/sequences_dash_with_negatives.csv", "data/interim/msas_af3_tcrmodel2/{name}_unpaired.a3m", "GLCTLVAML"),
        paired_msas=_get_msas_for_dash("data/interim/sequences_dash_with_negatives.csv", "data/interim/msas_af3_tcrmodel2/{name}_paired.a3m", "GLCTLVAML"),
        templates=_get_msas_for_dash("data/interim/sequences_dash_with_negatives.csv", "data/interim/templates_af3_tcrmodel2/{name}", "GLCTLVAML")

rule combine_msas_to_parquet:
    input:
        unpaired_msa="data/interim/msas_{msa_method}_{data_source}/{name}_unpaired.a3m",
        paired_msa="data/interim/msas_{msa_method}_{data_source}/{name}_paired.a3m"
    output: "data/interim/msas_{msa_method}_{data_source}/{name}.aligned.pqt"
    resources:
        runtime="5m",
        mem="5GB"
    script: '../scripts/combine_msas_to_parquet.py'

rule combine_msas_to_csv:
    input:
        unpaired_msa="data/interim/msas_{msa_method}_{data_source}/{name}_unpaired.a3m",
        paired_msa="data/interim/msas_{msa_method}_{data_source}/{name}_paired.a3m"
    output: "data/interim/msas_{msa_method}_{data_source}/{name}.csv"
    resources:
        runtime="5m",
        mem="5GB"
    script: '../scripts/combine_msas_to_csv.py'

rule convert_template_indices_to_structure:
    input: "data/interim/templates_{msa_method}_{data_source}/{name}"
    output: directory("data/interim/templates_{msa_method}_{data_source}_converted_indices/{name}")
    resources:
        runtime="1m",
        mem="1GB"
    shell:
        """
        mkdir {output}

        for template in {input}/*.cif; do
            template_name=$(basename $template .cif)
            indices="$(dirname $template)/${{template_name}}_indices.json"

            python -m tcr_pmhc_predicted_structures.apps.convert_index_mapping_from_sequence_to_structure \
                --log-level {config[log_level]} \
                -o {output}/${{template_name}}_indices.json \
                $template \
                $indices

            cp $template {output}/${{template_name}}.cif
        done
        """

rule create_m8_file_from_templates:
    input:
        msa_input="data/interim/msa-inputs_{msa_method}/{name}.json",
        templates="data/interim/templates_{msa_method}_{data_source}_converted_indices/{name}"
    output: "data/interim/templates_{msa_method}_{data_source}_converted_indices/{name}.m8"
    resources:
        runtime="5m",
        mem="1GB"
    script: '../scripts/create_m8_file_from_templates.py'
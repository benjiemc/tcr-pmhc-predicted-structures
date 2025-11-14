include: "utils.smk"

wildcard_constraints:
    name="[A-Za-z0-9_]+",
    seed="\d+",
    sample="\d+",
    model="[A-Za-z0-9]+",
    msa_method="[A-Za-z0-9]+",
    data_source="[A-Za-z0-9]+",
    sample_size="\d+"

rule extract_confidence_scores:
    input:
        summary_confidences="data/interim/structures_{model}_{msa_method}_{data_source}/{name}_seed-{seed}_sample-{sample}_summary_confidences.json",
        confidences="data/interim/structures_{model}_{msa_method}_{data_source}/{name}_seed-{seed}_sample-{sample}_confidences.json",
        sequences="data/interim/sequences/{name}.json"
    output: temp("data/interim/{model}_{msa_method}_{data_source}_{name}_{seed}_{sample}_confidence_scores.csv")
    resources:
        runtime="20m",
        mem="1GB"
    script: "../scripts/extract_{wildcards.model}_confidence_metrics.py"

rule aggregate_confidence_scores:
    input: lambda wildcards: get_predicted_structure_names(f"data/interim/sequences_sample_{wildcards.sample_size}.csv", f"data/interim/{wildcards.model}_{wildcards.msa_method}_{wildcards.data_source}_{{name}}_{{seed}}_{{sample}}_confidence_scores.csv", config['seed'])
    output: "data/interim/confidence_scores_{model}_{msa_method}_{data_source}_{sample_size}.csv"
    resources:
        runtime="10m",
        mem="1GB"
    run:
        import pandas as pd

        entries = []
        for path in input:
            confidence_scores = pd.read_csv(path)

            _, _, _, name, seed, sample = os.path.basename(path).replace('_confidence_scores.csv', '').rsplit('_', 5)
            confidence_scores['name'] = name
            confidence_scores['seed'] = seed
            confidence_scores['sample'] = sample

            entries.append(confidence_scores)

        df = pd.concat(entries).reset_index(drop=True)
        df.to_csv(output[0], index=False)

rule summarise_predicted_structural_features:
    input: "data/interim/structures_{model}_{msa_method}_{data_source}/{name}_seed-{seed}_sample-{sample}_model.cif"
    output: temp("data/interim/{model}_{msa_method}_{data_source}_{name}_{seed}_{sample}_structural_summary.json")
    resources:
        runtime="20m",
        mem="2GB"
    script: "../scripts/summarise_structural_features.py"

rule aggregate_predicted_structural_features:
    input: lambda wildcards: get_predicted_structure_names(f"data/interim/sequences_sample_{wildcards.sample_size}.csv", f"data/interim/{wildcards.model}_{wildcards.msa_method}_{wildcards.data_source}_{{name}}_{{seed}}_{{sample}}_structural_summary.json", config['seed'])
    output: "data/processed/predicted_{model}_{msa_method}_{data_source}_structural_features_{sample_size}.csv"
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

            _, _, _, name, seed, sample = os.path.basename(path).replace('_structural_summary.json', '').rsplit('_', 5)
            structural_features['name'] = name
            structural_features['seed'] = seed
            structural_features['sample'] = sample

            entries.append(structural_features)

        df = pd.DataFrame(entries)
        df.to_csv(output[0], index=False)

rule process_predicted_structure:
    input: "data/interim/structures_{model}_{msa_method}_{data_source}/{name}_seed-{seed}_sample-{sample}_model.cif"
    output: "data/processed/structures_{model}_{msa_method}_{data_source}/{name}_{seed}_{sample}.pdb"
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

rule create_predicted_contact_map:
    input: "data/processed/structures_{model}_{msa_method}_{data_source}/{name}_{seed}_{sample}.pdb"
    output: temp("data/interim/{model}_{msa_method}_{data_source}_{name}_{seed}_{sample}_contacts.csv")
    params:
        contact_distance=5.0
    resources:
        runtime="10m",
        mem="1GB"
    shell:
        """
        columns="cdr_type,chain_id_tcr,chain_type_tcr,residue_seq_id_tcr,residue_insert_code_tcr,resi_tcr,residue_name_tcr,"
        columns+="cdr_tcr,cdr_position,relative_pos_centre_tcr,chain_id_pmhc,residue_seq_id_pmhc,residue_insert_code_pmhc,"
        columns+="resi_pmhc,residue_name_pmhc,chain_type_pmhc,peptide_position,relative_pos_centre_pmhc,mhc_type"

        python -m tcr_pmhc_predicted_structures.apps.create_contact_map \
            --log-level {config[log_level]} \
            --cutoff-distance {params.contact_distance} \
            -o {output} \
            {input} || echo $columns > {output}
        """

rule aggregate_predicted_contact_maps:
    input: lambda wildcards: get_predicted_structure_names(f"data/interim/sequences_sample_{wildcards.sample_size}.csv", f"data/interim/{wildcards.model}_{wildcards.msa_method}_{wildcards.data_source}_{{name}}_{{seed}}_{{sample}}_contacts.csv", config['seed'])
    output: "data/processed/predicted_{model}_{msa_method}_{data_source}_contact_maps_{sample_size}.csv"
    resources:
        runtime="10m",
        mem="1GB"
    run:
        import pandas as pd

        contact_maps = []
        for path in input:
            contact_map = pd.read_csv(path)

            _, _, _, name, seed, sample = os.path.basename(path).replace('_contacts.csv', '').rsplit('_', 5)
            contact_map['name'] = name
            contact_map['seed'] = seed
            contact_map['sample'] = sample

            contact_maps.append(contact_map)

        contact_maps = pd.concat(contact_maps)
        contact_maps.to_csv(output[0], index=False)

rule profile_predicted_interactions:
    input: "data/processed/structures_{model}_{msa_method}_{data_source}/{name}_{seed}_{sample}.pdb"
    output: temp("data/interim/{model}_{msa_method}_{data_source}_{name}_{seed}_{sample}_interactions.csv")
    resources:
        runtime="10m",
        mem="2GB"
    script: "../scripts/profile_interactions.py"

rule aggregate_predicted_interactions:
    input: lambda wildcards: get_predicted_structure_names(f"data/interim/sequences_sample_{wildcards.sample_size}.csv", f"data/interim/{wildcards.model}_{wildcards.msa_method}_{wildcards.data_source}_{{name}}_{{seed}}_{{sample}}_interactions.csv", config['seed'])
    output: "data/processed/predicted_{model}_{msa_method}_{data_source}_interactions_{sample_size}.csv"
    resources:
        runtime="10m",
        mem="1GB"
    run:
        import pandas as pd

        interactions = []
        for path in input:
            interaction = pd.read_csv(path)

            _, _, _, name, seed, sample = os.path.basename(path).replace('_interactions.csv', '').rsplit('_', 5)
            interaction['name'] = name
            interaction['seed'] = seed
            interaction['sample'] = sample

            interactions.append(interaction)

        interactions = pd.concat(interactions)
        interactions.to_csv(output[0], index=False)

rule calculate_dockq_score:
    input:
        prediction="data/interim/structures_{model}_{msa_method}_{data_source}/{name}_seed-{seed}_sample-{sample}_model.cif",
        structure="data/raw/structures_no_crop/{name}.pdb"
    output: temp("data/interim/{model}_{msa_method}_{data_source}_{name}_{seed}_{sample}_dockq.json")
    resources:
        runtime="5m",
        mem="1GB"
    script: "../scripts/calculate_dockq_score.py"

rule aggregate_dockq_scores:
    input: lambda wildcards: get_predicted_structure_names(f"data/interim/benchmark_structures.csv", f"data/interim/{wildcards.model}_{wildcards.msa_method}_{wildcards.data_source}_{{name}}_{{seed}}_{{sample}}_dockq.json", config['seed'], column_name='id')
    output: "data/processed/predicted_{model}_{msa_method}_{data_source}_dockq_benchmark.csv"
    resources:
        runtime="1m",
        mem="10MB"
    run:
        import json
        import pandas as pd

        dockq_scores = []
        for path in input:
            with open(path, 'r') as fh:
                dockq_score = json.load(fh)

            info = re.search(
                r'[A-Za-z0-9]+_[A-Za-z0-9]+_[A-Za-z0-9]+_(?P<name>[A-Za-z0-9_]+)_(?P<seed>\d+)_(?P<sample>\d+)'
                '_dockq.json',
                os.path.basename(path),
            ).groupdict()

            for key, value in info.items():
                dockq_score[key] = value

            dockq_scores.append(dockq_score)

        dockq_scores = pd.DataFrame(dockq_scores)
        dockq_scores = dockq_scores[[
            'name',
            'seed',
            'sample',
            'DockQ',
            'F1',
            'iRMSD',
            'LRMSD',
            'fnat',
            'nat_correct',
            'nat_total',
            'fnonnat',
            'nonnat_count',
            'model_total',
            'clashes',
            'len1',
            'len2',
        ]]
        dockq_scores.to_csv(output[0], index=False)

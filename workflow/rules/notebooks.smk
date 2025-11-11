include: "utils.smk"

rule run_compare_structural_features_to_confidence:
    input:
        "data/external/tcr_pmhc_contacts.csv",
        "data/interim/sequences_sample_1500.csv",
        "data/processed/confidence_scores_af3_tcrmodel2_1500.csv",
        "data/processed/predicted_af3_tcrmodel2_structural_features_1500.csv",
        "data/processed/predicted_af3_tcrmodel2_contact_maps_1500.csv",
        "data/processed/predicted_af3_tcrmodel2_interactions_1500.csv",
        "data/processed/tcr_pmhc_structural_features.csv",
        "data/processed/tcr_pmhc_interactions.csv"
    output: "data/flags/compare_structural_features_to_confidence_done.flag"
    resources:
        runtime="10m",
        mem="1GB"
    shell:
        """
        jupyter nbconvert \
            --to notebook \
            --execute \
            --inplace \
            notebooks/compare_structural_features_to_confidence.ipynb
        touch {output}
        """

rule run_investigate_specificity_information_in_af3_scores:
    input:
        "data/processed/confidence_scores_af3_tcrmodel2_1500.csv",
        "data/interim/sequences_sample_1500.csv"
    output: "data/flags/investigate_specificity_information_in_af3_scores_done.flag"
    resources:
        runtime="5m",
        mem="1GB"
    shell:
        """
        jupyter nbconvert \
            --to notebook \
            --execute \
            --inplace \
            notebooks/investigate_specificity_information_in_af3_scores.ipynb
        touch {output}
        """

def _get_chai_names(path: str) -> list[str]:
    if not os.path.exists(path):
        print(path, 'does not exist!')
        return []

    with open(path, 'r') as fh:
        return [f"data/interim/structures-chai-restraints/{entry['index']}" for entry in csv.DictReader(fh)]

rule run_investigate_specificity_information_in_chai_scores:
    input: _get_chai_names("data/interim/sequences_sample_50.csv")
    output: "data/flags/investigate_specificity_information_in_chai_scores_done.flag"
    resources:
        runtime="10m",
        mem="1GB"
    shell:
        """
        jupyter nbconvert \
            --to notebook \
            --execute \
            --inplace \
            notebooks/investigate_specificity_information_in_chai_scores.ipynb
        touch {output}
        """

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

rule run_compare_msa_generation:
    input:
        'data/interim/sequences_sample_50.csv',
        'data/interim/sequence_key.txt',
        'data/interim/sequences_full.csv',
        _get_msas_for_sample("data/interim/sequences_sample_50.csv", "data/interim/sequences_full.csv", "data/interim/sequence_key.txt", "data/benchmarks/msa-generation_af3_af3_{name}.tsv"),
        _get_msas_for_sample("data/interim/sequences_sample_50.csv", "data/interim/sequences_full.csv", "data/interim/sequence_key.txt", "data/benchmarks/msa-generation_af3_tcrmodel2_{name}.tsv"),
        _get_msas_for_sample("data/interim/sequences_sample_50.csv", "data/interim/sequences_full.csv", "data/interim/sequence_key.txt", "data/benchmarks/msa-generation_mmseqs2_colabfold_{name}.tsv"),
        'data/processed/confidence_scores_af3_af3_50.csv',
        'data/processed/confidence_scores_af3_tcrmodel2_50.csv',
        'data/processed/confidence_scores_mmseqs2_colabfold_50.csv',
        'data/processed/predicted_af3_af3_dockq_benchmark.csv',
        'data/processed/predicted_af3_tcrmodel2_dockq_benchmark.csv',
        get_predicted_structure_names(f"data/interim/benchmark_structures.csv", "data/interim/structures_af3_af3/{name}_seed-{seed}_sample-{sample}_summary_confidences.json", config['seed'], column_name='id'),
        get_predicted_structure_names(f"data/interim/benchmark_structures.csv", "data/interim/structures_af3_tcrmodel2/{name}_seed-{seed}_sample-{sample}_summary_confidences.json", config['seed'], column_name='id')
    output: "data/flags/compare_msa_generation_done.flag"
    resources:
        runtime="5m",
        mem="1GB"
    shell:
        """
        jupyter nbconvert \
            --to notebook \
            --execute \
            --inplace \
            notebooks/compare_msa_generation.ipynb
        touch {output}
        """

rule run_specificity_predictor:
    input:
        "data/processed/confidence_scores_af3_tcrmodel2_1500.csv",
        "data/processed/predicted_af3_tcrmodel2_structural_features_1500.csv",
        "data/interim/sequences_sample_1500.csv"
    output: "data/flags/specificity_predictor_done.flag"
    resources:
        runtime="15m",
        mem="5GB"
    shell:
        """
        jupyter nbconvert \
            --to notebook \
            --execute \
            --inplace \
            notebooks/specificity_predictor.ipynb
        touch {output}
        """
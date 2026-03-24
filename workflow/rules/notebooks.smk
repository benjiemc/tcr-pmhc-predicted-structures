include: "utils.smk"

rule run_exploring_dash_dataset:
    input:
        "data/external/dash_human.csv",
        "data/interim/confidence_scores_af3_af3_tcrmodel2_dash_bmlf.csv",
        "data/processed/predicted_af3_af3_tcrmodel2_structural_features_dash_bmlf.csv",
        "data/processed/predicted_af3_af3_tcrmodel2_contact_maps_dash_bmlf.csv",
        "data/processed/predicted_af3_af3_tcrmodel2_interactions_dash_bmlf.csv"
    output:
        "data/flags/exploring_dash_dataset_done.flag"
    resources:
        runtime="10m",
        mem="1GB"
    shell:
        """
        jupyter nbconvert \
            --to notebook \
            --execute \
            --inplace \
            notebooks/exploring_dash_dataset.ipynb
        touch {output}
        """

rule run_exploring_coeliac_tcrs:
    input:
        "data/interim/sequences_coeliac.csv",
        "data/processed/predicted_af3_af3_tcrmodel2_structural_features_sequences_coeliac.csv",
        "data/processed/predicted_af3_af3_tcrmodel2_contact_maps_coeliac.csv",
        "data/processed/predicted_af3_af3_tcrmodel2_interactions_coeliac.csv"
    output:
        "data/flags/exploring_coeliac_tcrs_done.flag"
    resources:
        runtime="10m",
        mem="1GB"
    shell:
        """
        jupyter nbconvert \
            --to notebook \
            --execute \
            --inplace \
            notebooks/exploring_coeliac_tcrs.ipynb
        touch {output}
        """

rule run_coeliac_antigen_specificity_with_alphafold3:
    input:
        "data/interim/sequences_coeliac_with_negatives.csv",
        "data/processed/confidence_scores_af3_af3_tcrmodel2_coeliac.csv"
    output:
        "data/flags/coeliac_antigen_specificity_with_alphafold3_done.flag"
    resources:
        runtime="10m",
        mem="1GB"
    shell:
        """
        jupyter nbconvert \
            --to notebook \
            --execute \
            --inplace \
            notebooks/coeliac_antigen_specificity_with_alphafold3.ipynb
        touch {output}
        """

rule run_coeliac_mutagenesis_analysis:
    input:
        "data/interim/sequences_coeliac_mutations.csv",
        "data/processed/predicted_af3_af3_tcrmodel2_dockq_coeliac.csv",
        "data/processed/confidence_scores_af3_af3_tcrmodel2_coeliac_mutations.csv"
    output:
        "data/flags/coeliac_mutagenesis_analysis_done.flag"
    resources:
        runtime="10m",
        mem="1GB"
    shell:
        """
        jupyter nbconvert \
            --to notebook \
            --execute \
            --inplace \
            notebooks/coeliac_mutagenesis_analysis.ipynb
        touch {output}
        """

rule run_compare_structural_features_to_confidence:
    input:
        "data/external/tcr_pmhc_contacts.csv",
        "data/interim/sequences_sample_1500.csv",
        "data/processed/confidence_scores_af3_af3_tcrmodel2_1500.csv",
        "data/processed/predicted_af3_af3_tcrmodel2_structural_features_1500.csv",
        "data/processed/predicted_af3_af3_tcrmodel2_contact_maps_1500.csv",
        "data/processed/predicted_af3_af3_tcrmodel2_interactions_1500.csv",
        "data/processed/tcr_pmhc_structural_features.csv",
        "data/processed/tcr_pmhc_interactions.csv"
    output:
        "report/figures/scanning_angle.svg",
        "report/figures/pitch_angle.svg",
        "report/figures/tilt_angle.svg",
        "report/figures/roll_angle.svg",
        "report/figures/x_distance.svg",
        "report/figures/y_distance.svg",
        "report/figures/z_distance.svg",
        "report/figures/interaction_numbers.svg",
        "report/figures/interaction_types.svg",
        "report/figures/structrual_feature_intra_correlations.svg",
        "report/figures/ranked_confidence_metric_correlations.svg",
        "report/figures/structural_feature_inter_correlations.svg",
        "report/figures/structural_feature_inter_correlations_no_outliers.svg",
        "report/figures/structural_feature_inter_correlations_no_outliers_bootstrap.svg",
        "data/flags/compare_structural_features_to_confidence_done.flag"
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
        "data/processed/confidence_scores_af3_af3_tcrmodel2_1500.csv",
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
        'data/processed/confidence_scores_af3_af3_af3_50.csv',
        'data/processed/confidence_scores_af3_af3_tcrmodel2_50.csv',
        'data/processed/confidence_scores_af3_mmseqs2_colabfold_50.csv',
        'data/processed/predicted_af3_af3_af3_dockq_benchmark.csv',
        'data/processed/predicted_af3_af3_tcrmodel2_dockq_benchmark.csv',
        get_predicted_structure_names(f"data/interim/benchmark_structures.csv", "data/interim/structures_af3_af3_af3/{name}_seed-{seed}_sample-{sample}_summary_confidences.json", config['seed'], column_name='id'),
        get_predicted_structure_names(f"data/interim/benchmark_structures.csv", "data/interim/structures_af3_af3_tcrmodel2/{name}_seed-{seed}_sample-{sample}_summary_confidences.json", config['seed'], column_name='id')
    output:
        "report/figures/msa_strategy_generation_time.svg",
        "report/figures/msa_strategy_specificity.svg",
        "report/figures/msa_strategy_dockq_boxplot.svg",
        "report/figures/msa_strategy_dockq_capri.svg",
        "data/flags/compare_msa_generation_done.flag"
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

rule run_evaluating_antigen_specificity_model:
    input:
        "data/interim/sequences_sample_1500.csv",
        "data/processed/confidence_scores_af3_af3_tcrmodel2_1500.csv",
        "data/processed/af3_af3_tcrmodel2_pae_aggregator_predictions_untrained.csv",
        "data/processed/af3_af3_tcrmodel2_pae_aggregator_predictions_trained.csv",
        "models/PAEAggregator_af3_af3_tcrmodel2"
    output: "data/flags/evaluating_antigen_specificity_model_done.flag"
    resources:
        runtime="10m",
        mem="5GB"
    shell:
        """
        jupyter nbconvert \
            --to notebook \
            --execute \
            --inplace \
            notebooks/evaluating_antigen_specificity_model.ipynb
        touch {output}
        """

rule run_exploring_specificity_prediction_normalisation:
    input:
        "data/interim/sequences_sample_1500_af_training.csv",
        "data/processed/confidence_scores_af3_af3_tcrmodel2_1500.csv",
        "data/interim/confidence_scores_af3_af3_tcrmodel2_af_training.csv",
        "data/processed/af3_af3_tcrmodel2_pae_aggregator_predictions_untrained.csv",
        "data/processed/af3_af3_tcrmodel2_pae_aggregator_predictions_trained.csv",
        "data/processed/af3_af3_tcrmodel2_pae.h5",
        "models/PAEAggregator_af3_af3_tcrmodel2",
    output: "data/flags/exploring_specificity_prediction_normalisation_done.flag"
    resources:
        runtime="15m",
        mem="10GB"
    shell:
        """
        jupyter nbconvert \
            --to notebook \
            --execute \
            --inplace \
            notebooks/exploring_specificity_prediction_normalisation.ipynb
        touch {output}
        """

rule run_evaluating_out_of_distribution_detector:
    input:
        "data/interim/sequences_sample_1500_af_training.csv",
        "data/interim/af3_af3_tcrmodel2_full_pae_pca.csv"
    output: "data/flags/evaluating_out_of_distribution_detector_done.flag"
    resources:
        runtime="15m",
        mem="10GB"
    shell:
        """
        jupyter nbconvert \
            --to notebook \
            --execute \
            --inplace \
            notebooks/evaluating_out_of_distribution_detector.ipynb
        touch {output}
        """

rule run_benchmark_structure_prediction_methods:
    input:
        "data/external/stcrdab/db_summary.dat",
        "data/interim/sequences_sample_1500.csv",
        "data/processed/confidence_scores_af3_af3_tcrmodel2_1500.csv",
        "data/interim/confidence_scores_boltz_af3_tcrmodel2_1500.csv",
        "data/interim/confidence_scores_chai_af3_tcrmodel2_1500.csv",
        "data/interim/confidence_scores_af3_af3_tcrmodel2_benchmark.csv",
        "data/interim/confidence_scores_boltz_af3_tcrmodel2_benchmark.csv",
        "data/interim/confidence_scores_chai_af3_tcrmodel2_benchmark.csv",
        "data/interim/confidence_scores_af2_af3_tcrmodel2_benchmark.csv",
        "data/processed/predicted_af3_af3_tcrmodel2_dockq_benchmark.csv",
        "data/processed/predicted_boltz_af3_tcrmodel2_dockq_benchmark.csv",
        "data/processed/predicted_chai_af3_tcrmodel2_dockq_benchmark.csv",
        "data/processed/predicted_af2_af3_tcrmodel2_dockq_benchmark.csv"
    output:
        "report/figures/method_benchmark_dockq_boxplot.svg",
        "report/figures/method_benchmark_dockq_capri.svg",
        "data/flags/benchmark_structure_prediction_methods_done.flag"
    resources:
        runtime="5m",
        mem="1GB"
    shell:
        """
        jupyter nbconvert \
            --to notebook \
            --execute \
            --inplace \
            notebooks/benchmark_structure_prediction_methods.ipynb
        touch {output}
        """
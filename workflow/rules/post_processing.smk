include: "utils.smk"

wildcard_constraints:
    name="[A-Za-z0-9_]+",
    coeliac_name="[A-Za-z0-9_]+",
    structure_name="[A-Za-z0-9_]+",
    seed="\d+",
    sample="\d+",
    model="[A-Za-z0-9]+",
    msa_method="[A-Za-z0-9]+",
    data_source="[A-Za-z0-9]+",
    ablations="(_no_unpaired_msas|_no_paired_msas|_no_templates)*",
    unbalance="(_unbalanced)?",
    sample_size="\d+"

rule extract_confidence_scores:
    input:
        summary_confidences="data/interim/structures_{model}_{msa_method}_{data_source}{ablations}/{name}_seed-{seed}_sample-{sample}_summary_confidences.json",
        confidences="data/interim/structures_{model}_{msa_method}_{data_source}{ablations}/{name}_seed-{seed}_sample-{sample}_confidences.json",
        sequences="data/interim/sequences/{name}.json"
    output: temp("data/interim/{model}_{msa_method}_{data_source}{ablations}_{name}_{seed}_{sample}_confidence_scores.csv")
    resources:
        runtime="20m",
        mem="1GB"
    script: "../scripts/extract_{wildcards.model}_confidence_metrics.py"

rule aggregate_confidence_scores:
    input: lambda wildcards: get_predicted_structure_names(f"data/interim/sequences_sample_{wildcards.sample_size}{wildcards.unbalance}.csv", f"data/interim/{wildcards.model}_{wildcards.msa_method}_{wildcards.data_source}{wildcards.ablations}_{{name}}_{{seed}}_{{sample}}_confidence_scores.csv", config['seed'])
    output: "data/interim/confidence_scores_{model}_{msa_method}_{data_source}{ablations}_{sample_size}{unbalance}.csv"
    resources:
        runtime="10m",
        mem="1GB"
    script: "../scripts/aggregate_post_processing.py"

rule aggregate_benchmark_confidence_scores:
    input: lambda wildcards: get_predicted_structure_names("data/interim/benchmark_structures.csv", f"data/interim/{wildcards.model}_{wildcards.msa_method}_{wildcards.data_source}{wildcards.ablations}_{{name}}_{{seed}}_{{sample}}_confidence_scores.csv", config['seed'], column_name='id')
    output: "data/interim/confidence_scores_{model}_{msa_method}_{data_source}{ablations}_benchmark.csv"
    resources:
        runtime="10m",
        mem="1GB"
    script: "../scripts/aggregate_post_processing.py"

rule aggregate_af_training_confidence_scores:
    input: lambda wildcards: get_predicted_structure_names("data/interim/af_training_structures.csv", f"data/interim/{wildcards.model}_{wildcards.msa_method}_{wildcards.data_source}{wildcards.ablations}_{{name}}_{{seed}}_{{sample}}_confidence_scores.csv", config['seed'], column_name='id')
    output: "data/interim/confidence_scores_{model}_{msa_method}_{data_source}{ablations}_af_training.csv"
    resources:
        runtime="10m",
        mem="1GB"
    script: "../scripts/aggregate_post_processing.py"

rule aggregate_dash_bmlf_confidence_scores:
    input: lambda wildcards: get_dash_names("data/interim/sequences_dash_with_negatives.csv", f"data/interim/{wildcards.model}_{wildcards.msa_method}_{wildcards.data_source}_{{name}}_{{seed}}_{{sample}}_confidence_scores.csv", 'GLCTLVAML', config['seed'])
    output: "data/interim/confidence_scores_{model}_{msa_method}_{data_source}_dash_bmlf.csv"
    resources:
        runtime="10m",
        mem="1GB"
    script: "../scripts/aggregate_post_processing.py"

rule aggregate_mutated_dash_bmlf_confidence_scores:
    input: lambda wildcards: get_predicted_structure_names("data/interim/sequences_dash_BMLF_mutated.csv", f"data/interim/{wildcards.model}_{wildcards.msa_method}_{wildcards.data_source}_{{name}}_{{seed}}_{{sample}}_confidence_scores.csv", config['seed'], column_name='id')
    output: "data/interim/confidence_scores_{model}_{msa_method}_{data_source}_mutated_dash_bmlf.csv"
    resources:
        runtime="10m",
        mem="1GB"
    script: "../scripts/aggregate_post_processing.py"

rule aggregate_coeliac_confidence_scores:
    input: lambda wildcards: get_predicted_structure_names("data/interim/sequences_coeliac_with_negatives.csv", f"data/interim/{wildcards.model}_{wildcards.msa_method}_{wildcards.data_source}{wildcards.ablations}_{{name}}_{{seed}}_{{sample}}_confidence_scores.csv", config['seed'], column_name='id')
    output: "data/interim/confidence_scores_{model}_{msa_method}_{data_source}{ablations}_coeliac.csv"
    resources:
        runtime="10m",
        mem="1GB"
    script: "../scripts/aggregate_post_processing.py"

rule aggregate_coeliac_mutations_confidence_scores:
    input: lambda wildcards: get_predicted_structure_names("data/interim/sequences_coeliac_mutations.csv", f"data/interim/{wildcards.model}_{wildcards.msa_method}_{wildcards.data_source}{wildcards.ablations}_{{name}}_{{seed}}_{{sample}}_confidence_scores.csv", config['seed'], column_name='id')
    output: "data/interim/confidence_scores_{model}_{msa_method}_{data_source}{ablations}_coeliac_mutations.csv"
    resources:
        runtime="10m",
        mem="1GB"
    script: "../scripts/aggregate_post_processing.py"

rule summarise_predicted_structural_features:
    input: "data/interim/structures_{model}_{msa_method}_{data_source}{ablations}/{name}_seed-{seed}_sample-{sample}_model.cif"
    output: temp("data/interim/{model}_{msa_method}_{data_source}{ablations}_{name}_{seed}_{sample}_structural_summary.json")
    resources:
        runtime="20m",
        mem="2GB"
    script: "../scripts/summarise_structural_features.py"

rule aggregate_predicted_structural_features:
    input: lambda wildcards: get_predicted_structure_names(f"data/interim/sequences_sample_{wildcards.sample_size}{wildcards.unbalance}.csv", f"data/interim/{wildcards.model}_{wildcards.msa_method}_{wildcards.data_source}{wildcards.ablations}_{{name}}_{{seed}}_{{sample}}_structural_summary.json", config['seed'])
    output: "data/processed/predicted_{model}_{msa_method}_{data_source}{ablations}_structural_features_{sample_size}{unbalance}.csv"
    resources:
        runtime="10m",
        mem="1GB"
    script: "../scripts/aggregate_post_processing.py"

rule aggregate_dash_bmlf_predicted_structural_features:
    input: lambda wildcards: get_dash_names("data/interim/sequences_dash_with_negatives.csv", f"data/interim/{wildcards.model}_{wildcards.msa_method}_{wildcards.data_source}_{{name}}_{{seed}}_{{sample}}_structural_summary.json", 'GLCTLVAML', config['seed'])
    output: "data/processed/predicted_{model}_{msa_method}_{data_source}_structural_features_dash_bmlf.csv"
    resources:
        runtime="10m",
        mem="1GB"
    script: "../scripts/aggregate_post_processing.py"

rule aggregate_coeliac_predicted_structural_features:
    input: lambda wildcards: get_predicted_structure_names(f"data/interim/sequences_coeliac_with_negatives.csv", f"data/interim/{wildcards.model}_{wildcards.msa_method}_{wildcards.data_source}{wildcards.ablations}_{{name}}_{{seed}}_{{sample}}_structural_summary.json", config['seed'], column_name='id')
    output: "data/processed/predicted_{model}_{msa_method}_{data_source}{ablations}_structural_features_coeliac.csv"
    resources:
        runtime="10m",
        mem="1GB"
    script: "../scripts/aggregate_post_processing.py"

rule aggregate_coeliac_mutations_predicted_structural_features:
    input: lambda wildcards: get_predicted_structure_names(f"data/interim/sequences_coeliac_mutations.csv", f"data/interim/{wildcards.model}_{wildcards.msa_method}_{wildcards.data_source}{wildcards.ablations}_{{name}}_{{seed}}_{{sample}}_structural_summary.json", config['seed'], column_name='id')
    output: "data/processed/predicted_{model}_{msa_method}_{data_source}{ablations}_structural_features_coeliac_mutations.csv"
    resources:
        runtime="10m",
        mem="1GB"
    script: "../scripts/aggregate_post_processing.py"

rule process_predicted_structure:
    input: "data/interim/structures_{model}_{msa_method}_{data_source}{ablations}/{name}_seed-{seed}_sample-{sample}_model.cif"
    output: "data/processed/structures_{model}_{msa_method}_{data_source}{ablations}/{name}_{seed}_{sample}.pdb"
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
    input: "data/processed/structures_{model}_{msa_method}_{data_source}{ablations}/{name}_{seed}_{sample}.pdb"
    output: temp("data/interim/{model}_{msa_method}_{data_source}{ablations}_{name}_{seed}_{sample}_contacts.csv")
    params:
        contact_distance=5.0
    resources:
        runtime="10m",
        mem="5GB"
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
    input: lambda wildcards: get_predicted_structure_names(f"data/interim/sequences_sample_{wildcards.sample_size}{wildcards.unbalance}.csv", f"data/interim/{wildcards.model}_{wildcards.msa_method}_{wildcards.data_source}{wildcards.ablations}_{{name}}_{{seed}}_{{sample}}_contacts.csv", config['seed'])
    output: "data/processed/predicted_{model}_{msa_method}_{data_source}{ablations}_contact_maps_{sample_size}{unbalance}.csv"
    resources:
        runtime="10m",
        mem="1GB"
    script: "../scripts/aggregate_post_processing.py"

rule aggregate_dash_bmlf_predicted_contact_maps:
    input: lambda wildcards: get_dash_names("data/interim/sequences_dash_with_negatives.csv", f"data/interim/{wildcards.model}_{wildcards.msa_method}_{wildcards.data_source}_{{name}}_{{seed}}_{{sample}}_contacts.csv", 'GLCTLVAML', config['seed'])
    output: "data/processed/predicted_{model}_{msa_method}_{data_source}_contact_maps_dash_bmlf.csv"
    resources:
        runtime="10m",
        mem="1GB"
    script: "../scripts/aggregate_post_processing.py"

rule aggregate_coeliac_predicted_contact_maps:
    input: lambda wildcards: get_predicted_structure_names(f"data/interim/sequences_coeliac_with_negatives.csv", f"data/interim/{wildcards.model}_{wildcards.msa_method}_{wildcards.data_source}{wildcards.ablations}_{{name}}_{{seed}}_{{sample}}_contacts.csv", config['seed'], column_name='id')
    output: "data/processed/predicted_{model}_{msa_method}_{data_source}{ablations}_contact_maps_coeliac.csv"
    resources:
        runtime="10m",
        mem="1GB"
    script: "../scripts/aggregate_post_processing.py"

rule aggregate_coeliac_mutations_predicted_contact_maps:
    input: lambda wildcards: get_predicted_structure_names(f"data/interim/sequences_coeliac_mutations.csv", f"data/interim/{wildcards.model}_{wildcards.msa_method}_{wildcards.data_source}{wildcards.ablations}_{{name}}_{{seed}}_{{sample}}_contacts.csv", config['seed'], column_name='id')
    output: "data/processed/predicted_{model}_{msa_method}_{data_source}{ablations}_contact_maps_coeliac_mutations.csv"
    resources:
        runtime="10m",
        mem="1GB"
    script: "../scripts/aggregate_post_processing.py"

rule profile_predicted_interactions:
    input: "data/processed/structures_{model}_{msa_method}_{data_source}{ablations}/{name}_{seed}_{sample}.pdb"
    output: temp("data/interim/{model}_{msa_method}_{data_source}{ablations}_{name}_{seed}_{sample}_interactions.csv")
    resources:
        runtime="10m",
        mem="2GB"
    script: "../scripts/profile_interactions.py"

rule aggregate_predicted_interactions:
    input: lambda wildcards: get_predicted_structure_names(f"data/interim/sequences_sample_{wildcards.sample_size}{wildcards.unbalance}.csv", f"data/interim/{wildcards.model}_{wildcards.msa_method}_{wildcards.data_source}{wildcards.ablations}_{{name}}_{{seed}}_{{sample}}_interactions.csv", config['seed'])
    output: "data/processed/predicted_{model}_{msa_method}_{data_source}{ablations}_interactions_{sample_size}{unbalance}.csv"
    resources:
        runtime="10m",
        mem="1GB"
    script: "../scripts/aggregate_post_processing.py"

rule aggregate_dash_bmlf_predicted_interactions:
    input: lambda wildcards: get_dash_names("data/interim/sequences_dash_with_negatives.csv", f"data/interim/{wildcards.model}_{wildcards.msa_method}_{wildcards.data_source}_{{name}}_{{seed}}_{{sample}}_interactions.csv", 'GLCTLVAML', config['seed'])
    output: "data/processed/predicted_{model}_{msa_method}_{data_source}_interactions_dash_bmlf.csv"
    resources:
        runtime="10m",
        mem="1GB"
    script: "../scripts/aggregate_post_processing.py"

rule aggregate_coeliac_predicted_interactions:
    input: lambda wildcards: get_predicted_structure_names(f"data/interim/sequences_coeliac_with_negatives.csv", f"data/interim/{wildcards.model}_{wildcards.msa_method}_{wildcards.data_source}{wildcards.ablations}_{{name}}_{{seed}}_{{sample}}_interactions.csv", config['seed'], column_name='id')
    output: "data/processed/predicted_{model}_{msa_method}_{data_source}{ablations}_interactions_coeliac.csv"
    resources:
        runtime="10m",
        mem="1GB"
    script: "../scripts/aggregate_post_processing.py"

rule aggregate_coeliac_mutations_predicted_interactions:
    input: lambda wildcards: get_predicted_structure_names(f"data/interim/sequences_coeliac_mutations.csv", f"data/interim/{wildcards.model}_{wildcards.msa_method}_{wildcards.data_source}{wildcards.ablations}_{{name}}_{{seed}}_{{sample}}_interactions.csv", config['seed'], column_name='id')
    output: "data/processed/predicted_{model}_{msa_method}_{data_source}{ablations}_interactions_coeliac_mutations.csv"
    resources:
        runtime="10m",
        mem="1GB"
    script: "../scripts/aggregate_post_processing.py"

rule calculate_predicted_sc_and_bsa:
    input: "data/processed/structures_{model}_{msa_method}_{data_source}{ablations}/{name}_{seed}_{sample}.pdb"
    output: temp("data/interim/{model}_{msa_method}_{data_source}{ablations}_{name}_{seed}_{sample}_sc_and_bsa.json")
    params:
        seed=123
    resources:
        runtime="10m",
        mem="2GB"
    script: "../scripts/calculate_sc_and_bsa.py"

rule aggregate_predicted_sc_and_bsa:
    input: lambda wildcards: get_predicted_structure_names(f"data/interim/sequences_sample_{wildcards.sample_size}{wildcards.unbalance}.csv", f"data/interim/{wildcards.model}_{wildcards.msa_method}_{wildcards.data_source}{wildcards.ablations}_{{name}}_{{seed}}_{{sample}}_sc_and_bsa.json", config['seed'])
    output: "data/processed/predicted_{model}_{msa_method}_{data_source}{ablations}_sc_and_bsa_{sample_size}{unbalance}.csv"
    resources:
        runtime="10m",
        mem="1GB"
    script: "../scripts/aggregate_post_processing.py"

rule aggregate_dash_bmlf_predicted_sc_and_bsa:
    input: lambda wildcards: get_dash_names("data/interim/sequences_dash_with_negatives.csv", f"data/interim/{wildcards.model}_{wildcards.msa_method}_{wildcards.data_source}{wildcards.ablations}_{{name}}_{{seed}}_{{sample}}_sc_and_bsa.json", 'GLCTLVAML', config['seed'])
    output: "data/processed/predicted_{model}_{msa_method}_{data_source}{ablations}_sc_and_bsa_dash_bmlf.csv"

rule aggregate_coeliac_predicted_sc_and_bsa:
    input: lambda wildcards: get_predicted_structure_names(f"data/interim/sequences_coeliac.csv", f"data/interim/{wildcards.model}_{wildcards.msa_method}_{wildcards.data_source}{wildcards.ablations}_{{name}}_{{seed}}_{{sample}}_sc_and_bsa.json", config['seed'], column_name='id')
    output: "data/processed/predicted_{model}_{msa_method}_{data_source}{ablations}_sc_and_bsa_coeliac.csv"
    resources:
        runtime="10m",
        mem="1GB"
    script: "../scripts/aggregate_post_processing.py"

rule aggregate_coeliac_mutations_predicted_sc_and_bsa:
    input: lambda wildcards: get_predicted_structure_names(f"data/interim/sequences_coeliac_mutations.csv", f"data/interim/{wildcards.model}_{wildcards.msa_method}_{wildcards.data_source}{wildcards.ablations}_{{name}}_{{seed}}_{{sample}}_sc_and_bsa.json", config['seed'], column_name='id')
    output: "data/processed/predicted_{model}_{msa_method}_{data_source}{ablations}_sc_and_bsa_coeliac_mutations.csv"
    resources:
        runtime="10m",
        mem="1GB"
    script: "../scripts/aggregate_post_processing.py"

rule calculate_dockq_score:
    input:
        prediction="data/interim/structures_{model}_{msa_method}_{data_source}{ablations}/{name}_seed-{seed}_sample-{sample}_model.cif",
        structure="data/raw/structures_no_crop/{name}.pdb"
    output: temp("data/interim/{model}_{msa_method}_{data_source}{ablations}_{name}_{seed}_{sample}_dockq.json")
    resources:
        runtime="10m",
        mem="1GB"
    script: "../scripts/calculate_dockq_score.py"

rule calculate_coeliac_dockq_score:
    input:
        prediction="data/interim/structures_{model}_{msa_method}_{data_source}{ablations}/{coeliac_name}_seed-{seed}_sample-{sample}_model.cif",
        structure="data/raw/external_structures/{structure_name}.pdb"
    output: temp("data/interim/{model}_{msa_method}_{data_source}{ablations}_{coeliac_name}-{structure_name}_{seed}_{sample}_dockq.json")
    resources:
        runtime="10m",
        mem="1GB"
    script: "../scripts/calculate_dockq_score.py"

rule aggregate_dockq_scores:
    input: lambda wildcards: get_predicted_structure_names(f"data/interim/benchmark_structures.csv", f"data/interim/{wildcards.model}_{wildcards.msa_method}_{wildcards.data_source}{wildcards.ablations}_{{name}}_{{seed}}_{{sample}}_dockq.json", config['seed'], column_name='id')
    output: "data/processed/predicted_{model}_{msa_method}_{data_source}{ablations}_dockq_benchmark.csv"
    resources:
        runtime="10m",
        mem="500MB"
    script: "../scripts/aggregate_dockq_scores.py"

def _get_coeliac_comparisons(template: str, seeds: list[int]):
    name_pairs = (
        ('coeliac_a1_mut_original', 'DQ2_5a1_572G2_DECAB'), ('coeliac_w1_mut_original', 'DQ2_5w1_572G2_DECAB')
    )

    return [
        template.format(coeliac_name=coeliac_name, structure_name=structure_name, seed=seed, sample=sample)
        for coeliac_name, structure_name in name_pairs
        for seed in seeds
        for sample in range(5)
    ]

rule aggregate_coeliac_dockq_scores:
    input: lambda wildcards: _get_coeliac_comparisons(f"data/interim/{wildcards.model}_{wildcards.msa_method}_{wildcards.data_source}{wildcards.ablations}_{{coeliac_name}}-{{structure_name}}_{{seed}}_{{sample}}_dockq.json", config['seed'])
    output: "data/processed/predicted_{model}_{msa_method}_{data_source}{ablations}_dockq_coeliac.csv"
    resources:
        runtime="10m",
        mem="500MB"
    script: "../scripts/aggregate_dockq_scores.py"

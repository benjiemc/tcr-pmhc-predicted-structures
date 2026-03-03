include: "utils.smk"

wildcard_constraints:
    name="[A-Za-z0-9_]+",
    seed="\d+",
    sample="\d+",
    model="[A-Za-z0-9]+",
    msa_method="[A-Za-z0-9]+",
    data_source="[A-Za-z0-9]+",
    sample_size="\d+"

rule extract_padded_pae_matrices:
    input:
        confidences="data/interim/structures_{model}_{msa_method}_{data_source}/{name}_seed-{seed}_sample-{sample}_confidences.json",
        sequences="data/interim/sequences/{name}.json"
    output: "data/interim/pae_padded_{model}_{msa_method}_{data_source}/{name}_{seed}_{sample}_pae.txt"
    resources:
        runtime="10m",
        mem="1GB"
    shell:
        """
        file_name=$(basename {input.confidences}  _confidences.json)
        tmp_file_name=$TMPDIR/${{file_name}}_pae.txt
        jq -r '.pae[] | map(tostring) | join(" ")' {input.confidences} > $tmp_file_name


        output_dir=$(dirname {output})
        mkdir -p $output_dir

        python -m tcr_pmhc_predicted_structures.apps.extract_tcr_pmhc_confidence \
            --log-level {config[log_level]} \
            --score pae \
            --binding-only \
            --pad \
            --pad-direction centre \
            -o {output} \
            $tmp_file_name \
            {input.sequences}
        """

rule create_pae_dataset:
    input:
        sequences="data/interim/sequences_sample_1500_af_training.csv",
        confidences=lambda wildcards: get_predicted_structure_names("data/interim/sequences_sample_1500_af_training.csv", f"data/interim/pae_padded_{wildcards.model}_{wildcards.msa_method}_{wildcards.data_source}/{{name}}_{{seed}}_{{sample}}_pae.txt", config['seed'], column_name='id')
    output: "data/processed/{model}_{msa_method}_{data_source}_pae.h5"
    resources:
        runtime="10m",
        mem="5GB"
    shell:
        """
        python -m tcr_pmhc_predicted_structures.apps.create_pae_dataset \
            --log-level {config[log_level]} \
            --sequences {input.sequences} \
            --confidences {input.confidences} \
            -o {output}
        """

rule train_pae_aggregator:
    input:
        paes="data/processed/{model}_{msa_method}_{data_source}_pae.h5",
        pcas="data/interim/{model}_{msa_method}_{data_source}_full_pae_pca.csv"
    output: directory("models/PAEAggregator_{model}_{msa_method}_{data_source}")
    log: "data/logs/train_PAEAggregator_{model}_{msa_method}_{data_source}.log"
    resources:
        runtime="1h",
        mem="10GB",
        gpu=1
    params:
        num_epochs=15,
        batch_size=32,
        learning_rate=0.00001,
        eval_batch_size=300,
        eval_interval=50
    shell:
        """
        python -m tcr_pmhc_predicted_structures.apps.train_pae_aggregator \
            --log-level {config[log_level]} \
            --log-file {log} \
            --seed {config[seed]} \
            --num-epochs {params.num_epochs} \
            --batch-size {params.batch_size} \
            --learning-rate {params.learning_rate} \
            --eval-batch-size {params.eval_batch_size} \
            --eval-interval {params.eval_interval} \
            --scaling-data {input.pcas} \
            -o {output} \
            {input.paes}
        """

rule run_pae_aggregator:
    input:
        paes="data/processed/{model}_{msa_method}_{data_source}_pae.h5",
        pcas="data/interim/{model}_{msa_method}_{data_source}_full_pae_pca.csv",
        model="models/PAEAggregator_{model}_{msa_method}_{data_source}"
    output: "data/processed/{model}_{msa_method}_{data_source}_pae_aggregator_predictions_{training_status}.csv"
    wildcard_constraints:
        training_status="(untrained)|(trained)"
    params:
        model_size=3040
    resources:
        runtime="10m",
        mem="5GB"
    script: "../scripts/run_pae_aggregator.py"

rule create_full_pae_dataset:
    input:
        sequences="data/interim/sequences_sample_1500_af_training.csv",
        confidences=lambda wildcards: get_predicted_structure_names("data/interim/sequences_sample_1500_af_training.csv", f"data/interim/structures_{wildcards.model}_{wildcards.msa_method}_{wildcards.data_source}/{{name}}_seed-{{seed}}_sample-{{sample}}_confidences.json", config['seed'], column_name='id')
    output: "data/processed/{model}_{msa_method}_{data_source}_full_pae.h5"
    resources:
        runtime="3h",
        mem="10GB"
    script: "../scripts/create_pae_dataset_full.py"

rule decompose_pae_dataset:
    input: "data/processed/{model}_{msa_method}_{data_source}{portion}_pae.h5"
    output: "data/interim/{model}_{msa_method}_{data_source}{portion}_pae_pca.csv"
    wildcard_constraints:
        portion="(_full)?"
    params:
        num_components=2
    resources:
        runtime="1h",
        mem="80GB"
    script: "../scripts/decompose_pae_dataset.py"

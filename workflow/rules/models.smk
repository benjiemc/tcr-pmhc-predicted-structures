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
        sequences="data/interim/sequences_sample_1500.csv",
        confidences=lambda wildcards: get_predicted_structure_names("data/interim/sequences_sample_1500.csv", f"data/interim/pae_padded_{wildcards.model}_{wildcards.msa_method}_{wildcards.data_source}/{{name}}_{{seed}}_{{sample}}_pae.txt", config['seed'])
    output: "data/processed/{model}_{msa_method}_{data_source}_pae.h5"
    resources:
        runtime="10m",
        mem="5GB"
    script: "../scripts/create_pae_dataset.py"

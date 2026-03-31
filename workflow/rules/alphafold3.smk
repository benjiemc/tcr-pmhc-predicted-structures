import csv
import os

include: "utils.smk"

envvars:
    "ALPHAFOLD3_WEIGHTS_DIR"

wildcard_constraints:
    sample_size="\d+",
    unbalance="(_unbalanced)?",
    seed="\d+",
    sample="\d+",
    method="[A-Za-z0-9]+",
    data="[A-Za-z0-9]+",
    name="[A-Za-z0-9_]+"

rule run_create_af3_jsons:
    input:
        sequences="data/interim/sequences/{name}.json",
        unpaired_msas=lambda wildcards: get_sequence_ids(f"data/interim/sequences/{wildcards.name}.json", f'data/interim/msas_{wildcards.method}_{wildcards.data}/{{name}}_unpaired.a3m'),
        paired_msas=lambda wildcards: get_sequence_ids(f"data/interim/sequences/{wildcards.name}.json", f'data/interim/msas_{wildcards.method}_{wildcards.data}/{{name}}_paired.a3m'),
        templates=lambda wildcards: get_sequence_ids(f"data/interim/sequences/{wildcards.name}.json", f'data/interim/templates_{wildcards.method}_{wildcards.data}/{{name}}')
    output: "data/interim/sequences_{method}_{data}/{name}_{seed}.json"
    resources:
        runtime="5m",
        mem="1GB"
    shell:
        """
        out_dir_name=$(dirname {output[0]})

        python -m tcr_pmhc_predicted_structures.apps.convert_to_input \
            --log-level {config[log_level]} \
            --format alphafold3 \
            --alphafold3-model-seeds {wildcards.seed} \
            --add-paths \
            --relative-paths \
            --unpaired-msa-paths {input.unpaired_msas} \
            --paired-msa-paths {input.paired_msas} \
            --template-paths {input.templates} \
            -o $out_dir_name \
            {input.sequences}
        """

rule run_create_af3_jsons_no_templates:
    input:
        sequences="data/interim/sequences/{name}.json",
        unpaired_msas=lambda wildcards: get_sequence_ids(f"data/interim/sequences/{wildcards.name}.json", f'data/interim/msas_{wildcards.method}_{wildcards.data}/{{name}}_unpaired.a3m'),
        paired_msas=lambda wildcards: get_sequence_ids(f"data/interim/sequences/{wildcards.name}.json", f'data/interim/msas_{wildcards.method}_{wildcards.data}/{{name}}_paired.a3m'),
    output: "data/interim/sequences_{method}_{data}_no_templates/{name}_{seed}.json"
    resources:
        runtime="5m",
        mem="1GB"
    shell:
        """
        out_dir_name=$(dirname {output[0]})

        python -m tcr_pmhc_predicted_structures.apps.convert_to_input \
            --log-level {config[log_level]} \
            --format alphafold3 \
            --alphafold3-model-seeds {wildcards.seed} \
            --add-paths \
            --relative-paths \
            --unpaired-msa-paths {input.unpaired_msas} \
            --paired-msa-paths {input.paired_msas} \
            -o $out_dir_name \
            {input.sequences}
        """

rule run_create_af3_jsons_no_paired_msas:
    input:
        sequences="data/interim/sequences/{name}.json",
        unpaired_msas=lambda wildcards: get_sequence_ids(f"data/interim/sequences/{wildcards.name}.json", f'data/interim/msas_{wildcards.method}_{wildcards.data}/{{name}}_unpaired.a3m'),
        templates=lambda wildcards: get_sequence_ids(f"data/interim/sequences/{wildcards.name}.json", f'data/interim/templates_{wildcards.method}_{wildcards.data}/{{name}}')
    output: "data/interim/sequences_{method}_{data}_no_paired_msas/{name}_{seed}.json"
    resources:
        runtime="5m",
        mem="1GB"
    shell:
        """
        out_dir_name=$(dirname {output[0]})

        python -m tcr_pmhc_predicted_structures.apps.convert_to_input \
            --log-level {config[log_level]} \
            --format alphafold3 \
            --alphafold3-model-seeds {wildcards.seed} \
            --add-paths \
            --relative-paths \
            --unpaired-msa-paths {input.unpaired_msas} \
            --template-paths {input.templates} \
            -o $out_dir_name \
            {input.sequences}
        """

rule run_create_af3_jsons_no_unpaired_msas:
    input:
        sequences="data/interim/sequences/{name}.json",
        paired_msas=lambda wildcards: get_sequence_ids(f"data/interim/sequences/{wildcards.name}.json", f'data/interim/msas_{wildcards.method}_{wildcards.data}/{{name}}_paired.a3m'),
        templates=lambda wildcards: get_sequence_ids(f"data/interim/sequences/{wildcards.name}.json", f'data/interim/templates_{wildcards.method}_{wildcards.data}/{{name}}')
    output: "data/interim/sequences_{method}_{data}_no_unpaired_msas/{name}_{seed}.json"
    resources:
        runtime="5m",
        mem="1GB"
    shell:
        """
        out_dir_name=$(dirname {output[0]})

        python -m tcr_pmhc_predicted_structures.apps.convert_to_input \
            --log-level {config[log_level]} \
            --format alphafold3 \
            --alphafold3-model-seeds {wildcards.seed} \
            --add-paths \
            --relative-paths \
            --paired-msa-paths {input.paired_msas} \
            --template-paths {input.templates} \
            -o $out_dir_name \
            {input.sequences}
        """

rule run_create_af3_jsons_no_paired_msas_no_templates:
    input:
        sequences="data/interim/sequences/{name}.json",
        unpaired_msas=lambda wildcards: get_sequence_ids(f"data/interim/sequences/{wildcards.name}.json", f'data/interim/msas_{wildcards.method}_{wildcards.data}/{{name}}_unpaired.a3m'),
    output: "data/interim/sequences_{method}_{data}_no_paired_msas_no_templates/{name}_{seed}.json"
    resources:
        runtime="5m",
        mem="1GB"
    shell:
        """
        out_dir_name=$(dirname {output[0]})

        python -m tcr_pmhc_predicted_structures.apps.convert_to_input \
            --log-level {config[log_level]} \
            --format alphafold3 \
            --alphafold3-model-seeds {wildcards.seed} \
            --add-paths \
            --relative-paths \
            --unpaired-msa-paths {input.unpaired_msas} \
            -o $out_dir_name \
            {input.sequences}
        """

rule run_create_af3_jsons_no_unpaired_msas_no_templates:
    input:
        sequences="data/interim/sequences/{name}.json",
        paired_msas=lambda wildcards: get_sequence_ids(f"data/interim/sequences/{wildcards.name}.json", f'data/interim/msas_{wildcards.method}_{wildcards.data}/{{name}}_paired.a3m')
    output: "data/interim/sequences_{method}_{data}_no_unpaired_msas_no_templates/{name}_{seed}.json"
    resources:
        runtime="5m",
        mem="1GB"
    shell:
        """
        out_dir_name=$(dirname {output[0]})

        python -m tcr_pmhc_predicted_structures.apps.convert_to_input \
            --log-level {config[log_level]} \
            --format alphafold3 \
            --alphafold3-model-seeds {wildcards.seed} \
            --add-paths \
            --relative-paths \
            --paired-msa-paths {input.paired_msas} \
            -o $out_dir_name \
            {input.sequences}
        """

rule run_create_af3_jsons_no_unpaired_msas_no_paired_msas:
    input:
        sequences="data/interim/sequences/{name}.json",
        templates=lambda wildcards: get_sequence_ids(f"data/interim/sequences/{wildcards.name}.json", f'data/interim/templates_{wildcards.method}_{wildcards.data}/{{name}}')
    output: "data/interim/sequences_{method}_{data}_no_unpaired_msas_no_paired_msas/{name}_{seed}.json"
    resources:
        runtime="5m",
        mem="1GB"
    shell:
        """
        out_dir_name=$(dirname {output[0]})

        python -m tcr_pmhc_predicted_structures.apps.convert_to_input \
            --log-level {config[log_level]} \
            --format alphafold3 \
            --alphafold3-model-seeds {wildcards.seed} \
            --add-paths \
            --relative-paths \
            --template-paths {input.templates} \
            -o $out_dir_name \
            {input.sequences}
        """

rule run_create_af3_jsons_no_unpaired_msas_no_paired_msas_no_templates:
    input:
        sequences="data/interim/sequences/{name}.json"
    output: "data/interim/sequences_{method}_{data}_no_unpaired_msas_no_paired_msas_no_templates/{name}_{seed}.json"
    resources:
        runtime="5m",
        mem="1GB"
    shell:
        """
        out_dir_name=$(dirname {output[0]})

        python -m tcr_pmhc_predicted_structures.apps.convert_to_input \
            --log-level {config[log_level]} \
            --format alphafold3 \
            --alphafold3-model-seeds {wildcards.seed} \
            --add-paths \
            --relative-paths \
            -o $out_dir_name \
            {input.sequences}
        """

rule run_alphafold_predictions:
    input:
        sequence="data/interim/sequences_{method}_{data}/{name}_{seed}.json",
        unpaired_msas=lambda wildcards: get_sequence_ids(f"data/interim/sequences/{wildcards.name}.json", f'data/interim/msas_{wildcards.method}_{wildcards.data}/{{name}}_unpaired.a3m'),
        paired_msas=lambda wildcards: get_sequence_ids(f"data/interim/sequences/{wildcards.name}.json", f'data/interim/msas_{wildcards.method}_{wildcards.data}/{{name}}_paired.a3m'),
        templates=lambda wildcards: get_sequence_ids(f"data/interim/sequences/{wildcards.name}.json", f'data/interim/templates_{wildcards.method}_{wildcards.data}/{{name}}')
    output:
        structures=protected(expand("data/interim/structures_af3_{{method}}_{{data}}/{{name}}_seed-{{seed}}_sample-{sample}_model.cif", sample=range(5))),
        summary_confidences=protected(expand("data/interim/structures_af3_{{method}}_{{data}}/{{name}}_seed-{{seed}}_sample-{sample}_summary_confidences.json", sample=range(5))),
        confidences=protected(expand("data/interim/structures_af3_{{method}}_{{data}}/{{name}}_seed-{{seed}}_sample-{sample}_confidences.json", sample=range(5)))
    params:
        alphafold3_weights_dir=os.environ['ALPHAFOLD3_WEIGHTS_DIR']
    resources:
        runtime="10m",
        mem="40GB",
        gpu=1
    script: "../scripts/predict_structure_with_alphafold3.py"

rule run_alphafold_predictions_no_templates:
    input:
        sequence="data/interim/sequences_{method}_{data}_no_templates/{name}_{seed}.json",
        unpaired_msas=lambda wildcards: get_sequence_ids(f"data/interim/sequences/{wildcards.name}.json", f'data/interim/msas_{wildcards.method}_{wildcards.data}/{{name}}_unpaired.a3m'),
        paired_msas=lambda wildcards: get_sequence_ids(f"data/interim/sequences/{wildcards.name}.json", f'data/interim/msas_{wildcards.method}_{wildcards.data}/{{name}}_paired.a3m'),
    output:
        structures=protected(expand("data/interim/structures_af3_{{method}}_{{data}}_no_templates/{{name}}_seed-{{seed}}_sample-{sample}_model.cif", sample=range(5))),
        summary_confidences=protected(expand("data/interim/structures_af3_{{method}}_{{data}}_no_templates/{{name}}_seed-{{seed}}_sample-{sample}_summary_confidences.json", sample=range(5))),
        confidences=protected(expand("data/interim/structures_af3_{{method}}_{{data}}_no_templates/{{name}}_seed-{{seed}}_sample-{sample}_confidences.json", sample=range(5)))
    params:
        alphafold3_weights_dir=os.environ['ALPHAFOLD3_WEIGHTS_DIR']
    resources:
        runtime="10m",
        mem="40GB",
        gpu=1
    script: "../scripts/predict_structure_with_alphafold3.py"

rule run_alphafold_predictions_no_paired_msas:
    input:
        sequence="data/interim/sequences_{method}_{data}_no_paired_msas/{name}_{seed}.json",
        unpaired_msas=lambda wildcards: get_sequence_ids(f"data/interim/sequences/{wildcards.name}.json", f'data/interim/msas_{wildcards.method}_{wildcards.data}/{{name}}_unpaired.a3m'),
        templates=lambda wildcards: get_sequence_ids(f"data/interim/sequences/{wildcards.name}.json", f'data/interim/templates_{wildcards.method}_{wildcards.data}/{{name}}')
    output:
        structures=protected(expand("data/interim/structures_af3_{{method}}_{{data}}_no_paired_msas/{{name}}_seed-{{seed}}_sample-{sample}_model.cif", sample=range(5))),
        summary_confidences=protected(expand("data/interim/structures_af3_{{method}}_{{data}}_no_paired_msas/{{name}}_seed-{{seed}}_sample-{sample}_summary_confidences.json", sample=range(5))),
        confidences=protected(expand("data/interim/structures_af3_{{method}}_{{data}}_no_paired_msas/{{name}}_seed-{{seed}}_sample-{sample}_confidences.json", sample=range(5)))
    params:
        alphafold3_weights_dir=os.environ['ALPHAFOLD3_WEIGHTS_DIR']
    resources:
        runtime="10m",
        mem="40GB",
        gpu=1
    script: "../scripts/predict_structure_with_alphafold3.py"

rule run_alphafold_predictions_no_unpaired_msas:
    input:
        sequence="data/interim/sequences_{method}_{data}_no_unpaired_msas/{name}_{seed}.json",
        paired_msas=lambda wildcards: get_sequence_ids(f"data/interim/sequences/{wildcards.name}.json", f'data/interim/msas_{wildcards.method}_{wildcards.data}/{{name}}_paired.a3m'),
        templates=lambda wildcards: get_sequence_ids(f"data/interim/sequences/{wildcards.name}.json", f'data/interim/templates_{wildcards.method}_{wildcards.data}/{{name}}')
    output:
        structures=protected(expand("data/interim/structures_af3_{{method}}_{{data}}_no_unpaired_msas/{{name}}_seed-{{seed}}_sample-{sample}_model.cif", sample=range(5))),
        summary_confidences=protected(expand("data/interim/structures_af3_{{method}}_{{data}}_no_unpaired_msas/{{name}}_seed-{{seed}}_sample-{sample}_summary_confidences.json", sample=range(5))),
        confidences=protected(expand("data/interim/structures_af3_{{method}}_{{data}}_no_unpaired_msas/{{name}}_seed-{{seed}}_sample-{sample}_confidences.json", sample=range(5)))
    params:
        alphafold3_weights_dir=os.environ['ALPHAFOLD3_WEIGHTS_DIR']
    resources:
        runtime="10m",
        mem="40GB",
        gpu=1
    script: "../scripts/predict_structure_with_alphafold3.py"

rule run_alphafold_predictions_no_paired_msas_no_templates:
    input:
        sequence="data/interim/sequences_{method}_{data}_no_paired_msas_no_templates/{name}_{seed}.json",
        unpaired_msas=lambda wildcards: get_sequence_ids(f"data/interim/sequences/{wildcards.name}.json", f'data/interim/msas_{wildcards.method}_{wildcards.data}/{{name}}_unpaired.a3m'),
    output:
        structures=protected(expand("data/interim/structures_af3_{{method}}_{{data}}_no_paired_msas_no_templates/{{name}}_seed-{{seed}}_sample-{sample}_model.cif", sample=range(5))),
        summary_confidences=protected(expand("data/interim/structures_af3_{{method}}_{{data}}_no_paired_msas_no_templates/{{name}}_seed-{{seed}}_sample-{sample}_summary_confidences.json", sample=range(5))),
        confidences=protected(expand("data/interim/structures_af3_{{method}}_{{data}}_no_paired_msas_no_templates/{{name}}_seed-{{seed}}_sample-{sample}_confidences.json", sample=range(5)))
    params:
        alphafold3_weights_dir=os.environ['ALPHAFOLD3_WEIGHTS_DIR']
    resources:
        runtime="10m",
        mem="40GB",
        gpu=1
    script: "../scripts/predict_structure_with_alphafold3.py"

rule run_alphafold_predictions_no_unpaired_msas_no_templates:
    input:
        sequence="data/interim/sequences_{method}_{data}_no_unpaired_msas_no_templates/{name}_{seed}.json",
        paired_msas=lambda wildcards: get_sequence_ids(f"data/interim/sequences/{wildcards.name}.json", f'data/interim/msas_{wildcards.method}_{wildcards.data}/{{name}}_paired.a3m')
    output:
        structures=protected(expand("data/interim/structures_af3_{{method}}_{{data}}_no_unpaired_msas_no_templates/{{name}}_seed-{{seed}}_sample-{sample}_model.cif", sample=range(5))),
        summary_confidences=protected(expand("data/interim/structures_af3_{{method}}_{{data}}_no_unpaired_msas_no_templates/{{name}}_seed-{{seed}}_sample-{sample}_summary_confidences.json", sample=range(5))),
        confidences=protected(expand("data/interim/structures_af3_{{method}}_{{data}}_no_unpaired_msas_no_templates/{{name}}_seed-{{seed}}_sample-{sample}_confidences.json", sample=range(5)))
    params:
        alphafold3_weights_dir=os.environ['ALPHAFOLD3_WEIGHTS_DIR']
    resources:
        runtime="10m",
        mem="40GB",
        gpu=1
    script: "../scripts/predict_structure_with_alphafold3.py"

rule run_alphafold_predictions_no_unpaired_msas_no_paired_msas:
    input:
        sequence="data/interim/sequences_{method}_{data}_no_unpaired_msas_no_paired_msas/{name}_{seed}.json",
        templates=lambda wildcards: get_sequence_ids(f"data/interim/sequences/{wildcards.name}.json", f'data/interim/templates_{wildcards.method}_{wildcards.data}/{{name}}')
    output:
        structures=protected(expand("data/interim/structures_af3_{{method}}_{{data}}_no_unpaired_msas_no_paired_msas/{{name}}_seed-{{seed}}_sample-{sample}_model.cif", sample=range(5))),
        summary_confidences=protected(expand("data/interim/structures_af3_{{method}}_{{data}}_no_unpaired_msas_no_paired_msas/{{name}}_seed-{{seed}}_sample-{sample}_summary_confidences.json", sample=range(5))),
        confidences=protected(expand("data/interim/structures_af3_{{method}}_{{data}}_no_unpaired_msas_no_paired_msas/{{name}}_seed-{{seed}}_sample-{sample}_confidences.json", sample=range(5)))
    params:
        alphafold3_weights_dir=os.environ['ALPHAFOLD3_WEIGHTS_DIR']
    resources:
        runtime="10m",
        mem="40GB",
        gpu=1
    script: "../scripts/predict_structure_with_alphafold3.py"

rule run_alphafold_predictions_no_unpaired_msas_no_paired_msas_no_templates:
    input:
        sequence="data/interim/sequences_{method}_{data}_no_unpaired_msas_no_paired_msas_no_templates/{name}_{seed}.json"
    output:
        structures=protected(expand("data/interim/structures_af3_{{method}}_{{data}}_no_unpaired_msas_no_paired_msas_no_templates/{{name}}_seed-{{seed}}_sample-{sample}_model.cif", sample=range(5))),
        summary_confidences=protected(expand("data/interim/structures_af3_{{method}}_{{data}}_no_unpaired_msas_no_paired_msas_no_templates/{{name}}_seed-{{seed}}_sample-{sample}_summary_confidences.json", sample=range(5))),
        confidences=protected(expand("data/interim/structures_af3_{{method}}_{{data}}_no_unpaired_msas_no_paired_msas_no_templates/{{name}}_seed-{{seed}}_sample-{sample}_confidences.json", sample=range(5)))
    params:
        alphafold3_weights_dir=os.environ['ALPHAFOLD3_WEIGHTS_DIR']
    resources:
        runtime="10m",
        mem="40GB",
        gpu=1
    script: "../scripts/predict_structure_with_alphafold3.py"

rule create_af3_structures_sample_size_50:
    input: get_predicted_structure_names("data/interim/sequences_sample_50.csv", "data/interim/structures_af3_af3_af3/{name}_seed-{seed}_sample-{sample}_model.cif", config['seed'])

rule create_af3_structures_benchmark:
    input: get_predicted_structure_names("data/interim/benchmark_structures.csv", "data/interim/structures_af3_af3_af3/{name}_seed-{seed}_sample-{sample}_model.cif", config['seed'], column_name='id')

rule create_colabfold_structures_sample_size_50:
    input: get_predicted_structure_names("data/interim/sequences_sample_50.csv", "data/interim/structures_af3_mmseqs2_colabfold/{name}_seed-{seed}_sample-{sample}_model.cif", config['seed'])

rule create_colabfold_structures_benchmark:
    input: get_predicted_structure_names("data/interim/benchmark_structures.csv", "data/interim/structures_af3_mmseqs2_colabfold/{name}_seed-{seed}_sample-{sample}_model.cif", config['seed'], column_name='id')

rule create_af3_structures_with_tcrmodel2_data_sample_size_50:
    input: get_predicted_structure_names("data/interim/sequences_sample_50.csv", "data/interim/structures_af3_af3_tcrmodel2/{name}_seed-{seed}_sample-{sample}_model.cif", config['seed'])

rule create_af3_structures_with_tcrmodel2_data_sample_size_100_unbalanced:
    input: get_predicted_structure_names("data/interim/sequences_sample_100_unbalanced.csv", "data/interim/structures_af3_af3_tcrmodel2/{name}_seed-{seed}_sample-{sample}_model.cif", config['seed'])

rule create_af3_structures_with_tcrmodel2_data_sample_size_1500:
    input: get_predicted_structure_names("data/interim/sequences_sample_1500.csv", "data/interim/structures_af3_af3_tcrmodel2/{name}_seed-{seed}_sample-{sample}_model.cif", config['seed'])

rule create_af3_structures_with_tcrmodel2_data_benchmark:
    input: get_predicted_structure_names("data/interim/benchmark_structures.csv", "data/interim/structures_af3_af3_tcrmodel2/{name}_seed-{seed}_sample-{sample}_model.cif", config['seed'], column_name='id')

rule create_af3_structures_with_tcrmodel2_data_af_training:
    input: get_predicted_structure_names("data/interim/af_training_structures.csv", "data/interim/structures_af3_af3_tcrmodel2/{name}_seed-{seed}_sample-{sample}_model.cif", config['seed'], column_name='id')

rule create_af3_structures_with_tcrmodel2_data_no_templates_sample_size_1500:
    input: get_predicted_structure_names("data/interim/sequences_sample_1500.csv", "data/interim/structures_af3_af3_tcrmodel2_no_templates/{name}_seed-{seed}_sample-{sample}_model.cif", config['seed'])

rule create_af3_structures_with_tcrmodel2_data_no_templates_benchmark:
    input: get_predicted_structure_names("data/interim/benchmark_structures.csv", "data/interim/structures_af3_af3_tcrmodel2_no_templates/{name}_seed-{seed}_sample-{sample}_model.cif", config['seed'], column_name='id')

rule create_af3_structures_with_tcrmodel2_data_no_paired_msas_sample_size_1500:
    input: get_predicted_structure_names("data/interim/sequences_sample_1500.csv", "data/interim/structures_af3_af3_tcrmodel2_no_paired_msas/{name}_seed-{seed}_sample-{sample}_model.cif", config['seed'])

rule create_af3_structures_with_tcrmodel2_data_no_paired_msas_benchmark:
    input: get_predicted_structure_names("data/interim/benchmark_structures.csv", "data/interim/structures_af3_af3_tcrmodel2_no_paired_msas/{name}_seed-{seed}_sample-{sample}_model.cif", config['seed'], column_name='id')

rule create_af3_structures_with_tcrmodel2_data_no_unpaired_msas_sample_size_1500:
    input: get_predicted_structure_names("data/interim/sequences_sample_1500.csv", "data/interim/structures_af3_af3_tcrmodel2_no_unpaired_msas/{name}_seed-{seed}_sample-{sample}_model.cif", config['seed'])

rule create_af3_structures_with_tcrmodel2_data_no_unpaired_msas_benchmark:
    input: get_predicted_structure_names("data/interim/benchmark_structures.csv", "data/interim/structures_af3_af3_tcrmodel2_no_unpaired_msas/{name}_seed-{seed}_sample-{sample}_model.cif", config['seed'], column_name='id')

rule create_af3_structures_with_tcrmodel2_data_no_paired_msas_no_templates_sample_size_1500:
    input: get_predicted_structure_names("data/interim/sequences_sample_1500.csv", "data/interim/structures_af3_af3_tcrmodel2_no_paired_msas_no_templates/{name}_seed-{seed}_sample-{sample}_model.cif", config['seed'])

rule create_af3_structures_with_tcrmodel2_data_no_paired_msas_no_templates_benchmark:
    input: get_predicted_structure_names("data/interim/benchmark_structures.csv", "data/interim/structures_af3_af3_tcrmodel2_no_paired_msas_no_templates/{name}_seed-{seed}_sample-{sample}_model.cif", config['seed'], column_name='id')

rule create_af3_structures_with_tcrmodel2_data_no_unpaired_msas_no_templates_sample_size_1500:
    input: get_predicted_structure_names("data/interim/sequences_sample_1500.csv", "data/interim/structures_af3_af3_tcrmodel2_no_unpaired_msas_no_templates/{name}_seed-{seed}_sample-{sample}_model.cif", config['seed'])

rule create_af3_structures_with_tcrmodel2_data_no_unpaired_msas_no_templates_benchmark:
    input: get_predicted_structure_names("data/interim/benchmark_structures.csv", "data/interim/structures_af3_af3_tcrmodel2_no_unpaired_msas_no_templates/{name}_seed-{seed}_sample-{sample}_model.cif", config['seed'], column_name='id')

rule create_af3_structures_with_tcrmodel2_data_no_unpaired_msas_no_paired_msas_sample_size_1500:
    input: get_predicted_structure_names("data/interim/sequences_sample_1500.csv", "data/interim/structures_af3_af3_tcrmodel2_no_unpaired_msas_no_paired_msas/{name}_seed-{seed}_sample-{sample}_model.cif", config['seed'])

rule create_af3_structures_with_tcrmodel2_data_no_unpaired_msas_no_paired_msas_benchmark:
    input: get_predicted_structure_names("data/interim/benchmark_structures.csv", "data/interim/structures_af3_af3_tcrmodel2_no_unpaired_msas_no_paired_msas/{name}_seed-{seed}_sample-{sample}_model.cif", config['seed'], column_name='id')

rule create_af3_structures_with_tcrmodel2_data_no_unpaired_msas_no_paired_msas_no_templates_sample_size_1500:
    input: get_predicted_structure_names("data/interim/sequences_sample_1500.csv", "data/interim/structures_af3_af3_tcrmodel2_no_unpaired_msas_no_paired_msas_no_templates/{name}_seed-{seed}_sample-{sample}_model.cif", config['seed'])

rule create_af3_structures_with_tcrmodel2_data_no_unpaired_msas_no_paired_msas_no_templates_benchmark:
    input: get_predicted_structure_names("data/interim/benchmark_structures.csv", "data/interim/structures_af3_af3_tcrmodel2_no_unpaired_msas_no_paired_msas_no_templates/{name}_seed-{seed}_sample-{sample}_model.cif", config['seed'], column_name='id')

rule create_af3_structures_with_tcrmodel2_data_dash_bmlf:
    input: get_dash_names("data/interim/sequences_dash_with_negatives.csv", "data/interim/structures_af3_af3_tcrmodel2/{name}_seed-{seed}_sample-{sample}_model.cif", "GLCTLVAML", config['seed'])

rule create_af3_structures_with_tcrmodel2_data_mutated_dash_bmlf:
    input: get_predicted_structure_names("data/interim/sequences_dash_BMLF_mutated.csv", "data/interim/structures_af3_af3_tcrmodel2/{name}_seed-{seed}_sample-{sample}_model.cif", config['seed'], column_name='id')

rule create_af3_structures_with_tcrmodel2_data_coeliac:
    input: get_predicted_structure_names("data/interim/sequences_coeliac_with_negatives.csv", "data/interim/structures_af3_af3_tcrmodel2/{name}_seed-{seed}_sample-{sample}_model.cif", config['seed'], column_name='id')

rule create_af3_structures_with_tcrmodel2_data_coeliac_mutations:
    input: get_predicted_structure_names("data/interim/sequences_coeliac_mutations.csv", "data/interim/structures_af3_af3_tcrmodel2/{name}_seed-{seed}_sample-{sample}_model.cif", config['seed'], column_name='id')

rule create_af3_structures_with_tcrmodel2_data_c259_affinity_mutations:
    input: get_predicted_structure_names("data/interim/sequences_c259_mutational_scan_affinity.csv", "data/interim/structures_af3_af3_tcrmodel2/{name}_seed-{seed}_sample-{sample}_model.cif", config['seed'], column_name='id')

rule compute_ipsae_on_af3_structure:
    input:
        confidences="data/interim/structures_af3_{method}_{data}/{name}_seed-{seed}_sample-{sample}_confidences.json",
        summary_confidences="data/interim/structures_af3_{method}_{data}/{name}_seed-{seed}_sample-{sample}_summary_confidences.json",
        structure="data/interim/structures_af3_{method}_{data}/{name}_seed-{seed}_sample-{sample}_model.cif"
    output:
        "data/interim/ipsaes_af3_{method}_{data}/{name}_seed-{seed}_sample-{sample}_model_{distance_cutoff}_{pae_cutoff}_byres.txt",
        "data/interim/ipsaes_af3_{method}_{data}/{name}_seed-{seed}_sample-{sample}_model_{distance_cutoff}_{pae_cutoff}.pml",
        "data/interim/ipsaes_af3_{method}_{data}/{name}_seed-{seed}_sample-{sample}_model_{distance_cutoff}_{pae_cutoff}.txt"
    resources:
        runtime="10m",
        mem="2GB"
    shell:
        """
        output_dir=$(dirname {output[0]})
        mkdir -p $output_dir

        cp {input.confidences} $TMPDIR/
        cp {input.summary_confidences} $TMPDIR/
        cp {input.structure} $TMPDIR/

        python third_party/IPSAE/ipsae.py \
            $TMPDIR/$(basename {input.confidences}) \
            $TMPDIR/$(basename {input.structure}) \
            {wildcards.distance_cutoff} \
            {wildcards.pae_cutoff}

        mv $TMPDIR/$(basename {output[0]}) $output_dir/
        mv $TMPDIR/$(basename {output[1]}) $output_dir/
        mv $TMPDIR/$(basename {output[2]}) $output_dir/
        """

def _get_ipsaes_names(path: str, seeds: list[int], distance_cutoff: float, pae_cutoff: float) -> list[str]:
    if not os.path.exists(path):
        print(path, 'does not exist!')
        return []

    template = "data/interim/ipsaes_af3_{{{{method}}}}_{{{{data}}}}/{{name}}_seed-{{seed}}_sample-{{sample}}_model_{distance_cutoff:02}_{pae_cutoff:02}.txt"

    with open(path, 'r') as fh:
        indices = [entry['index'] for entry in csv.DictReader(fh)]

    return expand(template.format(distance_cutoff=distance_cutoff, pae_cutoff=pae_cutoff), name=indices, seed=seeds, sample=range(5))

rule add_ipsaes_to_confidence_scores:
    input:
        confidence_scores="data/interim/confidence_scores_af3_{method}_{data}_{sample_size}{unbalance}.csv",
        ipaes=lambda wildcards: _get_ipsaes_names(f"data/interim/sequences_sample_{wildcards.sample_size}{wildcards.unbalance}.csv", seeds=config['seed'], distance_cutoff=10, pae_cutoff=10)
    output: "data/processed/confidence_scores_af3_{method}_{data}_{sample_size}{unbalance}.csv"
    resources:
        runtime="10m",
        mem="2GB"
    run:
        import pandas as pd

        chain_map = {'A': 'mhc_chain1', 'B': 'mhc_chain2', 'C': 'antigen_chain', 'D': 'alpha_chain', 'E': 'beta_chain'}

        ipsaes = []
        for path in input.ipaes:
            name, seed, sample, _, _, _ = os.path.basename(path).replace('.txt', '').rsplit('_', 5)
            ipsae = pd.read_csv(path, sep='\s+')

            ipsae = ipsae[ipsae['Type'] != 'max'].copy()

            ipsae['chain_type1'] = ipsae['Chn1'].map(chain_map)
            ipsae['chain_type2'] = ipsae['Chn2'].map(chain_map)
            index = ipsae['chain_type1'] + '_' + ipsae['chain_type2'] + '_ipsae'
            ipsae = ipsae.set_index(index)

            ipsae = ipsae[['ipSAE']].T

            ipsae['name'] = name
            ipsae['seed'] = int(seed.replace('seed-', ''))
            ipsae['sample'] = int(sample.replace('sample-', ''))

            ipsaes.append(ipsae)

        ipsaes = pd.concat(ipsaes).reset_index(drop=True)

        ipsaes['combined_ipsaes'] = ipsaes.filter(
            regex=(
                '((antigen_chain|mhc_chain_1)_(alpha_chain|beta_chain))'
                '|((alpha_chain|beta_chain)_(antigen_chain|mhc_chain_1))_ipsae'
            )
        ).apply('mean', axis=1)

        confidence_scores = pd.read_csv(input.confidence_scores)
        confidence_scores['name'] = confidence_scores['name'].astype(str)
        confidence_scores = confidence_scores.merge(ipsaes, how='left', on=['name', 'seed', 'sample'])
        confidence_scores.to_csv(output[0], index=False)

rule process_colabfold_structures_sample_size_50:
    input: get_predicted_structure_names("data/interim/sequences_sample_50.csv", "data/processed/structures_af3_mmseqs2_colabfold/{name}_{seed}_{sample}.pdb", config['seed'])

rule process_af3_structures_sample_size_50:
    input: get_predicted_structure_names("data/interim/sequences_sample_50.csv", "data/processed/structures_af3_af3_af3/{name}_{seed}_{sample}.pdb", config['seed'])

rule process_af3_structures_with_tcrmodel2_data_sample_size_50:
    input: get_predicted_structure_names("data/interim/sequences_sample_50.csv", "data/processed/structures_af3_af3_tcrmodel2/{name}_{seed}_{sample}.pdb", config['seed'])

rule process_af3_structures_with_tcrmodel2_data_sample_size_1500:
    input: get_predicted_structure_names("data/interim/sequences_sample_1500.csv", "data/processed/structures_af3_af3_tcrmodel2/{name}_{seed}_{sample}.pdb", config['seed'])

rule process_af3_structures_with_tcrmodel2_data_dash_bmlf:
    input: get_dash_names("data/interim/sequences_dash_with_negatives.csv", "data/processed/structures_af3_af3_tcrmodel2/{name}_{seed}_{sample}.pdb", "GLCTLVAML", config['seed'])

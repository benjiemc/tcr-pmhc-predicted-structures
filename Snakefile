import os

envvars:
    "ALPHAFOLD_DATA_DIR",
    "TCRMODEL2_ALPHAFOLD_DATA_DIR",
    "ALPHAFOLD3_WEIGHTS_DIR"

rule lint:
    shell:
        """
        ruff check
        ruff format --check
        """

rule test:
    shell: "pytest --cov=src tests"

rule run_alphafold:
    input:
        json="data/external/af_input"
    output:
        directory("data/interim/af_output")
    params:
        alphafold_data_dir=os.environ['ALPHAFOLD_DATA_DIR'],
        alphafold3_weights_dir=os.environ['ALPHAFOLD3_WEIGHTS_DIR']
    shell:
        """
        mkdir -p {output}

        input_abs="$(pwd)/{input}"
        output_abs="$(pwd)/{output}"

        docker run -it \
            --volume $output_abs:/tmp/af_output \
            --entrypoint bash \
            alphafold3 \
            -c "chown -R $(id -u):$(id -g) /tmp/af_output"


        docker run -it \
            --user $(id -u):$(id -g) \
            --volume $input_abs:/tmp/af_input \
            --volume $output_abs:/tmp/af_output \
            --volume {params.alphafold3_weights_dir}:/models \
            --volume {params.alphafold_data_dir}:/public_databases \
            --gpus all \
            alphafold3 \
            python run_alphafold.py \
                --json_path=/tmp/af_input/fold_input.json \
                --model_dir=/models \
                --output_dir=/tmp/af_output
        """


rule predict_tcr_structures_using_tcrmodel2:
    input: "data/interim/sampled_vdjdb_norm.csv"
    output: directory("data/interim/tcrmodel2_output/tcrs")
    script: "scripts/predict_tcr_structures.py"

rule predict_tcr_pmhc_I_structures_using_tcrmodel2:
    input: "data/interim/sampled_vdjdb_norm.csv"
    output: directory("data/interim/tcrmodel2_output/tcr-pmhc-Is")
    script: "scripts/predict_tcr_mhcI_structures.py"

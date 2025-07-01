envvars:
    "TCRMODEL2_ALPHAFOLD_DATA_DIR"

rule predict_tcr_structures_using_tcrmodel2:
    input: "data/interim/sampled_vdjdb_norm.csv"
    output: directory("data/interim/tcrmodel2_output/tcrs")
    script: "../scripts/predict_tcr_structures.py"

rule predict_tcr_pmhc_I_structures_using_tcrmodel2:
    input: "data/interim/sampled_vdjdb_norm.csv"
    output: directory("data/interim/tcrmodel2_output/tcr-pmhc-Is")
    script: "../scripts/predict_tcr_mhcI_structures.py"

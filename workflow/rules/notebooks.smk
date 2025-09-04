rule run_compare_structural_features_to_confidence:
    input:
        "data/external/tcr_pmhc_contacts.csv",
        "data/interim/sequences_sample_50.csv",
        "data/processed/confidence_scores_af3_af3_50.csv",
        "data/processed/predicted_af3_af3_structural_features_50.csv",
        "data/processed/predicted_af3_af3_contact_maps_50.csv",
        "data/processed/predicted_af3_af3_interactions_50.csv",
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
        "data/processed/confidence_scores_af3_af3_50.csv",
        "data/interim/sequences_sample_50.csv"
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
import hashlib
import json
import os


def get_sequence_ids(sequence_path: str, template: str) -> list[str]:
    if not os.path.exists(sequence_path):
        print(sequence_path, 'does not exist!')
        return []

    with open(sequence_path, 'r') as fh:
        sequences = json.load(fh)

    ids = [
        hashlib.sha256(sequence.encode()).hexdigest()[:8]
        for sequence_type, sequence in sequences.items()
        if sequence_type in ('tcr_alpha_chain', 'tcr_beta_chain', 'mhc_chain_1', 'mhc_chain_2')
    ]

    return expand(template, name=ids)


def get_predicted_structure_names(path: str, template: str, seeds: list, column_name: str = 'index') -> list[str]:
    if not os.path.exists(path):
        print(path, 'does not exist!')
        return []

    with open(path, 'r') as fh:
        indices = [entry[column_name] for entry in csv.DictReader(fh)]

    return expand(template, name=indices, seed=seeds, sample=list(range(5)))

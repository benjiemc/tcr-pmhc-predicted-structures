#!/usr/bin/env python
"""Extract confidence metrics from AF2 predictions."""

import json

import numpy as np
from snakemake.script import snakemake

from tcr_pmhc_predicted_structures.confidence_scores import (
    compute_alpha_bridge_scores,
    extract_confidence_score,
    get_sequence_information,
)

with open(snakemake.input.sequences, 'r') as fh:
    sequences = json.load(fh)

sequences_flat, chain_ids, sequence_ids, insert_codes = get_sequence_information(sequences)

scores = {}

with open(snakemake.input.summary_confidences, 'r') as fh:
    summary_confidences = json.load(fh)

scores = {**scores, **summary_confidences}

with open(snakemake.input.confidences, 'r') as fh:
    confidences = json.load(fh)

pae = np.array(confidences['pae'])
scores['tcr_pmhc_pae_score'] = np.mean(
    extract_confidence_score(
        sequences_flat,
        chain_ids,
        sequence_ids,
        insert_codes,
        pae,
        mhc_type=sequences['mhc_type'],
        score_type='pae',
        portion='all',
    )
)
scores['cdr_pmhc_pae_score'] = np.mean(
    extract_confidence_score(
        sequences_flat,
        chain_ids,
        sequence_ids,
        insert_codes,
        pae,
        mhc_type=sequences['mhc_type'],
        score_type='pae',
        portion='binding',
    )
)

plddt = np.array(confidences['plddt'])
scores['tcr_pmhc_plddt_score'] = np.mean(
    extract_confidence_score(
        sequences_flat,
        chain_ids,
        sequence_ids,
        insert_codes,
        plddt,
        mhc_type=sequences['mhc_type'],
        score_type='plddt',
        portion='all',
    )
)
scores['cdr_pmhc_plddt_score'] = np.mean(
    extract_confidence_score(
        sequences_flat,
        chain_ids,
        sequence_ids,
        insert_codes,
        plddt,
        mhc_type=sequences['mhc_type'],
        score_type='plddt',
        portion='binding',
    )
)

scores['tcr_pmhc_alphabridge_score'] = np.mean(
    compute_alpha_bridge_scores(
        sequences_flat,
        chain_ids,
        sequence_ids,
        insert_codes,
        plddt,
        pae,
        mhc_type=sequences['mhc_type'],
        atom_level_plddt=False,
        portion='all',
    )
)
scores['cdr_pmhc_alphabridge_score'] = np.mean(
    compute_alpha_bridge_scores(
        sequences_flat,
        chain_ids,
        sequence_ids,
        insert_codes,
        plddt,
        pae,
        mhc_type=sequences['mhc_type'],
        atom_level_plddt=False,
        portion='binding',
    )
)

with open(snakemake.output[0], 'w') as fh:
    fh.write(','.join(scores.keys()))
    fh.write('\n')
    fh.write(','.join(str(value) for value in scores.values()))
    fh.write('\n')

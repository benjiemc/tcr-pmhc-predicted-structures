#!/usr/bin/env python
"""Extract confidence metrics from Chai predictions."""

import json
import re

import numpy as np
from snakemake.script import snakemake

from tcr_pmhc_predicted_structures.confidence_scores import (
    compute_alpha_bridge_scores,
    extract_confidence_score,
    get_sequence_information,
)
from tcr_pmhc_predicted_structures.predictions import CHAIN_ORDER


def collect_summary_scores(summary_confidences: dict) -> dict:
    """Collect all values from summary score dictionary."""
    scores = {key: val for key, val in summary_confidences.items() if not isinstance(val, list)}

    for chain_name, chain_ptm in zip(CHAIN_ORDER, summary_confidences['per_chain_ptm'], strict=False):
        scores[chain_name + '_ptm'] = chain_ptm

    for i, (chain_name_1, chain_iptms) in enumerate(
        zip(CHAIN_ORDER, summary_confidences['per_chain_pair_iptm'], strict=False)
    ):
        for chain_name_2, iptm in zip(CHAIN_ORDER[i + 1 :], chain_iptms[i + 1 :], strict=False):
            if chain_name_1 != chain_name_2:
                scores[chain_name_1 + '_' + chain_name_2 + '_iptm'] = iptm

    for chain_name_1, chain_clashes in zip(CHAIN_ORDER, summary_confidences['chain_chain_clashes'], strict=False):
        for chain_name_2, clash in zip(CHAIN_ORDER, chain_clashes, strict=False):
            scores[chain_name_1 + '_' + chain_name_2 + '_clash'] = clash

    return scores


with open(snakemake.input.sequences, 'r') as fh:
    sequences = json.load(fh)

sequences_flat, chain_ids, sequence_ids, insert_codes = get_sequence_information(sequences)

scores = {}

with open(snakemake.input.summary_confidences, 'r') as fh:
    summary_confidences = json.load(fh)

summary_confidences = {key: val[0] for key, val in summary_confidences.items()}

scores = {**scores, **collect_summary_scores(summary_confidences)}

if sequences['mhc_type'] == 'MH1':
    scores['combined_iptm'] = np.mean(
        [score for key, score in scores.items() if re.match('(peptide|mhc_chain_1)_tcr_(alpha|beta)_chain_iptm', key)]
    )

else:
    scores['combined_iptm'] = np.mean(
        [
            score
            for key, score in scores.items()
            if re.match('(peptide|mhc_chain_1|mhc_chain_2)_tcr_(alpha|beta)_chain_iptm', key)
        ]
    )

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

pde = np.array(confidences['pde'])
scores['tcr_pmhc_pde_score'] = np.mean(
    extract_confidence_score(
        sequences_flat,
        chain_ids,
        sequence_ids,
        insert_codes,
        pde,
        mhc_type=sequences['mhc_type'],
        score_type='pde',
        portion='all',
    )
)
scores['cdr_pmhc_pde_score'] = np.mean(
    extract_confidence_score(
        sequences_flat,
        chain_ids,
        sequence_ids,
        insert_codes,
        pde,
        mhc_type=sequences['mhc_type'],
        score_type='pde',
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

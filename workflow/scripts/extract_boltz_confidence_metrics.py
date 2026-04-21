#!/usr/bin/env python
"""Extract confidence metrics from Boltz predictions."""

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
    scores = {key: val for key, val in summary_confidences.items() if not isinstance(val, dict)}

    chain_map = {
        chain_id: chain_name
        for chain_name, chain_id in zip(
            CHAIN_ORDER, sorted(summary_confidences['chains_ptm'].keys(), key=int), strict=False
        )
    }

    for chain_id, chain_ptm in summary_confidences['chains_ptm'].items():
        scores[chain_map[chain_id] + '_ptm'] = chain_ptm

    for chain_id_1, chain_iptms in summary_confidences['pair_chains_iptm'].items():
        for chain_id_2, iptm in chain_iptms.items():
            if chain_id_1 != chain_id_2:
                scores[chain_map[chain_id_1] + '_' + chain_map[chain_id_2] + '_iptm'] = iptm

    return scores


with open(snakemake.input.sequences, 'r') as fh:
    sequences = json.load(fh)

sequences_flat, chain_ids, sequence_ids, insert_codes = get_sequence_information(sequences)

scores = {}

with open(snakemake.input.summary_confidences, 'r') as fh:
    summary_confidences = json.load(fh)

scores = {**scores, **collect_summary_scores(summary_confidences)}

if sequences['mhc_type'] == 'MH1':
    scores['combined_iptm'] = np.mean(
        [
            score
            for key, score in scores.items()
            if re.match(
                '((peptide|mhc_chain_1)_tcr_(alpha|beta)_chain)|(tcr_(alpha|beta)_chain_(peptide|mhc_chain_1))_iptm',
                key,
            )
        ]
    )

else:
    scores['combined_iptm'] = np.mean(
        [
            score
            for key, score in scores.items()
            if re.match(
                '((peptide|mhc_chain_1|mhc_chain_2)_tcr_(alpha|beta)_chain)'
                '|(tcr_(alpha|beta)_chain_(peptide|mhc_chain_1|mhc_chain_2))_iptm',
                key,
            )
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

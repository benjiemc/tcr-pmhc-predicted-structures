"""Extract confidence metrics from alphafold3 structures."""

import argparse
import json
import logging
import os
import re
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
from Bio.SeqUtils import IUPACData

from tcr_pmhc_predicted_structures.apps._log import add_logging_arguments, setup_logger
from tcr_pmhc_predicted_structures.utils import (
    IMGT_CDR,
    IMGT_MHC_TCR_CONTACT_POSITIONS,
    MHC_I_IMGT_BETA_HELIX_START,
)

CHAIN_ORDER = ('mhc_chain_1', 'mhc_chain_2', 'antigen_chain', 'alpha_chain', 'beta_chain')

logger = logging.getLogger()

try:
    import anarci

except ImportError:
    logger.warning('Anarci is unavailable and is required by some functions')

parser = argparse.ArgumentParser(
    prog=f'python -m {sys.modules[__name__].__spec__.name}',
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)

parser.add_argument('--summary-confidence-paths', required=True, nargs='+', help='paths to the af3 summary confidences')
parser.add_argument('--confidence-paths', required=True, nargs='+', help='paths to the af3 confidences')
parser.add_argument('--sequence-paths', required=True, nargs='+', help='paths to the af3 input sequence jsons')
parser.add_argument('--output', '-o', help='paths to the output csv')

add_logging_arguments(parser)


RESIDUE_NUMBER_OF_ATOMS = {
    'ALA': 5,
    'ARG': 11,
    'ASN': 8,
    'ASP': 8,
    'CYS': 6,
    'GLN': 9,
    'GLU': 9,
    'GLY': 4,
    'HIS': 10,
    'ILE': 8,
    'LEU': 8,
    'LYS': 9,
    'MET': 8,
    'PHE': 11,
    'PRO': 7,
    'SER': 6,
    'THR': 7,
    'TRP': 14,
    'TYR': 12,
    'VAL': 7,
}
'''Number of atoms in each amino acid residue.'''


def collect_summary_scores(paths: list[str]) -> pd.DataFrame:
    af3_scores = []

    for path in paths:
        with open(path) as fh:
            summary_confidences = json.load(fh)

        name = re.match(
            r'(?P<structure_name>\d+(?:_\d+)?)_seed-(?P<seed>\d+)_sample-(?P<sample>\d+)_summary_confidences.json',
            os.path.basename(path),
        )

        score = {
            **name.groupdict(),
            **{key: val for key, val in summary_confidences.items() if not isinstance(val, list)},
        }

        for chain_name, chain_ptm in zip(CHAIN_ORDER, summary_confidences['chain_ptm'], strict=True):
            score[chain_name + '_ptm'] = chain_ptm

        for chain_name, chain_ptm in zip(CHAIN_ORDER, summary_confidences['chain_iptm'], strict=True):
            score[chain_name + '_iptm'] = chain_ptm

        for i, (chain_name_1, chain_iptms) in enumerate(
            zip(CHAIN_ORDER, summary_confidences['chain_pair_iptm'], strict=True)
        ):
            for chain_name_2, iptm in zip(CHAIN_ORDER[i + 1 :], chain_iptms[i + 1 :], strict=True):
                if chain_name_1 != chain_name_2:
                    score[chain_name_1 + '_' + chain_name_2 + '_iptm'] = iptm

        for chain_name_1, chain_pae_min in zip(CHAIN_ORDER, summary_confidences['chain_pair_pae_min'], strict=True):
            for chain_name_2, pae_min in zip(CHAIN_ORDER, chain_pae_min, strict=True):
                score[chain_name_1 + '_' + chain_name_2 + '_pae_min'] = pae_min

        af3_scores.append(score)

    return pd.DataFrame(af3_scores)


def collect_pae_scores(confidence: dict) -> np.array:
    pae_numbering = [
        (idx, (chain_id, res_seq_id))
        for idx, (chain_id, res_seq_id) in enumerate(
            zip(confidence['token_chain_ids'], confidence['token_res_ids'], strict=True)
        )
    ]

    tcr_idxs = np.array([idx for idx, (chain_id, _) in pae_numbering if chain_id in ('D', 'E')])
    pmhc_idxs = np.array([idx for idx, (chain_id, _) in pae_numbering if chain_id in ('A', 'C')])

    return np.array(confidence['pae'])[np.ix_(tcr_idxs, pmhc_idxs)]


def number_sequences(sequences: dict[str, str]) -> dict[str, list[tuple[tuple[str, str], str]]]:
    """Numbers sequences with ANARCI."""
    numbered_sequences = {}
    for chain_id, sequence in sequences.items():
        try:
            numbered_sequence, _ = anarci.number(sequence)

        except:  # noqa: E722 Anarci throws base exception on B2M
            numbered_sequence = None

        if numbered_sequence:
            if len(numbered_sequence) == 2:  # noqa: PLR2004 Sometimes anarci will find two valid numberings
                numbered_sequence = numbered_sequence[0]

            numbered_sequences[chain_id] = [(seq_id, olc) for seq_id, olc in numbered_sequence if olc != '-']

        else:
            numbered_sequences[chain_id] = [((res_id, ' '), olc) for res_id, olc in enumerate(sequence, 1)]

    return numbered_sequences


def collect_specific_pae_scores(confidence: dict, input_sequences: dict) -> np.array:
    pae_numbering = [
        (idx, (chain_id, res_seq_id))
        for idx, (chain_id, res_seq_id) in enumerate(
            zip(confidence['token_chain_ids'], confidence['token_res_ids'], strict=True)
        )
    ]

    numbered_sequences = number_sequences(input_sequences)

    # Calculate offset between anarci sequence and input sequences
    offsets = {
        chain_id: sequence.find(''.join([olc for _, olc in numbered_sequences[chain_id]]))
        for chain_id, sequence in input_sequences.items()
    }

    # Get relevant ranges for input sequences
    relative_ranges = {
        chain_id: tuple(
            [
                idx + offsets[chain_id]
                for idx, ((seq_id, _), _) in enumerate(numbered_sequences[chain_id])
                if seq_id in IMGT_CDR
            ]
        )
        for chain_id in ('D', 'E')
    }

    if max([seq_id for (seq_id, _), _ in numbered_sequences['A']]) > MHC_I_IMGT_BETA_HELIX_START:
        tcr_contact_positions = []

        for helix in IMGT_MHC_TCR_CONTACT_POSITIONS:
            offset = MHC_I_IMGT_BETA_HELIX_START if helix == 'beta' else 0

            for res_id in IMGT_MHC_TCR_CONTACT_POSITIONS[helix]:
                seq_id = int(''.join([char for char in res_id if char.isnumeric()])) + offset
                insert_code = (
                    ''.join(
                        [char for char in res_id if not char.isnumeric()],
                    )
                    if len([char for char in res_id if not char.isnumeric()]) > 0
                    else ' '
                )

                tcr_contact_positions.append((seq_id, insert_code))

        tcr_contact_positions = set(tcr_contact_positions)

        relative_ranges['A'] = tuple(
            [
                idx + offsets['A']
                for idx, (res_id, _) in enumerate(numbered_sequences['A'])
                if res_id in tcr_contact_positions
            ]
        )
        relative_ranges['B'] = None

    else:
        tcr_contact_positions = {'A': [], 'B': []}

        for helix in IMGT_MHC_TCR_CONTACT_POSITIONS:
            for res_id in IMGT_MHC_TCR_CONTACT_POSITIONS[helix]:
                seq_id = int(''.join([char for char in res_id if char.isnumeric()]))
                insert_code = (
                    ''.join(
                        [char for char in res_id if not char.isnumeric()],
                    )
                    if len([char for char in res_id if not char.isnumeric()]) > 0
                    else ' '
                )

                tcr_contact_positions['A' if helix == 'alpha' else 'B'].append((seq_id, insert_code))

        relative_ranges['A'] = tuple(
            [
                idx + offsets['A']
                for idx, (res_id, _) in enumerate(numbered_sequences['A'])
                if res_id in tcr_contact_positions['A']
            ]
        )

        relative_ranges['B'] = tuple(
            [
                idx + offsets['B']
                for idx, (res_id, _) in enumerate(numbered_sequences['B'])
                if res_id in tcr_contact_positions['B']
            ]
        )

    relative_ranges['C'] = [idx + offsets['C'] for idx, _ in enumerate(numbered_sequences['C'])]

    # Translate ranges to concatenated PAE space
    current_chain = None
    relative_offset = {}

    for pae_idx, (chain_id, _) in pae_numbering:
        if chain_id != current_chain:
            relative_offset[chain_id] = pae_idx
            current_chain = chain_id

    ranges = {
        chain_id: [idx + relative_offset[chain_id] for idx in relative_range]
        for chain_id, relative_range in relative_ranges.items()
        if relative_range is not None
    }

    # Get PAE values and average
    tcr_idxs = np.array([idx for chain_id, range_ in ranges.items() for idx in range_ if chain_id in ('D', 'E')])
    pmhc_idxs = np.array([idx for chain_id, range_ in ranges.items() for idx in range_ if chain_id not in ('D', 'E')])

    return np.array(confidence['pae'])[np.ix_(tcr_idxs, pmhc_idxs)]


def collect_plddt_scores(confidence: dict) -> float:
    tcr_pmhc_plddts = np.array(
        [
            plddt
            for chain_id, plddt in zip(confidence['atom_chain_ids'], confidence['atom_plddts'], strict=True)
            if chain_id in {'A', 'C', 'D', 'E'}
        ]
    )

    return tcr_pmhc_plddts


def collect_specific_plddt_scores(confidence: dict, input_sequences: dict) -> np.array:
    res_id_grouped_by_chain = defaultdict(list)
    for chain_id, res_id in zip(confidence['token_chain_ids'], confidence['token_res_ids'], strict=True):
        res_id_grouped_by_chain[chain_id].append(res_id)

    chain_order = []
    for chain_id in confidence['token_chain_ids']:
        if chain_id not in chain_order:
            chain_order.append(chain_id)

    atom_res_ids = [
        res_id
        for chain_id in chain_order
        for i, (amino_acid, res_id) in enumerate(
            zip(input_sequences[chain_id], res_id_grouped_by_chain[chain_id], strict=True)
        )
        for _ in range(
            RESIDUE_NUMBER_OF_ATOMS[IUPACData.protein_letters_1to3[amino_acid].upper()]
            if i != len(input_sequences[chain_id]) - 1
            else RESIDUE_NUMBER_OF_ATOMS[IUPACData.protein_letters_1to3[amino_acid].upper()] + 1
        )
    ]

    plddt_numbering = [
        (idx, (chain_id, res_seq_id))
        for idx, (chain_id, res_seq_id) in enumerate(zip(confidence['atom_chain_ids'], atom_res_ids, strict=True))
    ]

    numbered_sequences = number_sequences(input_sequences)

    # Calculate offset between anarci sequence and input sequences
    offsets = {
        chain_id: sequence.find(''.join([olc for _, olc in numbered_sequences[chain_id]]))
        for chain_id, sequence in input_sequences.items()
    }

    # Get relevant ranges for input sequences
    relative_ranges = {
        chain_id: tuple(
            [
                idx + offsets[chain_id]
                for idx, ((seq_id, _), _) in enumerate(numbered_sequences[chain_id])
                if seq_id in IMGT_CDR
            ]
        )
        for chain_id in ('D', 'E')
    }

    if max([seq_id for (seq_id, _), _ in numbered_sequences['A']]) > MHC_I_IMGT_BETA_HELIX_START:
        tcr_contact_positions = []

        for helix in IMGT_MHC_TCR_CONTACT_POSITIONS:
            offset = MHC_I_IMGT_BETA_HELIX_START if helix == 'beta' else 0

            for res_id in IMGT_MHC_TCR_CONTACT_POSITIONS[helix]:
                seq_id = int(''.join([char for char in res_id if char.isnumeric()])) + offset
                insert_code = (
                    ''.join(
                        [char for char in res_id if not char.isnumeric()],
                    )
                    if len([char for char in res_id if not char.isnumeric()]) > 0
                    else ' '
                )

                tcr_contact_positions.append((seq_id, insert_code))

        tcr_contact_positions = set(tcr_contact_positions)

        relative_ranges['A'] = tuple(
            [
                idx + offsets['A']
                for idx, (res_id, _) in enumerate(numbered_sequences['A'])
                if res_id in tcr_contact_positions
            ]
        )
        relative_ranges['B'] = None

    else:
        tcr_contact_positions = {'A': [], 'B': []}

        for helix in IMGT_MHC_TCR_CONTACT_POSITIONS:
            for res_id in IMGT_MHC_TCR_CONTACT_POSITIONS[helix]:
                seq_id = int(''.join([char for char in res_id if char.isnumeric()]))
                insert_code = (
                    ''.join(
                        [char for char in res_id if not char.isnumeric()],
                    )
                    if len([char for char in res_id if not char.isnumeric()]) > 0
                    else ' '
                )

                tcr_contact_positions['A' if helix == 'alpha' else 'B'].append((seq_id, insert_code))

        relative_ranges['A'] = tuple(
            [
                idx + offsets['A']
                for idx, (res_id, _) in enumerate(numbered_sequences['A'])
                if res_id in tcr_contact_positions['A']
            ]
        )

        relative_ranges['B'] = tuple(
            [
                idx + offsets['B']
                for idx, (res_id, _) in enumerate(numbered_sequences['B'])
                if res_id in tcr_contact_positions['B']
            ]
        )

    relative_ranges['C'] = [idx + offsets['C'] for idx, _ in enumerate(numbered_sequences['C'])]

    indices = np.array(
        [
            idx
            for idx, (chain_id, res_id) in plddt_numbering
            if res_id - 1 in (relative_ranges[chain_id] if relative_ranges[chain_id] else [])
        ]
    )

    return (
        np.array(confidence['atom_plddts'])[indices],
        np.array(confidence['atom_chain_ids'])[indices],
        np.array(atom_res_ids)[indices],
    )


def compute_alpha_bridge_scores(pae_scores, plddt_scores, atom_chain_ids, atom_res_ids) -> pd.DataFrame:
    idx = 0
    plddt_indices = []

    ids = list(zip(atom_chain_ids, atom_res_ids, strict=True))
    previous_id = ids[0]

    for current_id in ids:
        if current_id != previous_id:
            idx += 1
            previous_id = current_id

        plddt_indices.append(idx)

    plddt_indices = np.array(plddt_indices)

    (tcr_plddt_indices,) = np.where((atom_chain_ids == 'D') | (atom_chain_ids == 'E'))
    (pmhc_plddt_indices,) = np.where((atom_chain_ids != 'D') & (atom_chain_ids != 'E'))

    tcr_indices = np.unique(plddt_indices[tcr_plddt_indices])
    pmhc_indices = np.unique(plddt_indices[pmhc_plddt_indices])

    plddt_agg = np.zeros((len(tcr_indices), len(pmhc_indices)), dtype=np.float64)
    for i, tcr_idx in enumerate(tcr_indices):
        for j, pmhc_idx in enumerate(pmhc_indices):
            tcr_part = plddt_scores[np.where(plddt_indices == tcr_idx)[0]]
            pmhc_part = plddt_scores[np.where(plddt_indices == pmhc_idx)[0]]
            plddt_agg[i, j] = 1 / (len(tcr_part) + len(pmhc_part)) * (np.sum(tcr_part) + np.sum(pmhc_part)) - 100

    scores = pae_scores - 1 / 3 * plddt_agg

    return scores


def collect_raw_scores(confidence_paths: list[str], sequence_paths: list[str]) -> pd.DataFrame:
    raw_scores = []
    for path in confidence_paths:
        logger.debug('Getting scores from %s', path)

        name = re.match(
            r'(?P<structure_name>\d+(?:_\d+)?)_seed-(?P<seed>\d+)_sample-(?P<sample>\d+)_confidences.json',
            os.path.basename(path),
        )

        scores = {**name.groupdict()}

        with open(path, 'r') as fh:
            confidence = json.load(fh)

        with open(
            next(
                sequence_path
                for sequence_path in sequence_paths
                if name['structure_name'] == os.path.basename(sequence_path).split('.')[0].split('_')[0]
            ),
            'r',
        ) as fh:
            input_data = json.load(fh)

        input_sequences = {
            sequence_info['protein']['id']: sequence_info['protein']['sequence']
            for sequence_info in input_data['sequences']
        }

        logger.debug('Collecting PAE scores')
        pae_scores = collect_pae_scores(confidence)
        scores['tcr_pmhc_pae_score'] = np.mean(pae_scores)

        logger.debug('Collecting specific PAE scores')
        specific_pae_scores = collect_specific_pae_scores(confidence, input_sequences)
        scores['cdr_pmhc_pae_score'] = np.mean(specific_pae_scores)

        logger.debug('Collecting pLDDT scores')
        plddt_scores = collect_plddt_scores(confidence)
        scores['tcr_pmhc_plddt_score'] = np.mean(plddt_scores)

        logger.debug('Collecting specific pLDDT scores')
        specifc_plddt_scores, atom_chain_ids, atom_res_ids = collect_specific_plddt_scores(confidence, input_sequences)
        scores['cdr_pmhc_plddt_score'] = np.mean(specifc_plddt_scores)

        logger.debug('Computing AlphaBridge scores')
        alpha_bridge_scores = compute_alpha_bridge_scores(
            specific_pae_scores, specifc_plddt_scores, atom_chain_ids, atom_res_ids
        )
        scores['alphabridge_score'] = np.mean(alpha_bridge_scores)

        raw_scores.append(scores)

    return pd.DataFrame(raw_scores)


def main() -> None:
    args = parser.parse_args()
    setup_logger(logger, args.log_level)

    logger.info('Collecting summary scores')
    af3_scores = collect_summary_scores(args.summary_confidence_paths)

    logger.debug('Identifying top model')
    af3_scores['top_model'] = False
    af3_scores.loc[af3_scores.groupby('structure_name')['ranking_score'].idxmax(), 'top_model'] = True

    logger.info('Creating composite summary scores')

    logger.debug('Creating composite PAE-min scores')
    af3_scores['combined_pae_min'] = af3_scores.filter(
        regex=(
            '((antigen_chain|mhc_chain_1)_(alpha_chain|beta_chain))'
            '|((alpha_chain|beta_chain)_(antigen_chain|mhc_chain_1))_pae_min'
        )
    ).apply('mean', axis=1)

    logger.debug('Creating composite iPTM scores')
    af3_scores['combined_iptm'] = af3_scores.filter(regex='(antigen_chain|mhc_chain_1)_(alpha|beta)_chain_iptm').apply(
        'mean', axis=1
    )

    logger.info('Collecting raw confidence scores')
    raw_scores = collect_raw_scores(args.confidence_paths, args.sequence_paths)
    af3_scores = af3_scores.merge(raw_scores, how='inner', on=['structure_name', 'seed', 'sample'])

    logger.info('Outputting to %s', args.output)
    af3_scores.to_csv(args.output, index=False)


if __name__ == '__main__':
    main()

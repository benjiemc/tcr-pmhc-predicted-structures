"""Functions for working with confidence scores from structure predictions."""

import logging
from collections.abc import Callable
from functools import cache

import numpy as np
from Bio.SeqUtils import IUPACData

from tcr_pmhc_predicted_structures.imgt_numbering import (
    IMGT_CDR,
    IMGT_CDR1,
    IMGT_CDR2,
    IMGT_CDR3,
    IMGT_MH1_ABD,
    IMGT_MH2_ABD,
    IMGT_MHC_TCR_CONTACT_POSITIONS,
    IMGT_VARIABLE_DOMAIN,
    MHC_I_IMGT_BETA_HELIX_START,
    number,
)
from tcr_pmhc_predicted_structures.padding import centre_pad, left_pad, right_pad
from tcr_pmhc_predicted_structures.predictions import CHAIN_ORDER, PDB_CHAIN_IDS, RESIDUE_NUMBER_OF_ATOMS

PADDING_FUNCTIONS = {'centre': centre_pad, 'left': left_pad, 'right': right_pad}

logger = logging.getLogger(__name__)


def get_sequence_information(sequences: dict[str, str]) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Create arrays from sequences dictionary with information about the sequences, chain IDs, and IMGT numbering."""
    sequences_flat = []
    chain_ids = []
    sequence_ids = []
    insert_codes = []

    for chain_id, chain_name in zip(PDB_CHAIN_IDS, CHAIN_ORDER, strict=False):
        if chain_name not in sequences:
            continue

        for amino_acid in sequences[chain_name]:
            sequences_flat.append(amino_acid)
            chain_ids.append(chain_id)

        if chain_name == 'peptide':
            sequence_ids.extend(list(range(1, len(sequences[chain_name]) + 1)))
            insert_codes.extend([''] * len(sequences[chain_name]))

        else:
            numbering = number(sequences[chain_name])
            sequence_ids.extend([seq_id for (seq_id, _), _ in numbering])
            insert_codes.extend([insert_code if insert_code != ' ' else '' for (_, insert_code), _ in numbering])

    return np.array(sequences_flat), np.array(chain_ids), np.array(sequence_ids), np.array(insert_codes)


@cache
def _get_mhc_tcr_contact_positions(
    mhc_type: str,
) -> tuple[tuple[np.ndarray, np.ndarray], tuple[np.ndarray | None, np.ndarray | None]]:
    match mhc_type:
        case 'MH1':
            seq_ids = []
            insert_codes = []
            for helix in IMGT_MHC_TCR_CONTACT_POSITIONS:
                offset = MHC_I_IMGT_BETA_HELIX_START if helix == 'beta' else 0

                for res_id in IMGT_MHC_TCR_CONTACT_POSITIONS[helix]:
                    seq_ids.append(int(''.join([char for char in res_id if char.isnumeric()])) + offset)
                    insert_codes.append(''.join([char for char in res_id if not char.isnumeric()]))

            return (np.array(seq_ids), np.array(insert_codes)), (None, None)

        case 'MH2':
            contact_positions = {
                'alpha': {'seq_ids': [], 'insert_codes': []},
                'beta': {'seq_ids': [], 'insert_codes': []},
            }
            for helix in IMGT_MHC_TCR_CONTACT_POSITIONS:
                for res_id in IMGT_MHC_TCR_CONTACT_POSITIONS[helix]:
                    contact_positions[helix]['seq_ids'].append(
                        int(''.join([char for char in res_id if char.isnumeric()]))
                    )
                    contact_positions[helix]['insert_codes'].append(
                        ''.join([char for char in res_id if not char.isnumeric()])
                    )

            return (
                (np.array(contact_positions['alpha']['seq_ids']), np.array(contact_positions['alpha']['insert_codes'])),
                (np.array(contact_positions['beta']['seq_ids']), np.array(contact_positions['beta']['insert_codes'])),
            )

        case _:
            msg = f"Invalid MHC type: {mhc_type}. Please enter 'MH1' or 'MH2'"
            raise ValueError(msg)


def _get_extension(sequence: np.ndarray, chain_ids: np.ndarray) -> np.ndarray:
    extension = np.array(
        [RESIDUE_NUMBER_OF_ATOMS[IUPACData.protein_letters_1to3[amino_acid].upper()] for amino_acid in sequence]
    )
    end_of_chains = []
    prev_chain_id = chain_ids[0]

    for i, chain_id in enumerate(chain_ids[1:], 1):
        if chain_id != prev_chain_id:
            end_of_chains.append(i - 1)
            prev_chain_id = chain_id

    end_of_chains.append(len(chain_ids) - 1)

    end_of_chains = np.array(end_of_chains)
    extension[end_of_chains] += 1

    return extension


def _extend_array(array: np.ndarray, extensions: np.ndarray) -> np.ndarray:
    new_array = []
    for element, repeats in zip(array, extensions, strict=True):
        new_array.extend([element] * repeats)

    return np.array(new_array)


def _aggregate_from_extension(array: np.ndarray, extensions: np.ndarray) -> np.ndarray:
    new_array = []
    idx = 0
    for repeats in extensions:
        new_array.append(np.mean(array[idx : idx + repeats]))
        idx += repeats

    return np.array(new_array)


def _create_masks(
    chain_ids: np.ndarray,
    sequence_ids: np.ndarray,
    insert_codes: np.ndarray,
    *,
    mhc_type: str,
    portion: bool = True,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    tcr_mask = np.isin(chain_ids, np.array(['D', 'E']))
    peptide_mask = chain_ids == 'C'
    mhc_1_mask = chain_ids == 'A'

    if portion == 'cropped':
        tcr_mask = tcr_mask & np.isin(sequence_ids, np.array(sorted(IMGT_VARIABLE_DOMAIN)))
        mhc_1_mask = mhc_1_mask & np.isin(
            sequence_ids, np.array(sorted(IMGT_MH1_ABD if mhc_type == 'MH1' else IMGT_MH2_ABD))
        )

    elif portion == 'binding':
        tcr_mask = tcr_mask & np.isin(sequence_ids, np.array(sorted(IMGT_CDR)))
        ((mhc_alpha_seq_ids, mhc_alpha_insert_codes), (mhc_beta_seq_ids, mhc_beta_insert_codes)) = (
            _get_mhc_tcr_contact_positions(mhc_type)
        )

        combined_ids = np.array(
            [f'{seq_id}{insert_code}' for seq_id, insert_code in zip(sequence_ids, insert_codes, strict=True)]
        )
        mhc_alpha = np.array(
            [
                f'{seq_id}{insert_code}'
                for seq_id, insert_code in zip(mhc_alpha_seq_ids, mhc_alpha_insert_codes, strict=True)
            ]
        )

        mhc_1_mask = mhc_1_mask & np.isin(combined_ids, mhc_alpha)

    if mhc_type == 'MH2':
        mhc_2_mask = chain_ids == 'B'

        if portion == 'cropped':
            mhc_2_mask = mhc_2_mask & np.isin(sequence_ids, np.array(sorted(IMGT_MH2_ABD)))

        elif portion == 'binding':
            mhc_beta = np.array(
                [
                    f'{seq_id}{insert_code}'
                    for seq_id, insert_code in zip(mhc_beta_seq_ids, mhc_beta_insert_codes, strict=True)
                ]
            )
            mhc_2_mask = mhc_2_mask & np.isin(combined_ids, mhc_beta)
    else:
        mhc_2_mask = np.array([False for _ in sequence_ids])

    mhc_mask = mhc_1_mask | mhc_2_mask

    return tcr_mask, peptide_mask, mhc_mask


def _add_padding(
    confidence: np.ndarray,
    indices: np.ndarray,
    pad_length: int,
    padding_function: Callable,
    adjustment: int = 0,
    axis: int | None = None,
) -> tuple[np.ndarray, int]:
    start = np.min(indices)
    end = np.max(indices)

    action, left_num, right_num = padding_function(len(indices), pad_length)

    match action:
        case 'pad':
            for _ in range(right_num):
                confidence = np.insert(confidence, end + adjustment + 1, np.nan, axis=axis)

            for _ in range(left_num):
                confidence = np.insert(confidence, start + adjustment, np.nan, axis=axis)

            adjustment += right_num + left_num

        case 'crop':
            for _ in range(right_num):
                confidence = np.delete(confidence, end + adjustment, axis=axis)

            for _ in range(left_num):
                confidence = np.delete(confidence, start + adjustment, axis=axis)

            adjustment -= right_num + left_num

    return confidence, adjustment


def pad_2d_confidence_scores(
    confidence: np.ndarray,
    row_chain_ids: np.ndarray,
    row_sequence_ids: np.ndarray,
    column_chain_ids: np.ndarray,
    padding: dict,
) -> np.ndarray:
    """Pad 2D confidence scores (PAE, PDE, etc) following the config provided."""
    padding_function = PADDING_FUNCTIONS[padding['direction']]

    peptide_indices = np.where(column_chain_ids == 'C')[0]
    confidence, _ = _add_padding(confidence, peptide_indices, padding['peptide_length'], padding_function, axis=1)

    adjustment = 0
    for chain_id, imgt_numbers, pad_length in (
        ('D', IMGT_CDR1, padding['cdr1_alpha_length']),
        ('D', IMGT_CDR2, padding['cdr2_alpha_length']),
        ('D', IMGT_CDR3, padding['cdr3_alpha_length']),
        ('E', IMGT_CDR1, padding['cdr1_beta_length']),
        ('E', IMGT_CDR2, padding['cdr2_beta_length']),
        ('E', IMGT_CDR3, padding['cdr3_beta_length']),
    ):
        cdr_indices = np.where((row_chain_ids == chain_id) & (np.isin(row_sequence_ids, list(imgt_numbers))))[0]
        confidence, adjustment = _add_padding(
            confidence, cdr_indices, pad_length, padding_function, adjustment=adjustment, axis=0
        )

    return confidence


def pad_1d_confidence_scores(
    confidence: np.ndarray, chain_ids: np.ndarray, sequence_ids: np.ndarray, padding: dict
) -> np.ndarray:
    """Pad 1D confidence scores (plddt) following the config provided."""
    padding_function = PADDING_FUNCTIONS[padding['direction']]

    peptide_indices = np.where(chain_ids == 'C')[0]
    confidence, adjustment = _add_padding(confidence, peptide_indices, padding['peptide_length'], padding_function)

    for chain_id, imgt_numbers, pad_length in (
        ('D', IMGT_CDR1, padding['cdr1_alpha_length']),
        ('D', IMGT_CDR2, padding['cdr2_alpha_length']),
        ('D', IMGT_CDR3, padding['cdr3_alpha_length']),
        ('E', IMGT_CDR1, padding['cdr1_beta_length']),
        ('E', IMGT_CDR2, padding['cdr2_beta_length']),
        ('E', IMGT_CDR3, padding['cdr3_beta_length']),
    ):
        cdr_indices = np.where((chain_ids == chain_id) & (np.isin(sequence_ids, list(imgt_numbers))))[0]
        confidence, adjustment = _add_padding(
            confidence, cdr_indices, pad_length, padding_function, adjustment=adjustment
        )

    return confidence


def extract_confidence_score(
    sequences: np.ndarray,
    chain_ids: np.ndarray,
    sequence_ids: np.ndarray,
    insert_codes: np.ndarray,
    confidence: np.ndarray,
    *,
    score_type: str,
    mhc_type: str,
    portion: str = 'binding',
    avg_atoms_for_residue: bool = False,
    padding: None | dict = None,
) -> np.ndarray:
    """Extract TCR:pMHC specific confidence score from given confidence."""
    if score_type.startswith('atom'):
        logger.debug('Extending look up arrays from residue level to atom level')
        extension = _get_extension(sequences, chain_ids)

        chain_ids_processed = _extend_array(chain_ids, extension)
        sequence_ids_processed = _extend_array(sequence_ids, extension)
        insert_codes_processed = _extend_array(insert_codes, extension)

    else:
        chain_ids_processed = chain_ids
        sequence_ids_processed = sequence_ids
        insert_codes_processed = insert_codes

    logger.debug('Selecting relevant portion of the matrix')
    logger.debug('Creating masks')
    tcr_mask, peptide_mask, mhc_mask = _create_masks(
        chain_ids_processed,
        sequence_ids_processed,
        insert_codes_processed,
        mhc_type=mhc_type,
        portion=portion,
    )
    pmhc_mask = peptide_mask | mhc_mask

    match score_type:
        case 'pae' | 'pde':
            selected_confidence = confidence[np.ix_(tcr_mask, pmhc_mask)]  # TODO get square version of this

            if padding:
                row_chain_ids = chain_ids[tcr_mask]
                row_sequence_ids = sequence_ids[tcr_mask]
                column_chain_ids = chain_ids[pmhc_mask]

                selected_confidence = pad_2d_confidence_scores(
                    selected_confidence, row_chain_ids, row_sequence_ids, column_chain_ids, padding
                )

        case 'atom_plddt' | 'plddt':
            selected_confidence = confidence[tcr_mask | pmhc_mask]

            if avg_atoms_for_residue:
                tcr_mask, peptide_mask, mhc_mask = _create_masks(
                    chain_ids,
                    sequence_ids,
                    insert_codes,
                    mhc_type=mhc_type,
                    portion=portion,
                )
                pmhc_mask = peptide_mask | mhc_mask

                selected_confidence = _aggregate_from_extension(selected_confidence, extension[tcr_mask | pmhc_mask])

            if padding:
                if score_type.startswith('atom') and not avg_atoms_for_residue:
                    msg = 'Padding is not supported for atom level confidence scores'
                    raise ValueError(msg)

                selected_confidence = pad_1d_confidence_scores(
                    selected_confidence, chain_ids[tcr_mask | pmhc_mask], sequence_ids[tcr_mask | pmhc_mask], padding
                )

    return selected_confidence


def compute_alpha_bridge_scores(
    sequences: np.ndarray,
    chain_ids: np.ndarray,
    sequence_ids: np.ndarray,
    insert_codes: np.ndarray,
    plddt: np.ndarray,
    pae: np.ndarray,
    *,
    mhc_type: str,
    atom_level_plddt: bool = True,
    portion: str = 'all',
) -> np.ndarray:
    """Compute AlphaBridge score.

    AlphaBridge is described in `AlphaBridge: tools for the analysis of predicted macromolecular complexes
    <https://www.biorxiv.org/content/10.1101/2024.10.23.619601v1>`__.

    Args:
        sequences: 1D array of length num_residues of one letter amino acid codes for MHC chain 1, MHC chain 2, peptide,
            TCR alpha chain, and TCR beta chain
        chain_ids: 1D array of length num_residues of chain ids. The MHC chain 1 should be 'A', the MHC chain 2 should
            be 'B', the peptide should be 'C', the TCR alpha chain should be 'D', and the TCR beta chain should be 'E'.
        sequence_ids: 1D array of length num_residues of sequence ids for each residue in IMGT numbering format for the
            TCR:pMHC chains
        insert_codes: 1D array of length num_residues of insert codes for each residue in IMGT numbering format for the
            TCR:pMHC chains. Blank insert codes are represented as empty strings ('').
        plddt: 1D array of length num_atoms if atom_level_plddt is True, or else of length num_residues with the plddt
            scores for each token
        pae: 2D array of shape num_residues x num_residues of the PAE score for each pairwise residue combination
        mhc_type: either 'MH1' for class I TCR:pMHC complexes or 'MH2' for class II TCR:pMHC complexes
        atom_level_plddt: flag to say if the resolution of the plddt array is at atom level (True) or residue level
            (False)
        portion: either 'all', 'cropped', or 'binding' to consider different portions of the structure

    Returns:
        a 2D array of shape num_selected_residues x num_selected_residues with the AlphaBridge score computed for each
        pairwise residue combination

    """
    tcr_mask, peptide_mask, mhc_mask = _create_masks(
        chain_ids, sequence_ids, insert_codes, mhc_type=mhc_type, portion=portion
    )
    pmhc_mask = peptide_mask | mhc_mask

    tcr_indices = np.arange(len(sequences))[tcr_mask]
    pmhc_indices = np.arange(len(sequences))[pmhc_mask]

    plddt_agg = np.zeros((len(tcr_indices), len(pmhc_indices)), dtype=np.float64)

    if atom_level_plddt:
        extension = _get_extension(sequences, chain_ids)

        for i, tcr_idx in enumerate(tcr_indices):
            tcr_atom_start_idx = np.sum(extension[:tcr_idx], dtype=int)
            tcr_atom_end_idx = tcr_atom_start_idx + extension[tcr_idx]

            tcr_plddt = plddt[tcr_atom_start_idx:tcr_atom_end_idx]

            for j, pmhc_idx in enumerate(pmhc_indices):
                pmhc_atom_start_idx = np.sum(extension[:pmhc_idx], dtype=int)
                pmhc_atom_end_idx = pmhc_atom_start_idx + extension[pmhc_idx]

                pmhc_plddt = plddt[pmhc_atom_start_idx:pmhc_atom_end_idx]

                plddt_agg[i, j] = (
                    1 / (len(tcr_plddt) + len(pmhc_plddt)) * (np.sum(tcr_plddt) + np.sum(pmhc_plddt)) - 100
                )

    else:
        for i, tcr_idx in enumerate(tcr_indices):
            for j, pmhc_idx in enumerate(pmhc_indices):
                plddt_agg[i, j] = 1 / 2 * (plddt[tcr_idx] + plddt[pmhc_idx]) - 100

    scores = pae[np.ix_(tcr_mask, pmhc_mask)] - 1 / 3 * plddt_agg

    return scores

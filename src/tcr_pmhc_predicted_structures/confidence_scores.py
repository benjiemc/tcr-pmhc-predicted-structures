"""Functions for working with confidence scores from structure predictions."""

import logging
from functools import cache

import numpy as np
from Bio.SeqUtils import IUPACData

from tcr_pmhc_predicted_structures.imgt_numbering import (
    IMGT_CDR,
    IMGT_MHC_TCR_CONTACT_POSITIONS,
    MHC_I_IMGT_BETA_HELIX_START,
    number,
)
from tcr_pmhc_predicted_structures.predictions import CHAIN_ORDER, PDB_CHAIN_IDS, RESIDUE_NUMBER_OF_ATOMS

logger = logging.getLogger(__name__)


def get_sequence_information(sequences: dict[str, str]) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Create arrays from sequences dictionary with information about the sequences, chain IDs, and IMGT numbering."""
    sequences_flat = []
    chain_ids = []
    sequence_ids = []
    insert_codes = []

    for chain_id, chain_name in zip(PDB_CHAIN_IDS, CHAIN_ORDER, strict=False):
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
    binding_only: bool = True,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    tcr_mask = np.isin(chain_ids, np.array(['D', 'E']))
    peptide_mask = chain_ids == 'C'
    mhc_mask = chain_ids == 'A'

    if binding_only:
        tcr_mask = tcr_mask & np.isin(sequence_ids, np.array(sorted(IMGT_CDR)))
        ((mhc_alpha_seq_ids, mhc_alpha_insert_codes), (mhc_beta_seq_ids, mhc_beta_insert_codes)) = (
            _get_mhc_tcr_contact_positions(mhc_type)
        )

        mhc_mask = mhc_mask & np.isin(sequence_ids, mhc_alpha_seq_ids) & np.isin(insert_codes, mhc_alpha_insert_codes)

    if mhc_type == 'MH2':
        mhc_mask = mhc_mask | (chain_ids == 'B')

        if binding_only:
            mhc_mask = mhc_mask & np.isin(sequence_ids, mhc_beta_seq_ids) & np.isin(insert_codes, mhc_beta_insert_codes)

    return tcr_mask, peptide_mask, mhc_mask


def extract_confidence_score(
    sequences: np.ndarray,
    chain_ids: np.ndarray,
    sequence_ids: np.ndarray,
    insert_codes: np.ndarray,
    confidence: np.ndarray,
    *,
    score_type: str,
    mhc_type: str,
    binding_only: bool = True,
    avg_atoms_for_residue: bool = False,
) -> np.ndarray:
    """Extract TCR:pMHC specific confidence score from given confidence."""
    if score_type.startswith('atom'):
        logger.debug('Extending look up arrays from residue level to atom level')
        extension = _get_extension(sequences, chain_ids)

        chain_ids = _extend_array(chain_ids, extension)
        sequence_ids = _extend_array(sequence_ids, extension)
        insert_codes = _extend_array(insert_codes, extension)

    logger.debug('Selecting relevant portion of the matrix')
    logger.debug('Creating masks')
    tcr_mask, peptide_mask, mhc_mask = _create_masks(
        chain_ids, sequence_ids, insert_codes, mhc_type=mhc_type, binding_only=binding_only
    )
    pmhc_mask = peptide_mask | mhc_mask

    match score_type:
        case 'pae' | 'pde':
            selected_confidence = confidence[np.ix_(tcr_mask, pmhc_mask)]  # TODO get square version of this

        case 'atom_plddt' | 'plddt':
            selected_confidence = confidence[tcr_mask | pmhc_mask]

            if avg_atoms_for_residue:
                selected_confidence = _aggregate_from_extension(selected_confidence, extension)

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
    binding_only: bool = True,
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
        binding_only: Whether to only to consider the binding interface of the TCR:pMHC complex in the score or not. The
            binding interface is formed of the CDR loops, peptide, and MHC-TCR contact positions.

    Returns:
        a 2D array of shape num_selected_residues x num_selected_residues with the AlphaBridge score computed for each
        pairwise residue combination

    """
    tcr_mask, peptide_mask, mhc_mask = _create_masks(
        chain_ids, sequence_ids, insert_codes, mhc_type=mhc_type, binding_only=binding_only
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

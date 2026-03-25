"""Functions for adding missing atom and residue information to sequences."""

import logging
import re

import numpy as np
import pandas as pd
from Bio.PDB import Structure
from Bio.SeqUtils import IUPACData

from tcr_pmhc_predicted_structures.aligners import align_sequences
from tcr_pmhc_predicted_structures.chemistry import RESIDUE_ATOMS

logger: logging.Logger = logging.getLogger(__name__)


def get_missing_residues(header: str) -> list[dict]:
    """Extract all missing residues from the header of a pdb file.

    Returns residue name, chain id, residue/sequence id, and insert id if available.
    """
    lines = re.findall(r'^REMARK 465\s+\w?\s+(\w{3}) (\w)\s+(\d+)(\w)?', header, flags=re.MULTILINE)
    return [
        {
            'residue_name': res_name,
            'chain_id': chain,
            'residue_seq_id': int(seq_id),
            'residue_insert_code': insert_id or None,
        }
        for res_name, chain, seq_id, insert_id in lines
    ]


def get_missing_atoms(header: str) -> list[dict]:
    """Extract residues with missing atoms from the header of a pdb file."""
    lines = re.findall(r'^REMARK 470\s+\w?\s+(\w{3}) (\w)\s+(\d+)(\w)?\s+(\w+(?: +\w+)*)', header, re.MULTILINE)
    return [
        {
            'residue_name': res_name,
            'chain_id': chain,
            'residue_seq_id': int(seq_id),
            'residue_insert_code': insert_id or None,
            'atoms': atoms.split(),
        }
        for res_name, chain, seq_id, insert_id, atoms in lines
    ]


def annotate_with_missing_residues(structure_df: pd.DataFrame, missing_residues: list[dict]) -> pd.DataFrame:
    """Add missing residue information to the structure dataframe."""
    structure_df = structure_df.copy()

    # Move any HETATMS to the end of the chain (sometimes these are separated to the end of the file)
    structure_df = structure_df.sort_values(['chain_id', 'atom_number']).reset_index(drop=True)
    smallest_atom_number = structure_df['atom_number'].min()
    structure_df['atom_number'] = range(smallest_atom_number, smallest_atom_number + len(structure_df))

    chains = structure_df['chain_id'].unique()

    for missing_residue in missing_residues:
        chain_before = structure_df[
            (structure_df['chain_id'] == missing_residue['chain_id'])
            & (structure_df['residue_seq_id'] < missing_residue['residue_seq_id'])
        ]

        if len(chain_before) > 0:
            atom_number_before = chain_before.iloc[-1]['atom_number']

        else:
            chain_idx = np.where(chains == missing_residue['chain_id'])[0][0]
            if chain_idx == 0:
                atom_number_before = 0

            else:
                atom_number_before = structure_df[structure_df['chain_id'] == chains[chain_idx - 1]].iloc[-1][
                    'atom_number'
                ]

        df_before = structure_df[structure_df['atom_number'] <= atom_number_before].copy()

        chain_after = structure_df[
            (structure_df['chain_id'] == missing_residue['chain_id'])
            & (structure_df['residue_seq_id'] > missing_residue['residue_seq_id'])
        ]

        if len(chain_after) > 0:
            atom_number_after = chain_after.iloc[0]['atom_number']

        else:
            chain_idx = np.where(chains == missing_residue['chain_id'])[0][0]
            if chain_idx == len(chains) - 1:
                atom_number_after = np.inf

            else:
                atom_number_after = structure_df[structure_df['chain_id'] == chains[chain_idx + 1]].iloc[0][
                    'atom_number'
                ]

        df_after = structure_df[structure_df['atom_number'] >= atom_number_after].copy()

        missing_residue_df = pd.DataFrame(
            [
                {
                    'record_type': 'ATOM',
                    'atom_number': i + (atom_number_before + 1),
                    'atom_name': atom_name,
                    'element': atom_name[0],
                    'missing': True,
                    **missing_residue,
                }
                for i, atom_name in enumerate(sorted(RESIDUE_ATOMS[missing_residue['residue_name']]))
            ],
        )

        offset = (missing_residue_df.iloc[-1]['atom_number'] + 1) - atom_number_after
        if offset > 0:
            df_after['atom_number'] += offset

        structure_df = pd.concat([df_before, missing_residue_df, df_after]).reset_index(drop=True)

    structure_df['missing'] = structure_df['missing'].fillna(value=False)
    return structure_df


def annotate_with_missing_atoms(structure_df: pd.DataFrame, missing_atoms: list[dict]) -> pd.DataFrame:
    """Add missing atom information information to the structure data frame."""
    structure_df = structure_df.copy()

    for residue_missing_atoms in missing_atoms:
        if residue_missing_atoms['residue_insert_code'] is not None:
            atom_number_before = structure_df[
                (structure_df['chain_id'] == residue_missing_atoms['chain_id'])
                & (structure_df['residue_seq_id'] == residue_missing_atoms['residue_seq_id'])
                & (structure_df['residue_insert_code'] == residue_missing_atoms['residue_insert_code'])
            ].iloc[-1]['atom_number']

        else:
            atom_number_before = structure_df[
                (structure_df['chain_id'] == residue_missing_atoms['chain_id'])
                & (structure_df['residue_seq_id'] == residue_missing_atoms['residue_seq_id'])
                & pd.isna(structure_df['residue_insert_code'])
            ].iloc[-1]['atom_number']

        df_before = structure_df[structure_df['atom_number'] <= atom_number_before].copy()
        df_after = structure_df[structure_df['atom_number'] > atom_number_before].copy()

        atom_number_after = df_after.iloc[0]['atom_number'] if len(df_after) > 0 else np.inf

        missing_atoms_df = pd.DataFrame(
            [
                {
                    'record_type': 'ATOM',
                    'atom_number': i + (atom_number_before + 1),
                    'atom_name': atom_name,
                    'alt_loc': None,
                    'element': atom_name[0],
                    'missing': True,
                    **{key: val for key, val in residue_missing_atoms.items() if key != 'atoms'},
                }
                for i, atom_name in enumerate(residue_missing_atoms['atoms'])
            ],
        )

        offset = (missing_atoms_df.iloc[-1]['atom_number'] + 1) - atom_number_after
        if offset > 0:
            df_after['atom_number'] += offset

        structure_df = pd.concat([df_before, missing_atoms_df, df_after]).reset_index(drop=True)

    structure_df['missing'] = structure_df['missing'].fillna(value=False)

    return structure_df


def annotate_with_missing_entities(
    structure_df: pd.DataFrame,
    missing_residues: list[dict],
    missing_atoms: list[dict],
) -> pd.DataFrame:
    """Add the missing entities information to the structure data frame.

    This function does not modify the original data frame. There is a new boolean column called 'missing' created in the
    output data frame. The atom numbering may change if space was not left for the missing atoms in the original PDB
    structure.

    Args:
        structure_df: pandas data frame with structure information
        missing_residues: missing residue information taken from PDB header
        missing_atoms: missing atom information taken from PBD header

    Returns:
        new data frame with filled out missing residues and atoms with an additional column called 'missing'

    """
    structure_df = structure_df.copy()
    structure_df['missing'] = False

    structure_df = annotate_with_missing_residues(structure_df, missing_residues)
    structure_df = annotate_with_missing_atoms(structure_df, missing_atoms)

    return structure_df


def compare_sequences_for_missing(
    sequence: list[list], annotated_raw_sequence: list[list], numbering: set[int]
) -> bool:
    """Check if there are missing residues in the numbered regions."""
    seq = ''.join([res[0] for res in sequence])
    raw_seq = ''.join([res[0] for res in annotated_raw_sequence])

    alignment, _ = align_sequences(raw_seq, seq)

    raw_index = 0
    index = 0
    current_seq_id = 0

    for raw_res, res in alignment:
        if raw_res == '-':
            index += 1
            continue

        if res != '-':
            current_seq_id = sequence[index][1]

        if annotated_raw_sequence[raw_index][3] == 'missing' and current_seq_id in numbering:
            return False

        if res == '-':
            raw_index += 1
            continue

        index += 1
        raw_index += 1

    return True


def get_alignment(
    structure: Structure.Structure, raw_structure: Structure.Structure, chain_id: str, missing_residues: list[dict]
) -> list[tuple[str, str]]:
    """Get alignment between original sequence and full sequence missing residues."""
    missing_residues_on_chain = [
        (
            IUPACData.protein_letters_3to1[entity['residue_name'].title()],
            entity['residue_seq_id'],
            entity['residue_insert_code'],
        )
        for entity in missing_residues
        if entity['chain_id'] == chain_id
    ]

    raw_sequence_info = [
        (IUPACData.protein_letters_3to1[res.get_resname().title()], res.id[1], res.id[2].strip())
        for res in raw_structure[0][chain_id].get_residues()
        if res.id[0] == ' '
    ]

    combined_raw_residues = sorted(missing_residues_on_chain + raw_sequence_info, key=lambda res: (res[1], res[2]))
    raw_seq = ''.join([res[0] for res in combined_raw_residues])

    imgt_sequence = ''.join(
        [
            IUPACData.protein_letters_3to1[res.get_resname().title()]
            for res in structure[0][chain_id].get_residues()
            if res.id[0] == ' '
        ]
    )

    alignment, _ = align_sequences(raw_seq, imgt_sequence)

    return alignment


def trim_start_and_end(annotated_raw_sequence: list[list]) -> list[list]:
    """Remove missing residues at start and end of the sequence."""
    start_index = 0
    for res in annotated_raw_sequence:
        if res[3] != 'missing':
            break

        start_index += 1

    end_index = len(annotated_raw_sequence)
    for res in reversed(annotated_raw_sequence):
        if res[3] != 'missing':
            break

        end_index -= 1

    return annotated_raw_sequence[start_index:end_index]


def screen_chain(
    structure: Structure.Structure,
    raw_structure: Structure.Structure,
    chain_id: str,
    missing_residues: list[dict],
    missing_atoms: list[dict],
    numbering: set[int],
    *,
    trim_ends: bool = True,
) -> bool:
    """Screen chain for missing residues or atoms."""
    missing_residues_on_chain = [
        [
            IUPACData.protein_letters_3to1[entity['residue_name'].title()],
            entity['residue_seq_id'],
            entity['residue_insert_code'],
            'missing',
        ]
        for entity in missing_residues
        if entity['chain_id'] == chain_id
    ]

    missing_atoms_on_chain = [
        [
            IUPACData.protein_letters_3to1[entity['residue_name'].title()],
            entity['residue_seq_id'],
            entity['residue_insert_code'],
        ]
        for entity in missing_atoms
        if entity['chain_id'] == chain_id
    ]

    raw_sequence_info = [
        [
            IUPACData.protein_letters_3to1[res.get_resname().title()],
            res.id[1],
            res.id[2].strip() or None,
            'not-missing',
        ]
        for res in raw_structure[0][chain_id].get_residues()
        if res.id[0] == ' '
    ]
    imgt_sequence_info = [
        [
            IUPACData.protein_letters_3to1[res.get_resname().title()],
            res.id[1],
            res.id[2].strip() or None,
        ]
        for res in structure[0][chain_id].get_residues()
        if res.id[0] == ' '
    ]

    combined_raw_residues = sorted(
        missing_residues_on_chain + raw_sequence_info,
        key=lambda res: (res[1], res[2] if res[2] is not None else ''),
    )

    if trim_ends:
        combined_raw_residues = trim_start_and_end(combined_raw_residues)

    for res in combined_raw_residues:
        if res[:3] in missing_atoms_on_chain:
            res[3] = 'missing'

    return compare_sequences_for_missing(imgt_sequence_info, combined_raw_residues, numbering)

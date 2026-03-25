"""Functions and classes for interacting with PDB structures."""

import logging
from collections.abc import Iterable

import pandas as pd
from Bio.PDB import Chain, Model, Structure
from Bio.SeqUtils import IUPACData

from tcr_pmhc_predicted_structures.imgt_numbering import IMGT_MH1_ABD, IMGT_MH2_ABD, IMGT_VARIABLE_DOMAIN

logger = logging.getLogger(__name__)


def get_sequence(structure: Structure.Structure, chain_id: str, residue_range: set[int] | None = None) -> str:
    """Get the sequence of amino acids from a biopython strcture.

    Args:
        structure: Biopython structure.
        chain_id: ID of the chain to get the sequence from.
        residue_range: Range of residues to include in the sequence (Optional). Default is to include them all.

    Returns:
        Amino acid sequence as one-letter codes.

    """
    chain = structure[0][chain_id]
    sequence = []

    for residue in chain:
        if residue.id[0] != ' ':
            continue

        if residue_range is not None and residue.id[1] not in residue_range:
            continue

        sequence.append(residue.get_resname())

    sequence = [IUPACData.protein_letters_3to1[res.title()] for res in sequence]

    return ''.join(sequence)


def get_header(pdb_contents: str) -> str:
    """Get the header lines from the contents of a pdb file."""
    header = []
    for line in pdb_contents.split('\n'):
        record_type = line[0:6]
        if record_type in ('ATOM  ', 'HETATM', 'MODEL '):
            break

        header.append(line)

    return '\n'.join(header)


def extract_chains(structure: Structure.Structure, chains: Iterable[str]) -> Structure.Structure:
    """Get only the selected chain from a PDB structure."""
    new_structure = Structure.Structure(structure.id)

    for model in structure:
        new_model = Model.Model(model.id)
        for chain in model:
            if chain.id in chains:
                new_model.add(chain.copy())

        new_structure.add(new_model)

    return new_structure


def bio_to_pandas(structure: Structure.Structure) -> pd.DataFrame:
    """Convert a biopython structure to a pandas dataframe."""
    records = []

    multiple_models = len(structure) > 1

    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    records.append(
                        {
                            'record_type': 'ATOM' if residue.id[0] == ' ' else 'HETATM',
                            'atom_number': atom.serial_number,
                            'atom_name': atom.id,
                            'alt_loc': atom.altloc if atom.altloc != ' ' else None,
                            'residue_name': residue.resname,
                            'chain_id': chain.id,
                            'residue_seq_id': residue.id[1],
                            'residue_insert_code': residue.id[2] if residue.id[2] != ' ' else None,
                            'pos_x': atom.coord[0],
                            'pos_y': atom.coord[1],
                            'pos_z': atom.coord[2],
                            'occupancy': atom.occupancy,
                            'b_factor': atom.bfactor,
                            'element': atom.element,
                            'charge': atom.get_charge(),
                        }
                    )

                    if multiple_models:
                        records[-1]['model_index'] = model.id

        column_names = [
            'record_type',
            'atom_number',
            'atom_name',
            'alt_loc',
            'residue_name',
            'chain_id',
            'residue_seq_id',
            'residue_insert_code',
            'pos_x',
            'pos_y',
            'pos_z',
            'occupancy',
            'b_factor',
            'element',
            'charge',
        ]
    if multiple_models:
        column_names.append('model_index')

    return pd.DataFrame(records, columns=column_names)


def replace_chain(structure: Structure.Structure, new_chain: Chain.Chain) -> Structure.Structure:
    """Replace chain with a new chain in a PDB structure (does not modify original structure)."""
    output_structure = structure.copy()
    chain_id = new_chain.id

    for model in output_structure:
        model.detach_child(chain_id)
        model.add(new_chain)

    return output_structure


def crop_chain(chain: Chain.Chain, numbering: set[int]) -> Chain.Chain:
    """Crop chain based on numbering."""
    new_chain = Chain.Chain(chain.id)

    for residue in chain:
        if residue.id[1] in numbering or residue.id[0] != ' ':
            new_chain.add(residue.copy())

    return new_chain


def crop_structure(  # noqa: C901
    structure: Structure.Structure,
    tcr_chains: Iterable[str] | None = None,
    mhc_chains: Iterable[str] | None = None,
    mhc_type: str | None = None,
) -> Structure.Structure:
    """Crop TCR and MHC structures to the TCR variable domain and MHC antigen binding domain.

    Args:
        structure: structure to crop
        tcr_chains: TCR chains to crop (optional)
        mhc_chains: MHC chains to crop (optional)
        mhc_type: either 'MH1' or 'MH2', required if `mhc_chains` is specified

    Returns:
        structure with chains cropped to variable domain or antigen binding domain

    Raises:
        ValueError: invalid `mhc_type` is entered or `mhc_type` is not specified and `mhc_chains` are input

    """
    if tcr_chains is None:
        tcr_chains = []

    if mhc_chains is None:
        mhc_chains = []

    mhc_numbering = None
    if len(mhc_chains) > 0:
        match mhc_type:
            case 'MH1':
                mhc_numbering = IMGT_MH1_ABD

            case 'MH2':
                mhc_numbering = IMGT_MH2_ABD

            case None:
                msg = 'MHC chains inputted but no MHC type specified'
                raise ValueError(msg)

            case _:
                msg = f'Invalid MHC type: {mhc_type}'
                raise ValueError(msg)

    cropped_structure = Structure.Structure(structure.id)

    for model in structure:
        new_model = Model.Model(model.id)

        for chain in model:
            if chain.id in tcr_chains:
                logger.debug('Cropping TCR chain %s', chain.id)
                new_chain = crop_chain(chain, IMGT_VARIABLE_DOMAIN)

            elif chain.id in mhc_chains:
                logger.debug('Cropping MHC chain %s', chain.id)
                new_chain = crop_chain(chain, mhc_numbering)

            else:
                new_chain = chain.copy()

            new_model.add(new_chain)

        cropped_structure.add(new_model)

    return cropped_structure


def remove_het_atoms(structure: Structure.Structure) -> Structure.Structure:
    """Remove hetero atoms from structure."""
    output_structure = Structure.Structure(structure.id)

    for model in structure:
        new_model = Model.Model(model.id)

        for chain in model:
            new_chain = Chain.Chain(chain.id)

            for res in chain:
                if res.id[0] == ' ':
                    new_chain.add(res.copy())

            new_model.add(new_chain)

        output_structure.add(new_model)

    return output_structure

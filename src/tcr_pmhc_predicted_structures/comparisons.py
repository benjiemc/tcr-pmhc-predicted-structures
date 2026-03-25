"""Functions for comparing data points."""

import numpy as np
from Bio.PDB import Atom, Structure, Superimposer
from Bio.SeqUtils import IUPACData

from tcr_pmhc_predicted_structures.aligners import align_sequences
from tcr_pmhc_predicted_structures.chemistry import HEAVY_ATOMS
from tcr_pmhc_predicted_structures.imgt_numbering import IMGT_MH1_ABD, IMGT_MH2_ABD, IMGT_VARIABLE_DOMAIN


def rmsd(coords1: np.ndarray, coords2: np.ndarray) -> float:
    """Calculate Root Mean Squard Deviation between two sets of coordinates."""
    diff = coords2 - coords1
    distance = np.sqrt(np.sum(diff * diff, axis=1))

    return np.sqrt(np.sum(distance * distance) / coords1.shape[0])


def find_equivalent_sequences(
    struct1: Structure.Structure, chain_map1: dict[str, str], struct2: Structure.Structure, chain_map2: dict[str, str]
) -> list[tuple]:
    """Find equivalent residues between two structures."""
    equivalent_residues = []

    for chain_type in set(chain_map1.keys()) & set(chain_map1.keys()):
        residues1 = [res for res in struct1[0][chain_map1[chain_type]] if res.id[0] == ' ']
        residues2 = [res for res in struct2[0][chain_map2[chain_type]] if res.id[0] == ' ']

        seq1 = ''.join([IUPACData.protein_letters_3to1[res.get_resname().title()] for res in residues1])
        seq2 = ''.join([IUPACData.protein_letters_3to1[res.get_resname().title()] for res in residues2])

        alignment, _ = align_sequences(seq1, seq2)

        index1 = 0
        index2 = 0

        for res1, res2 in alignment:
            if res1 == '-' or res2 == '-':
                if res1 == '-':
                    index2 += 1

                if res2 == '-':
                    index1 += 1

            else:
                if res1 == res2:
                    equivalent_residues.append(
                        (
                            (
                                chain_map1[chain_type],
                                residues1[index1].get_resname(),
                                residues1[index1].id[1],
                                residues1[index1].id[2].rstrip(),
                            ),
                            (
                                chain_map2[chain_type],
                                residues2[index2].get_resname(),
                                residues2[index2].id[1],
                                residues2[index2].id[2].rstrip(),
                            ),
                        )
                    )

                index1 += 1
                index2 += 1

    return equivalent_residues


def get_relevant_atoms(  # noqa: C901, PLR0912
    struct1: Structure.Structure,
    chain_map1: dict[str, str],
    struct2: Structure.Structure,
    chain_map2: dict[str, str],
    equivalent_residues: list[tuple],
    mhc_type: str,
) -> tuple[list[Atom.Atom], list[Atom.Atom]]:
    """Get the relevant atoms for the TCR variable domain, antigen, and MHC binding domain.

    TODO: Refactor if statements

    Args:
        struct1: biopython TCR:pMHC structure
        chain_map1: dictionary mapping the alpha_chain, beta_chain, antigen_chain, mhc_chain1, and mhc_chain2 to
                    their specific PDB chain ids in the first structure
        struct2: biopython TCR:pMHC structure
        chain_map2: dictionary mapping the alpha_chain, beta_chain, antigen_chain, mhc_chain1, and mhc_chain2 to
                    their specific PDB chain ids in the second structure
        equivalent_residues: equivalent residues between the two structures
        mhc_type: type of MHC molecule ('MH1' or 'MH2')

    Returns:
        two lists of relevant atoms in the structures to align or compute RMSDs from

    """
    atoms1 = []
    atoms2 = []

    chain_map1_inv = {chain_id: chain_type for chain_type, chain_id in chain_map1.items()}
    chain_map2_inv = {chain_id: chain_type for chain_type, chain_id in chain_map2.items()}

    for res1, res2 in equivalent_residues:
        res_fixed_id = (' ', res1[2], res1[3] if res1[3] != '' else ' ')
        res_mobile_id = (' ', res2[2], res2[3] if res2[3] != '' else ' ')

        relevant_atoms = False

        if (chain_map1_inv[res1[0]] == chain_map2_inv[res2[0]] == 'alpha_chain') or (
            chain_map1_inv[res1[0]] == chain_map2_inv[res2[0]] == 'beta_chain'
        ):
            if res1[2] in IMGT_VARIABLE_DOMAIN and res2[2] in IMGT_VARIABLE_DOMAIN:
                relevant_atoms = True

        elif chain_map1_inv[res1[0]] == chain_map2_inv[res2[0]] == 'antigen_chain':
            relevant_atoms = True

        else:  # noqa: PLR5501
            if mhc_type == 'MH1':
                if chain_map1_inv[res1[0]] == chain_map2_inv[res2[0]] == 'mhc_chain1':  # noqa: SIM102
                    if res1[2] in IMGT_MH1_ABD and res2[2] in IMGT_MH1_ABD:
                        relevant_atoms = True

            elif mhc_type == 'MH2':
                if (  # noqa: SIM102
                    chain_map1_inv[res1[0]] == chain_map2_inv[res2[0]] == 'mhc_chain1'
                    or chain_map1_inv[res1[0]] == chain_map2_inv[res2[0]] == 'mhc_chain2'
                ):
                    if res1[1] in IMGT_MH2_ABD and res2[1] in IMGT_MH2_ABD:
                        relevant_atoms = True

            else:
                msg = f'Incorrect MHC type: {mhc_type}. mhc_type should be MH1 or MH2'
                raise ValueError(msg)

        if relevant_atoms:
            res_atoms1 = [atom for atom in struct1[0][res1[0]][res_fixed_id].get_atoms() if atom.element in HEAVY_ATOMS]
            res_atoms2 = [
                atom for atom in struct2[0][res2[0]][res_mobile_id].get_atoms() if atom.element in HEAVY_ATOMS
            ]

            if len(res_atoms1) == len(res_atoms2):
                atoms1 += res_atoms1
                atoms2 += res_atoms2

    return atoms1, atoms2


def compute_structural_distances(
    structures: list[Structure.Structure], chain_maps: list[dict], mhc_type: str
) -> np.array:
    """Create a distance matrix of RMSD between the input TCR:pMHC structures.

    The comparison is done between the TCR variable domain, peptide, and MHC antigen binding domain.

    Args:
        structures: biopython structures to compare
        chain_maps: dictionaries mapping 'alpha_chain', 'beta_chain', 'antigen_chain', 'mhc_chain1', and 'mhc_chain2' to
                    the PDB chain ids
        mhc_type: either 'MH1' or 'MH2' to select the antigen binding domain of the MHC molecules

    Returns:
        distance matrix of the RMSD between equivalent atoms in the structures

    """
    distance_matrix = np.zeros((len(structures), len(structures)))

    for i, (struct1, chain_map1) in enumerate(zip(structures[:-1], chain_maps[:-1], strict=False)):
        for j, (struct2, chain_map2) in enumerate(zip(structures[i + 1 :], chain_maps[i + 1 :], strict=False), i + 1):
            equivalent_residues = find_equivalent_sequences(struct1, chain_map1, struct2, chain_map2)
            atoms1, atoms2 = get_relevant_atoms(struct1, chain_map1, struct2, chain_map2, equivalent_residues, mhc_type)

            struct2_aligned = struct2.copy()
            super_imposer = Superimposer()
            super_imposer.set_atoms(atoms1, atoms2)
            super_imposer.apply(struct2_aligned.get_atoms())

            atoms1, atoms2 = get_relevant_atoms(
                struct1, chain_map1, struct2_aligned, chain_map2, equivalent_residues, mhc_type
            )

            distance = rmsd(
                np.array([atom.get_coord() for atom in atoms1]), np.array([atom.get_coord() for atom in atoms2])
            )

            distance_matrix[i, j] = distance

    return np.maximum(distance_matrix, distance_matrix.T)

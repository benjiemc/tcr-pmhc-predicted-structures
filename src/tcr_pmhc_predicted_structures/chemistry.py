"""Chemistry constants."""

RESIDUE_ATOMS: dict[set[str]] = {
    'ALA': {'C', 'CA', 'CB', 'N', 'O'},
    'ARG': {'C', 'CA', 'CB', 'CD', 'CG', 'CZ', 'N', 'NE', 'NH1', 'NH2', 'O'},
    'ASN': {'C', 'CA', 'CB', 'CG', 'N', 'ND2', 'O', 'OD1'},
    'ASP': {'C', 'CA', 'CB', 'CG', 'N', 'O', 'OD1', 'OD2'},
    'CYS': {'C', 'CA', 'CB', 'N', 'O', 'SG'},
    'GLN': {'C', 'CA', 'CB', 'CD', 'CG', 'N', 'NE2', 'O', 'OE1'},
    'GLU': {'C', 'CA', 'CB', 'CD', 'CG', 'N', 'O', 'OE1', 'OE2'},
    'GLY': {'C', 'CA', 'N', 'O'},
    'HIS': {'C', 'CA', 'CB', 'CD2', 'CE1', 'CG', 'N', 'ND1', 'NE2', 'O'},
    'ILE': {'C', 'CA', 'CB', 'CD1', 'CG1', 'CG2', 'N', 'O'},
    'LEU': {'C', 'CA', 'CB', 'CD1', 'CD2', 'CG', 'N', 'O'},
    'LYS': {'C', 'CA', 'CB', 'CD', 'CE', 'CG', 'N', 'NZ', 'O'},
    'MET': {'C', 'CA', 'CB', 'CE', 'CG', 'N', 'O', 'SD'},
    'PHE': {'C', 'CA', 'CB', 'CD1', 'CD2', 'CE1', 'CE2', 'CG', 'CZ', 'N', 'O'},
    'PRO': {'C', 'CA', 'CB', 'CD', 'CG', 'N', 'O'},
    'SER': {'C', 'CA', 'CB', 'N', 'O', 'OG'},
    'THR': {'C', 'CA', 'CB', 'CG2', 'N', 'O', 'OG1'},
    'TRP': {'C', 'CA', 'CB', 'CD1', 'CD2', 'CE2', 'CE3', 'CG', 'CH2', 'CZ2', 'CZ3', 'N', 'NE1', 'O'},
    'TYR': {'C', 'CA', 'CB', 'CD1', 'CD2', 'CE1', 'CE2', 'CG', 'CZ', 'N', 'O', 'OH'},
    'VAL': {'C', 'CA', 'CB', 'CG1', 'CG2', 'N', 'O'},
}
'''Atoms in each amino acid residue.'''

HEAVY_ATOMS: set[str] = {'C', 'N', 'O', 'S'}
'''Heavy atoms that make up amino acids.'''

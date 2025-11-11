"""Constants used to work with predicted structures."""

CHAIN_ORDER = ('mhc_chain_1', 'mhc_chain_2', 'peptide', 'tcr_alpha_chain', 'tcr_beta_chain')
'''Order the chains appear in output structure files.'''

PDB_CHAIN_IDS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789'
'''Possible chain IDS for PDB structures.'''

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

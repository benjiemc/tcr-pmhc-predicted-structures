"""Constants for annotating sequences as CDR domains in T cell receptors."""

IMGT_CDR1: set[int] = set(range(27, 38 + 1))
'''IMGT residue numbers corresponding to CDR 1 domains.'''
IMGT_CDR2: set[int] = set(range(56, 65 + 1))
'''IMGT residue numbers corresponding to CDR 2 domains.'''
IMGT_CDR3: set[int] = set(range(105, 117 + 1))
'''IMGT residue numbers corresponding to CDR 3 domains.'''
IMGT_CDR: set[int] = IMGT_CDR1.union(IMGT_CDR2).union(IMGT_CDR3)
'''IMGT residue numbers corresponding to all CDR domains.'''

IMGT_VARIABLE_DOMAIN: set[int] = set(range(1, 128 + 1))
'''Variable domain range for IMGT numbered TCR structures.'''

IMGT_FRAMEWORK_REGION: set[int] = IMGT_VARIABLE_DOMAIN - IMGT_CDR
'''Framework (Fw) region range for IMGT numbered TCR structures.'''

IMGT_MH1_ABD: set[int] = set(range(1, 92)) | set(range(1001, 1092))
'''IMGT ranges of the antigen binding domain of MHC class I molecules.'''

IMGT_MH2_ABD: set[int] = set(range(1, 92))
'''IMGT ranges of the antigen binding domain of MHC class II molecules.'''

MHC_I_IMGT_BETA_HELIX_START: int = 1000
'''Start of the beta helix of class I molecules.'''


def assign_cdr_number(seq_id: int) -> int | None:
    """Assign CDR number for a sequence ID or return None if the sequence ID if the ID is not in CDR range."""
    if seq_id in IMGT_CDR1:
        return 1

    if seq_id in IMGT_CDR2:
        return 2

    if seq_id in IMGT_CDR3:
        return 3

    return None

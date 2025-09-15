#!/usr/bin/env python
"""Profile TCR:peptide interactions."""

import stcrpy
from snakemake.script import snakemake

tcr = stcrpy.load_TCR(snakemake.input[0])
tcr.standardise_chain_names()

try:
    interactions = tcr.profile_peptide_interactions()

except Exception as err:
    if type(err).__name__ == 'ArgumentError' and type(err).__module__ == 'Boost.Python':
        print('Interaction profiling raised an error:', err)
        interactions = None

    else:
        raise

if interactions is not None:
    interactions = interactions.query("protein_chain == 'D' or protein_chain == 'E'")

    interactions.to_csv(snakemake.output[0], index=False)

else:
    with open(snakemake.output[0], 'w') as fh:
        fh.write(
            'type,protein_atom,protein_chain,protein_residue,protein_number,ligand_atom,distance,angle,plip_id,'
            'original_numbering,ligand_residue,ligand_number,domain\n'
        )

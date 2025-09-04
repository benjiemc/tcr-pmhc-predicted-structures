#!/usr/bin/env python
"""Profile TCR:peptide interactions."""

import stcrpy
from snakemake.script import snakemake

tcr = stcrpy.load_TCR(snakemake.input[0])
tcr.standardise_chain_names()

interactions = tcr.profile_peptide_interactions()

if interactions is not None:
    interactions = interactions.query("protein_chain == 'D' or protein_chain == 'E'")

    interactions.to_csv(snakemake.output[0], index=False)

else:
    with open(snakemake.output[0], 'w') as fh:
        fh.write(
            'type,protein_atom,protein_chain,protein_residue,protein_number,ligand_atom,distance,angle,plip_id,'
            'original_numbering,ligand_residue,ligand_number,domain\n'
        )

#!/usr/bin/env python
"""Profile TCR:peptide interactions."""

import stcrpy
from snakemake.script import snakemake

tcr = stcrpy.load_TCR(snakemake.input[0])
tcr.standardise_chain_names()

interactions = tcr.profile_peptide_interactions()
interactions = interactions.query("protein_chain == 'D' or protein_chain == 'E'")

interactions.to_csv(snakemake.output[0], index=False)

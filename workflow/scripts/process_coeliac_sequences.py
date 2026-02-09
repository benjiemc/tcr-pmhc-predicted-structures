#!/usr/bin/env python
"""Process Coeliac TCR sequences and create full length sequences."""

import logging
import re

import pandas as pd
import tidytcells
from snakemake.script import snakemake

from tcr_pmhc_predicted_structures.apps._log import setup_logger
from tcr_pmhc_predicted_structures.utils import stitch_sequence

logger = logging.getLogger()
setup_logger(logger, snakemake.config['log_level'])

logger.info('Loading data from %s', snakemake.input[0])
coeliac_tcrs = pd.read_excel(snakemake.input[0])

logger.info('Separate out multiple designations')

logger.debug('Separating multiple TCR alpha chains')
for original, second in (('va_gene', 'va2_gene'), ('ja_gene', 'ja2_gene'), ('cdr3a', 'cdr3a2')):
    coeliac_tcrs[original] = coeliac_tcrs[[original, second]].apply(
        lambda row: [entry for entry in row if not pd.isna(entry)], axis=1
    )

coeliac_tcrs = coeliac_tcrs.explode(['va_gene', 'ja_gene', 'cdr3a']).reset_index(drop=True)

logger.debug('Separating multiple TCR beta chains')
extra_cdr3bs = pd.Series([set() for _ in range(len(coeliac_tcrs))])

double_beta_mask = coeliac_tcrs['clone_comment'].map(
    lambda comment: not pd.isna(comment) and '2 beta chains' in comment
)

extra_cdr3bs[double_beta_mask] = coeliac_tcrs[double_beta_mask]['clone_comment'].map(
    lambda comment: {entry.split(':')[1] for entry in re.findall('TR[A-B]:[A-Z]+', comment) if entry.startswith('TRB')}
)

coeliac_tcrs['cdr3b'] = (
    coeliac_tcrs['cdr3b'].map(lambda cdr3: {cdr3}).combine(extra_cdr3bs, lambda s1, s2: s1 | s2).map(list)
)

coeliac_tcrs = coeliac_tcrs.explode('cdr3b').reset_index(drop=True)

logger.debug('Separating multiple antigens')
coeliac_tcrs['pMHC_calls'] = coeliac_tcrs['pMHC_calls'].replace('a1w1', 'a1,w1').str.split(',')
coeliac_tcrs = coeliac_tcrs.explode('pMHC_calls').reset_index(drop=True)

logger.info('Adding epitope sequence information')

logger.debug('Adding peptide information')
# Jones, Lim et al 2026, From Gliadin-Specific CD4⁺ T Cells Exhibit Functional Diversity Independent of TCR Usage or
# Epitope Specificity
peptide_map = {
    'a1': 'QPFPQPELPYP',
    'a2': 'FPQPELPYPQP',
    'w1': 'QPFPQPEQPFP',
    'w2': 'PQPEQPFPWQ',
}
coeliac_tcrs['peptide'] = coeliac_tcrs['pMHC_calls'].map(peptide_map)

logger.debug('Adding MHC information')
coeliac_tcrs['mhc_type'] = 'MH2'
coeliac_tcrs['mhc1'] = 'HLA-DQA1*05:01'
coeliac_tcrs['mhc2'] = 'HLA-DQB1*02:01'

# HLA-DQA1*05:01: https://www.ebi.ac.uk/ipd/imgt/hla/alleles/allele/?accession=HLA00613
# HLA-DQB1*02:01: https://www.ebi.ac.uk/ipd/imgt/hla/alleles/allele/?accession=HLA00622
mhc_sequences = {
    'HLA-DQA1*05:01': (
        'IVADHVASYGVNLYQSYGPSGQYTHEFDGDEQFYVDLGRKETVWCLPVLRQFRFDPQFALTNIAVLKHNLNSLIKRSNSTAATNEVPEVTVFSKSPVTLGQPNILICLVD'
        'NIFPPVVNITWLSNGHSVTEGVSETSFLSKSDHSFFKISYLTLLPSAEESYDCKVEHWGLDKPLLKHWEP'
    ),
    'HLA-DQB1*02:01': (
        'EDFVYQFKGMCYFTNGTERVRLVSRSIYNREEIVRFDSDVGEFRAVTLLGLPAAEYWNSQKDILERKRAAVDRVCRHNYQLELRTTLQRRVEPTVTISPSRTEALNHHNL'
        'LVCSVTDFYPAQIKVRWFRNDQEETAGVVSTPLIRNGDWTFQILVMLEMTPQRGDVYTCHVEHPSLQSPITVEWRAQSES'
    ),
}

coeliac_tcrs['mhc_chain_1'] = mhc_sequences['HLA-DQA1*05:01']
coeliac_tcrs['mhc_chain_2'] = mhc_sequences['HLA-DQB1*02:01']

logger.info('Standardising gene assignments')
coeliac_tcrs[['va_gene', 'ja_gene', 'vb_gene', 'jb_gene']] = coeliac_tcrs[
    ['va_gene', 'ja_gene', 'vb_gene', 'jb_gene']
].map(tidytcells.tr.standardise, species='Homo Sapiens')

logger.info('Standardising CDR3 assignments')
coeliac_tcrs[['cdr3a', 'cdr3b']] = coeliac_tcrs[['cdr3a', 'cdr3b']].map(
    tidytcells.junction.standardise, species='Homo Sapiens'
)

logger.info('Creating full length TCR sequences')

logger.debug('Creating full length alpha chains')
coeliac_tcrs['tcr_alpha_chain'] = coeliac_tcrs[['va_gene', 'ja_gene', 'cdr3a']].apply(
    lambda row: stitch_sequence(*row, 'Human'), axis=1
)

logger.debug('Creating full length beta chains')
coeliac_tcrs['tcr_beta_chain'] = coeliac_tcrs[['vb_gene', 'jb_gene', 'cdr3b']].apply(
    lambda row: stitch_sequence(*row, 'Human'), axis=1
)

logger.debug('Adding new IDs')
coeliac_tcrs['clone_id_repeat'] = coeliac_tcrs.groupby('clone_id').cumcount() + 1
coeliac_tcrs['id'] = (
    'coeliac_'
    + coeliac_tcrs['clone_id'].str.replace('clonotype', '')
    + '_'
    + coeliac_tcrs['clone_id_repeat'].astype(str)
)
coeliac_tcrs = coeliac_tcrs.set_index('id')

logger.info('Outputting to %s', snakemake.output[0])
output = coeliac_tcrs[
    [
        'clone_id',
        'Donor_ID',
        'pMHC_calls',
        'va_gene',
        'cdr3a',
        'ja_gene',
        'vb_gene',
        'cdr3b',
        'jb_gene',
        'mhc1',
        'mhc2',
        'mhc_type',
        'tcr_alpha_chain',
        'tcr_beta_chain',
        'peptide',
        'mhc_chain_1',
        'mhc_chain_2',
    ]
]
output.columns = [column_name.lower() for column_name in output.columns]

output.to_csv(snakemake.output[0])

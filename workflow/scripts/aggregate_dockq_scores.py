#!/usr/bin/env python
"""Aggregate DockQ score Jsons into one CSV."""

import json
import os
import re

import pandas as pd
from snakemake.script import snakemake

DOCKQ_COLUMNS = [
    'name',
    'seed',
    'sample',
    'DockQ',
    'F1',
    'iRMSD',
    'LRMSD',
    'fnat',
    'nat_correct',
    'nat_total',
    'fnonnat',
    'nonnat_count',
    'model_total',
    'clashes',
    'len1',
    'len2',
]

dockq_scores = []
for path in snakemake.input:
    with open(path, 'r') as fh:
        dockq_score = json.load(fh)

    info = re.search(
        r'[A-Za-z0-9]+'
        r'_[A-Za-z0-9]+'
        r'_[A-Za-z0-9]+'
        r'(?:_no_unpaired_msas|_no_paired_msas|_no_templates)*'
        r'_(?P<name>[A-Za-z0-9_-]+)_(?P<seed>\d+)_(?P<sample>\d+)'
        '_dockq.json',
        os.path.basename(path),
    ).groupdict()

    for key, value in info.items():
        dockq_score[key] = value

    dockq_scores.append(dockq_score)

dockq_scores = pd.DataFrame(dockq_scores)

for column_name in DOCKQ_COLUMNS:
    if column_name not in dockq_scores.columns:
        dockq_scores[column_name] = None

dockq_scores = dockq_scores[DOCKQ_COLUMNS]
dockq_scores.to_csv(snakemake.output[0], index=False)

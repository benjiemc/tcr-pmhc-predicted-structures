#!/usr/bin/env python
"""Aggregate results from post-processing steps."""

import os
import re

import pandas as pd
from snakemake.script import snakemake

results = []
for path in snakemake.input:
    basename = os.path.basename(path)
    filetype = basename.split('.', 1)[-1]

    match filetype:
        case 'csv':
            result = pd.read_csv(path)

        case 'json':
            result = pd.read_json(path, typ='series').to_frame().T

        case _:
            msg = f'{filetype} not supported.'
            raise ValueError(msg)

    info = re.search(
        r'^[A-Za-z0-9]+_[A-Za-z0-9]+_[A-Za-z0-9]+_(?P<name>[A-Za-z0-9_]+)_(?P<seed>\d+)_(?P<sample>\d+)', basename
    ).groupdict()

    for key, value in info.items():
        result[key] = value

    results.append(result)

df = pd.concat(results).reset_index(drop=True)
df.to_csv(snakemake.output[0], index=False)

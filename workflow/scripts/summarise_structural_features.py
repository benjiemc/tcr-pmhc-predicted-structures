#!/usr/bin/env python
"""Calculate structural properties of TCR:pMHC structures."""

import json
import warnings

import numpy as np
import stcrpy
from snakemake.script import snakemake


def proj(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Projection of vector a onto b."""
    return (np.dot(a, b) * b) / np.linalg.norm(b) ** 2


X_AXIS = np.array([1, 0, 0])
Y_AXIS = np.array([0, 1, 0])
Z_AXIS = np.array([0, 0, 1])

with warnings.catch_warnings():
    warnings.simplefilter('ignore', UserWarning)
    tcr = stcrpy.load_TCR(snakemake.input[0])

tcr_info = {}

docking_geometry = tcr.calculate_docking_geometry(mode='cys')

if docking_geometry is not None:
    tcr_info['pitch_angle'] = float(docking_geometry['pitch_angle'])
    tcr_info['scanning_angle'] = float(docking_geometry['scanning_angle'])
    tcr_info['x_score'] = float(docking_geometry['tcr_com'][0][0])
    tcr_info['y_score'] = float(docking_geometry['tcr_com'][0][1])
    tcr_info['z_score'] = float(docking_geometry['tcr_com'][0][2])

    tcr_vec = np.array(docking_geometry['tcr_VB_cys_centroid'][0] - docking_geometry['tcr_VA_cys_centroid'][0])
    tcr_vec /= np.linalg.norm(tcr_vec)

    tcr_vec_xz = tcr_vec - proj(tcr_vec, Y_AXIS)
    tcr_vec_xz /= np.linalg.norm(tcr_vec_xz)
    tcr_info['roll_angle'] = 90 - np.degrees(np.arccos(np.dot(tcr_vec_xz, Z_AXIS)))

    tcr_vec_yz = tcr_vec - proj(tcr_vec, X_AXIS)
    tcr_vec_yz /= np.linalg.norm(tcr_vec_yz)
    tcr_info['tilt_angle'] = 90 - np.degrees(np.arccos(np.dot(tcr_vec_yz, Z_AXIS)))

with open(snakemake.output[0], 'w') as fh:
    json.dump(tcr_info, fh)

"""Functions for calculating various diversity metrics."""

import numpy as np


def shannon_entropy(probabilities: np.ndarray) -> float:
    """Calculate the Shannon Entropy of probabilities."""
    return -np.sum(probabilities * np.log(probabilities))


def pielou_entropy(probabilities: np.ndarray) -> float:
    """Calculate the Pielou Entropy (normalised Shannon Entropy) of probabilities."""
    return shannon_entropy(probabilities) / np.log(len(probabilities))


def simpson_index(probabilities: np.ndarray) -> float:
    """Calculate the Simpson Diversity Index of probabilities."""
    return np.sum(probabilities * probabilities)

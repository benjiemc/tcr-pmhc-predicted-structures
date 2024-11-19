import numpy as np


def shannon_entropy(probabilities: np.ndarray) -> float:
    return -np.sum(probabilities * np.log(probabilities))

def pielou_entropy(probabilities: np.ndarray) -> float:
    return shannon_entropy(probabilities) / np.log(len(probabilities))

def simpson_index(probabilities: np.ndarray) -> float:
    return np.sum(probabilities * probabilities)
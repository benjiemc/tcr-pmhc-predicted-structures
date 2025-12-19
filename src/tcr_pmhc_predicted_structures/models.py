"""Models to predict T cell antigen specificity using predicted structures."""

import torch
from torch import nn


class PAEAggregator(nn.Module):
    """Model to learn the best weighting to aggregate PAE matrices.

    Starts with equal weighting which is the same as taking an average of the PAE matrix.

    Args:
        size: length of the flattened PAE matrix

    Attributes:
        size: expected length of the input PAE matrix
        weights: weights to use when summing the elements of the PAE matrix

    """

    def __init__(self, size: int) -> None:
        super().__init__()
        self.size = size

        initial_weight = 1 / self.size
        self.weights = nn.Parameter(torch.Tensor([initial_weight for _ in range(size)]))

    def forward(self, pae: torch.Tensor) -> torch.Tensor:
        """Forward pass of the model."""
        x = torch.mul(self.weights, pae)
        return torch.sum(x, dim=-1) / (torch.count_nonzero(pae, dim=-1) / self.size)

"""Models to predict T cell antigen specificity using predicted structures."""

import torch
from torch import nn


class PAEAggregator(nn.Module):
    """Model to learn the best weighting to aggregate PAE matrices.

    Starts with equal weighting which is the same as taking an average of the PAE matrix.

    Args:
        size: length of the flattened PAE matrix
        scaling_data_size: length of the scaling data (Default: 2)

    Attributes:
        size: expected length of the input PAE matrix
        weights: weights to use when summing the elements of the PAE matrix
        scaler: linear layer for assessing scaling factor for average

    """

    def __init__(self, size: int, scaling_data_size: int = 2) -> None:
        super().__init__()
        self.size = size

        initial_weight = 1 / self.size
        self.weights = nn.Parameter(torch.Tensor([initial_weight for _ in range(size)]))

        self.scaler = nn.Linear(scaling_data_size, 1)
        self.scaler.weight = nn.Parameter(torch.zeros(1, scaling_data_size))
        self.scaler.bias = nn.Parameter(torch.Tensor([1]))

    def forward(self, pae: torch.Tensor, scaling_data: torch.Tensor) -> torch.Tensor:
        """Forward pass of the model."""
        # Average
        x = torch.mul(self.weights, pae)
        x = torch.sum(x, dim=-1) / (torch.count_nonzero(pae, dim=-1) / self.size)

        # Scale
        scale = self.scaler(scaling_data)
        scale = torch.squeeze(scale, dim=-1)

        x = torch.mul(scale, x)

        # ReLu1
        x = torch.where(torch.le(x, 0), 0, x)
        x = torch.where(torch.ge(x, 1), 1, x)

        return x

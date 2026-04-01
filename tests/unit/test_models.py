from unittest import TestCase

import numpy as np
import torch

from tcr_pmhc_predicted_structures.models import PAEAggregator


class TestPAEAggregator(TestCase):
    def test(self):
        model = PAEAggregator(5)
        output = model.forward(torch.Tensor([[0.16666666666, 0.3333333333, 0.0, 0.3333333333, 1.0]]))

        np.testing.assert_array_almost_equal(output.cpu().detach().numpy(), np.array([0.45833333333]))

    def test_one(self):
        model = PAEAggregator(5)
        output = model.forward(torch.Tensor([[1.0, 0.0, 0.0, 1.0, 1.0]]))

        np.testing.assert_array_almost_equal(output.cpu().detach().numpy(), np.array([1.0]))

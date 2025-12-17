from unittest import TestCase

import numpy as np
import torch

from tcr_pmhc_predicted_structures.models import PAEAggregator


class TestPAEAggregator(TestCase):
    def test(self):
        model = PAEAggregator(5)
        output = model.forward(torch.Tensor([[1.0, 2.0, np.nan, 2.0, 6.0]]))

        np.testing.assert_array_almost_equal(output.cpu().detach().numpy(), np.array([2.75]))

    def test_one(self):
        model = PAEAggregator(5)
        output = model.forward(torch.Tensor([[1.0, np.nan, np.nan, 1.0, 1.0]]))

        np.testing.assert_array_almost_equal(output.cpu().detach().numpy(), np.array([1.0]))

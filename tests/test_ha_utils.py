import pytest
import sys
import types

import numpy as np
torch = pytest.importorskip("torch")

from phylofoundry.methods.ha_sites import _build_combined_attention_tensor, _compute_for_sequence, _received_attention_by_layer


def _install_mock_pwlf(break_frac=0.5, slopes=(-1.0, 1.0)):
    module = types.ModuleType("pwlf")

    class PiecewiseLinFit:
        def __init__(self, x, y):
            self.slopes = slopes

        def fit(self, n_segments):
            return [0.0, break_frac, 1.0]

    module.PiecewiseLinFit = PiecewiseLinFit
    sys.modules["pwlf"] = module


def test_combined_attention_shape_and_default_mean():
    _install_mock_pwlf(break_frac=0.5)
    att = torch.zeros((1, 2, 2, 6, 6), dtype=torch.float32)
    att[:, :, 0, :, :] = 1.0
    att[:, :, 1, :, :] = 5.0
    combined = _build_combined_attention_tensor(att)
    assert tuple(combined.shape) == (2, 1, 2, 6, 6)

    received = _received_attention_by_layer(combined, pooling_idx=0)
    assert received.shape == (2, 4)  # BOS/EOS removed

    heat, mask, info = _compute_for_sequence(
        "AAAA",
        combined,
        {"mode": "loc", "call_mode": "loc_break", "pooling_used": "mean", "loc_break_adjust": -1, "loc_theta_target_deg": 90},
    )
    assert heat.shape == (2, 4)
    assert info["n_ha"] == int(np.floor(0.5 * 4) - 1 if int(np.floor(0.5 * 4) - 1) > 0 else 1)
    assert int(mask.sum()) == info["n_ha"]


def test_bos_eos_slicing_respected():
    # Inner region low attention, BOS/EOS very high. Correct slicing should ignore BOS/EOS.
    att = torch.zeros((1, 1, 1, 7, 7), dtype=torch.float32)
    att[..., 0, :] = 100.0
    att[..., -1, :] = 100.0
    att[..., :, 0] = 100.0
    att[..., :, -1] = 100.0
    att[..., 1:-1, 1:-1] = 1.0
    combined = _build_combined_attention_tensor(att)
    received = _received_attention_by_layer(combined, pooling_idx=0)
    assert np.allclose(received[0], np.ones(5))

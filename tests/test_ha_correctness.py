import sys
import types

import numpy as np
import pandas as pd
import pytest

torch = pytest.importorskip("torch")

from phylofoundry.methods.ha_sites import _build_combined_attention_tensor, _compute_for_sequence, compute_ha_sites_for_hmm


def _install_mock_pwlf(break_frac=0.6, slopes=(-1.0, 1.0)):
    module = types.ModuleType("pwlf")

    class PiecewiseLinFit:
        def __init__(self, x, y):
            self.slopes = slopes

        def fit(self, n_segments):
            return [0.0, break_frac, 1.0]

    module.PiecewiseLinFit = PiecewiseLinFit
    sys.modules["pwlf"] = module


def test_ha_attention_processing_and_layer_ids(tmp_path):
    _install_mock_pwlf()
    # raw attentions: [batch, layers, heads, T, T]
    raw = torch.zeros((1, 3, 1, 6, 6), dtype=torch.float32)
    raw[:, :, :, 1:-1, 1:-1] = 1.0
    raw[:, 2, :, 1:-1, 1:-1] = 5.0
    combined = _build_combined_attention_tensor(raw)

    heatmap, mask, info = _compute_for_sequence(
        "ACDE",
        combined,
        {
            "mode": "loc",
            "call_mode": "loc_break",
            "pooling_used": "mean",
            "loc_break_adjust": 0,
            "loc_theta_target_deg": 90,
            "layers": [1, 2],
        },
    )

    assert heatmap.shape == (2, 4)  # BOS/EOS stripped with 1:-1
    assert np.allclose(heatmap[0], np.ones(4))
    assert info["loc_layer_idx"] in (0, 1)
    assert info["loc_layer_id"] == info["layers_used"][info["loc_layer_idx"]]
    assert mask.sum() > 0

    fasta = tmp_path / "HMMX.faa"
    fasta.write_text(">s1\nACDE\n")
    aln = tmp_path / "HMMX.aln.faa"
    aln.write_text(">s1\nA-CD-E\n")

    def artifacts(_seq_id, _seq):
        return combined

    cfg = {
        "ha": {"mode": "loc", "call_mode": "loc_break", "pooling_used": "mean", "loc_break_adjust": 0, "loc_theta_target_deg": 90, "layers": [1, 2]},
        "embeddings": {"model": "mock", "device": "cpu"},
    }
    out = tmp_path / "ha" / "HMMX"
    res = compute_ha_sites_for_hmm("HMMX", fasta, out, cfg, artifacts, alignment_path=aln)
    df = pd.read_csv(res["ha_sites"], sep="\t")
    assert not df.empty
    assert list(df["ungapped_pos"]) == [1, 2, 3, 4]
    assert list(df["msa_col"]) == [1, 3, 4, 6]

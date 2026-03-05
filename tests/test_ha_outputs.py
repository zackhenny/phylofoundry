import pytest
from pathlib import Path

import pandas as pd
torch = pytest.importorskip("torch")

from phylofoundry.methods.ha_sites import _build_combined_attention_tensor, compute_ha_sites_for_hmm
from phylofoundry.tasks.discover import discover_motifs


def test_ha_stage_outputs_and_discover_consumption(tmp_path):
    fasta = tmp_path / "HMM1.faa"
    fasta.write_text(">s1\nACDE\n>s2\nACDF\n")

    def artifacts(seq_id, seq):
        del seq_id, seq
        raw = torch.ones((1, 2, 2, 1 + 4 + 1, 1 + 4 + 1))
        return _build_combined_attention_tensor(raw)

    cfg = {
        "ha": {"mode": "middle", "pooling_used": "mean", "call_mode": "topk", "topk": 2},
        "embeddings": {"model": "mock", "device": "cpu"},
        "workflow": {"max_failure_rate": 1.0},
    }
    out = tmp_path / "ha" / "HMM1"
    res = compute_ha_sites_for_hmm("HMM1", fasta, out, cfg, artifacts)

    assert res["ha_sites"].exists()
    assert res["ha_counts"].exists()
    assert res["loc_layers"].exists()

    summary_dir = tmp_path / "summary"
    summary_dir.mkdir()
    pd.read_csv(res["ha_sites"], sep="\t").to_csv(summary_dir / "ha_sites.tsv", sep="\t", index=False)

    disc_cfg = {"discover": {"enabled": True, "use_ha": True}}
    out_df = discover_motifs(disc_cfg, str(tmp_path), str(summary_dir), hmm_keep=None, force=True)
    assert not out_df.empty

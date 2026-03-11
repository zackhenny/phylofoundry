import numpy as np
import pandas as pd
import pytest

from phylofoundry.utils.ha import (
    call_ha_sites,
    compute_loc_and_ha,
    compute_ha_from_layer_profiles,
    select_layer_range,
    ungapped_to_msa_column,
)
from phylofoundry.tasks.discover_motifs import _attention_received_by_layer, _compute_candidate_residues_for_hmm, _compute_ha_enrichment_for_hmm


def test_select_layer_range_middle_and_range():
    s, e = select_layer_range(12, mode="middle")
    assert (s, e) == (4, 8)
    s2, e2 = select_layer_range(12, mode="range", start=2, end=6)
    assert (s2, e2) == (2, 6)


def test_call_ha_sites_percentile_topk_and_caps():
    scores = np.array([0.9, 0.9, 0.3, 0.1, 0.05])
    m = call_ha_sites(scores, call_mode="percentile", percentile=0.2, min_sites=2, max_sites=3)
    assert m.sum() == 2
    assert m[0] == 1 and m[1] == 1  # stable tie-break lower index first

    m2 = call_ha_sites(scores, call_mode="topk", topk=4, min_sites=1, max_sites=3)
    assert m2.sum() == 3


def test_compute_ha_from_layer_profiles_middle_mode():
    # 6 layers x 5 residues
    layer = np.array(
        [
            [1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1],
            [10, 2, 1, 1, 1],
            [9, 2, 1, 1, 1],
            [1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1],
        ],
        dtype=float,
    )
    cfg = {"layer_mode": "middle", "call_mode": "topk", "topk": 1, "min_sites": 1, "max_sites": 2}
    scores, mask, s, e, loc_layer = compute_ha_from_layer_profiles(layer, cfg)
    assert (s, e) == (2, 4)
    assert loc_layer is None
    assert int(mask.sum()) == 1
    assert int(np.argmax(scores)) == 0


def test_loc_mode_selects_layer_and_mask_size():
    pytest.importorskip("pwlf")
    layer = np.array(
        [
            [0.22, 0.2, 0.19, 0.18, 0.17, 0.16],
            [1.0, 0.95, 0.9, 0.12, 0.1, 0.09],
            [0.21, 0.2, 0.2, 0.19, 0.18, 0.17],
        ],
        dtype=float,
    )
    loc_layer, scores, mask = compute_loc_and_ha(layer, {"loc_norm_mode": "max", "loc_break_adjust": 0})
    assert loc_layer == 1
    assert np.argmax(scores) == 0
    assert int(mask.sum()) >= 2


def test_ungapped_to_msa_column_mapping():
    m = ungapped_to_msa_column("A-CDE-")
    assert m == {1: 0, 2: 2, 3: 3, 4: 4}


def test_ha_enrichment_outputs(tmp_path):
    attention_dir = tmp_path / "attention"
    discover_dir = tmp_path / "discover"
    attention_dir.mkdir(parents=True)

    pd.DataFrame(
        [
            {"seq_id": "s1", "pos_ungapped": 1, "is_ha": 1},
            {"seq_id": "s2", "pos_ungapped": 1, "is_ha": 1},
            {"seq_id": "s3", "pos_ungapped": 3, "is_ha": 1},
            {"seq_id": "s4", "pos_ungapped": 3, "is_ha": 1},
        ]
    ).to_csv(attention_dir / "HMMX.ha_sites.tsv", sep="\t", index=False)

    aln = {"s1": "ACD-", "s2": "ACD-", "s3": "A-DG", "s4": "A-DG"}
    clade_df = pd.DataFrame(
        [
            {"protein": "s1", "cluster_id": 0},
            {"protein": "s2", "cluster_id": 0},
            {"protein": "s3", "cluster_id": 1},
            {"protein": "s4", "cluster_id": 1},
        ]
    )
    _compute_ha_enrichment_for_hmm(
        "HMMX",
        aln,
        clade_df,
        ha_cfg={},
        disc_cfg={"ha_window": 3, "ha_delta_min": 0.0, "ha_gap_frac_max": 1.0},
        attention_dir=str(attention_dir),
        discover_dir=str(discover_dir),
    )

    assert (discover_dir / "HMMX.ha_enrichment.tsv").exists()
    assert (discover_dir / "HMMX.ha_hubs.tsv").exists()


def test_attention_received_by_layer_shape():
    torch = pytest.importorskip("torch")
    attentions = torch.rand((4, 8, 7, 7))
    out = _attention_received_by_layer(attentions, seq_len=5)
    assert out.shape == (4, 5)


def test_candidate_residue_scoring(tmp_path):
    attention_dir = tmp_path / "attention"
    discover_dir = tmp_path / "discover"
    attention_dir.mkdir(parents=True)

    rows = []
    for sid in ["a1", "a2", "a3"]:
        for pos in [2, 3]:
            rows.append({"seq_id": sid, "pos_ungapped": pos, "is_ha": 1})
    for sid in ["b1", "b2", "b3"]:
        rows.append({"seq_id": sid, "pos_ungapped": 1, "is_ha": 1})
    pd.DataFrame(rows).to_csv(attention_dir / "HMMY.ha_sites.tsv", sep="\t", index=False)

    aln = {
        "a1": "ACDE",
        "a2": "ACDE",
        "a3": "ACDE",
        "b1": "ASDE",
        "b2": "ASDE",
        "b3": "ASDE",
    }
    clade_df = pd.DataFrame(
        [
            {"protein": "a1", "cluster_id": 0},
            {"protein": "a2", "cluster_id": 0},
            {"protein": "a3", "cluster_id": 0},
            {"protein": "b1", "cluster_id": 1},
            {"protein": "b2", "cluster_id": 1},
            {"protein": "b3", "cluster_id": 1},
        ]
    )
    _compute_candidate_residues_for_hmm(
        "HMMY",
        aln,
        clade_df,
        {"candidates": {"enabled": True, "min_delta_ha": 0.1, "min_cons_frac": 0.6, "max_gap_frac": 1.0, "top_n": 20}},
        str(attention_dir),
        str(discover_dir),
    )

    out = pd.read_csv(discover_dir / "HMMY.candidate_residues.tsv", sep="\t")
    top = out.sort_values("candidate_score", ascending=False).iloc[0]
    assert int(top["msa_col"]) == 1
    assert int(top["aa_shift"]) == 1

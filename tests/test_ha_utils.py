import numpy as np
import pandas as pd

from phylofoundry.utils.ha import (
    call_ha_sites,
    compute_ha_from_layer_profiles,
    select_layer_range,
    ungapped_to_msa_column,
)
from phylofoundry.tasks.discover import _compute_ha_enrichment_for_hmm


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
    scores, mask, s, e = compute_ha_from_layer_profiles(layer, cfg)
    assert (s, e) == (2, 4)
    assert int(mask.sum()) == 1
    assert int(np.argmax(scores)) == 0


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

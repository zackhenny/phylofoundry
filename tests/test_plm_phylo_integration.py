import json
from pathlib import Path

import numpy as np
import pandas as pd

from phylofoundry.methods.evidence_join import run_evidence_join
from phylofoundry.methods.ha_sites import run_ha_sites
from phylofoundry.methods.regime_shift import run_regime_shift
from phylofoundry.config import validate_config


def test_regime_shift_and_ha_and_evidence(tmp_path: Path):
    out = tmp_path / "out"
    fasta = out / "fasta_per_hmm"
    emb = out / "embeddings"
    trees = out / "trees_iqtree"
    summary = out / "summary"
    qc = out / "qc"
    for d in [fasta, emb, trees, summary, qc]:
        d.mkdir(parents=True, exist_ok=True)

    hmm = "toy"
    seq_ids = [f"g{i}|p{i}" for i in range(6)]
    seqs = ["AAAAAMMMMM", "AAAAAMMMMM", "AAAAAMMMMM", "VVVVVMMMMM", "VVVVVMMMMM", "VVVVVMMMMM"]
    with open(fasta / f"{hmm}.faa", "w") as f:
        for sid, s in zip(seq_ids, seqs):
            f.write(f">{sid}\n{s}\n")

    pca = pd.DataFrame({"hmm": hmm, "tip": seq_ids, "PC1": [0, 0.1, -0.1, 4, 4.1, 3.9], "PC2": [0, 0, 0, 0, 0, 0]})
    pca.to_csv(emb / f"{hmm}.pca.tsv", sep="\t", index=False)
    (trees / f"{hmm}.treefile").write_text("((g0|p0,g1|p1,g2|p2),(g3|p3,g4|p4,g5|p5));")

    cfg = {"regime_shift": {"metric": "centroid", "min_size": 2, "n_permutations": 10}, "ha": {"mode": "middle", "call_mode": "topk", "topk": 2}, "qc": {"dpi": 80}, "evidence_join": {"classification_thresholds": {"delta_ha": 0.1, "js": 0.1, "meme_p": 0.1}}}

    rs = run_regime_shift(cfg, str(trees), str(emb), str(summary))
    assert not rs.empty

    run_ha_sites(cfg, str(fasta), str(emb), str(summary), str(qc))
    ha = pd.read_csv(summary / "ha_sites.tsv", sep="\t")
    assert ha["is_ha"].sum() > 0

    pd.DataFrame([{"hmm": hmm, "msa_col": 1, "delta_ha": 0.5, "js_divergence": 0.4}]).to_csv(summary / "discover_candidates.tsv", sep="\t", index=False)
    pd.DataFrame([{"hmm_id": hmm, "msa_col": 1, "MEME_p": 0.01, "MEME_omega": 2.0}]).to_csv(summary / "hyphy_results_summary.tsv", sep="\t", index=False)
    ev = run_evidence_join(cfg, str(summary))
    assert "classification" in ev.columns
    assert (ev["classification"] == "adaptive_shift").any()


def test_config_validation_catches_errors():
    bad = {"inputs": {"faa_dir": None, "hmm_input": None}, "output": {"outdir": None}}
    try:
        validate_config(bad)
    except SystemExit as e:
        assert "must specify" in str(e)
    else:
        raise AssertionError("validate_config should fail")

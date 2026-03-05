import pandas as pd
import pytest

pytest.importorskip("matplotlib")
from phylofoundry.pipeline import run_pipeline


def test_pipeline_runs_ha_before_discover_with_canonical_artifacts(tmp_path, monkeypatch):
    data_dir = tmp_path / "inputs"
    data_dir.mkdir()
    (data_dir / "g1.faa").write_text(">g1|s1\nACDE\n")
    (data_dir / "HMM1.hmm").write_text("HMMER3/f\n")

    outdir = tmp_path / "out"
    order = []

    monkeypatch.setattr("phylofoundry.pipeline.prep.run_prep", lambda *a, **k: order.append("prep"))
    monkeypatch.setattr(
        "phylofoundry.pipeline.hmmer.run_hmmer",
        lambda *a, **k: (
            pd.DataFrame([]),
            pd.DataFrame([]),
            pd.DataFrame([{"protein": "g1|s1", "hmm": "HMM1"}]),
        ),
    )

    def _extract(*_args, **kwargs):
        order.append("extract")
        fasta_dir = outdir / "fasta_per_hmm"
        fasta_dir.mkdir(parents=True, exist_ok=True)
        (fasta_dir / "HMM1.faa").write_text(">g1|s1\nACDE\n")
        return {"HMM1": ["g1|s1"]}

    monkeypatch.setattr("phylofoundry.pipeline.extract.run_extract", _extract)
    monkeypatch.setattr("phylofoundry.pipeline.embed.run_embed", lambda *a, **k: order.append("embed"))
    monkeypatch.setattr("phylofoundry.pipeline.phylo.run_phylo", lambda *a, **k: order.append("phylo"))

    def _ha(*_args, **_kwargs):
        order.append("ha_sites")
        hmm_ha = outdir / "ha" / "HMM1"
        hmm_ha.mkdir(parents=True, exist_ok=True)
        df = pd.DataFrame(
            [
                {"hmm_id": "HMM1", "seq_id": "g1|s1", "msa_col": 1, "ha_score": 0.8, "is_ha": 1},
                {"hmm_id": "HMM1", "seq_id": "g1|s1", "msa_col": 2, "ha_score": 0.2, "is_ha": 0},
            ]
        )
        df.to_csv(hmm_ha / "ha_sites.tsv", sep="\t", index=False)
        summary_dir = outdir / "summary"
        summary_dir.mkdir(parents=True, exist_ok=True)
        df.to_csv(summary_dir / "ha_sites.tsv", sep="\t", index=False)

    monkeypatch.setattr("phylofoundry.pipeline.run_ha_sites", _ha)

    cfg = {
        "inputs": {"faa_dir": str(data_dir), "hmm_input": str(data_dir)},
        "output": {"outdir": str(outdir)},
        "workflow": {"start_at": None, "stop_after": None, "force": True, "hmm_manifest": None},
        "phylo": {"no_asr": True},
        "embeddings": {"enabled": True},
        "post": {"enabled": False},
        "synteny": {"enabled": False},
        "codon": {"enabled": False},
        "hyphy": {"enabled": False},
        "motifs": {"enabled": False},
        "discover": {"enabled": True, "use_ha": True},
        "regime_shift": {"enable": False},
        "ha": {"enabled": True},
        "evidence_join": {"enable": False},
        "qc": {"enable": False},
        "prep": {},
    }

    run_pipeline(cfg)

    status = pd.read_csv(outdir / "summary" / "step_status.tsv", sep="\t")
    steps = status["step"].tolist()
    assert steps.index("ha_sites") < steps.index("discover_motifs")
    candidates = pd.read_csv(outdir / "summary" / "discover_candidates.tsv", sep="\t")
    assert not candidates.empty
    assert (outdir / "ha" / "HMM1" / "ha_sites.tsv").exists()

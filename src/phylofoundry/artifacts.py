from __future__ import annotations

from pathlib import Path


def ha_hmm_dir(outdir: str | Path, hmm_id: str) -> Path:
    return Path(outdir) / "ha" / hmm_id


def ha_hmm_path(outdir: str | Path, hmm_id: str) -> Path:
    """Canonical per-HMM HA table path."""
    return ha_hmm_dir(outdir, hmm_id) / "ha_sites.tsv"


def ha_global_path(outdir: str | Path) -> Path:
    """Canonical global HA table path."""
    return Path(outdir) / "summary" / "ha_sites.tsv"


def alignment_path(outdir: str | Path, hmm_id: str) -> Path:
    """Canonical alignment FASTA path for one HMM (ClipKIT output)."""
    return Path(outdir) / "alignments_clipkit" / f"{hmm_id}.faa"


def discover_out_paths(outdir: str | Path, hmm_id: str | None = None) -> dict[str, Path]:
    summary_dir = Path(outdir) / "summary"
    out = {
        "candidates": summary_dir / "discover_candidates.tsv",
        "status": summary_dir / "discover_status.tsv",
    }
    if hmm_id is not None:
        out["per_hmm_candidates"] = Path(outdir) / "discover" / hmm_id / "candidate_residues.tsv"
    return out

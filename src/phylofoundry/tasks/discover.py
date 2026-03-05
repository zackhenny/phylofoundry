"""discover.py — motif discovery consuming precomputed optional HA artifacts."""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from ..artifacts import discover_out_paths, ha_global_path, ha_hmm_path

logger = logging.getLogger(__name__)


def _load_ha_for_hmm(outdir: str | Path, hmm_id: str) -> pd.DataFrame:
    canonical_fp = ha_hmm_path(outdir, hmm_id)
    if canonical_fp.exists():
        ha = pd.read_csv(canonical_fp, sep="\t")
        if "msa_col" not in ha.columns:
            raise RuntimeError(
                f"HA artifact {canonical_fp} missing msa_col. Re-run `phylofoundry ha --hmm {hmm_id}` with alignments available."
            )
        return ha

    global_fp = ha_global_path(outdir)
    if global_fp.exists():
        ha = pd.read_csv(global_fp, sep="\t")
        if {"hmm_id", "msa_col"}.issubset(ha.columns):
            sub = ha[ha["hmm_id"] == hmm_id].copy()
            if not sub.empty:
                return sub

    legacy_fp = Path(outdir) / "attention" / f"{hmm_id}.ha_sites.tsv"
    if legacy_fp.exists():
        logger.warning("Using deprecated HA artifact path: %s. Please migrate to %s", legacy_fp, canonical_fp)
        ha = pd.read_csv(legacy_fp, sep="\t")
        if "msa_col" not in ha.columns:
            raise RuntimeError(
                f"Legacy HA artifact {legacy_fp} missing msa_col; discovery requires MSA-mapped HA tables. Run `phylofoundry ha --hmm {hmm_id}`."
            )
        if "hmm_id" not in ha.columns:
            ha["hmm_id"] = hmm_id
        return ha

    raise RuntimeError(
        "discover.use_ha=true but HA artifacts are missing for "
        f"{hmm_id}. Expected {canonical_fp} (preferred), {global_fp} (filtered by hmm_id), "
        f"or legacy {legacy_fp}. Run the HA stage first: `phylofoundry ha --hmm {hmm_id}` "
        "or run the full pipeline with ha.enabled=true starting at or before ha_sites."
    )


def discover_motifs(cfg, fasta_dir, summary_dir, hmm_keep, force=False, clade_assign_dir=None):
    del clade_assign_dir
    disc_cfg = cfg.get("discover", {})
    if not disc_cfg.get("enabled", False):
        return None

    outdir = Path(summary_dir).parent
    out_fp = discover_out_paths(outdir)["candidates"]
    out_fp.parent.mkdir(parents=True, exist_ok=True)
    if out_fp.exists() and not force:
        return pd.read_csv(out_fp, sep="\t")

    fasta_hmms = {Path(fp).stem for fp in Path(fasta_dir).glob("*.faa")}
    target_hmms = sorted(hmm_keep & fasta_hmms if hmm_keep else fasta_hmms)

    use_ha = bool(disc_cfg.get("use_ha", False))
    if use_ha:
        if not target_hmms:
            raise RuntimeError(
                "discover.use_ha=true but no HMM FASTA inputs were found in "
                f"{Path(fasta_dir)}. Run extract/phylo first, then run HA and discovery."
            )
        rows = []
        per_hmm_status = []
        for hmm in target_hmms:
            ha = _load_ha_for_hmm(outdir, hmm)
            ha = ha.copy()
            ha["msa_col"] = pd.to_numeric(ha["msa_col"], errors="coerce")
            ha = ha.dropna(subset=["msa_col"])
            if ha.empty:
                raise RuntimeError(
                    f"HA table for {hmm} has no valid msa_col values. Re-run HA with alignments so msa_col can be mapped."
                )
            base = float(ha["ha_score"].mean())
            for col, sub in ha.groupby("msa_col"):
                freq = float(sub["is_ha"].mean())
                score = float(sub["ha_score"].mean())
                rows.append(
                    {
                        "hmm": hmm,
                        "msa_col": int(col),
                        "delta_ha": freq - 0.5,
                        "js_divergence": abs(score - base),
                        "candidate_score": (freq * 0.7) + (score * 0.3),
                    }
                )
            per_hmm_status.append({"hmm_id": hmm, "used_ha": True, "ha_rows": int(len(ha))})
        out = pd.DataFrame(rows).sort_values("candidate_score", ascending=False)
        pd.DataFrame(per_hmm_status).to_csv(discover_out_paths(outdir)["status"], sep="\t", index=False)
    else:
        out = pd.DataFrame(columns=["hmm", "msa_col", "delta_ha", "js_divergence", "candidate_score"])

    out.to_csv(out_fp, sep="\t", index=False)
    return out

"""discover.py — motif discovery consuming precomputed optional HA artifacts."""

from __future__ import annotations

import os
import pandas as pd


def discover_motifs(cfg, fasta_dir, summary_dir, hmm_keep, force=False, clade_assign_dir=None):
    del fasta_dir, hmm_keep, clade_assign_dir
    disc_cfg = cfg.get("discover", {})
    if not disc_cfg.get("enabled", False):
        return None

    out_fp = os.path.join(summary_dir, "discover_candidates.tsv")
    if os.path.exists(out_fp) and not force:
        return pd.read_csv(out_fp, sep="\t")

    use_ha = bool(disc_cfg.get("use_ha", False))
    if use_ha:
        ha_fp = os.path.join(summary_dir, "ha_sites.tsv")
        if not os.path.exists(ha_fp):
            raise SystemExit("HA artifacts missing; run `phylofoundry ha` or enable ha.enabled.")
        ha = pd.read_csv(ha_fp, sep="\t")
        if ha.empty:
            raise SystemExit("HA artifacts missing; run `phylofoundry ha` or enable ha.enabled.")
        base = float(ha["ha_score"].mean())
        rows = []
        for (hmm, col), sub in ha.groupby(["hmm_id", "msa_col"]):
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
        out = pd.DataFrame(rows).sort_values("candidate_score", ascending=False)
    else:
        # HA-independent placeholder output for decoupled discovery runs.
        out = pd.DataFrame(columns=["hmm", "msa_col", "delta_ha", "js_divergence", "candidate_score"])

    out.to_csv(out_fp, sep="\t", index=False)
    return out

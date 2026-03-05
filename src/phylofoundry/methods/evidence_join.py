from __future__ import annotations

import json
import os

import numpy as np
import pandas as pd


def _classify(row, th):
    meme = row.get("MEME_p", np.nan)
    delta = row.get("delta_ha", 0.0)
    js = row.get("js_divergence", 0.0)
    if pd.notna(meme) and meme <= th.get("meme_p", 0.1) and delta >= th.get("delta_ha", 0.2):
        return "adaptive_shift"
    if delta >= th.get("delta_ha", 0.2) or js >= th.get("js", 0.15):
        return "functional_candidate"
    if pd.notna(meme) and meme > 0.1:
        return "drift"
    return "unsupported"


def enrichment_test(selected_mask: np.ndarray, ha_mask: np.ndarray):
    try:
        from scipy.stats import fisher_exact

        a = int(np.sum((selected_mask == 1) & (ha_mask == 1)))
        b = int(np.sum((selected_mask == 1) & (ha_mask == 0)))
        c = int(np.sum((selected_mask == 0) & (ha_mask == 1)))
        d = int(np.sum((selected_mask == 0) & (ha_mask == 0)))
        _, p = fisher_exact([[a, b], [c, d]], alternative="greater")
        return p
    except Exception:
        return np.nan


def correlation_ha_selection(ha: np.ndarray, stat: np.ndarray):
    try:
        from scipy.stats import spearmanr

        return float(spearmanr(ha, stat, nan_policy="omit").correlation)
    except Exception:
        return np.nan


def run_evidence_join(cfg: dict, summary_dir: str):
    if not os.path.exists(os.path.join(summary_dir, "ha_sites.tsv")):
        raise SystemExit(f"EvidenceJoin missing required file: {os.path.join(summary_dir, 'ha_sites.tsv')}")

    ha = pd.read_csv(os.path.join(summary_dir, "ha_sites.tsv"), sep="\t")
    by_col = ha.groupby(["hmm_id", "msa_col"]).agg(
        ha_freq_family=("is_ha", "mean"), ha_score_family_mean=("ha_score", "mean")
    ).reset_index()

    cand_fp = os.path.join(summary_dir, "discover_candidates.tsv")
    cand = pd.read_csv(cand_fp, sep="\t") if os.path.exists(cand_fp) else pd.DataFrame(columns=["hmm", "msa_col", "delta_ha", "js_divergence"])
    cand = cand.rename(columns={"hmm": "hmm_id"})

    meme_fp = os.path.join(summary_dir, "hyphy_results_summary.tsv")
    hy = pd.read_csv(meme_fp, sep="\t") if os.path.exists(meme_fp) else pd.DataFrame(columns=["hmm_id", "msa_col", "MEME_p", "MEME_omega"])

    out = by_col.merge(cand, on=["hmm_id", "msa_col"], how="left").merge(hy, on=["hmm_id", "msa_col"], how="left")
    if "delta_ha" not in out.columns:
        out["delta_ha"] = 0.0
    else:
        out["delta_ha"] = out["delta_ha"].fillna(0.0)
    if "js_divergence" not in out.columns:
        out["js_divergence"] = 0.0
    else:
        out["js_divergence"] = out["js_divergence"].fillna(0.0)

    th = cfg.get("evidence_join", {}).get("classification_thresholds", {})
    out["classification"] = out.apply(lambda r: _classify(r, th), axis=1)
    out["confidence_grade"] = np.where(out["classification"] == "adaptive_shift", "A", np.where(out["classification"] == "functional_candidate", "B", "C"))
    out["notes"] = "PLM-derived signals are hypotheses; phylogenetic tests provide mechanistic evidence"
    out.to_csv(os.path.join(summary_dir, "site_evidence.tsv"), sep="\t", index=False)

    summary = {
        "n_sites": int(len(out)),
        "n_adaptive_shift": int((out["classification"] == "adaptive_shift").sum()),
        "n_functional_candidate": int((out["classification"] == "functional_candidate").sum()),
        "top_sites": out.sort_values(["classification", "ha_score_family_mean"], ascending=[True, False]).head(20).to_dict("records"),
    }
    with open(os.path.join(summary_dir, "evidence_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    return out

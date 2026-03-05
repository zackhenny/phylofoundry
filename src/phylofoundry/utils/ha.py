"""High-Attention (HA) site utilities.

Paper-inspired HA scoring collapses transformer attention into per-residue
"attention received" vectors and calls top residues as HA sites.
"""

from __future__ import annotations

import math
import os
from typing import Dict, Iterable

import numpy as np
import pandas as pd


def select_layer_range(n_layers: int, mode: str = "middle", start=None, end=None) -> tuple[int, int]:
    """Select [start, end) layer interval deterministically."""
    if n_layers <= 0:
        raise ValueError("n_layers must be > 0")

    mode = str(mode or "middle").strip().lower()
    if mode == "range":
        if start is None or end is None:
            raise ValueError("ha.layer_mode=range requires ha.layer_start and ha.layer_end")
        s = int(start)
        e = int(end)
        if s < 0 or e > n_layers or s >= e:
            raise ValueError(f"Invalid layer range [{s}, {e}) for model with {n_layers} layers")
        return s, e

    # middle-third [floor(N/3), floor(2N/3)); ensure at least one layer
    s = int(math.floor(n_layers / 3))
    e = int(math.floor(2 * n_layers / 3))
    if e <= s:
        e = min(n_layers, s + 1)
    return s, e


def _normalize_layer_vector(v: np.ndarray) -> np.ndarray:
    tot = float(np.sum(v))
    if tot <= 0:
        return np.zeros_like(v, dtype=float)
    return v / tot


def normalize_layer_vector(v: np.ndarray, mode: str = "sum") -> np.ndarray:
    """Normalize a per-position attention vector using sum or max scaling."""
    vals = np.asarray(v, dtype=float)
    mode = str(mode or "sum").strip().lower()
    if vals.size == 0:
        return vals.astype(float)
    if mode == "max":
        den = float(np.max(vals))
    else:
        den = float(np.sum(vals))
    if den <= 0:
        return np.zeros_like(vals, dtype=float)
    return vals / den


def loc_breakpoint_pwlf(sorted_scores: np.ndarray, segments: int = 2):
    """Fit pwlf on sorted scores and return breakpoint (rank index, 1-based) and model."""
    y = np.asarray(sorted_scores, dtype=float)
    n = len(y)
    if n == 0:
        return 0, None
    if n == 1:
        return 1, None

    if segments != 2:
        raise ValueError("LoC currently supports exactly 2 PWLF segments")

    try:
        import pwlf
    except ImportError as e:
        raise ImportError("pwlf is required for ha.method='loc'. Install with extras [attention].") from e

    x = np.linspace(0.0, 1.0, n)
    model = pwlf.PiecewiseLinFit(x, y)
    breaks = model.fit(segments)
    frac = float(breaks[1]) if len(breaks) > 1 else 1.0
    k = int(np.floor(frac * n))
    k = max(1, min(k, n))
    return k, model


def loc_theta_metric(theta_target_deg: float, pwlf_model) -> float:
    """Return |theta - target| where theta is computed from pwlf slopes."""
    if pwlf_model is None:
        return float("inf")
    slopes = np.asarray(getattr(pwlf_model, "slopes", []), dtype=float)
    if slopes.size < 2:
        return float("inf")
    m1, m2 = float(slopes[0]), float(slopes[1])
    denom = (1.0 + (m1 * m2))
    if abs(denom) < 1e-12:
        theta = 90.0
    else:
        theta = abs(np.degrees(np.arctan((m2 - m1) / denom)))
    return abs(theta - float(theta_target_deg))


def compute_loc_and_ha(layer_by_pos: np.ndarray, cfg: dict) -> tuple[int, np.ndarray, np.ndarray]:
    """Compute convergence layer (LoC) and HA mask from per-layer attention received."""
    if layer_by_pos.ndim != 2:
        raise ValueError("layer_by_pos must be 2D: (n_layers, seq_len)")
    n_layers, seq_len = layer_by_pos.shape
    if seq_len == 0:
        return 0, np.zeros(0, dtype=float), np.zeros(0, dtype=np.uint8)

    norm_mode = cfg.get("loc_norm_mode", "max")
    theta_target = float(cfg.get("loc_theta_target_deg", 90))
    segments = int(cfg.get("loc_pwlf_segments", 2))
    break_adjust = int(cfg.get("loc_break_adjust", -1))

    best_idx = 0
    best_dist = float("inf")
    best_norm = normalize_layer_vector(layer_by_pos[0], mode=norm_mode)
    best_order = np.argsort(-best_norm)
    best_k = 1

    for layer_idx in range(n_layers):
        norm_v = normalize_layer_vector(layer_by_pos[layer_idx], mode=norm_mode)
        order = np.argsort(-norm_v)
        sorted_scores = norm_v[order]
        k, model = loc_breakpoint_pwlf(sorted_scores, segments=segments)
        dist = loc_theta_metric(theta_target, model)
        if dist < best_dist:
            best_idx = layer_idx
            best_dist = dist
            best_norm = norm_v
            best_order = order
            best_k = k

    cutoff = int(best_k + break_adjust)
    cutoff = max(1, min(cutoff, seq_len))
    mask = np.zeros(seq_len, dtype=np.uint8)
    mask[best_order[:cutoff]] = 1
    return best_idx, best_norm, mask


def aggregate_attention_received(layer_by_pos: np.ndarray, agg: str = "median") -> np.ndarray:
    """Aggregate per-layer per-position attention-received into one score vector.

    Parameters
    ----------
    layer_by_pos : np.ndarray
        Shape (n_layers, seq_len) values of attention-received.
    """
    if layer_by_pos.ndim != 2:
        raise ValueError("layer_by_pos must be 2D: (n_layers, seq_len)")
    normalized = np.vstack([_normalize_layer_vector(v.astype(float)) for v in layer_by_pos])
    agg = str(agg or "median").lower()
    if agg == "mean":
        return normalized.mean(axis=0)
    return np.median(normalized, axis=0)


def call_ha_sites(
    scores: np.ndarray,
    call_mode: str = "percentile",
    percentile: float = 0.05,
    topk: int = 20,
    min_sites: int = 8,
    max_sites: int = 60,
) -> np.ndarray:
    """Call HA positions with deterministic tie-breaking by residue index."""
    n = int(len(scores))
    if n == 0:
        return np.zeros(0, dtype=np.uint8)

    call_mode = str(call_mode or "percentile").lower()
    if call_mode == "topk":
        k = int(topk)
    else:
        frac = max(0.0, min(1.0, float(percentile)))
        k = int(math.ceil(n * frac))

    k = max(int(min_sites), k)
    if max_sites is not None:
        k = min(k, int(max_sites))
    k = max(1, min(k, n))

    # stable deterministic order: score desc, index asc
    order = sorted(range(n), key=lambda i: (-float(scores[i]), i))
    mask = np.zeros(n, dtype=np.uint8)
    for i in order[:k]:
        mask[i] = 1
    return mask


def compute_ha_from_layer_profiles(
    layer_by_pos: np.ndarray, cfg: dict
) -> tuple[np.ndarray, np.ndarray, int, int, int | None]:
    """Compute HA scores and mask from layer-by-position attention vectors."""
    if layer_by_pos.ndim != 2:
        raise ValueError("layer_by_pos must be 2D: (n_layers, seq_len)")
    n_layers = layer_by_pos.shape[0]
    method = str(cfg.get("method", "middle")).strip().lower()
    if method == "loc":
        loc_layer, scores, mask = compute_loc_and_ha(layer_by_pos, cfg)
        return scores, mask, loc_layer, loc_layer + 1, loc_layer

    s, e = select_layer_range(
        n_layers,
        mode=cfg.get("layer_mode", "middle"),
        start=cfg.get("layer_start"),
        end=cfg.get("layer_end"),
    )
    agg_scores = aggregate_attention_received(layer_by_pos[s:e, :], agg=cfg.get("agg", "median"))
    mask = call_ha_sites(
        agg_scores,
        call_mode=cfg.get("call_mode", "percentile"),
        percentile=cfg.get("percentile", 0.05),
        topk=cfg.get("topk", 20),
        min_sites=cfg.get("min_sites", 8),
        max_sites=cfg.get("max_sites", 60),
    )
    return agg_scores, mask, s, e, None


def ungapped_to_msa_column(aln_seq: str) -> dict[int, int]:
    """Map 1-based ungapped residue index to 0-based MSA column."""
    m = {}
    ungapped = 0
    for col, aa in enumerate(str(aln_seq)):
        if aa != "-":
            ungapped += 1
            m[ungapped] = col
    return m


def write_ha_outputs(
    attention_dir: str,
    summary_dir: str,
    hmm: str,
    seqs: Dict[str, str],
    per_seq_scores: Dict[str, np.ndarray],
    per_seq_mask: Dict[str, np.ndarray],
    layer_start: int,
    layer_end: int,
    cfg: dict,
    per_seq_loc_layer: Dict[str, int] | None = None,
):
    os.makedirs(attention_dir, exist_ok=True)

    rows = []
    summary_rows = []
    method = str(cfg.get("method", "middle")).strip().lower()
    norm_mode = cfg.get("loc_norm_mode") if method == "loc" else np.nan

    for seq_id, seq in seqs.items():
        scores = per_seq_scores.get(seq_id)
        mask = per_seq_mask.get(seq_id)
        if scores is None or mask is None:
            continue
        n_sites = int(mask.sum())
        loc_layer = np.nan
        if per_seq_loc_layer and seq_id in per_seq_loc_layer:
            loc_layer = int(per_seq_loc_layer[seq_id])
        summary_rows.append(
            {
                "hmm": hmm,
                "seq_id": seq_id,
                "n_sites": n_sites,
                "percentile_used": cfg.get("percentile", 0.05) if cfg.get("call_mode", "percentile") == "percentile" else np.nan,
                "layer_start": layer_start,
                "layer_end": layer_end,
                "method": method,
                "loc_layer": loc_layer,
                "norm_mode": norm_mode,
            }
        )
        for i, aa in enumerate(seq[: len(scores)], start=1):
            rows.append(
                {
                    "hmm": hmm,
                    "seq_id": seq_id,
                    "pos_ungapped": i,
                    "aa": aa,
                    "ha_score": float(scores[i - 1]),
                    "is_ha": int(mask[i - 1]),
                    "layer_start": layer_start,
                    "layer_end": layer_end,
                    "loc_layer": loc_layer,
                }
            )

    ha_fp = os.path.join(attention_dir, f"{hmm}.ha_sites.tsv")
    pd.DataFrame(rows).to_csv(ha_fp, sep="\t", index=False)

    # explicit summary outputs per HMM
    hmm_summary_dir = os.path.join(summary_dir, "summary") if os.path.basename(summary_dir) != "summary" else summary_dir
    os.makedirs(hmm_summary_dir, exist_ok=True)
    per_site_fp = os.path.join(hmm_summary_dir, f"{hmm}.ha_sites.tsv")
    per_count_fp = os.path.join(hmm_summary_dir, f"{hmm}.ha_counts.tsv")
    pd.DataFrame(rows).to_csv(per_site_fp, sep="\t", index=False)
    pd.DataFrame(summary_rows)[["hmm","seq_id","n_sites"]].rename(columns={"n_sites":"ha_sites"}).to_csv(per_count_fp, sep="\t", index=False)

    summary_fp = os.path.join(summary_dir, "ha_summary.tsv")
    new_df = pd.DataFrame(summary_rows)
    if os.path.exists(summary_fp):
        old = pd.read_csv(summary_fp, sep="\t")
        old = old[old["hmm"] != hmm]
        out = pd.concat([old, new_df], ignore_index=True)
    else:
        out = new_df
    out.to_csv(summary_fp, sep="\t", index=False)

    return ha_fp, summary_fp


def bh_fdr(pvals: Iterable[float]) -> np.ndarray:
    pvals = np.asarray(list(pvals), dtype=float)
    n = len(pvals)
    if n == 0:
        return np.array([], dtype=float)
    order = np.argsort(np.nan_to_num(pvals, nan=1.0))
    ranked = pvals[order]
    q = np.full(n, np.nan)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        p = ranked[i]
        rank = i + 1
        if np.isnan(p):
            val = np.nan
        else:
            val = min(prev, p * n / rank)
            prev = val
        q[i] = val
    out = np.full(n, np.nan)
    out[order] = q
    return out

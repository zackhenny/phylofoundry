from __future__ import annotations

import glob
import logging
import math
import os
import traceback
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..utils.bio import read_fasta


logger = logging.getLogger(__name__)
_ESM_CACHE: dict[tuple[str, str], tuple[Any, Any, Any]] = {}
_RESULT_KEYS_LOGGED = False


def _build_combined_attention_tensor(attentions):
    """Return tensor shaped like reference [2, B, layers, T, T]."""
    import torch

    if attentions.ndim != 5:
        raise ValueError(f"Expected attention tensor with ndim=5; got shape={tuple(attentions.shape)}")
    mean_pool = attentions.mean(dim=2, keepdim=False)
    max_pool = attentions.max(dim=2).values
    return torch.stack([mean_pool, max_pool], dim=0)


def _received_attention_by_layer(combined_attn, pooling_idx: int = 0) -> np.ndarray:
    """Reference behavior: data[pool,0,layer,1:-1,1:-1], sum over rows(dim=0), max-normalize."""
    import torch

    data = combined_attn
    if isinstance(data, np.ndarray):
        data = torch.from_numpy(data)
    if data.ndim != 5:
        raise ValueError(f"combined_attn_tensor must be 5D [2,B,L,T,T], got {tuple(data.shape)}")

    n_layers = int(data.shape[2])
    out = []
    for layer in range(n_layers):
        attn_matrix = data[pooling_idx, 0, layer, 1:-1, 1:-1]
        col_list = torch.sum(attn_matrix, dim=0)
        maxv = float(torch.max(col_list).item()) if col_list.numel() else 0.0
        if maxv > 0:
            col_norm = col_list / maxv
        else:
            col_norm = torch.zeros_like(col_list)
        out.append(col_norm.detach().cpu().numpy())
    return np.stack(out, axis=0) if out else np.zeros((0, 0), dtype=float)


def _theta_from_slopes(m1: float, m2: float) -> float:
    denom = 1.0 + (m1 * m2)
    if abs(denom) < 1e-12:
        return 90.0
    return float(np.degrees(np.arctan((m2 - m1) / denom)))


def _fit_loc_break(sorted_scores: np.ndarray) -> tuple[float, float]:
    import pwlf

    n = len(sorted_scores)
    x = np.linspace(0, 1, n)
    model = pwlf.PiecewiseLinFit(x, sorted_scores)
    breaks = model.fit(2)
    break_frac = float(breaks[1])
    slopes = np.asarray(model.slopes, dtype=float)
    theta = abs(_theta_from_slopes(float(slopes[0]), float(slopes[1]))) if slopes.size >= 2 else 0.0
    return break_frac, theta


def _select_loc_layer(received_layer_pos: np.ndarray, theta_target: float, loc_break_adjust: int) -> dict[str, Any]:
    l_count, seq_len = received_layer_pos.shape
    best = None
    for layer in range(l_count):
        v = received_layer_pos[layer]
        order = np.argsort(-v)
        sorted_vals = v[order]
        break_frac, theta = _fit_loc_break(sorted_vals)
        theta_delta = abs(theta - theta_target)
        n_ha = int(math.floor(break_frac * seq_len) + int(loc_break_adjust))
        n_ha = max(1, min(n_ha, seq_len))
        cand = {
            "layer": int(layer),
            "theta_deg": float(theta),
            "theta_delta": float(theta_delta),
            "break_frac": float(break_frac),
            "n_ha": int(n_ha),
            "scores": v,
            "order": order,
        }
        if best is None or cand["theta_delta"] < best["theta_delta"]:
            best = cand
    assert best is not None
    return best


def _call_sites(scores: np.ndarray, mode: str, threshold: float | None = None, topk: int | None = None) -> np.ndarray:
    n = len(scores)
    order = np.argsort(-scores)
    mask = np.zeros(n, dtype=np.uint8)
    if mode == "percentile":
        thr = float(threshold if threshold is not None else 0.95)
        q = float(np.quantile(scores, thr)) if n else 0.0
        mask[scores >= q] = 1
    elif mode == "topk":
        k = max(1, min(int(topk or 1), n))
        mask[order[:k]] = 1
    else:
        raise ValueError(f"Unsupported call_mode='{mode}'")
    return mask


def _compute_for_sequence(seq: str, combined_attn, cfg: dict[str, Any]) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    pool = str(cfg.get("pooling_used", "mean")).lower()
    pool_idx = 0 if pool == "mean" else 1
    received = _received_attention_by_layer(combined_attn, pooling_idx=pool_idx)

    layer_filter = cfg.get("layers")
    if layer_filter:
        received = received[layer_filter, :]
        layer_lookup = list(layer_filter)
    else:
        layer_lookup = list(range(received.shape[0]))

    mode = str(cfg.get("mode", "loc")).lower()
    call_mode = str(cfg.get("call_mode", "loc_break")).lower()

    if mode == "loc":
        best = _select_loc_layer(
            received,
            theta_target=float(cfg.get("loc_theta_target_deg", 90)),
            loc_break_adjust=int(cfg.get("loc_break_adjust", -1)),
        )
        loc_layer_idx = int(best["layer"])
        loc_layer_id = int(layer_lookup[loc_layer_idx])
        scores = best["scores"]
        if call_mode == "loc_break":
            mask = np.zeros(len(scores), dtype=np.uint8)
            mask[best["order"][: best["n_ha"]]] = 1
            threshold = best["break_frac"]
        elif call_mode == "percentile":
            threshold = float(cfg.get("percentile", 0.95))
            mask = _call_sites(scores, "percentile", threshold=threshold)
        elif call_mode == "topk":
            k = int(cfg.get("topk", max(1, best["n_ha"])))
            threshold = float(k)
            mask = _call_sites(scores, "topk", topk=k)
        else:
            raise ValueError(f"Unsupported ha.call_mode={call_mode}")
        info = {
            "loc_layer_idx": loc_layer_idx,
            "loc_layer_id": loc_layer_id,
            "layers_used": layer_lookup,
            "theta_deg": best["theta_deg"],
            "break_frac": best["break_frac"],
            "n_ha": int(mask.sum()),
            "threshold": threshold,
        }
        return received, mask, info

    mid = received.shape[0] // 2
    scores = received[mid]
    if call_mode == "topk":
        k = int(cfg.get("topk", 20))
        mask = _call_sites(scores, "topk", topk=k)
        threshold = float(k)
    else:
        q = float(cfg.get("percentile", 0.95))
        mask = _call_sites(scores, "percentile", threshold=q)
        threshold = q
    info = {
        "loc_layer_idx": int(mid),
        "loc_layer_id": int(layer_lookup[mid]),
        "layers_used": layer_lookup,
        "theta_deg": np.nan,
        "break_frac": np.nan,
        "n_ha": int(mask.sum()),
        "threshold": threshold,
    }
    return received, mask, info


def _build_msa_pos_map(aligned_seq: str) -> dict[int, int]:
    pos_map: dict[int, int] = {}
    ungapped_pos = 0
    for msa_col, residue in enumerate(aligned_seq, start=1):
        if residue != "-":
            ungapped_pos += 1
            pos_map[ungapped_pos] = msa_col
    return pos_map


def compute_ha_sites_for_hmm(hmm_id, fasta_path, outdir, config, embeddings_artifacts, alignment_path: str | Path | None = None):
    """Compute HA artifacts for one HMM and write summary, tensor, and QC outputs."""
    ha_cfg = config.get("ha", {})
    seqs = read_fasta(str(fasta_path))
    hmm_out = Path(outdir)
    summary_dir = hmm_out / "summary"
    attn_dir = hmm_out / "ha" / "attn"
    heatmap_dir = hmm_out / "ha" / "heatmaps"
    summary_dir.mkdir(parents=True, exist_ok=True)
    attn_dir.mkdir(parents=True, exist_ok=True)
    heatmap_dir.mkdir(parents=True, exist_ok=True)

    msa_map_by_seq: dict[str, dict[int, int]] = {}
    if alignment_path and Path(alignment_path).exists():
        msa_map_by_seq = {sid: _build_msa_pos_map(s) for sid, s in read_fasta(str(alignment_path)).items()}
    elif alignment_path:
        logger.warning("HA alignment missing for hmm_id=%s at %s; msa_col will be blank", hmm_id, alignment_path)

    qc_rows, site_rows, count_rows, loc_rows = [], [], [], []
    failures: list[dict[str, str]] = []
    max_logged_failures = int(ha_cfg.get("max_logged_failures", 5))
    fail_threshold = float(ha_cfg.get("fail_fast_threshold", config.get("workflow", {}).get("max_failure_rate", 0.3)))

    for seq_id, seq in seqs.items():
        try:
            combined_attn = embeddings_artifacts(seq_id, seq)
            import torch

            attn_fp = attn_dir / f"{seq_id}.pt"
            torch.save(combined_attn, attn_fp)

            heatmap, mask, info = _compute_for_sequence(seq, combined_attn, ha_cfg)
            np.save(heatmap_dir / f"{seq_id}.npy", heatmap)

            seq_len = len(seq)
            count_rows.append({"seq_id": seq_id, "seq_len": seq_len, "n_ha": int(mask.sum()), "frac_ha": float(mask.mean()) if seq_len else 0.0})
            loc_rows.append(
                {
                    "seq_id": seq_id,
                    "loc_layer_idx": info["loc_layer_idx"],
                    "loc_layer_id": info["loc_layer_id"],
                    "theta_deg": info["theta_deg"],
                    "break_frac": info["break_frac"],
                    "n_ha": int(mask.sum()),
                }
            )
            msa_pos_map = msa_map_by_seq.get(seq_id, {})
            if alignment_path and not msa_pos_map:
                logger.warning("No MSA mapping for hmm_id=%s seq_id=%s; msa_col will be blank", hmm_id, seq_id)
            for idx, _aa in enumerate(seq, start=1):
                site_rows.append(
                    {
                        "hmm_id": hmm_id,
                        "seq_id": seq_id,
                        "loc_layer_id": info["loc_layer_id"],
                        "loc_layer_idx": info["loc_layer_idx"],
                        "msa_col": msa_pos_map.get(idx, ""),
                        "ungapped_pos": idx,
                        "ha_score": float(heatmap[info["loc_layer_idx"]][idx - 1]) if info["loc_layer_idx"] < heatmap.shape[0] else np.nan,
                        "is_ha": bool(mask[idx - 1]),
                        "call_mode": ha_cfg.get("call_mode", "loc_break"),
                        "threshold": info["threshold"],
                        "pooling_used": ha_cfg.get("pooling_used", "mean"),
                        "model": config.get("embeddings", {}).get("model", "esm2_t33_650M_UR50D"),
                        "device": config.get("embeddings", {}).get("device", "cpu"),
                    }
                )
            qc_rows.append({"seq_id": seq_id, "heatmap": heatmap, "loc_layer_id": info["loc_layer_id"], "loc_layer_idx": info["loc_layer_idx"], "n_ha": int(mask.sum())})
        except Exception as exc:  # no silent swallowing
            failures.append({"seq_id": seq_id, "error": repr(exc), "traceback": traceback.format_exc()})
            if len(failures) <= max_logged_failures:
                logger.error("[ha] ERROR seq_id=%s\n%s", seq_id, failures[-1]["traceback"])
            failure_rate = len(failures) / max(1, len(seqs))
            if failure_rate > fail_threshold:
                raise RuntimeError(f"HA failure rate {failure_rate:.2%} exceeded threshold {fail_threshold:.2%} for {hmm_id}") from exc

    n_total = len(seqs)
    n_fail = len(failures)
    n_ok = n_total - n_fail
    fail_frac = (n_fail / n_total) if n_total else 0.0
    qc_df = pd.DataFrame([{"hmm_id": hmm_id, "n_total": n_total, "n_ok": n_ok, "n_fail": n_fail, "fail_frac": fail_frac}])
    qc_df.to_csv(summary_dir / "ha_run_qc.tsv", sep="\t", index=False)

    pd.DataFrame(site_rows, columns=["hmm_id", "seq_id", "loc_layer_id", "loc_layer_idx", "msa_col", "ungapped_pos", "ha_score", "is_ha", "call_mode", "threshold", "pooling_used", "model", "device"]).to_csv(summary_dir / "ha_sites.tsv", sep="\t", index=False)
    pd.DataFrame(count_rows, columns=["seq_id", "seq_len", "n_ha", "frac_ha"]).to_csv(summary_dir / "ha_counts.tsv", sep="\t", index=False)
    pd.DataFrame(loc_rows, columns=["seq_id", "loc_layer_idx", "loc_layer_id", "theta_deg", "break_frac", "n_ha"]).to_csv(summary_dir / "loc_layers.tsv", sep="\t", index=False)

    return {
        "ha_sites": summary_dir / "ha_sites.tsv",
        "ha_counts": summary_dir / "ha_counts.tsv",
        "loc_layers": summary_dir / "loc_layers.tsv",
        "ha_run_qc": summary_dir / "ha_run_qc.tsv",
        "qc_rows": qc_rows,
    }


def _make_embedding_loader(cfg: dict[str, Any]):
    emb_cfg = cfg.get("embeddings", {})

    def _loader(seq_id: str, seq: str):
        import esm
        import torch

        global _RESULT_KEYS_LOGGED

        model_name = emb_cfg.get("model", "esm2_t33_650M_UR50D")
        device = emb_cfg.get("device", "cpu")
        cache_key = (model_name, device)
        if cache_key not in _ESM_CACHE:
            model, alphabet = esm.pretrained.load_model_and_alphabet(model_name)
            model = model.eval().to(device)
            batch_converter = alphabet.get_batch_converter()
            _ESM_CACHE[cache_key] = (model, alphabet, batch_converter)
        model, _alphabet, batch_converter = _ESM_CACHE[cache_key]
        _, _, toks = batch_converter([("seq", seq)])
        toks = toks.to(device)
        with torch.no_grad():
            try:
                results = model(toks, repr_layers=[], return_contacts=True)
            except TypeError:
                results = model(toks, repr_layers=[], need_head_weights=True)
        if "attentions" not in results:
            if not _RESULT_KEYS_LOGGED:
                logger.error("ESM forward output keys: %s", sorted(results.keys()))
                _RESULT_KEYS_LOGGED = True
            raise RuntimeError(
                "ESM forward pass did not return 'attentions'. Ensure return_contacts=True or need_head_weights=True are supported by the model."
            )
        return _build_combined_attention_tensor(results["attentions"])

    return _loader


def _write_qc_plots(hmm_id: str, qc_dir: Path, ha_sites: pd.DataFrame, ha_counts: pd.DataFrame, loc_layers: pd.DataFrame):
    qc_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(6, 3))
    if "n_ha" in ha_counts.columns and not ha_counts.empty:
        ax.hist(ha_counts["n_ha"], bins=min(20, max(5, len(ha_counts))))
    ax.set_title("HA counts per sequence")
    fig.tight_layout()
    fig.savefig(qc_dir / "ha_counts_per_seq.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 3))
    layer_col = "loc_layer_id" if "loc_layer_id" in loc_layers.columns else "loc_layer"
    if layer_col in loc_layers.columns and not loc_layers.empty:
        loc_layers[layer_col].value_counts().sort_index().plot(kind="bar", ax=ax)
    ax.set_title("LoC layer distribution")
    fig.tight_layout()
    fig.savefig(qc_dir / "loc_layer_distribution.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 3))
    if not ha_sites.empty and {"ungapped_pos", "is_ha"}.issubset(ha_sites.columns):
        prof = ha_sites.groupby("ungapped_pos")["is_ha"].mean()
        ax.plot(prof.index, prof.values)
    ax.set_title(f"{hmm_id} HA density track")
    fig.tight_layout()
    fig.savefig(qc_dir / "ha_density_track.png", dpi=150)
    plt.close(fig)


def run_ha_sites(cfg: dict, fasta_dir: str, emb_dir: str, summary_dir: str, qc_dir: str, hmm_keep=None, alignment_dir: str | None = None):
    """Stage wrapper: compute HA for each HMM independently."""
    del emb_dir
    all_sites, all_counts, all_loc, all_run_qc = [], [], [], []
    loader = _make_embedding_loader(cfg)

    for faa in sorted(glob.glob(os.path.join(fasta_dir, "*.faa"))):
        hmm_id = Path(faa).stem
        if hmm_keep and hmm_id not in hmm_keep:
            continue
        hmm_out = Path(summary_dir).parent / "ha" / hmm_id
        alignment_path = None
        if alignment_dir:
            candidate = Path(alignment_dir) / f"{hmm_id}.faa"
            alignment_path = candidate if candidate.exists() else candidate
        result = compute_ha_sites_for_hmm(hmm_id, faa, hmm_out, cfg, loader, alignment_path=alignment_path)

        ha_sites = pd.read_csv(result["ha_sites"], sep="\t")
        ha_counts = pd.read_csv(result["ha_counts"], sep="\t")
        loc_layers = pd.read_csv(result["loc_layers"], sep="\t")
        run_qc = pd.read_csv(result["ha_run_qc"], sep="\t")
        _write_qc_plots(hmm_id, Path(qc_dir) / hmm_id / "ha", ha_sites, ha_counts, loc_layers)

        all_sites.append(ha_sites)
        all_counts.append(ha_counts.assign(hmm_id=hmm_id))
        all_loc.append(loc_layers.assign(hmm_id=hmm_id))
        all_run_qc.append(run_qc)

    if all_sites:
        pd.concat(all_sites, ignore_index=True).to_csv(Path(summary_dir) / "ha_sites.tsv", sep="\t", index=False)
        pd.concat(all_counts, ignore_index=True).to_csv(Path(summary_dir) / "ha_counts.tsv", sep="\t", index=False)
        pd.concat(all_loc, ignore_index=True).to_csv(Path(summary_dir) / "loc_layers.tsv", sep="\t", index=False)
        pd.concat(all_run_qc, ignore_index=True).to_csv(Path(summary_dir) / "ha_run_qc.tsv", sep="\t", index=False)
    else:
        pd.DataFrame([]).to_csv(Path(summary_dir) / "ha_sites.tsv", sep="\t", index=False)
        pd.DataFrame([]).to_csv(Path(summary_dir) / "ha_counts.tsv", sep="\t", index=False)
        pd.DataFrame([]).to_csv(Path(summary_dir) / "loc_layers.tsv", sep="\t", index=False)
        pd.DataFrame([]).to_csv(Path(summary_dir) / "ha_run_qc.tsv", sep="\t", index=False)

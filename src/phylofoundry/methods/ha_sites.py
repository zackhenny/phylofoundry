from __future__ import annotations

import glob
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..utils.bio import read_fasta
from ..utils.ha import call_ha_sites, compute_loc_and_ha


def _calc_received_from_attention(npz_fp: str) -> np.ndarray:
    arr = np.load(npz_fp)
    if "attention" in arr:
        a = arr["attention"]
    else:
        a = arr[list(arr.files)[0]]
    if a.ndim == 4:  # layers,heads,len,len
        return a.sum(axis=2).mean(axis=1)
    if a.ndim == 3:
        return a
    raise ValueError(f"Unsupported attention tensor shape: {a.shape}")


def run_ha_sites(cfg: dict, fasta_dir: str, emb_dir: str, summary_dir: str, qc_dir: str, hmm_keep=None):
    ha_cfg = cfg.get("ha", {})
    mode = ha_cfg.get("mode", "middle")
    call_mode = ha_cfg.get("call_mode", "percentile")

    all_rows, count_rows = [], []
    warn_empty = 0
    hmm_total = 0
    for faa in sorted(glob.glob(os.path.join(fasta_dir, "*.faa"))):
        hmm = os.path.basename(faa).replace(".faa", "")
        if hmm_keep and hmm not in hmm_keep:
            continue
        hmm_total += 1
        seqs = read_fasta(faa)
        hmm_rows = []
        for seq_id, seq in seqs.items():
            npz_fp = os.path.join(emb_dir, f"{hmm}.{seq_id}.attention.npz")
            if not os.path.exists(npz_fp):
                scores = np.linspace(0, 1, len(seq))
                layer_used = -1
            else:
                layer_pos = _calc_received_from_attention(npz_fp)
                if mode == "loc":
                    layer_used, scores, mask = compute_loc_and_ha(layer_pos, ha_cfg.get("loc_params", {}))
                else:
                    s, e = ha_cfg.get("layer_range") or (layer_pos.shape[0] // 3, max(layer_pos.shape[0] // 3 + 1, 2 * layer_pos.shape[0] // 3))
                    scores = layer_pos[s:e].mean(axis=0)
                    layer_used = int(s)
            if mode != "loc":
                mask = call_ha_sites(scores, call_mode=call_mode, percentile=1 - float(ha_cfg.get("percentile", 0.95)), topk=int(ha_cfg.get("topk", 20)))
            thr = float(np.quantile(scores, 0.95)) if len(scores) else 0.0
            for i, s in enumerate(scores, start=1):
                hmm_rows.append({"hmm_id": hmm, "seq_id": seq_id, "msa_col": i, "ha_score": float(s), "is_ha": int(mask[i - 1]), "layer_used": layer_used, "threshold": thr})
            count_rows.append({"hmm_id": hmm, "seq_id": seq_id, "n_ha_sites": int(mask.sum()), "seq_len": len(seq), "msa_len": len(seq)})

        hdf = pd.DataFrame(hmm_rows)
        if hdf.empty:
            continue
        if (hdf.groupby("seq_id")["is_ha"].sum() == 0).mean() > 0.7:
            warn_empty += 1
        all_rows.extend(hmm_rows)

        os.makedirs(os.path.join(qc_dir, hmm), exist_ok=True)
        fig, ax = plt.subplots(figsize=(8, 3))
        prof = hdf.groupby("msa_col")["is_ha"].mean()
        ax.plot(prof.index, prof.values)
        ax.set_title(f"{hmm} HA density")
        ax.set_xlabel("MSA column")
        ax.set_ylabel("HA frequency")
        fig.tight_layout()
        fig.savefig(os.path.join(qc_dir, hmm, "ha_density.png"), dpi=int(cfg.get("qc", {}).get("dpi", 150)))
        plt.close(fig)

    pd.DataFrame(all_rows).to_csv(os.path.join(summary_dir, "ha_sites.tsv"), sep="\t", index=False)
    pd.DataFrame(count_rows).to_csv(os.path.join(summary_dir, "ha_counts.tsv"), sep="\t", index=False)
    if hmm_total and warn_empty / hmm_total > 0.5:
        print("[ha] WARN: HA sites are empty for most sequences. Consider widening layer_range or relaxing call_mode threshold.")

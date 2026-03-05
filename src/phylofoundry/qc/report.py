from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def _hist(vals, fp, title, x):
    if len(vals) == 0:
        return
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(vals, bins=25)
    ax.set_title(title)
    ax.set_xlabel(x)
    ax.set_ylabel("count")
    fig.tight_layout(); fig.savefig(fp); plt.close(fig)


def _line(x, y, fp, title, yl):
    fig, ax = plt.subplots(figsize=(8, 3))
    ax.plot(x, y)
    ax.set_title(title)
    ax.set_ylabel(yl)
    ax.set_xlabel("MSA column")
    fig.tight_layout(); fig.savefig(fp); plt.close(fig)


def generate_qc_summary(outdir: str, summary_dir: str, cfg: dict | None = None):
    cfg = cfg or {}
    qc_dir = os.path.join(outdir, "qc")
    os.makedirs(os.path.join(qc_dir, "combined"), exist_ok=True)

    manifest_rows = []
    scan_fp = os.path.join(summary_dir, "all_scan_hits.tsv")
    best_fp = os.path.join(summary_dir, "best_hits.tsv")
    if os.path.exists(scan_fp):
        df = pd.read_csv(scan_fp, sep="\t")
        if "bitscore" in df.columns: _hist(df["bitscore"].dropna(), os.path.join(qc_dir, "combined", "hmmer_bitscore_pre.png"), "HMMER bitscore", "bitscore")
        if "qcov" in df.columns: _hist(df["qcov"].dropna(), os.path.join(qc_dir, "combined", "hmmer_coverage_pre.png"), "HMMER coverage", "coverage")
    if os.path.exists(best_fp):
        df = pd.read_csv(best_fp, sep="\t")
        if "bitscore" in df.columns: _hist(df["bitscore"].dropna(), os.path.join(qc_dir, "combined", "hmmer_bitscore_post.png"), "HMMER bitscore (post)", "bitscore")
        if "qcov" in df.columns: _hist(df["qcov"].dropna(), os.path.join(qc_dir, "combined", "hmmer_coverage_post.png"), "HMMER coverage (post)", "coverage")
        for hmm, sub in df.groupby("hmm"):
            manifest_rows.append({"hmm": hmm, "n_hits": len(sub)})

    ha_fp = os.path.join(summary_dir, "ha_sites.tsv")
    if os.path.exists(ha_fp):
        ha = pd.read_csv(ha_fp, sep="\t")
        for hmm, sub in ha.groupby("hmm_id"):
            d = os.path.join(qc_dir, hmm); os.makedirs(d, exist_ok=True)
            _hist(sub.groupby("seq_id")["is_ha"].sum(), os.path.join(d, "ha_site_counts_per_sequence.png"), f"{hmm} HA counts", "n HA")
            freq = sub.groupby("msa_col")["is_ha"].mean()
            _line(freq.index, freq.values, os.path.join(d, "ha_frequency_msa.png"), f"{hmm} HA frequency", "HA freq")

    rs_fp = os.path.join(summary_dir, "regime_shifts.tsv")
    if os.path.exists(rs_fp):
        rs = pd.read_csv(rs_fp, sep="\t")
        _hist(rs.get("shift_score", pd.Series(dtype=float)).dropna(), os.path.join(qc_dir, "combined", "regime_shift_scores.png"), "Regime shift scores", "score")
        if "p_perm" in rs.columns:
            _hist(rs["p_perm"].dropna(), os.path.join(qc_dir, "combined", "regime_shift_pvalues.png"), "Regime shift p-values", "p")

    qc_manifest = os.path.join(summary_dir, "qc_manifest.tsv")
    pd.DataFrame(manifest_rows).to_csv(qc_manifest, sep="\t", index=False)


def write_run_manifest(outdir: str, summary_dir: str, step_status_rows: list[dict]):
    pd.DataFrame(step_status_rows).to_csv(os.path.join(summary_dir, "step_status.tsv"), sep="\t", index=False)
    manifest = {"outdir": outdir, "summary_dir": summary_dir, "n_steps": len(step_status_rows)}
    with open(os.path.join(summary_dir, "run_manifest.json"), "w") as f:
        json.dump(manifest, f, indent=2)

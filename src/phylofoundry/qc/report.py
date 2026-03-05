from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def _safe_hist(series, out_png: str, title: str, xlabel: str):
    if series is None or len(series) == 0:
        return
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(series, bins=25)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("count")
    fig.tight_layout()
    fig.savefig(out_png, dpi=120)
    plt.close(fig)


def generate_qc_summary(outdir: str, summary_dir: str):
    qc_dir = os.path.join(outdir, "qc")
    os.makedirs(os.path.join(qc_dir, "combined"), exist_ok=True)

    manifest_rows = []
    scan_fp = os.path.join(summary_dir, "all_scan_hits.tsv")
    if os.path.exists(scan_fp):
        df = pd.read_csv(scan_fp, sep="\t")
        if "bitscore" in df.columns:
            _safe_hist(df["bitscore"].dropna(), os.path.join(qc_dir, "combined", "hits_bitscore.png"), "bitscore distribution", "bitscore")
        if "qcov" in df.columns:
            _safe_hist(df["qcov"].dropna(), os.path.join(qc_dir, "combined", "hits_qcov.png"), "coverage distribution", "coverage")
        for hmm, sub in df.groupby("hmm"):
            manifest_rows.append({"hmm": hmm, "n_hits": int(len(sub)), "mean_bitscore": float(sub["bitscore"].mean()) if "bitscore" in sub.columns else None})

    qc_manifest = os.path.join(summary_dir, "qc_manifest.tsv")
    pd.DataFrame(manifest_rows).to_csv(qc_manifest, sep="\t", index=False)


def write_run_manifest(outdir: str, summary_dir: str, step_status_rows: list[dict]):
    manifest = {
        "outdir": outdir,
        "summary_dir": summary_dir,
        "artifacts": {},
    }
    for rel in ["resolved_config.json", "qc_manifest.tsv", "step_status.tsv"]:
        fp = os.path.join(summary_dir, rel)
        if os.path.exists(fp):
            manifest["artifacts"][rel] = fp
    with open(os.path.join(summary_dir, "run_manifest.json"), "w") as f:
        json.dump(manifest, f, indent=2)

    pd.DataFrame(step_status_rows).to_csv(os.path.join(summary_dir, "step_status.tsv"), sep="\t", index=False)

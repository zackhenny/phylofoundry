"""Amino acid composition analysis for selected gene hits.

Runs after the phylo (and optionally curate / detect_clades) stage to compute
per-gene amino acid composition statistics and biochemical metrics derived from
the approach in *genomeSPOT* (Barnum et al.), then optionally performs
clade-level comparative statistics and generates visualisation plots.

Outputs
-------
summary/aa_composition/
    aa_comp_per_gene.tsv   – gene × metric matrix
                             rows = individual selected-gene sequences
                             cols = gene_id, hmm, len, tip_or_clade,
                                    aa_A_comp … aa_Y_comp,
                                    Zc, gravy, nH2O, pi, S_content,
                                    N_content, thermo_freq
    aa_comp_stats.tsv      – per-AA/metric Kruskal-Wallis + BH-FDR table
                             (written only when ≥2 clade groups present)
    plots/                 – boxplots and heatmap PNGs
                             (written only when matplotlib/seaborn available)

Usage
-----
Standalone::

    from phylofoundry.tasks.aa_composition import run_aa_composition
    run_aa_composition(cfg, fasta_dir, summary_dir, hmm_keep=None,
                       clade_assign_dir=None, force=False)

Via pipeline::

    phylofoundry aa-composition --outdir ./results

Configuration keys (``aa_composition`` section)
------------------------------------------------
enabled : bool
    Whether this step runs.  Default ``false``.
remove_initial_met : bool
    Strip leading methionine before computing composition.  Default ``true``.
compute_pi : bool
    Compute isoelectric point per gene (requires BioPython).  Default ``true``.
generate_plots : bool
    Generate boxplot and heatmap PNGs.  Default ``true``.
top_n_aas_heatmap : int
    Number of AAs (by variance) shown in the heatmap.  Default ``10``.
"""

from __future__ import annotations

import glob
import os
import sys
import warnings
from collections import Counter
from typing import Optional

import numpy as np
import pandas as pd

from ..constants import AA_ALPHABET
from ..utils.bio import read_fasta
from ..utils.helpers import safe_mkdir

# ---------------------------------------------------------------------------
# Biochemical constants (adapted from genomeSPOT / Barnum et al.)
# ---------------------------------------------------------------------------

STANDARD_AAS: list[str] = sorted(AA_ALPHABET)  # 20 standard amino acids

# Average oxidation state of carbon (Zc) — per-residue weighted Zc weights
WEIGHTED_ZC: dict[str, float] = {
    "A":  0.0,  "C":  2.0,  "D":  4.0,  "E":  2.0,  "F": -4.0,
    "G":  2.0,  "H":  4.0,  "I": -6.0,  "K":  4.0,  "L": -6.0,
    "M": -1.6,  "N":  4.0,  "P": -2.0,  "Q":  2.0,  "R":  2.0,
    "S":  1.98, "T":  0.0,  "V": -4.0,  "W": -2.0,  "Y": -2.0,
}

# Number of carbon atoms per residue (used as denominator for Zc)
CARBON_NUMBER: dict[str, int] = {
    "A": 3,  "C": 3,  "D": 4,  "E": 5,  "F": 9,
    "G": 2,  "H": 6,  "I": 6,  "K": 6,  "L": 6,
    "M": 5,  "N": 4,  "P": 5,  "Q": 5,  "R": 6,
    "S": 3,  "T": 4,  "V": 5,  "W": 11, "Y": 9,
}

# Water demand (nH2O) per residue — QEC thermodynamic scale
NH2O_QEC: dict[str, float] = {
    "A": -0.4,  "C": -1.0,  "D": -1.2,  "E": -1.0,  "F": -3.2,
    "G": -0.6,  "H": -2.8,  "I":  0.2,  "K":  0.2,  "L":  0.2,
    "M": -0.6,  "N": -1.2,  "P": -1.0,  "Q": -1.0,  "R": -0.8,
    "S": -0.4,  "T": -0.2,  "V":  0.0,  "W": -4.8,  "Y": -3.2,
}

# Kyte-Doolittle hydrophobicity (GRAVY metric)
HYDROPHOBICITY: dict[str, float] = {
    "A":  1.8,  "R": -4.5,  "N": -3.5,  "D": -3.5,  "C":  2.5,
    "Q": -3.5,  "E": -3.5,  "G": -0.4,  "H": -3.2,  "I":  4.5,
    "L":  3.8,  "K": -3.9,  "M":  1.9,  "F":  2.8,  "P": -1.6,
    "S": -0.8,  "T": -0.7,  "W": -0.9,  "Y": -1.3,  "V":  4.2,
}

# Sulfur-containing residues
SULFUR_CONTENT: dict[str, int] = {aa: (1 if aa in {"C", "M"} else 0)
                                   for aa in STANDARD_AAS}

# Nitrogen atoms per residue
NITROGEN_CONTENT: dict[str, int] = {
    "A": 1, "C": 1, "D": 1, "E": 1, "F": 1,
    "G": 1, "H": 3, "I": 1, "K": 2, "L": 1,
    "M": 1, "N": 2, "P": 1, "Q": 2, "R": 4,
    "S": 1, "T": 1, "V": 1, "W": 2, "Y": 1,
}

# Residues associated with thermostability
THERMOSTABLE_RESIDUES: frozenset[str] = frozenset({"I", "V", "Y", "W", "R", "E", "L"})


# ---------------------------------------------------------------------------
# Core per-sequence computation
# ---------------------------------------------------------------------------

def _clean_sequence(seq: str, remove_initial_met: bool = True) -> list[str]:
    """Return a list of standard amino acids from *seq*.

    Non-standard characters (gaps, X, *, etc.) are silently dropped.
    If *remove_initial_met* is ``True`` and the first residue is ``M``, it is
    removed before further analysis.
    """
    seq = seq.upper().replace("-", "").replace(".", "").replace("*", "")
    if remove_initial_met and seq.startswith("M"):
        seq = seq[1:]
    return [aa for aa in seq if aa in AA_ALPHABET]


def compute_gene_metrics(
    gene_id: str,
    raw_seq: str,
    hmm: str,
    tip_or_clade: str = "unassigned",
    remove_initial_met: bool = True,
    compute_pi: bool = True,
) -> dict:
    """Compute per-gene amino acid composition and biochemical metrics.

    Parameters
    ----------
    gene_id:
        Identifier for this gene/sequence (FASTA header).
    raw_seq:
        Raw (possibly gapped) protein sequence string.
    hmm:
        HMM / gene family this sequence belongs to.
    tip_or_clade:
        Tree tip label or clade assignment for this gene.
    remove_initial_met:
        Whether to remove the leading methionine before analysis.
    compute_pi:
        Whether to compute the isoelectric point (requires BioPython).

    Returns
    -------
    dict
        Mapping of column name → value for all metrics.
    """
    seq = _clean_sequence(raw_seq, remove_initial_met=remove_initial_met)
    length = len(seq)

    row: dict = {
        "gene_id": gene_id,
        "hmm": hmm,
        "len": length,
        "tip_or_clade": tip_or_clade,
    }

    if length == 0:
        # Degenerate sequence — fill NaN for all metrics
        for aa in STANDARD_AAS:
            row[f"aa_{aa}_comp"] = float("nan")
        for col in ("Zc", "gravy", "nH2O", "S_content", "N_content", "thermo_freq"):
            row[col] = float("nan")
        if compute_pi:
            row["pi"] = float("nan")
        return row

    counts = Counter(seq)

    # AA fractional composition (fraction of each residue in the gene, 0–1)
    for aa in STANDARD_AAS:
        row[f"aa_{aa}_comp"] = counts.get(aa, 0) / length

    # Zc — average oxidation state of carbon
    total_zc_weighted = sum(WEIGHTED_ZC[a] for a in seq)
    total_carbon = sum(CARBON_NUMBER[a] for a in seq)
    row["Zc"] = (total_zc_weighted / total_carbon) if total_carbon > 0 else float("nan")

    # GRAVY — grand average of hydropathicity (Kyte-Doolittle)
    row["gravy"] = sum(HYDROPHOBICITY[a] for a in seq) / length

    # nH2O — water demand per residue
    row["nH2O"] = sum(NH2O_QEC[a] for a in seq) / length

    # Sulfur content fraction
    row["S_content"] = sum(SULFUR_CONTENT[a] for a in seq) / length

    # Nitrogen content fraction (atoms per residue)
    row["N_content"] = sum(NITROGEN_CONTENT[a] for a in seq) / length

    # Thermostable residue frequency
    row["thermo_freq"] = sum(counts.get(a, 0) for a in THERMOSTABLE_RESIDUES) / length

    # Isoelectric point (optional — requires BioPython)
    if compute_pi:
        try:
            from Bio.SeqUtils.ProtParam import ProteinAnalysis
            analysis = ProteinAnalysis("".join(seq))
            row["pi"] = analysis.isoelectric_point()
        except Exception:
            row["pi"] = float("nan")

    return row


# ---------------------------------------------------------------------------
# Comparative statistics helpers
# ---------------------------------------------------------------------------

def _benjamini_hochberg(pvalues: list[float]) -> list[float]:
    """Return BH-corrected FDR q-values for *pvalues* (1-D list)."""
    n = len(pvalues)
    if n == 0:
        return []
    indexed = sorted(enumerate(pvalues), key=lambda x: x[1])
    qvalues = [1.0] * n
    prev_qval = 1.0
    for rank, (orig_idx, pval) in enumerate(reversed(indexed), start=1):
        q = min(prev_qval, pval * n / (n - rank + 1))
        qvalues[orig_idx] = q
        prev_qval = q
    return qvalues


def run_comparative_stats(
    df: pd.DataFrame,
    metric_cols: list[str],
    group_col: str = "tip_or_clade",
) -> pd.DataFrame:
    """Run Kruskal-Wallis tests and apply BH-FDR correction.

    For each metric column in *metric_cols*, test whether the distribution
    differs significantly across the groups defined by *group_col*.

    Parameters
    ----------
    df:
        Gene-level metrics DataFrame (output of :func:`run_aa_composition`).
    metric_cols:
        Column names to test (e.g. ``['aa_A_comp', 'aa_C_comp', ..., 'Zc']``).
    group_col:
        Column whose unique values define the groups.

    Returns
    -------
    pd.DataFrame
        Table with columns ``metric``, ``H_stat``, ``p_value``, ``q_value``
        (BH-corrected FDR), ``n_groups``.
    """
    try:
        from scipy.stats import kruskal
    except ImportError:
        warnings.warn(
            "[aa_composition] scipy not available; skipping comparative statistics.",
            stacklevel=2,
        )
        return pd.DataFrame()

    groups = df[group_col].dropna().unique()
    if len(groups) < 2:
        return pd.DataFrame()

    rows = []
    pvalues = []
    for col in metric_cols:
        group_data = [
            df.loc[df[group_col] == g, col].dropna().values
            for g in groups
        ]
        # Kruskal-Wallis requires ≥1 observation per group and ≥2 groups
        group_data = [g for g in group_data if len(g) > 0]
        if len(group_data) < 2:
            rows.append({"metric": col, "H_stat": float("nan"),
                         "p_value": float("nan"), "n_groups": len(group_data)})
            pvalues.append(float("nan"))
            continue
        try:
            stat, p = kruskal(*group_data)
        except Exception:
            stat, p = float("nan"), float("nan")
        rows.append({"metric": col, "H_stat": stat,
                     "p_value": p, "n_groups": len(group_data)})
        pvalues.append(p)

    if not rows:
        return pd.DataFrame()

    # BH-FDR correction — skip NaN entries
    valid_idx = [i for i, p in enumerate(pvalues) if not np.isnan(p)]
    valid_p = [pvalues[i] for i in valid_idx]
    q_all = [float("nan")] * len(pvalues)
    if valid_p:
        bh_q = _benjamini_hochberg(valid_p)
        for rank_i, orig_i in enumerate(valid_idx):
            q_all[orig_i] = bh_q[rank_i]

    for i, row in enumerate(rows):
        row["q_value"] = q_all[i]

    return pd.DataFrame(rows)[["metric", "H_stat", "p_value", "q_value", "n_groups"]]


# ---------------------------------------------------------------------------
# Plotting helpers (optional — matplotlib / seaborn)
# ---------------------------------------------------------------------------

def _make_plots(df: pd.DataFrame, plots_dir: str, top_n: int = 10) -> None:
    """Generate boxplot and heatmap PNGs.

    Silently skips if matplotlib or seaborn are unavailable.

    Parameters
    ----------
    df:
        Gene-level metrics DataFrame.
    plots_dir:
        Directory where plot files are written.
    top_n:
        Number of amino acids shown in the heatmap (chosen by highest
        inter-clade variance).
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import seaborn as sns
    except ImportError:
        print(
            "[aa_composition] matplotlib/seaborn not available; skipping plots.",
            file=sys.stderr,
        )
        return

    safe_mkdir(plots_dir)
    aa_cols = [f"aa_{aa}_comp" for aa in STANDARD_AAS]
    present_aa_cols = [c for c in aa_cols if c in df.columns]

    group_col = "tip_or_clade"
    n_groups = df[group_col].nunique()

    # ── Boxplots for each AA ────────────────────────────────────────────────
    for col in present_aa_cols:
        aa = col.split("_")[1]
        fig, ax = plt.subplots(figsize=(max(6, n_groups * 0.8 + 2), 4))
        try:
            sns.boxplot(data=df, x=group_col, y=col, ax=ax, palette="tab20")
        except Exception:
            df.boxplot(column=col, by=group_col, ax=ax)
        ax.set_title(f"AA {aa} composition by clade/tip")
        ax.set_xlabel("Clade / tip")
        ax.set_ylabel(f"Fraction ({aa})")
        ax.tick_params(axis="x", rotation=45)
        fig.tight_layout()
        fig.savefig(os.path.join(plots_dir, f"boxplot_aa_{aa}.png"), dpi=150)
        plt.close(fig)

    # ── Heatmap: mean composition per clade ────────────────────────────────
    # Select top-N AAs by inter-clade variance of means
    clade_means = df.groupby(group_col)[present_aa_cols].mean()
    if len(present_aa_cols) > 0 and len(clade_means) > 1:
        try:
            variances = clade_means.var(axis=0)
            top_cols = variances.nlargest(min(top_n, len(present_aa_cols))).index.tolist()

            fig, ax = plt.subplots(figsize=(max(8, len(top_cols) * 0.7 + 2),
                                            max(4, n_groups * 0.5 + 2)))
            sns.heatmap(
                clade_means[top_cols],
                annot=True, fmt=".3f",
                cmap="viridis", ax=ax,
                linewidths=0.3,
            )
            ax.set_title("Mean AA composition by clade/tip (top variable AAs)")
            ax.set_xlabel("Amino acid")
            ax.set_ylabel("Clade / tip")
            fig.tight_layout()
            fig.savefig(os.path.join(plots_dir, "heatmap_aa_composition.png"), dpi=150)
            plt.close(fig)
        except Exception as exc:
            print(f"[aa_composition] Heatmap failed: {exc}", file=sys.stderr)

    # ── Biochemical metrics boxplots ───────────────────────────────────────
    metric_cols = [c for c in ("Zc", "gravy", "nH2O", "S_content",
                               "N_content", "thermo_freq", "pi")
                   if c in df.columns]
    for col in metric_cols:
        fig, ax = plt.subplots(figsize=(max(6, n_groups * 0.8 + 2), 4))
        try:
            sns.boxplot(data=df, x=group_col, y=col, ax=ax, palette="tab20")
        except Exception:
            df.boxplot(column=col, by=group_col, ax=ax)
        ax.set_title(f"{col} by clade/tip")
        ax.set_xlabel("Clade / tip")
        ax.set_ylabel(col)
        ax.tick_params(axis="x", rotation=45)
        fig.tight_layout()
        fig.savefig(os.path.join(plots_dir, f"boxplot_{col}.png"), dpi=150)
        plt.close(fig)


# ---------------------------------------------------------------------------
# Clade-assignment loading helper
# ---------------------------------------------------------------------------

def _load_clade_assignments(clade_assign_dir: str, hmm: str) -> dict[str, str]:
    """Return a ``{tip_label: clade_name}`` dict for *hmm*.

    Tries ``<clade_assign_dir>/<hmm>.tsv`` first; then looks for
    ``summary/detected_clades.tsv`` one level up.  Returns an empty dict if
    no file is found.
    """
    per_hmm_tsv = os.path.join(clade_assign_dir, f"{hmm}.tsv")
    if os.path.exists(per_hmm_tsv):
        try:
            df = pd.read_csv(per_hmm_tsv, sep="\t")
            # Expect columns: tip, clade_name (or clade)
            tip_col = next((c for c in df.columns if c.lower() in ("tip", "tip_label")), None)
            clade_col = next(
                (c for c in df.columns if c.lower() in ("clade_name", "clade")), None
            )
            if tip_col and clade_col:
                return dict(zip(df[tip_col], df[clade_col]))
        except Exception:
            pass

    # Fall back: detected_clades.tsv with columns [clade_name (or clade), tip (or tip_label)]
    detected = os.path.join(os.path.dirname(clade_assign_dir), "detected_clades.tsv")
    if not os.path.exists(detected):
        detected = os.path.join(clade_assign_dir, "..", "detected_clades.tsv")
    if os.path.exists(detected):
        try:
            df = pd.read_csv(detected, sep="\t")
            tip_col = next((c for c in df.columns if c.lower() in ("tip", "tip_label")), None)
            clade_col = next(
                (c for c in df.columns if c.lower() in ("clade_name", "clade")), None
            )
            if tip_col and clade_col:
                return dict(zip(df[tip_col], df[clade_col]))
        except Exception:
            pass

    return {}


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def run_aa_composition(
    cfg: dict,
    fasta_dir: str,
    summary_dir: str,
    hmm_keep: Optional[set] = None,
    clade_assign_dir: Optional[str] = None,
    force: bool = False,
) -> pd.DataFrame:
    """Compute per-gene AA composition and comparative statistics.

    Parameters
    ----------
    cfg:
        Resolved pipeline configuration dict.
    fasta_dir:
        Directory containing per-HMM ``.faa`` files of selected gene hits.
        When the curate step has run, this should already point to the
        curated overlay (the pipeline handles the redirect).
    summary_dir:
        Top-level ``summary/`` directory; outputs are written under
        ``summary/aa_composition/``.
    hmm_keep:
        Optional set of HMM names to restrict processing to.
    clade_assign_dir:
        Directory containing per-HMM clade assignment TSV files
        (``clade_assignments/<hmm>.tsv``).  When ``None``, clade labels are
        omitted from the output.
    force:
        Re-compute and overwrite existing outputs when ``True``.

    Returns
    -------
    pd.DataFrame
        The gene-level metrics table (also written to disk).
    """
    print("\n[aa_composition] Computing per-gene amino acid composition metrics...")

    aa_cfg = cfg.get("aa_composition", {})
    remove_initial_met: bool = aa_cfg.get("remove_initial_met", True)
    compute_pi: bool = aa_cfg.get("compute_pi", True)
    generate_plots: bool = aa_cfg.get("generate_plots", True)
    top_n_heatmap: int = int(aa_cfg.get("top_n_aas_heatmap", 10))

    out_dir = os.path.join(summary_dir, "aa_composition")
    safe_mkdir(out_dir)

    per_gene_tsv = os.path.join(out_dir, "aa_comp_per_gene.tsv")
    stats_tsv = os.path.join(out_dir, "aa_comp_stats.tsv")

    # ── Discover per-HMM FASTA files ──────────────────────────────────────
    faa_files = sorted(glob.glob(os.path.join(fasta_dir, "*.faa")))
    if not faa_files:
        print(
            f"[aa_composition] No .faa files found in {fasta_dir}; nothing to do.",
            file=sys.stderr,
        )
        return pd.DataFrame()

    # ── Compute metrics for each gene ─────────────────────────────────────
    all_rows: list[dict] = []
    for faa_fp in faa_files:
        hmm = os.path.basename(faa_fp).rsplit(".", 1)[0]
        if hmm_keep is not None and hmm not in hmm_keep:
            continue

        seqs = read_fasta(faa_fp)
        if not seqs:
            continue

        # Load clade assignments for this HMM if available
        clade_map: dict[str, str] = {}
        if clade_assign_dir and os.path.isdir(clade_assign_dir):
            clade_map = _load_clade_assignments(clade_assign_dir, hmm)

        for gene_id, raw_seq in seqs.items():
            tip_label = clade_map.get(gene_id, "unassigned")
            row = compute_gene_metrics(
                gene_id=gene_id,
                raw_seq=raw_seq,
                hmm=hmm,
                tip_or_clade=tip_label,
                remove_initial_met=remove_initial_met,
                compute_pi=compute_pi,
            )
            all_rows.append(row)

    if not all_rows:
        print("[aa_composition] No sequences processed.", file=sys.stderr)
        return pd.DataFrame()

    df = pd.DataFrame(all_rows)

    # ── Canonical column ordering ──────────────────────────────────────────
    front_cols = ["gene_id", "hmm", "len", "tip_or_clade"]
    aa_cols = [f"aa_{aa}_comp" for aa in STANDARD_AAS if f"aa_{aa}_comp" in df.columns]
    metric_cols = [c for c in ("Zc", "gravy", "nH2O", "S_content",
                               "N_content", "thermo_freq", "pi")
                   if c in df.columns]
    df = df[front_cols + aa_cols + metric_cols]

    # ── Write gene-level matrix ────────────────────────────────────────────
    df.to_csv(per_gene_tsv, sep="\t", index=False)
    print(f"[aa_composition] Wrote {len(df):,} rows → {per_gene_tsv}")

    # ── Comparative statistics (Kruskal-Wallis + BH-FDR) ──────────────────
    test_cols = aa_cols + [c for c in metric_cols if c != "len"]
    n_clades = df["tip_or_clade"].nunique()
    if n_clades >= 2:
        stats_df = run_comparative_stats(df, test_cols, group_col="tip_or_clade")
        if not stats_df.empty:
            stats_df.to_csv(stats_tsv, sep="\t", index=False)
            n_sig = (stats_df["q_value"].fillna(1.0) < 0.05).sum()
            print(
                f"[aa_composition] Comparative statistics written → {stats_tsv} "
                f"({n_sig} metrics significant at FDR<0.05)"
            )
    else:
        print(
            "[aa_composition] <2 clade groups; skipping comparative statistics."
        )

    # ── Plots ──────────────────────────────────────────────────────────────
    if generate_plots and n_clades >= 2:
        plots_dir = os.path.join(out_dir, "plots")
        _make_plots(df, plots_dir, top_n=top_n_heatmap)
        print(f"[aa_composition] Plots written → {plots_dir}/")

    print("[aa_composition] Done.")
    return df

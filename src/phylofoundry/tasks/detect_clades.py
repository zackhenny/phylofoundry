"""Clade detection step.

Responsibilities:
- Generate ``summary/detected_clades.tsv``.
- Support configured clade detection strategies:
  ``taxonomy``, ``treecluster``, ``tree_embed``, or a pre-built TSV.
- Emit node-level QC outputs such as ``node_scores.embedtree.tsv``.

This module contains logic extracted from the monolithic ``post`` step so that
clade detection can run independently of conservation metrics and taxonomy
integration.

HyPhy can depend on this step's outputs (``clade_assignments/``) rather than
relying on the monolithic ``post`` step.
"""

import os

import pandas as pd

from .post import (
    _detect_taxonomy_clades,
    _detect_treecluster_clades,
    _detect_tree_embed_clades,
    _load_taxonomy,
    _write_detected_clades,
    load_clades_tsv,
)
from ..utils.helpers import safe_mkdir


def run_detect_clades(cfg, tree_dir, emb_dir, summary_dir, hmm_keep, clade_assign_dir):
    """Detect clades and write assignment files.

    Parameters
    ----------
    cfg:
        Resolved pipeline configuration dict.
    tree_dir:
        Directory containing ``*.treefile`` files.
    emb_dir:
        Directory containing embedding PCA/vector TSV files (used only for
        the ``tree_embed`` detection method).
    summary_dir:
        Path to the ``summary/`` output directory where
        ``detected_clades.tsv`` and ``node_scores.embedtree.tsv`` are written.
    hmm_keep:
        Optional set of HMM names to restrict processing to, or ``None`` to
        process all HMMs.
    clade_assign_dir:
        Directory where per-HMM clade assignment TSV files are written
        (``clade_assignments/``).

    Returns
    -------
    dict
        Mapping ``{clade_name: [tip, ...]}`` for all detected clades.  Empty
        dict when no clades are detected.
    """
    print("\n[detect_clades] Detecting clades...")
    safe_mkdir(clade_assign_dir)

    dc_cfg = cfg.get("detect_clades", {})
    # Fall back to post.* config keys for backward compatibility.
    post_cfg = cfg.get("post", {})

    clades_tsv = dc_cfg.get("clades_tsv", post_cfg.get("clades_tsv", None))
    # "detect_method" is the canonical key; fall back to "detect_clades_method"
    # from the legacy post config.
    detect_method = str(
        dc_cfg.get("detect_method", post_cfg.get("detect_clades_method", ""))
    ).strip().lower()

    clades = None
    embed_node_scores = []

    if clades_tsv:
        clades = load_clades_tsv(clades_tsv)
        print(f"[detect_clades] Loaded pre-built clades from {clades_tsv}")
    elif detect_method == "taxonomy":
        gtdb_dir = cfg["inputs"].get("gtdb_dir")
        tax_file = cfg["inputs"].get("taxonomy_file")
        tax_map = {}
        if gtdb_dir or tax_file:
            tax_map = _load_taxonomy(gtdb_dir, tax_file)
        taxonomy_clade_level = dc_cfg.get(
            "taxonomy_clade_level",
            post_cfg.get("taxonomy_clade_level", "genus"),
        )
        clades = _detect_taxonomy_clades(summary_dir, tax_map, taxonomy_clade_level)
    elif detect_method == "treecluster":
        treecluster_threshold = dc_cfg.get(
            "treecluster_threshold",
            post_cfg.get("treecluster_threshold", 0.045),
        )
        treecluster_method = dc_cfg.get(
            "treecluster_method",
            post_cfg.get("treecluster_method", "max_clade"),
        )
        clades = _detect_treecluster_clades(
            tree_dir, treecluster_threshold, treecluster_method
        )
    elif detect_method == "tree_embed":
        # Merge detect_clades config on top of post config so that
        # embedtree_* keys from either section are honoured.
        combined_cfg = {**post_cfg, **dc_cfg}
        clades, embed_node_scores = _detect_tree_embed_clades(
            tree_dir,
            emb_dir,
            combined_cfg,
            hmm_keep=hmm_keep,
        )

    if clades:
        _write_detected_clades(summary_dir, clades, clade_assign_dir=clade_assign_dir)

    if embed_node_scores:
        pd.DataFrame(embed_node_scores).to_csv(
            os.path.join(summary_dir, "node_scores.embedtree.tsv"),
            sep="\t",
            index=False,
        )

    n = len(clades) if clades else 0
    print(f"[detect_clades] Done. Detected {n} clades.")
    return clades or {}

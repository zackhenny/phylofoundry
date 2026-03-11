"""Conservation metrics step using scikit-bio.

Responsibilities:
- Compute per-site conservation scores (inverse Shannon uncertainty or other
  scikit-bio metrics).
- Compute KL divergence between clade amino-acid distributions.
- Perform MRCA sanity checks.
- Write output files under ``summary/post_scikitbio/``.

This module contains logic extracted from the monolithic ``post`` step so that
conservation metrics can run independently of taxonomy integration and clade
detection.
"""

import glob
import os
import sys

import pandas as pd

from .post import (
    consurf_like_scores,
    kl_divergence,
    load_clades_tsv,
    site_counts_from_subset,
    tree_load,
    validate_tips_in_tree,
)
from ..utils.bio import read_fasta
from ..utils.helpers import safe_mkdir


def run_conservation_metrics(cfg, tree_dir, clipkit_dir, aln_dir, post_dir, hmm_keep):
    """Compute conservation metrics and write output TSVs.

    Parameters
    ----------
    cfg:
        Resolved pipeline configuration dict.
    tree_dir:
        Directory containing ``*.treefile`` files.
    clipkit_dir:
        Directory containing ``*.clipkit.faa`` trimmed alignments.
    aln_dir:
        Directory containing raw ``*.afa`` / ``*.mafft.fasta`` alignments.
    post_dir:
        Output directory for scikit-bio metric files
        (``summary/post_scikitbio/``).
    hmm_keep:
        Optional set of HMM names to restrict processing to, or ``None`` to
        process all HMMs.
    """
    print("\n[conservation_metrics] Computing conservation metrics...")
    safe_mkdir(post_dir)

    cm_cfg = cfg.get("conservation_metrics", {})
    # Fall back to post.* config keys for backward compatibility.
    post_cfg = cfg.get("post", {})

    compute_conservation = cm_cfg.get(
        "compute_conservation", post_cfg.get("compute_conservation", False)
    )
    conservation_metric = cm_cfg.get(
        "conservation_metric",
        post_cfg.get("conservation_metric", "inverse_shannon_uncertainty"),
    )
    compute_kl = cm_cfg.get("compute_kl", post_cfg.get("compute_kl", False))
    kl_pairs_str = cm_cfg.get("kl_pairs", post_cfg.get("kl_pairs", None))
    clades_tsv = cm_cfg.get("clades_tsv", post_cfg.get("clades_tsv", None))

    # Resolve clade assignments
    clades = None
    if clades_tsv:
        clades = load_clades_tsv(clades_tsv)
    else:
        # Fall back to detected_clades.tsv written by detect_clades or post
        summary_dir = os.path.dirname(post_dir)
        detected_fp = os.path.join(summary_dir, "detected_clades.tsv")
        if os.path.exists(detected_fp):
            try:
                clades = load_clades_tsv(detected_fp)
                print(f"[conservation_metrics] Loaded clades from {detected_fp}")
            except Exception:
                clades = None

    if compute_kl and not clades:
        raise SystemExit(
            "conservation_metrics.compute_kl requires clades from "
            "conservation_metrics.clades_tsv or summary/detected_clades.tsv"
        )

    kl_pairs = []
    has_manual_kl_pairs = bool(kl_pairs_str)
    if kl_pairs_str:
        for pair in str(kl_pairs_str).split(","):
            pair = pair.strip()
            if not pair:
                continue
            if ":" not in pair:
                raise SystemExit(
                    "conservation_metrics.kl_pairs must look like 'A:B' "
                    "or 'A:background'"
                )
            a, b = pair.split(":", 1)
            kl_pairs.append((a.strip(), b.strip()))

    cons_rows = []
    kl_rows = []
    mrca_rows = []

    hmm_names = sorted(
        [
            os.path.basename(x).replace(".treefile", "")
            for x in glob.glob(os.path.join(tree_dir, "*.treefile"))
        ]
    )
    if hmm_keep is not None:
        hmm_names = [h for h in hmm_names if h in hmm_keep]

    for hmm in hmm_names:
        tree_fp = os.path.join(tree_dir, f"{hmm}.treefile")

        clip_aln = os.path.join(clipkit_dir, f"{hmm}.clipkit.faa")
        if os.path.exists(clip_aln):
            aln_fp = clip_aln
        else:
            aln_fp = None
            for candidate in [
                os.path.join(aln_dir, f"{hmm}.mafft.fasta"),
                os.path.join(aln_dir, f"{hmm}.afa"),
            ]:
                if os.path.exists(candidate):
                    aln_fp = candidate
                    break

        if not aln_fp or not os.path.exists(tree_fp):
            continue

        aln_seqs = read_fasta(aln_fp)

        if compute_conservation:
            try:
                cons = consurf_like_scores(aln_seqs, metric=conservation_metric)
                for i, v in enumerate(cons, start=1):
                    cons_rows.append(
                        {
                            "hmm": hmm,
                            "scope": "global",
                            "site_1based": i,
                            "conservation": float(v),
                        }
                    )
            except Exception as e:
                print(
                    f"[conservation_metrics] Conservation scoring failed for "
                    f"{hmm}: {e}",
                    file=sys.stderr,
                )

        if compute_kl and clades:
            from collections import defaultdict

            counts = {}
            counts["background"] = site_counts_from_subset(aln_seqs, subset_tips=None)
            dynamic_pairs = []
            for cname, tips in clades.items():
                tips_in_aln = [x for x in tips if x in aln_seqs]
                counts[cname] = site_counts_from_subset(
                    aln_seqs, subset_tips=tips_in_aln
                )
                if not has_manual_kl_pairs and tips_in_aln:
                    others_key = f"others::{cname}"
                    others_tips = [x for x in aln_seqs if x not in set(tips_in_aln)]
                    if others_tips:
                        counts[others_key] = site_counts_from_subset(
                            aln_seqs, subset_tips=others_tips
                        )
                        dynamic_pairs.append((cname, others_key))

            L = len(counts["background"])
            pairs_to_use = kl_pairs if has_manual_kl_pairs else dynamic_pairs
            for a, b in pairs_to_use:
                if a not in counts or b not in counts:
                    continue
                if len(counts[a]) != L or len(counts[b]) != L:
                    continue
                for i in range(L):
                    kl = kl_divergence(counts[a][i], counts[b][i])
                    kl_rows.append(
                        {
                            "hmm": hmm,
                            "pair": f"{a}:{b}",
                            "site_1based": i + 1,
                            "kl_bits": kl,
                            "nA": sum(counts[a][i].values()),
                            "nB": sum(counts[b][i].values()),
                        }
                    )

        # MRCA sanity checks
        if clades:
            try:
                t = tree_load(tree_fp)
                for cname, tips in clades.items():
                    tips2 = [x for x in tips if x in aln_seqs]
                    missing = validate_tips_in_tree(t, tips2)
                    if missing:
                        mrca_rows.append(
                            {
                                "hmm": hmm,
                                "clade": cname,
                                "status": "missing_tips",
                                "detail": ",".join(missing[:10]),
                            }
                        )
                        continue
                    mrca = t.lca(tips2)
                    mrca_rows.append(
                        {
                            "hmm": hmm,
                            "clade": cname,
                            "status": "ok",
                            "mrca_id": mrca.id,
                            "mrca_name": mrca.name or "",
                        }
                    )
            except Exception as e:
                print(
                    f"[conservation_metrics] MRCA sanity checks failed for "
                    f"{hmm}: {e}",
                    file=sys.stderr,
                )

    if cons_rows:
        pd.DataFrame(cons_rows).to_csv(
            os.path.join(post_dir, "conservation.tsv"), sep="\t", index=False
        )
    if kl_rows:
        pd.DataFrame(kl_rows).to_csv(
            os.path.join(post_dir, "kl_divergence.tsv"), sep="\t", index=False
        )
    if mrca_rows:
        pd.DataFrame(mrca_rows).to_csv(
            os.path.join(post_dir, "mrca_sanity.tsv"), sep="\t", index=False
        )

    print("[conservation_metrics] Done.")

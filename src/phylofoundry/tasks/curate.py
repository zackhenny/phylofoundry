import os
import glob
import json
import shutil
import datetime
import pandas as pd
from Bio import Phylo, SeqIO
from ..utils.helpers import run_cmd, safe_mkdir

def _prune_fasta(in_fp, out_fp, keep_ids):
    """Write a pruned FASTA to *out_fp* keeping only sequences in *keep_ids*.

    The source *in_fp* is never modified.
    """
    if not os.path.exists(in_fp):
        return
    recs = [r for r in SeqIO.parse(in_fp, "fasta") if r.id in keep_ids]
    SeqIO.write(recs, out_fp, "fasta")

def run_treeshrink(tree_fp, hmm, tmp_dir):
    """Run TreeShrink to automatically detect and prune abnormally long branches."""
    out_tree = os.path.join(tmp_dir, f"{hmm}.treeshrink.treefile")
    cmd = ["run_treeshrink.py", "-t", tree_fp, "-O", hmm, "-o", tmp_dir]
    try:
        # treeshrink creates a directory named `hmm` inside `tmp_dir`
        run_cmd(cmd, quiet=True, shell=False)
        shrunk_fp = os.path.join(tmp_dir, hmm, "output.tree")
        if os.path.exists(shrunk_fp):
            shutil.copy2(shrunk_fp, out_tree)
            return out_tree
    except Exception as e:
        import sys
        print(f"[curate] TreeShrink warning for {hmm}: {e}", file=sys.stderr)
    return None

def run_esm_filter(tree_fp, emb_fp, out_tree_fp, min_samples=5):
    """Use DBSCAN on ESM embeddings to drop structural outliers from the tree."""
    try:
        from sklearn.cluster import DBSCAN
        import numpy as np
    except ImportError:
        print("[curate] scikit-learn not found. Skipping ESM filter.")
        return None
        
    if not os.path.exists(emb_fp):
        return None
        
    df = pd.read_parquet(emb_fp)
    if len(df) < min_samples * 2:
        return None # Too small to safely cluster
        
    # Extract structural embedding array
    X = np.stack(df["embedding"].values)
    
    # DBSCAN clustering (cosine distance is often better for high-dim embeddings)
    # eps=0.3 is a heuristic for cosine distance embeddings, might need tuning
    clusterer = DBSCAN(eps=0.3, min_samples=min_samples, metric="cosine")
    labels = clusterer.fit_predict(X)
    
    # Identify the 'core' cluster (usually the largest non-noise cluster)
    unique_labels, counts = np.unique(labels[labels != -1], return_counts=True)
    if len(unique_labels) == 0:
        return None # Everything is noise? Skip filtering to be safe
        
    core_label = unique_labels[np.argmax(counts)]
    
    # Identify valid IDs
    valid_ids = set(df[labels == core_label]["id"].values)
    
    # Prune from tree
    tree = Phylo.read(tree_fp, "newick")
    tree_terminals = [t.name for t in tree.get_terminals()]
    dropped = 0
    for t in tree_terminals:
        if t not in valid_ids:
            try:
                tree.prune(t)
                dropped += 1
            except ValueError:
                pass
                
    if dropped > 0:
        Phylo.write(tree, out_tree_fp, "newick")
        return out_tree_fp
    return None

def _write_manifest(manifest_fp, data):
    """Write a per-HMM curation manifest as JSON."""
    with open(manifest_fp, "w") as fh:
        json.dump(data, fh, indent=2)


def _append_audit(audit_fp, record):
    """Append a single JSON record (one line) to the audit JSONL log."""
    with open(audit_fp, "a") as fh:
        fh.write(json.dumps(record) + "\n")


def run_curate(cfg, tree_dir, fasta_dir, clipkit_dir, emb_dir, summary_dir,
               hmm_keep, force=False, curated_dir=None):
    """Run automated curation writing outputs to an explicit overlay directory.

    Raw files in *tree_dir*, *fasta_dir*, and *clipkit_dir* are never modified.
    All curated artifacts are written under *curated_dir* (default:
    ``<outdir>/curated/``):

    curated/
    ├── trees/              curated .treefile per HMM
    ├── fasta_per_hmm/      pruned raw FASTA per HMM
    ├── alignments_clipkit/ pruned clipkit alignment per HMM
    ├── manifests/          per-HMM JSON manifest
    └── audit/              append-only JSONL audit log
    """
    print("\n[curate] Running automated tree curation (Outlier detection)...")

    cur_cfg = cfg.get("curate", {})
    if not cur_cfg.get("enabled", False):
        print("[curate] Curate step is disabled.")
        return

    use_treeshrink = cur_cfg.get("use_treeshrink", True)
    use_esm_filter = cur_cfg.get("use_esm_filter", True)

    # ── Overlay layout ──────────────────────────────────────────────────────
    if curated_dir is None:
        # Derive from summary_dir's parent (outdir)
        curated_dir = os.path.join(os.path.dirname(summary_dir), "curated")

    cur_tree_dir    = os.path.join(curated_dir, "trees")
    cur_fasta_dir   = os.path.join(curated_dir, "fasta_per_hmm")
    cur_clipkit_dir = os.path.join(curated_dir, "alignments_clipkit")
    cur_manifest_dir = os.path.join(curated_dir, "manifests")
    cur_audit_dir   = os.path.join(curated_dir, "audit")

    for d in [cur_tree_dir, cur_fasta_dir, cur_clipkit_dir,
              cur_manifest_dir, cur_audit_dir]:
        safe_mkdir(d)

    tmp_dir = os.path.join(summary_dir, "curate_tmp")
    safe_mkdir(tmp_dir)

    audit_fp = os.path.join(
        cur_audit_dir,
        f"curation_run_{datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S')}.jsonl",
    )

    # ── Discover raw trees ──────────────────────────────────────────────────
    trees = sorted(glob.glob(os.path.join(tree_dir, "*.treefile")))
    # Exclude legacy .raw.treefile / .treeshrink.treefile backups that may
    # still exist in the raw directory from previous pipeline runs.
    trees = [t for t in trees
             if not t.endswith(".raw.treefile")
             and not t.endswith(".treeshrink.treefile")]

    if hmm_keep is not None:
        trees = [t for t in trees
                 if os.path.basename(t).replace(".treefile", "") in hmm_keep]

    for raw_tree_fp in trees:
        hmm = os.path.basename(raw_tree_fp).replace(".treefile", "")
        print(f"  Curating {hmm}...")

        current_tree_fp = raw_tree_fp
        filters_applied = []
        dropped_taxa = []

        # ── ESM Filter ──────────────────────────────────────────────────────
        if use_esm_filter:
            emb_fp = os.path.join(emb_dir, f"{hmm}.parquet")
            esm_out_fp = os.path.join(tmp_dir, f"{hmm}.esm_filtered.treefile")
            res = run_esm_filter(current_tree_fp, emb_fp, esm_out_fp)
            if res:
                # Track which taxa were dropped by this filter
                before = {t.name for t in Phylo.read(current_tree_fp, "newick").get_terminals()}
                after  = {t.name for t in Phylo.read(res, "newick").get_terminals()}
                dropped_taxa.extend(sorted(before - after))
                filters_applied.append("esm_filter")
                current_tree_fp = res

        # ── TreeShrink ──────────────────────────────────────────────────────
        if use_treeshrink:
            res = run_treeshrink(current_tree_fp, hmm, tmp_dir)
            if res:
                before = {t.name for t in Phylo.read(current_tree_fp, "newick").get_terminals()}
                after  = {t.name for t in Phylo.read(res, "newick").get_terminals()}
                dropped_taxa.extend(sorted(before - after))
                filters_applied.append("treeshrink")
                current_tree_fp = res

        # ── Write curated overlay outputs ───────────────────────────────────
        out_tree_fp    = os.path.join(cur_tree_dir,    f"{hmm}.treefile")
        out_fasta_fp   = os.path.join(cur_fasta_dir,   f"{hmm}.fasta")
        out_clipkit_fp = os.path.join(cur_clipkit_dir, f"{hmm}.clipkit.faa")

        curated_tree = Phylo.read(current_tree_fp, "newick")
        valid_terminals = {t.name for t in curated_tree.get_terminals()}

        Phylo.write(curated_tree, out_tree_fp, "newick")
        _prune_fasta(os.path.join(fasta_dir,   f"{hmm}.fasta"),        out_fasta_fp,   valid_terminals)
        _prune_fasta(os.path.join(clipkit_dir, f"{hmm}.clipkit.faa"),  out_clipkit_fp, valid_terminals)

        was_modified = bool(filters_applied)

        # ── Manifest ────────────────────────────────────────────────────────
        timestamp = datetime.datetime.now(datetime.timezone.utc).isoformat()
        manifest = {
            "hmm": hmm,
            "timestamp": timestamp,
            "source": "automated",
            "raw_tree": raw_tree_fp,
            "curated_tree": out_tree_fp,
            "filters_applied": filters_applied,
            "dropped_taxa": dropped_taxa,
            "retained_taxa_count": len(valid_terminals),
            "was_modified": was_modified,
        }
        _write_manifest(os.path.join(cur_manifest_dir, f"{hmm}.json"), manifest)
        _append_audit(audit_fp, manifest)

        if was_modified:
            print(f"    -> Pruned {len(dropped_taxa)} outlier(s) from {hmm} "
                  f"(filters: {', '.join(filters_applied)}).")
        else:
            print(f"    -> No outliers detected for {hmm}.")

    print(f"[curate] Curated outputs written to: {curated_dir}")
    return curated_dir

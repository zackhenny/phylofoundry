import os
import glob
import shutil
import pandas as pd
from Bio import Phylo, SeqIO
from ..utils.helpers import run_cmd, safe_mkdir

def _prune_fasta(in_fp, out_fp, keep_ids):
    if not os.path.exists(in_fp): return
    # Also backup original if not backed up
    raw_fp = in_fp.replace(os.path.basename(in_fp), os.path.basename(in_fp).replace(".", ".raw.", 1))
    if not os.path.exists(raw_fp):
        shutil.copy2(in_fp, raw_fp)
    
    recs = []
    for r in SeqIO.parse(raw_fp, "fasta"):
        if r.id in keep_ids:
            recs.append(r)
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

def run_curate(cfg, tree_dir, fasta_dir, clipkit_dir, emb_dir, summary_dir, hmm_keep, force=False):
    print("\n[curate] Running automated tree curation (Outlier detection)...")
    
    cur_cfg = cfg.get("curate", {})
    if not cur_cfg.get("enabled", False):
        print("[curate] Curate step is disabled.")
        return

    use_treeshrink = cur_cfg.get("use_treeshrink", True)
    use_esm_filter = cur_cfg.get("use_esm_filter", True)
    
    tmp_dir = os.path.join(summary_dir, "curate_tmp")
    safe_mkdir(tmp_dir)
    
    trees = sorted(glob.glob(os.path.join(tree_dir, "*.treefile")))
    trees = [t for t in trees if not t.endswith(".raw.treefile") and not t.endswith(".treeshrink.treefile")]
    
    if hmm_keep is not None:
        trees = [t for t in trees if os.path.basename(t).replace(".treefile", "") in hmm_keep]
        
    for tree_fp in trees:
        hmm = os.path.basename(tree_fp).replace(".treefile", "")
        print(f"  Curating {hmm}...")
        
        # 1. Backup raw tree
        raw_tree_fp = tree_fp.replace(".treefile", ".raw.treefile")
        if not os.path.exists(raw_tree_fp):
            shutil.copy2(tree_fp, raw_tree_fp)
        
        current_tree_fp = raw_tree_fp
        
        # 2. ESM Filter
        if use_esm_filter:
            emb_fp = os.path.join(emb_dir, f"{hmm}.parquet")
            esm_out_fp = os.path.join(tmp_dir, f"{hmm}.esm_filtered.treefile")
            res = run_esm_filter(current_tree_fp, emb_fp, esm_out_fp)
            if res:
                current_tree_fp = res
                
        # 3. TreeShrink
        if use_treeshrink:
            res = run_treeshrink(current_tree_fp, hmm, tmp_dir)
            if res:
                current_tree_fp = res
        
        # Write final curated tree back to primary file
        if current_tree_fp != raw_tree_fp: # Modifications were made
            shutil.copy2(current_tree_fp, tree_fp)
            
            # Sync alignments/fastas
            curated_tree = Phylo.read(tree_fp, "newick")
            valid_terminals = set([t.name for t in curated_tree.get_terminals()])
            
            _prune_fasta(os.path.join(fasta_dir, f"{hmm}.fasta"), os.path.join(fasta_dir, f"{hmm}.fasta"), valid_terminals)
            _prune_fasta(os.path.join(clipkit_dir, f"{hmm}.clipkit.faa"), os.path.join(clipkit_dir, f"{hmm}.clipkit.faa"), valid_terminals)
            
            print(f"    -> Pruned outliers from {hmm}. Downstream seqs synced.")
        else:
            print(f"    -> No outliers detected for {hmm}.")

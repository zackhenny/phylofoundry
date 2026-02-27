"""Post-processing module using scikit-bio for conservation,
KL divergence, and MRCA sanity checks."""

import os
import glob
import math
import shutil
import subprocess
import tempfile
import numpy as np
import pandas as pd
from collections import defaultdict, Counter
from ..utils.bio import read_fasta
from ..utils.helpers import safe_mkdir

AA_ALPHABET = set(list("ACDEFGHIKLMNPQRSTVWY"))
GAP_CHARS = set(["-", ".", "X", "x", "?", "*"])


def _import_skbio():
    """Lazy import of scikit-bio to avoid crashing when not installed."""
    try:
        from skbio import TreeNode, TabularMSA
        from skbio.sequence import Protein
        return TreeNode, TabularMSA, Protein
    except ImportError:
        raise ImportError(
            "scikit-bio is required for the 'post' step. "
            "Install it with: conda install -c conda-forge scikit-bio"
        )


def msa_from_fasta_dict(aln_seqs: dict):
    _, TabularMSA, Protein = _import_skbio()
    seqs = []
    for tip, seq in aln_seqs.items():
        # Replace degenerate characters known to break skbio conservation metrics
        for c in "XBZJUOxbzjuo?*":
            seq = seq.replace(c, "-")
        seq = seq.upper()
        seqs.append(Protein(seq, metadata={"id": tip}))
    return TabularMSA(seqs)


def consurf_like_scores(aln_seqs: dict, metric="inverse_shannon_uncertainty"):
    msa = msa_from_fasta_dict(aln_seqs)
    scores = msa.conservation(metric=metric, gap_mode="ignore")
    return list(scores)


def site_counts_from_subset(aln_seqs: dict, subset_tips=None):
    if not aln_seqs:
        return []
    msa = msa_from_fasta_dict(aln_seqs)
    tip_set = set(subset_tips) if subset_tips is not None else None
    msa_ids = [seq.metadata["id"] for seq in msa]
    counts = []
    for pos in msa.iter_positions():
        c = Counter()
        s = str(pos)
        for idx, char in enumerate(s):
            tip = msa_ids[idx]
            if tip_set is not None and tip not in tip_set:
                continue
            aa = char
            if aa in GAP_CHARS:
                continue
            aa = aa.upper()
            if aa in AA_ALPHABET:
                c[aa] += 1
        counts.append(c)
    return counts


def kl_divergence(p_counts, q_counts, pseudocount=1e-6):
    alphabet = list(AA_ALPHABET)
    p_total = sum(p_counts.get(a, 0) for a in alphabet) + pseudocount * len(alphabet)
    q_total = sum(q_counts.get(a, 0) for a in alphabet) + pseudocount * len(alphabet)
    kl = 0.0
    for a in alphabet:
        p = (p_counts.get(a, 0) + pseudocount) / p_total
        q = (q_counts.get(a, 0) + pseudocount) / q_total
        kl += p * math.log(p / q, 2)
    return float(kl)


def load_clades_tsv(clades_tsv):
    df = pd.read_csv(clades_tsv, sep="\t", dtype=str)
    if "clade_name" not in df.columns or "tip" not in df.columns:
        raise SystemExit("clades_tsv must have columns: clade_name, tip")
    out = defaultdict(list)
    for _, r in df.iterrows():
        if pd.isna(r["clade_name"]) or pd.isna(r["tip"]):
            continue
        out[str(r["clade_name"])].append(str(r["tip"]))
    return dict(out)


def tree_load(tree_fp):
    TreeNode, _, _ = _import_skbio()
    t = TreeNode.read(tree_fp, format="newick")
    t.create_caches()
    t.assign_ids()
    return t


def validate_tips_in_tree(tree, tips):
    tip_names = set([n.name for n in tree.tips()])
    missing = [x for x in tips if x not in tip_names]
    return missing


def _load_taxonomy(gtdb_dir, tax_file):
    """Load taxonomy map from GTDB-Tk output or custom TSV."""
    from ..utils.helpers import normalize_genome_id
    tax_map = {}

    # 1. Load from GTDB-Tk output if provided
    if gtdb_dir and os.path.isdir(gtdb_dir):
        # GTDB-Tk outputs: gtdbtk.bac120.summary.tsv, gtdbtk.ar122.summary.tsv
        summary_files = glob.glob(os.path.join(gtdb_dir, "gtdbtk.*.summary.tsv"))
        for fp in summary_files:
            try:
                df = pd.read_csv(fp, sep="\t")
                # Expected cols: user_genome, classification
                if "user_genome" in df.columns and "classification" in df.columns:
                    for _, r in df.iterrows():
                        key = normalize_genome_id(str(r["user_genome"]))
                        tax_map[key] = str(r["classification"])
            except Exception as e:
                print(f"[post] Warning: Failed to parse GTDB summary {fp}: {e}")

    # 2. Load from custom taxonomy file if provided (overrides GTDB)
    if tax_file and os.path.exists(tax_file):
        try:
            df = pd.read_csv(tax_file, sep="\t", dtype=str)
            if "genome" in df.columns and "lineage" in df.columns:
                 for _, r in df.iterrows():
                    key = normalize_genome_id(str(r["genome"]))
                    tax_map[key] = str(r["lineage"])
        except Exception as e:
             print(f"[post] Warning: Failed to parse taxonomy file {tax_file}: {e}")
             
    return tax_map


def _extract_taxon_label(lineage, level):
    rank_map = {
        "domain": "d__", "phylum": "p__", "class": "c__", "order": "o__",
        "family": "f__", "genus": "g__", "species": "s__"
    }
    level = str(level or "genus").strip().lower()
    prefix = rank_map.get(level, level)
    if len(prefix) == 1:
        prefix = f"{prefix}__"
    if "__" not in prefix and len(prefix) <= 3:
        prefix = f"{prefix}__"
    for token in str(lineage or "").split(";"):
        token = token.strip()
        if token.lower().startswith(prefix.lower()) and token:
            return token
    return None


def _detect_taxonomy_clades(summary_dir, tax_map, level):
    from ..utils.helpers import normalize_genome_id
    best_hits_fp = os.path.join(summary_dir, "best_hits.competitive.tsv")
    if not os.path.exists(best_hits_fp):
        return {}
    df = pd.read_csv(best_hits_fp, sep="\t", dtype=str)
    if "genome" not in df.columns or "protein" not in df.columns:
        return {}
    clades = defaultdict(list)
    for _, row in df.iterrows():
        genome = str(row["genome"])
        protein = str(row["protein"])
        tip = f"{genome}|{protein}"
        lineage = tax_map.get(normalize_genome_id(genome))
        label = _extract_taxon_label(lineage, level)
        if label:
            clades[label].append(tip)
    return dict(clades)


def _detect_treecluster_clades(tree_dir, threshold, method):
    treecluster_bin = shutil.which("TreeCluster.py") or shutil.which("treecluster")
    if not treecluster_bin:
        print("[post] TreeCluster not found in PATH; skipping TreeCluster clade detection.")
        return {}

    clades = defaultdict(list)
    for tree_fp in sorted(glob.glob(os.path.join(tree_dir, "*.treefile"))):
        hmm = os.path.basename(tree_fp).split(".")[0]
        with tempfile.NamedTemporaryFile(suffix=".tsv", delete=False) as tmp:
            out_fp = tmp.name
        try:
            cmd = [
                treecluster_bin, "-i", tree_fp, "-o", out_fp,
                "-t", str(threshold), "-m", str(method)
            ]
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            df = pd.read_csv(out_fp, sep="\t", dtype=str)
            if len(df.columns) < 2:
                continue
            tip_col = "SequenceName" if "SequenceName" in df.columns else df.columns[0]
            cluster_col = "ClusterNumber" if "ClusterNumber" in df.columns else df.columns[1]
            for _, row in df.iterrows():
                cluster_id = str(row[cluster_col])
                if cluster_id == "-1":
                    continue
                clade_name = f"{hmm}|cluster_{cluster_id}"
                clades[clade_name].append(str(row[tip_col]))
        except Exception as e:
            print(f"[post] TreeCluster failed for {tree_fp}: {e}")
        finally:
            try:
                os.remove(out_fp)
            except OSError:
                pass
    return dict(clades)


def _write_detected_clades(summary_dir, clades):
    if not clades:
        return
    rows = []
    for cname, tips in clades.items():
        for tip in tips:
            rows.append({"clade_name": cname, "tip": tip})
    if rows:
        pd.DataFrame(rows).to_csv(os.path.join(summary_dir, "detected_clades.tsv"), sep="\t", index=False)


def _load_embedding_coords(emb_dir, hmm, n_pcs):
    """Load per-tip embedding coordinates for one HMM, preferring precomputed PCA."""
    pca_fp = os.path.join(emb_dir, f"{hmm}.pca.tsv")
    if os.path.exists(pca_fp):
        df = pd.read_csv(pca_fp, sep="\t", dtype={"tip": str})
        pc_cols = [c for c in df.columns if c.startswith("PC")]
        if "tip" not in df.columns or not pc_cols:
            return None
        use_cols = pc_cols[:max(1, int(n_pcs))]
        M = df[use_cols].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
        keep = ~np.isnan(M).any(axis=1)
        return dict(zip(df.loc[keep, "tip"].astype(str), M[keep]))

    vec_fp = os.path.join(emb_dir, f"{hmm}.vectors.tsv")
    if os.path.exists(vec_fp):
        from sklearn.decomposition import PCA

        df = pd.read_csv(vec_fp, sep="\t", dtype={"tip": str})
        vec_cols = [c for c in df.columns if c.startswith("d")]
        if "tip" not in df.columns or not vec_cols:
            return None
        X = df[vec_cols].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
        keep = ~np.isnan(X).any(axis=1)
        X = X[keep]
        tips = df.loc[keep, "tip"].astype(str).tolist()
        if X.shape[0] < 2:
            return None
        ncomp = min(max(2, int(n_pcs)), X.shape[0] - 1, X.shape[1])
        Z = PCA(n_components=ncomp, random_state=0).fit_transform(X)
        return dict(zip(tips, Z))

    return None


def _parse_support(node):
    v = getattr(node, "support", None)
    if v is not None:
        try:
            return float(v)
        except (TypeError, ValueError):
            pass
    if node.name is not None:
        try:
            return float(str(node.name))
        except (TypeError, ValueError):
            return None
    return None


def _distance_fn(metric):
    metric = str(metric or "euclidean").strip().lower()

    def euclidean(A, b):
        return np.sqrt(((A - b) ** 2).sum(axis=1))

    def cosine(A, b):
        A_norm = np.linalg.norm(A, axis=1)
        b_norm = np.linalg.norm(b)
        denom = np.clip(A_norm * max(b_norm, 1e-12), 1e-12, None)
        return 1.0 - np.clip((A @ b) / denom, -1.0, 1.0)

    return cosine if metric == "cosine" else euclidean


def _point_distance(metric, a, b):
    if metric == "cosine":
        na = np.linalg.norm(a)
        nb = np.linalg.norm(b)
        denom = max(na * nb, 1e-12)
        return float(1.0 - np.clip(np.dot(a, b) / denom, -1.0, 1.0))
    return float(np.linalg.norm(a - b))


def _approx_tree_diameter(node):
    tips = list(node.tips())
    if len(tips) < 2:
        return 0.0
    try:
        a = tips[0]
        b = max(tips, key=lambda t: a.distance(t))
        c = max(tips, key=lambda t: b.distance(t))
        return float(b.distance(c))
    except Exception:
        return None


def _detect_tree_embed_clades(tree_dir, emb_dir, post_cfg, hmm_keep=None):
    """Detect clades by embedding shift across internal tree nodes."""
    support_min = float(post_cfg.get("embedtree_support_min", 80))
    min_size = int(post_cfg.get("embedtree_min_size", 5))
    max_size_cfg = post_cfg.get("embedtree_max_size", 5000)
    max_size = None if max_size_cfg in [None, "", 0] else int(max_size_cfg)
    top_k = int(post_cfg.get("embedtree_top_k", 10))
    n_pcs = int(post_cfg.get("embedtree_pcs", 10))
    distance = str(post_cfg.get("embedtree_distance", "euclidean")).strip().lower()
    allow_nested = bool(post_cfg.get("embedtree_allow_nested", False))
    emit_all = bool(post_cfg.get("embedtree_emit_all", False))
    eps = 1e-8

    per_hmm_selected = []
    score_rows = []

    tree_files = sorted(glob.glob(os.path.join(tree_dir, "*.treefile")))
    if hmm_keep is not None:
        tree_files = [fp for fp in tree_files if os.path.basename(fp).split(".")[0] in hmm_keep]

    for tree_fp in tree_files:
        hmm = os.path.basename(tree_fp).split(".")[0]
        coords = _load_embedding_coords(emb_dir, hmm, n_pcs=n_pcs)
        if not coords:
            continue

        tree = tree_load(tree_fp)
        tip_names = {n.name for n in tree.tips() if n.name is not None}
        shared = tip_names.intersection(coords.keys())
        if not shared:
            continue

        dist_to_centroid = _distance_fn(distance)
        candidates = []
        traversal = 0

        for node in tree.preorder(include_self=True):
            if node.is_tip():
                continue
            traversal += 1
            tips = [t.name for t in node.tips() if t.name in shared]
            size = len(tips)
            if size < min_size:
                continue
            if max_size is not None and size > max_size:
                continue

            support = _parse_support(node)
            if support is not None and support < support_min:
                continue

            children = [c for c in node.children if not c.is_tip() or c.name in shared]
            if len(children) < 2:
                continue
            child_tips = []
            for child in children:
                ct = [t.name for t in child.tips() if t.name in shared]
                if ct:
                    child_tips.append(ct)
            if len(child_tips) < 2:
                continue
            child_tips.sort(key=len, reverse=True)
            left_tips, right_tips = child_tips[0], child_tips[1]

            M = np.vstack([coords[t] for t in tips])
            centroid = M.mean(axis=0)
            within_dispersion = float(np.median(dist_to_centroid(M, centroid)))

            L = np.vstack([coords[t] for t in left_tips])
            R = np.vstack([coords[t] for t in right_tips])
            lc = L.mean(axis=0)
            rc = R.mean(axis=0)
            child_sep = _point_distance(distance, lc, rc)
            left_within = float(np.median(dist_to_centroid(L, lc)))
            right_within = float(np.median(dist_to_centroid(R, rc)))
            pooled = (left_within + right_within) / 2.0
            effect = float(child_sep / (pooled + eps))
            diameter = _approx_tree_diameter(node)

            node_label = f"{hmm}|node_{node.id}"
            score_rows.append({
                "node_id": node_label,
                "clade_size": int(size),
                "support_min": support,
                "within_dispersion": within_dispersion,
                "child_separation": child_sep,
                "effect_size": effect,
                "tree_diameter": diameter,
                "hmm_name": hmm,
            })
            candidates.append({
                "node": node,
                "tips": tips,
                "effect_size": effect,
                "traversal": traversal,
                "hmm": hmm,
            })

        candidates.sort(key=lambda x: (-x["effect_size"], x["traversal"]))
        chosen = []
        chosen_ids = set()
        for cand in candidates:
            node = cand["node"]
            if not allow_nested:
                p = node.parent
                blocked = False
                while p is not None:
                    if p.id in chosen_ids:
                        blocked = True
                        break
                    p = p.parent
                if blocked:
                    continue
            chosen.append(cand)
            chosen_ids.add(node.id)
            if not emit_all and len(chosen) >= top_k:
                break

        per_hmm_selected.extend(chosen)

    per_hmm_selected.sort(key=lambda x: (-x["effect_size"], x["hmm"], x["traversal"]))
    clades = {}
    for i, sel in enumerate(per_hmm_selected, start=1):
        cname = f"embed__C{i:05d}"
        clades[cname] = sel["tips"]
    return clades, score_rows



def run_post(cfg, tree_dir, clipkit_dir, aln_dir, post_dir, summary_dir, hmm_keep, force=False):
    print("\n[post] scikit-bio post-processing...")

    # ── Taxonomy Integration ──────────────────────────────────────────────
    gtdb_dir = cfg["inputs"].get("gtdb_dir")
    tax_file = cfg["inputs"].get("taxonomy_file")
    
    if gtdb_dir or tax_file:
        tax_map = _load_taxonomy(gtdb_dir, tax_file)
        if tax_map:
            print(f"[post] Loaded taxonomy for {len(tax_map)} genomes.")
            # 1. Save genome->taxonomy map
            tax_rows = [{"genome": k, "taxonomy": v} for k, v in tax_map.items()]
            pd.DataFrame(tax_rows).to_csv(os.path.join(summary_dir, "genome_taxonomy.tsv"), sep="\t", index=False)
            
            # 2. Merge into best_hits
            best_hits_fp = os.path.join(summary_dir, "best_hits.competitive.tsv")
            if os.path.exists(best_hits_fp):
                from ..utils.helpers import normalize_genome_id
                df = pd.read_csv(best_hits_fp, sep="\t")
                # Normalize genome column to match taxonomy keys
                # We assume df["genome"] is filename like "genomeA.faa"
                # We apply normalize to look up
                
                def lookup_tax(g_filename):
                    norm = normalize_genome_id(str(g_filename))
                    return tax_map.get(norm, "Unknown")

                if "genome" in df.columns:
                    df["taxonomy"] = df["genome"].apply(lookup_tax)
                    out_fp = os.path.join(summary_dir, "best_hits.with_taxonomy.tsv")
                    df.to_csv(out_fp, sep="\t", index=False)
                    print(f"[post] Wrote {out_fp}")
    
    post_cfg = cfg.get("post", {})
    emb_dir = os.path.join(os.path.dirname(summary_dir), "embeddings")
    clades = None
    embed_node_scores = []
    if post_cfg.get("clades_tsv", None):
        clades = load_clades_tsv(post_cfg["clades_tsv"])
    else:
        detect_method = str(post_cfg.get("detect_clades_method", "")).strip().lower()
        if detect_method == "taxonomy":
            clades = _detect_taxonomy_clades(summary_dir, tax_map if (gtdb_dir or tax_file) else {}, post_cfg.get("taxonomy_clade_level", "genus"))
        elif detect_method == "treecluster":
            clades = _detect_treecluster_clades(
                tree_dir,
                post_cfg.get("treecluster_threshold", 0.045),
                post_cfg.get("treecluster_method", "max_clade"),
            )
        elif detect_method == "tree_embed":
            clades, embed_node_scores = _detect_tree_embed_clades(
                tree_dir,
                emb_dir,
                post_cfg,
                hmm_keep=hmm_keep,
            )
        if clades:
            _write_detected_clades(summary_dir, clades)
        if embed_node_scores:
            pd.DataFrame(embed_node_scores).to_csv(
                os.path.join(summary_dir, "node_scores.embedtree.tsv"),
                sep="\t",
                index=False,
            )

    if not clades:
        detected_clades_fp = os.path.join(summary_dir, "detected_clades.tsv")
        if os.path.exists(detected_clades_fp):
            try:
                clades = load_clades_tsv(detected_clades_fp)
                print(f"[post] Loaded clades from {detected_clades_fp}")
            except Exception:
                clades = None

    if post_cfg.get("compute_kl", False) and not clades:
        raise SystemExit(
            "post.compute_kl requires clades from post.clades_tsv, "
            "post.detect_clades_method, or summary/detected_clades.tsv"
        )

    kl_pairs = []
    has_manual_kl_pairs = bool(post_cfg.get("kl_pairs", None))
    if post_cfg.get("kl_pairs", None):
        for pair in str(post_cfg["kl_pairs"]).split(","):
            pair = pair.strip()
            if not pair:
                continue
            if ":" not in pair:
                raise SystemExit("post.kl_pairs must look like 'A:B' or 'A:background'")
            a, b = pair.split(":", 1)
            kl_pairs.append((a.strip(), b.strip()))

    cons_rows = []
    kl_rows = []
    mrca_rows = []

    hmm_names = sorted([os.path.basename(x).split(".")[0] for x in glob.glob(os.path.join(tree_dir, "*.treefile"))])
    if hmm_keep is not None:
        hmm_names = [h for h in hmm_names if h in hmm_keep]

    for hmm in hmm_names:
        tree_fp = os.path.join(tree_dir, f"{hmm}.treefile")

        clip_aln = os.path.join(clipkit_dir, f"{hmm}.clipkit.faa")
        if os.path.exists(clip_aln):
            aln_fp = clip_aln
        else:
            cand1 = os.path.join(aln_dir, f"{hmm}.afa")
            cand2 = os.path.join(aln_dir, f"{hmm}.mafft.fasta")
            aln_fp = cand2 if os.path.exists(cand2) else cand1 if os.path.exists(cand1) else None

        if not aln_fp or not os.path.exists(tree_fp):
            continue

        aln_seqs = read_fasta(aln_fp)

        if post_cfg.get("compute_conservation", False):
            try:
                cons = consurf_like_scores(aln_seqs, metric=post_cfg.get("conservation_metric", "inverse_shannon_uncertainty"))
                for i, v in enumerate(cons, start=1):
                    cons_rows.append({"hmm": hmm, "scope": "global", "site_1based": i, "conservation": float(v)})
            except Exception as e:
                import sys
                print(f"[post] Conservation scoring failed for {hmm}: {e}", file=sys.stderr)

        if post_cfg.get("compute_kl", False) and clades:
            counts = {}
            counts["background"] = site_counts_from_subset(aln_seqs, subset_tips=None)
            dynamic_pairs = []
            for cname, tips in clades.items():
                tips_in_aln = [x for x in tips if x in aln_seqs]
                counts[cname] = site_counts_from_subset(aln_seqs, subset_tips=tips_in_aln)
                if not has_manual_kl_pairs and tips_in_aln:
                    others_key = f"others::{cname}"
                    others_tips = [x for x in aln_seqs if x not in set(tips_in_aln)]
                    if others_tips:
                        counts[others_key] = site_counts_from_subset(aln_seqs, subset_tips=others_tips)
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
                    kl_rows.append({
                        "hmm": hmm,
                        "pair": f"{a}:{b}",
                        "site_1based": i + 1,
                        "kl_bits": kl,
                        "nA": sum(counts[a][i].values()),
                        "nB": sum(counts[b][i].values())
                    })

        # MRCA sanity only (no internal node mapping)
        if clades:
            try:
                t = tree_load(tree_fp)
                for cname, tips in clades.items():
                    tips2 = [x for x in tips if x in aln_seqs]
                    missing = validate_tips_in_tree(t, tips2)
                    if missing:
                        mrca_rows.append({"hmm": hmm, "clade": cname, "status": "missing_tips", "detail": ",".join(missing[:10])})
                        continue
                    mrca = t.lca(tips2)
                    mrca_rows.append({"hmm": hmm, "clade": cname, "status": "ok", "mrca_id": mrca.id, "mrca_name": mrca.name or ""})
            except Exception as e:
                import sys
                print(f"[post] MRCA sanity checks failed for {hmm}: {e}", file=sys.stderr)

    if cons_rows:
        pd.DataFrame(cons_rows).to_csv(os.path.join(post_dir, "conservation.tsv"), sep="\t", index=False)
    if kl_rows:
        pd.DataFrame(kl_rows).to_csv(os.path.join(post_dir, "kl_divergence.tsv"), sep="\t", index=False)
    if mrca_rows:
        pd.DataFrame(mrca_rows).to_csv(os.path.join(post_dir, "mrca_sanity.tsv"), sep="\t", index=False)

"""Post-processing module using scikit-bio for conservation,
KL divergence, and MRCA sanity checks."""

import os
import glob
import math
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
    seqs = [Protein(seq, metadata={"id": tip}) for tip, seq in aln_seqs.items()]
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
    clades = None
    if post_cfg.get("clades_tsv", None):
        clades = load_clades_tsv(post_cfg["clades_tsv"])

    if post_cfg.get("compute_kl", False) and not clades:
        raise SystemExit("post.compute_kl requires post.clades_tsv")

    kl_pairs = []
    if post_cfg.get("kl_pairs", None):
        for pair in str(post_cfg["kl_pairs"]).split(","):
            pair = pair.strip()
            if not pair:
                continue
            if ":" not in pair:
                raise SystemExit("post.kl_pairs must look like 'A:B' or 'A:background'")
            a, b = pair.split(":", 1)
            kl_pairs.append((a.strip(), b.strip()))
    else:
        if post_cfg.get("compute_kl", False) and clades and len(clades) >= 2:
            keys = list(clades.keys())
            for i in range(len(keys)):
                for j in range(i + 1, len(keys)):
                    kl_pairs.append((keys[i], keys[j]))
        if post_cfg.get("compute_kl", False) and clades and len(clades) == 1:
            k = list(clades.keys())[0]
            kl_pairs.append((k, "background"))

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
            for cname, tips in clades.items():
                counts[cname] = site_counts_from_subset(aln_seqs, subset_tips=tips)

            L = len(counts["background"])
            for a, b in kl_pairs:
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

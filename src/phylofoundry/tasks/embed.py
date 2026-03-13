import os
import glob
import math
import subprocess
import shutil
import warnings
import pandas as pd
import numpy as np
from ..utils.bio import read_fasta, write_fasta
from ..utils.helpers import safe_mkdir, write_json

def _pca_fit_transform(X, n_components=3):
    from sklearn.decomposition import PCA
    pca = PCA(n_components=n_components, random_state=0)
    Z = pca.fit_transform(X)
    return Z, pca.explained_variance_ratio_.tolist()

def _embed_esm(seqs: dict, model_name: str, device: str, batch_size: int, repr_layer, model_dir: str | None = None):
    import torch
    import esm

    if model_dir is not None:
        os.makedirs(model_dir, exist_ok=True)
        torch.hub.set_dir(model_dir)

    model, alphabet = esm.pretrained.load_model_and_alphabet(model_name)
    model.eval()
    model = model.to(device)

    if repr_layer is None:
        repr_layer = model.num_layers

    batch_converter = alphabet.get_batch_converter()

    items = list(seqs.items())
    vectors = []
    ids = []

    with torch.no_grad():
        for i in range(0, len(items), batch_size):
            chunk = items[i:i + batch_size]
            batch = [(k, v) for k, v in chunk]
            labels, strs, toks = batch_converter(batch)
            toks = toks.to(device)
            out = model(toks, repr_layers=[repr_layer], return_contacts=False)
            reps = out["representations"][repr_layer]  # (B, T, C)

            pad_idx = alphabet.padding_idx
            mask = (toks != pad_idx).float()
            mask[:, 0] = 0.0
            eos_idx = alphabet.eos_idx
            mask = mask * (toks != eos_idx).float()

            denom = mask.sum(dim=1).clamp(min=1.0).unsqueeze(-1)  # (B,1)
            pooled = (reps * mask.unsqueeze(-1)).sum(dim=1) / denom  # (B,C)

            pooled = pooled.detach().cpu().numpy()
            for (k, _), v in zip(chunk, pooled):
                ids.append(k)
                vectors.append(v)

    X = np.vstack(vectors)
    return ids, X

def _embed_transformers(seqs: dict, model_id_or_path: str, device: str, batch_size: int, model_dir: str | None = None):
    import torch
    from transformers import AutoTokenizer, AutoModel

    cache_dir = model_dir  # HuggingFace uses cache_dir; local paths are passed as model_id_or_path
    tok = AutoTokenizer.from_pretrained(model_id_or_path, do_lower_case=False, cache_dir=cache_dir)
    model = AutoModel.from_pretrained(model_id_or_path, cache_dir=cache_dir)
    model.eval()
    model = model.to(device)

    items = list(seqs.items())
    vectors = []
    ids = []

    with torch.no_grad():
        for i in range(0, len(items), batch_size):
            chunk = items[i:i + batch_size]
            labels = [k for k, _ in chunk]
            strings = [v for _, v in chunk]
            enc = tok(strings, return_tensors="pt", padding=True, truncation=True)
            enc = {k: v.to(device) for k, v in enc.items()}
            out = model(**enc)
            reps = out.last_hidden_state  # (B,T,C)

            attn = enc.get("attention_mask", None)
            if attn is None:
                raise RuntimeError("transformers embedding requires attention_mask")
            mask = attn.float()  # (B,T)
            denom = mask.sum(dim=1).clamp(min=1.0).unsqueeze(-1)
            pooled = (reps * mask.unsqueeze(-1)).sum(dim=1) / denom  # (B,C)

            pooled = pooled.detach().cpu().numpy()
            for k, v in zip(labels, pooled):
                ids.append(k)
                vectors.append(v)

    X = np.vstack(vectors)
    return ids, X

def _save_umap_plot(U, ids, hmm_name, out_png, clades=None, cluster_labels=None, title_suffix=""):
    """Generate and save a UMAP scatter plot as PNG. Supports 2D and 3D projections."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        n_dims = U.shape[1] if U.ndim > 1 else 2
        is_3d = (n_dims >= 3)

        if is_3d:
            fig = plt.figure(figsize=(9, 7))
            ax = fig.add_subplot(111, projection="3d")
        else:
            fig, ax = plt.subplots(figsize=(8, 6))

        def _scatter(pts, color, label=None):
            if is_3d:
                ax.scatter(pts[:, 0], pts[:, 1], pts[:, 2], c=[color], s=20, alpha=0.7,
                           label=label, edgecolors="none")
            else:
                ax.scatter(pts[:, 0], pts[:, 1], c=[color], s=20, alpha=0.7,
                           label=label, edgecolors="none")

        if cluster_labels is not None:
            unique_labels = sorted(set(cluster_labels))
            cmap = plt.cm.get_cmap("tab20", max(len(unique_labels), 1))
            for label in unique_labels:
                mask = [cl == label for cl in cluster_labels]
                pts = U[[i for i, m in enumerate(mask) if m]]
                color = "lightgrey" if label == -1 else cmap(unique_labels.index(label) % 20)
                lbl = "Noise" if label == -1 else f"Cluster {label}"
                _scatter(pts, color, label=lbl)
            ax.legend(fontsize=7, loc="best", framealpha=0.7, markerscale=1.5)
        elif clades:
            # Color by user-provided clades
            tip_to_clade = {}
            for cname, tips in clades.items():
                for t in tips:
                    tip_to_clade[t] = cname
            clade_names = sorted(set(tip_to_clade.values()))
            cmap = plt.cm.get_cmap("tab10", max(len(clade_names), 1))
            assigned = [tip_to_clade.get(t, None) for t in ids]
            for ci, cn in enumerate(clade_names):
                mask = [a == cn for a in assigned]
                pts = U[[i for i, m in enumerate(mask) if m]]
                if len(pts) > 0:
                    _scatter(pts, cmap(ci), label=cn)
            # Unassigned
            mask = [a is None for a in assigned]
            pts = U[[i for i, m in enumerate(mask) if m]]
            if len(pts) > 0:
                _scatter(pts, "lightgrey", label="unassigned")
            ax.legend(fontsize=7, loc="best", framealpha=0.7, markerscale=1.5)
        else:
            if is_3d:
                ax.scatter(U[:, 0], U[:, 1], U[:, 2], c="steelblue", s=20, alpha=0.7, edgecolors="none")
            else:
                ax.scatter(U[:, 0], U[:, 1], c="steelblue", s=20, alpha=0.7, edgecolors="none")

        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
        if is_3d:
            ax.set_zlabel("UMAP 3")
        ax.set_title(f"{hmm_name} — UMAP{title_suffix}")
        fig.tight_layout()
        fig.savefig(out_png, dpi=200)
        plt.close(fig)
        print(f"[embed] Saved UMAP plot: {out_png}")
    except ImportError:
        import sys
        print("[embed] matplotlib not installed — skipping UMAP plot.", file=sys.stderr)
    except Exception as e:
        import sys
        print(f"[embed] UMAP plot failed for {hmm_name}: {e}", file=sys.stderr)


def _run_hdbscan(X, min_cluster_size=5, metric="euclidean"):
    """Cluster embeddings with HDBSCAN. Returns list of integer labels (-1 = noise).

    Args:
        X: Feature matrix (n_samples x n_features).
        min_cluster_size: Minimum cluster size for HDBSCAN.
        metric: Distance metric (e.g. 'euclidean', 'cosine').
    """
    try:
        from sklearn.cluster import HDBSCAN
        clusterer = HDBSCAN(min_cluster_size=min_cluster_size, metric=metric)
        labels = clusterer.fit_predict(X)
        return labels.tolist()
    except ImportError:
        import sys
        print("[embed] HDBSCAN requires scikit-learn >= 1.3. Skipping clustering.", file=sys.stderr)
        return None
    except Exception as e:
        import sys
        print(f"[embed] HDBSCAN clustering failed: {e}", file=sys.stderr)
        return None


def _run_leiden(X, n_neighbors=15, metric="cosine", resolution=1.0):
    """Cluster embeddings with Leiden algorithm on a kNN graph.

    Builds a cosine (or other metric) k-nearest-neighbour graph, then runs
    the Leiden community-detection algorithm.  Suitable for large datasets.

    Args:
        X: Feature matrix (n_samples x n_features).
        n_neighbors: Number of neighbours for the kNN graph.
        metric: Distance metric for kNN graph construction (default: 'cosine').
        resolution: Resolution parameter for Leiden (higher → more clusters).

    Returns:
        List of integer cluster labels (0-based), or None on failure.
    """
    try:
        import igraph as ig
        import leidenalg
        from sklearn.neighbors import NearestNeighbors

        n_samples = X.shape[0]
        k = min(n_neighbors, n_samples - 1)

        nn = NearestNeighbors(n_neighbors=k, metric=metric, algorithm="auto")
        nn.fit(X)
        distances, indices = nn.kneighbors(X)

        # Build undirected kNN graph (edges weighted by 1 - distance)
        edges = []
        weights = []
        for i, (nbrs, dists) in enumerate(zip(indices, distances)):
            for j, d in zip(nbrs, dists):
                if j > i:
                    edges.append((i, j))
        weights.append(max(0.0, 1.0 - d))  # similarity weight (clamped to [0, 1])

        g = ig.Graph(n=n_samples, edges=edges)
        g.es["weight"] = weights

        partition = leidenalg.find_partition(
            g,
            leidenalg.RBConfigurationVertexPartition,
            weights="weight",
            resolution_parameter=resolution,
            seed=42,
        )
        labels = [0] * n_samples
        for cluster_id, members in enumerate(partition):
            for node_idx in members:
                labels[node_idx] = cluster_id
        return labels
    except ImportError:
        import sys
        print("[embed] Leiden clustering requires igraph and leidenalg. Skipping.", file=sys.stderr)
        return None
    except Exception as e:
        import sys
        print(f"[embed] Leiden clustering failed: {e}", file=sys.stderr)
        return None


_VALID_CLUSTER_ON = {"pca", "embeddings"}
_VALID_CLUSTER_METHODS = {"hdbscan", "leiden"}
_N_AMINO_ACIDS = 20  # Standard amino acids — used for information-content calculations.

# ── Amino-acid color palette for sequence logos ───────────────────────────────
_AA_COLORS = {
    # Hydrophobic (non-polar aliphatic / aromatic)
    "A": "#f0a030", "V": "#f0a030", "L": "#f0a030", "I": "#f0a030",
    "M": "#f0a030", "F": "#f0a030", "W": "#f0a030", "P": "#c080ff",
    # Polar uncharged
    "S": "#30a830", "T": "#30a830", "N": "#30a830", "Q": "#30a830",
    "C": "#e0c040", "Y": "#30a830",
    # Positively charged
    "R": "#3070d0", "K": "#3070d0", "H": "#8080d0",
    # Negatively charged
    "D": "#d03030", "E": "#d03030",
    # Special
    "G": "#bbbbbb",
    # Gap / unknown
    "-": "#eeeeee", "X": "#aaaaaa",
}


# ── Cluster subworkflow helpers ───────────────────────────────────────────────

def _compute_knn_metrics(Z, ids, labels, k=10):
    """Compute per-sequence kNN neighborhood metrics for cluster interpretation.

    Parameters
    ----------
    Z : np.ndarray
        PCA-reduced feature matrix (n_samples x n_components).
    ids : list
        Sequence identifiers aligned with rows of Z.
    labels : list
        Cluster labels (int; -1 = noise) aligned with rows of Z.
    k : int
        Number of nearest neighbours to consider.

    Returns
    -------
    pd.DataFrame with columns: protein_id, cluster_id, dominant_cluster,
    neighborhood_purity, dist_weighted_purity, mutual_knn_support,
    neighborhood_entropy, median_neighbor_distance
    """
    try:
        from sklearn.neighbors import NearestNeighbors
        from collections import Counter

        k_eff = min(k, len(ids) - 1)
        if k_eff < 1:
            return pd.DataFrame()

        nbrs = NearestNeighbors(n_neighbors=k_eff, metric="euclidean").fit(Z)
        distances, indices = nbrs.kneighbors(Z)

        rows = []
        for i, (tip, my_label) in enumerate(zip(ids, labels)):
            protein = tip.split("|", 1)[1] if "|" in tip else tip
            neighbor_labels = [labels[j] for j in indices[i]]
            neighbor_dists = distances[i]

            # Dominant cluster among all neighbours (noise counts as -1)
            cnt_all = Counter(neighbor_labels)
            nn_nonnoise = [lb for lb in neighbor_labels if lb != -1]
            if nn_nonnoise:
                cnt_nn = Counter(nn_nonnoise)
                dominant_cluster = cnt_nn.most_common(1)[0][0]
                purity = cnt_nn[dominant_cluster] / len(neighbor_labels)
            else:
                dominant_cluster = -1
                purity = 0.0

            # Distance-weighted purity
            weights = 1.0 / (neighbor_dists + 1e-8)
            total_weight = weights.sum()
            if total_weight > 0:
                per_label_weight: dict = {}
                for lb, w in zip(neighbor_labels, weights):
                    per_label_weight[lb] = per_label_weight.get(lb, 0.0) + w
                if per_label_weight:
                    best_lb = max(per_label_weight, key=lambda x: per_label_weight[x])
                    dist_weighted_purity = per_label_weight[best_lb] / total_weight
                else:
                    dist_weighted_purity = 0.0
            else:
                dist_weighted_purity = 0.0

            # Mutual kNN support: fraction of neighbours that also list i in their kNN
            mutual_count = sum(1 for j in indices[i] if i in indices[j])
            mutual_knn_support = mutual_count / k_eff if k_eff > 0 else 0.0

            # Shannon entropy of neighbourhood label distribution
            entropy = 0.0
            total_n = len(neighbor_labels)
            for cnt_val in cnt_all.values():
                p = cnt_val / total_n
                if p > 0:
                    entropy -= p * math.log2(p)

            rows.append({
                "protein_id": protein,
                "cluster_id": int(my_label),
                "dominant_cluster": int(dominant_cluster),
                "neighborhood_purity": round(float(purity), 4),
                "dist_weighted_purity": round(float(dist_weighted_purity), 4),
                "mutual_knn_support": round(float(mutual_knn_support), 4),
                "neighborhood_entropy": round(float(entropy), 4),
                "median_neighbor_distance": round(float(np.median(neighbor_dists)), 6),
            })

        return pd.DataFrame(rows)

    except Exception as e:
        import sys
        print(f"[embed] kNN metrics computation failed: {e}", file=sys.stderr)
        return pd.DataFrame()


def _assign_membership_tiers(ids, labels, knn_df):
    """Assign core/affiliate/bridge/outlier membership tiers to each sequence.

    Rules
    -----
    * **core** — cluster member with neighbourhood purity ≥ 0.8 and dominant
      cluster matching assigned cluster and mutual kNN support ≥ 0.3.
    * **affiliate** — cluster member with purity in [0.5, 0.8) or mutual < 0.3.
    * **bridge** — cluster member whose dominant neighbour cluster differs from
      its own assigned cluster.
    * **noise_peripheral** — noise point (label −1) with purity ≥ 0.5 toward a
      specific cluster.
    * **noise_bridge** — noise point with moderate purity (0.25–0.5).
    * **outlier** — noise point with no clear affiliation.

    Returns
    -------
    dict : {protein_id: tier_string}
    """
    tiers: dict = {}

    if knn_df.empty:
        for tip, label in zip(ids, labels):
            protein = tip.split("|", 1)[1] if "|" in tip else tip
            tiers[protein] = "core" if label >= 0 else "outlier"
        return tiers

    for _, row in knn_df.iterrows():
        pid = row["protein_id"]
        cl = int(row["cluster_id"])
        dom = int(row["dominant_cluster"])
        purity = float(row["neighborhood_purity"])
        mutual = float(row["mutual_knn_support"])

        if cl >= 0:
            if dom != cl:
                tiers[pid] = "bridge"
            elif purity >= 0.8 and mutual >= 0.3:
                tiers[pid] = "core"
            elif purity >= 0.5:
                tiers[pid] = "affiliate"
            else:
                tiers[pid] = "bridge"
        else:
            if dom >= 0 and purity >= 0.5:
                tiers[pid] = "noise_peripheral"
            elif dom >= 0 and purity >= 0.25:
                tiers[pid] = "noise_bridge"
            else:
                tiers[pid] = "outlier"

    return tiers


def _write_cluster_fastas(seqs, ids, labels, tiers, outdir_fasta, force=False):
    """Write per-cluster FASTA files separated by membership tier.

    Parameters
    ----------
    seqs : dict
        {tip_id: sequence} mapping (may use full ``genome|protein`` keys or
        bare protein IDs — both are tried).
    ids : list
        Sequence identifiers in the same order as *labels*.
    labels : list
        Integer cluster labels (−1 = noise).
    tiers : dict
        {protein_id: tier_string} from :func:`_assign_membership_tiers`.
    outdir_fasta : str
        Directory to write FASTA files into.
    force : bool
        Overwrite existing files.

    Returns
    -------
    dict : {cluster_id: {"core": path|None, "affiliate": path|None,
                          "n_core": int, "n_affiliate": int}}
    """
    safe_mkdir(outdir_fasta)

    cluster_core: dict = {}
    cluster_affiliate: dict = {}

    for tip, label in zip(ids, labels):
        protein = tip.split("|", 1)[1] if "|" in tip else tip
        tier = tiers.get(protein, "outlier")
        seq = seqs.get(tip) or seqs.get(protein)
        if seq is None:
            continue
        if label < 0:
            continue  # noise handled separately

        if tier == "core":
            cluster_core.setdefault(label, {})[protein] = seq
        elif tier in ("affiliate", "bridge"):
            cluster_affiliate.setdefault(label, {})[protein] = seq

    result: dict = {}
    all_clusters = sorted(set(list(cluster_core.keys()) + list(cluster_affiliate.keys())))
    for cl_id in all_clusters:
        core_seqs = cluster_core.get(cl_id, {})
        aff_seqs = cluster_affiliate.get(cl_id, {})

        core_fp = os.path.join(outdir_fasta, f"cluster_{cl_id}.core.faa")
        aff_fp = os.path.join(outdir_fasta, f"cluster_{cl_id}.affiliate.faa")

        if core_seqs and (not os.path.exists(core_fp) or force):
            write_fasta(core_fp, core_seqs)
        if aff_seqs and (not os.path.exists(aff_fp) or force):
            write_fasta(aff_fp, aff_seqs)

        result[cl_id] = {
            "core": core_fp if core_seqs else None,
            "affiliate": aff_fp if aff_seqs else None,
            "n_core": len(core_seqs),
            "n_affiliate": len(aff_seqs),
        }

    return result


def _build_cluster_msa(fasta_path, output_aln, threads=4, force=False):
    """Run MAFFT on a cluster FASTA to produce a seed MSA.

    Returns ``True`` on success or if the output already exists, ``False``
    otherwise (e.g. MAFFT not found, too few sequences, or an error).
    """
    if os.path.exists(output_aln) and not force:
        return True
    if not os.path.exists(fasta_path):
        return False
    if not shutil.which("mafft"):
        import sys
        print("[embed/cluster] mafft not found in PATH — skipping MSA.", file=sys.stderr)
        return False

    n_seqs = sum(1 for line in open(fasta_path) if line.startswith(">"))
    if n_seqs < 2:
        import sys
        print(
            f"[embed/cluster] {os.path.basename(fasta_path)}: only {n_seqs} "
            "sequence(s) — skipping MSA.",
            file=sys.stderr,
        )
        return False

    try:
        cmd = f"mafft --auto --thread {threads} --quiet {fasta_path} > {output_aln}"
        subprocess.run(cmd, shell=True, check=True, stderr=subprocess.DEVNULL)
        return True
    except subprocess.CalledProcessError as e:
        import sys
        print(f"[embed/cluster] mafft failed for {fasta_path}: {e}", file=sys.stderr)
        return False
    except Exception as e:
        import sys
        print(f"[embed/cluster] MSA error for {fasta_path}: {e}", file=sys.stderr)
        return False


def _generate_sequence_logo(msa_path, out_png, out_svg, hmm_name, cluster_id,
                             formats=("png", "svg")):
    """Generate a sequence logo from an MSA FASTA file using matplotlib.

    Computes per-position information content (bits) and draws stacked
    coloured bars where each bar height reflects the residue's contribution
    to that position's information content.

    Parameters
    ----------
    msa_path : str
        Path to an aligned FASTA file.
    out_png / out_svg : str
        Output file paths for PNG and SVG respectively.
    hmm_name : str
        HMM name used in the plot title.
    cluster_id : int | str
        Cluster identifier used in the plot title.
    formats : tuple[str]
        Which formats to save (``"png"`` and/or ``"svg"``).

    Returns
    -------
    tuple : (success: bool, n_seqs: int, alignment_length: int, n_effective_sites: int)
        *n_effective_sites* counts positions with IC > 0.5 bits.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        import sys
        print("[embed/cluster] matplotlib not available — skipping logo.", file=sys.stderr)
        return False, 0, 0, 0

    try:
        msa = read_fasta(msa_path)
        if not msa:
            return False, 0, 0, 0
        seqs = list(msa.values())
        if len(seqs) < 2:
            return False, 0, 0, 0

        aln_len = max(len(s) for s in seqs)
        seqs_padded = [s.ljust(aln_len, "-") for s in seqs]
        n_seqs = len(seqs_padded)
        AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")

        # Compute per-position residue frequencies and information content
        positions = []
        for pos in range(aln_len):
            col = [s[pos].upper() for s in seqs_padded if pos < len(s)]
            non_gap = [c for c in col if c in AMINO_ACIDS]
            if not non_gap:
                positions.append({"ic": 0.0, "freqs": {}})
                continue
            from collections import Counter
            cnt = Counter(non_gap)
            total_ng = len(non_gap)
            freqs = {aa: cnt[aa] / total_ng for aa in AMINO_ACIDS if cnt.get(aa, 0) > 0}
            H = -sum(f * math.log2(f) for f in freqs.values() if f > 0)
            ic = max(0.0, math.log2(_N_AMINO_ACIDS) - H) * (total_ng / len(col))
            positions.append({"ic": ic, "freqs": freqs})

        n_effective_sites = sum(1 for p in positions if p["ic"] > 0.5)

        # Select positions to display (IC > threshold, skip nearly all-gap cols)
        show_idx = [i for i, p in enumerate(positions) if p["ic"] > 0.05]
        if not show_idx:
            return False, n_seqs, aln_len, n_effective_sites

        # Cap display length for readability
        if len(show_idx) > 80:
            show_idx = sorted(show_idx, key=lambda i: -positions[i]["ic"])[:80]
            show_idx = sorted(show_idx)

        fig_w = max(8, len(show_idx) * 0.22)
        fig, ax = plt.subplots(figsize=(fig_w, 3))

        for x_pos, col_idx in enumerate(show_idx):
            p = positions[col_idx]
            ic = p["ic"]
            if ic <= 0 or not p["freqs"]:
                continue
            sorted_aas = sorted(p["freqs"].items(), key=lambda kv: kv[1])
            y_bot = 0.0
            for aa, freq in sorted_aas:
                h = freq * ic
                ax.bar(x_pos, h, bottom=y_bot, width=0.85,
                       color=_AA_COLORS.get(aa, "#888888"),
                       edgecolor="none", linewidth=0)
                y_bot += h

        ax.set_xlim(-0.5, len(show_idx) - 0.5)
        ax.set_ylim(0, math.log2(_N_AMINO_ACIDS))
        ax.set_xlabel("Alignment position")
        ax.set_ylabel("Information content (bits)")
        ax.set_title(f"{hmm_name} — Cluster {cluster_id} logo  (n={n_seqs})")

        step = max(1, len(show_idx) // 20)
        tick_x = list(range(0, len(show_idx), step))
        ax.set_xticks(tick_x)
        ax.set_xticklabels([str(show_idx[i] + 1) for i in tick_x], fontsize=6, rotation=90)

        fig.tight_layout()
        saved = False
        if "png" in formats and out_png:
            os.makedirs(os.path.dirname(out_png), exist_ok=True)
            fig.savefig(out_png, dpi=150, format="png")
            saved = True
        if "svg" in formats and out_svg:
            os.makedirs(os.path.dirname(out_svg), exist_ok=True)
            fig.savefig(out_svg, format="svg")
            saved = True
        plt.close(fig)
        return saved, n_seqs, aln_len, n_effective_sites

    except Exception as e:
        import sys
        print(
            f"[embed/cluster] Logo generation failed for {hmm_name} "
            f"cluster {cluster_id}: {e}",
            file=sys.stderr,
        )
        return False, 0, 0, 0


def _build_cluster_hmm(msa_path, output_hmm, cluster_name, threads=4, force=False):
    """Run ``hmmbuild`` to create a profile HMM from a cluster seed MSA.

    Returns ``True`` on success or if the HMM already exists, ``False``
    otherwise (e.g. ``hmmbuild`` not in PATH or an error).
    """
    if os.path.exists(output_hmm) and not force:
        return True
    if not os.path.exists(msa_path):
        return False
    if not shutil.which("hmmbuild"):
        import sys
        print("[embed/cluster] hmmbuild not found — skipping HMM build.", file=sys.stderr)
        return False
    try:
        cmd = ["hmmbuild", "--cpu", str(threads), "-n", cluster_name,
               output_hmm, msa_path]
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
        return True
    except subprocess.CalledProcessError as e:
        import sys
        print(f"[embed/cluster] hmmbuild failed for {cluster_name}: {e}", file=sys.stderr)
        return False
    except Exception as e:
        import sys
        print(f"[embed/cluster] hmmbuild error for {cluster_name}: {e}", file=sys.stderr)
        return False


def _classify_noise_sequences(knn_df, hmm_name):
    """Classify noise (label −1) sequences using kNN neighbourhood evidence.

    Noise sequences are categorised as one of:

    * **peripheral_homolog** — strong kNN affinity to a single cluster
    * **partial_homolog** — moderate kNN affinity with low entropy
    * **bridge_sequence** — mixed neighbourhood suggesting inter-cluster
      bridging or domain fusion / truncation
    * **fusion_or_extension** — very high entropy neighbourhood implying
      an unusual architecture
    * **outlier** — no clear cluster affiliation

    Parameters
    ----------
    knn_df : pd.DataFrame
        Output of :func:`_compute_knn_metrics`.
    hmm_name : str
        HMM name, stored in the result table.

    Returns
    -------
    pd.DataFrame  with columns: hmm_name, protein_id, classification,
    dominant_cluster, neighborhood_purity, neighborhood_entropy, notes
    """
    if knn_df.empty:
        return pd.DataFrame()

    noise_df = knn_df[knn_df["cluster_id"] == -1].copy()
    if noise_df.empty:
        return pd.DataFrame()

    rows = []
    for _, row in noise_df.iterrows():
        pid = row["protein_id"]
        dom = int(row["dominant_cluster"])
        purity = float(row["neighborhood_purity"])
        entropy = float(row["neighborhood_entropy"])

        if dom >= 0 and purity >= 0.6:
            classification = "peripheral_homolog"
            notes = f"Strong kNN affinity to cluster {dom} (purity={purity:.2f})"
        elif dom >= 0 and purity >= 0.3:
            if entropy > 1.5:
                classification = "bridge_sequence"
                notes = (
                    f"Mixed neighbourhood, partial affinity to cluster {dom} "
                    f"(purity={purity:.2f}, entropy={entropy:.2f})"
                )
            else:
                classification = "partial_homolog"
                notes = f"Moderate kNN affinity to cluster {dom} (purity={purity:.2f})"
        elif entropy > 2.0:
            classification = "fusion_or_extension"
            notes = (
                f"High-entropy neighbourhood (entropy={entropy:.2f}), "
                "possible domain fusion or extension"
            )
        else:
            classification = "outlier"
            notes = "No clear cluster affiliation"

        rows.append({
            "hmm_name": hmm_name,
            "protein_id": pid,
            "classification": classification,
            "dominant_cluster": dom,
            "neighborhood_purity": purity,
            "neighborhood_entropy": round(entropy, 4),
            "notes": notes,
        })

    return pd.DataFrame(rows)


def run_cluster_subworkflow(hmm_name, seqs, ids, Z, labels, emb_cfg,
                            outdir, summary_dir, threads=4, force=False):
    """Orchestrate the cluster-aware subworkflow for a single HMM hit set.

    For each embedding-defined cluster this function:

    1. Computes per-sequence kNN neighbourhood metrics (local graph support).
    2. Assigns membership tiers: core / affiliate / bridge / outlier.
    3. Writes per-cluster FASTA files (core and affiliate seeds).
    4. Builds a seed MSA for each cluster using MAFFT (core sequences only).
    5. Generates a sequence logo (PNG + SVG) for each cluster MSA.
    6. Optionally builds a profile HMM for each cluster MSA via ``hmmbuild``.
    7. Classifies noise (HDBSCAN −1) sequences using kNN evidence.
    8. Writes summary TSV tables and a logo manifest.

    Parameters
    ----------
    hmm_name : str
        Identifier for this HMM / hit-set.
    seqs : dict
        {tip_id: amino_acid_sequence}.
    ids : list
        Sequence IDs in the same order as *labels* and rows of *Z*.
    Z : np.ndarray
        PCA-reduced feature matrix (n_samples × n_pca_components).
    labels : list
        Integer cluster labels (−1 = noise).
    emb_cfg : dict
        The validated embeddings config dict.  The relevant sub-section is
        ``emb_cfg["cluster_subworkflow"]``.
    outdir : str
        Root output directory (parent of ``cluster_fasta/``, etc.).
    summary_dir : str
        Directory for summary TSV files.
    threads : int
        CPU threads for MAFFT / ``hmmbuild``.
    force : bool
        Overwrite existing output files.
    """
    import sys

    sub_cfg = emb_cfg.get("cluster_subworkflow", {})
    if not sub_cfg.get("enabled", False):
        return

    build_msas = sub_cfg.get("build_cluster_msas", True)
    seed_membership = sub_cfg.get("seed_membership", "core_only")
    build_hmms = sub_cfg.get("build_cluster_hmms", True)
    classify_noise = sub_cfg.get("classify_noise", True)
    gen_logos = sub_cfg.get("generate_sequence_logos", True)
    logo_formats = sub_cfg.get("logo_format", ["png", "svg"])

    unique_clusters = sorted(set(labels) - {-1})
    noise_count = sum(1 for lb in labels if lb == -1)
    print(
        f"[embed/cluster] Running cluster subworkflow for {hmm_name}: "
        f"{len(unique_clusters)} clusters, {noise_count} noise points"
    )

    if not unique_clusters:
        print(
            f"[embed/cluster] No clusters found for {hmm_name} — "
            "skipping cluster subworkflow.",
            file=sys.stderr,
        )
        return

    # ── 1. kNN metrics ────────────────────────────────────────────────────────
    k = int(emb_cfg.get("knn_neighbors", 10))
    knn_df = _compute_knn_metrics(Z, ids, labels, k=k)

    if not knn_df.empty and summary_dir:
        safe_mkdir(summary_dir)
        knn_metrics_fp = os.path.join(summary_dir, f"{hmm_name}.cluster_membership.tsv")
        knn_df.to_csv(knn_metrics_fp, sep="\t", index=False)
        print(f"[embed/cluster] Wrote kNN metrics: {knn_metrics_fp}")

    # ── 2. Membership tiers ───────────────────────────────────────────────────
    tiers = _assign_membership_tiers(ids, labels, knn_df)

    # ── 3. Per-cluster FASTAs ─────────────────────────────────────────────────
    fasta_dir = os.path.join(outdir, "cluster_fasta", hmm_name)
    cluster_info = _write_cluster_fastas(seqs, ids, labels, tiers, fasta_dir, force=force)
    print(
        f"[embed/cluster] Wrote cluster FASTAs for {hmm_name} "
        f"({len(cluster_info)} clusters)"
    )

    # ── 4–6. Per-cluster MSA, logo, and HMM ──────────────────────────────────
    aln_dir = os.path.join(outdir, "cluster_alignments", hmm_name)
    logo_dir = os.path.join(outdir, "cluster_logos", hmm_name)
    hmm_dir = os.path.join(outdir, "cluster_hmms", hmm_name)
    safe_mkdir(aln_dir)
    safe_mkdir(logo_dir)
    if build_hmms:
        safe_mkdir(hmm_dir)

    logo_manifest_rows = []

    for cl_id, info in cluster_info.items():
        # Choose the seed FASTA according to seed_membership config
        if seed_membership == "core_only":
            seed_fasta = info["core"]
            n_seed = info["n_core"]
        else:  # "core_and_affiliate"
            # Merge core + affiliate into a combined seed FASTA
            core_fp = info["core"]
            aff_fp = info["affiliate"]
            if core_fp and aff_fp:
                combined_fp = os.path.join(fasta_dir, f"cluster_{cl_id}.seed.faa")
                if not os.path.exists(combined_fp) or force:
                    core_s = read_fasta(core_fp) if core_fp and os.path.exists(core_fp) else {}
                    aff_s = read_fasta(aff_fp) if aff_fp and os.path.exists(aff_fp) else {}
                    write_fasta(combined_fp, {**core_s, **aff_s})
                seed_fasta = combined_fp
                n_seed = info["n_core"] + info["n_affiliate"]
            elif core_fp:
                seed_fasta = core_fp
                n_seed = info["n_core"]
            else:
                seed_fasta = aff_fp
                n_seed = info["n_affiliate"]

        if not seed_fasta or not os.path.exists(seed_fasta) or n_seed < 2:
            print(
                f"[embed/cluster] Cluster {cl_id} of {hmm_name}: "
                f"insufficient seed sequences ({n_seed}) — skipping.",
                file=sys.stderr,
            )
            continue

        aln_fp = os.path.join(aln_dir, f"cluster_{cl_id}.seed.aln.faa")
        logo_png = os.path.join(logo_dir, f"cluster_{cl_id}.logo.png")
        logo_svg = os.path.join(logo_dir, f"cluster_{cl_id}.logo.svg")
        hmm_fp = os.path.join(hmm_dir, f"cluster_{cl_id}.hmm")

        # Build MSA
        msa_ok = False
        if build_msas:
            msa_ok = _build_cluster_msa(seed_fasta, aln_fp, threads=threads, force=force)

        # Generate logo
        logo_ok = False
        n_logo_seqs = 0
        n_eff_sites = 0
        aln_len_logo = 0
        if gen_logos and msa_ok:
            logo_ok, n_logo_seqs, aln_len_logo, n_eff_sites = _generate_sequence_logo(
                aln_fp, logo_png, logo_svg,
                hmm_name, cl_id,
                formats=logo_formats,
            )
            if logo_ok:
                print(
                    f"[embed/cluster] Logo saved for {hmm_name}/cluster_{cl_id} "
                    f"(n={n_logo_seqs}, aln_len={aln_len_logo}, effective_sites={n_eff_sites})"
                )

        # Build cluster HMM
        hmm_ok = False
        if build_hmms and msa_ok:
            cluster_name = f"{hmm_name}_cluster_{cl_id}"
            hmm_ok = _build_cluster_hmm(aln_fp, hmm_fp, cluster_name,
                                         threads=threads, force=force)
            if hmm_ok:
                print(f"[embed/cluster] HMM built: {hmm_fp}")

        logo_manifest_rows.append({
            "hmm_name": hmm_name,
            "cluster_id": cl_id,
            "membership_tier_used": seed_membership,
            "n_sequences": n_seed,
            "alignment_length": aln_len_logo,
            "n_effective_sites": n_eff_sites,
            "logo_png": logo_png if logo_ok and "png" in logo_formats else "",
            "logo_svg": logo_svg if logo_ok and "svg" in logo_formats else "",
            "hmm": hmm_fp if hmm_ok else "",
            "notes": "",
        })

    # ── Write logo manifest ───────────────────────────────────────────────────
    if logo_manifest_rows and summary_dir:
        safe_mkdir(summary_dir)
        manifest_fp = os.path.join(summary_dir, f"{hmm_name}.cluster_logo_manifest.tsv")
        pd.DataFrame(logo_manifest_rows).to_csv(manifest_fp, sep="\t", index=False)
        print(f"[embed/cluster] Wrote logo manifest: {manifest_fp}")

    # ── 7. Noise classification ───────────────────────────────────────────────
    if classify_noise and noise_count > 0:
        noise_cls_df = _classify_noise_sequences(knn_df, hmm_name)
        if not noise_cls_df.empty and summary_dir:
            noise_fp = os.path.join(summary_dir, f"{hmm_name}.noise_classification.tsv")
            noise_cls_df.to_csv(noise_fp, sep="\t", index=False)
            print(f"[embed/cluster] Wrote noise classification: {noise_fp}")

        # Optionally score noise against cluster HMMs if they were built
        if build_hmms and shutil.which("hmmscan") and summary_dir:
            _score_noise_against_hmms(
                ids, labels, seqs, hmm_dir, hmm_name,
                summary_dir, threads, force,
            )

    # ── 8. KL divergence between cluster MSAs ────────────────────────────────
    kl_cfg = sub_cfg.get("kl_divergence", {})
    if kl_cfg.get("enabled", False) and build_msas and len(unique_clusters) >= 2:
        min_cl_size = int(kl_cfg.get("min_cluster_size", 5))
        pseudocount = float(kl_cfg.get("pseudocount", 1e-6))
        top_n_sites = int(kl_cfg.get("top_n_sites", 20))

        # Collect alignment paths for clusters that have a valid MSA
        cluster_aln_paths = {}
        for cl_id in unique_clusters:
            aln_fp = os.path.join(aln_dir, f"cluster_{cl_id}.seed.aln.faa")
            if os.path.exists(aln_fp):
                cluster_aln_paths[cl_id] = aln_fp

        if len(cluster_aln_paths) >= 2:
            kl_df, top_df = _compute_cluster_kl_divergence(
                hmm_name,
                cluster_aln_paths,
                min_cluster_size=min_cl_size,
                pseudocount=pseudocount,
                top_n_sites=top_n_sites,
            )
            if not kl_df.empty and summary_dir:
                safe_mkdir(summary_dir)
                kl_fp = os.path.join(
                    summary_dir, f"{hmm_name}.cluster_kl_divergence.tsv"
                )
                kl_df.to_csv(kl_fp, sep="\t", index=False)
                print(f"[embed/cluster] Wrote KL divergence table: {kl_fp}")
            if not top_df.empty and summary_dir:
                safe_mkdir(summary_dir)
                top_fp = os.path.join(
                    summary_dir, f"{hmm_name}.cluster_kl_top_sites.tsv"
                )
                top_df.to_csv(top_fp, sep="\t", index=False)
                print(f"[embed/cluster] Wrote top divergent sites: {top_fp}")

    print(f"[embed/cluster] Cluster subworkflow complete for {hmm_name}.")


def _compute_cluster_kl_divergence(
    hmm_name,
    cluster_aln_paths,
    min_cluster_size=5,
    pseudocount=1e-6,
    top_n_sites=20,
):
    """Compute per-position KL divergence between all pairs of cluster MSAs.

    For each pair of clusters whose MSAs share a common alignment coordinate
    system (i.e. were produced by the same MAFFT run on the whole hit set, or
    are otherwise column-aligned), this function:

    1. Reads each cluster MSA and tallies amino-acid counts at every column.
    2. Adds a pseudocount to avoid zero probabilities.
    3. Computes the asymmetric KL divergence KL(P || Q) **and** the
       symmetric Jensen-Shannon divergence (JSD) for every alignment position.
    4. Returns two DataFrames:

       * **per_position** – one row per (cluster_pair, alignment_position).
       * **top_sites** – the ``top_n_sites`` highest-JSD positions per pair.

    Parameters
    ----------
    hmm_name : str
        HMM / hit-set identifier used for output labelling.
    cluster_aln_paths : dict
        Mapping of ``{cluster_id: aligned_fasta_path}``.  All MSAs must have
        the same number of alignment columns.
    min_cluster_size : int
        Clusters with fewer aligned sequences than this value are skipped.
    pseudocount : float
        Laplace-like pseudocount added to every residue count before
        normalisation.  Prevents log(0) errors.
    top_n_sites : int
        Number of top-divergent sites to include in the summary table.

    Returns
    -------
    per_position_df : pd.DataFrame
        Columns: ``hmm_name``, ``cluster_A``, ``cluster_B``, ``pair``,
        ``aln_position``, ``kl_A_to_B``, ``kl_B_to_A``, ``js_divergence``,
        ``top_aa_A``, ``top_aa_B``, ``n_seqs_A``, ``n_seqs_B``.
    top_sites_df : pd.DataFrame
        Same schema as *per_position_df*, restricted to the highest-JSD
        positions per cluster pair.
    """
    import sys
    import math
    from itertools import combinations
    from ..utils.bio import read_fasta

    AA_ALPHABET = list("ACDEFGHIKLMNPQRSTVWY")
    GAP_CHARS = set("-. X?*")

    def _site_counts(aln_seqs):
        """Return a list of Counter objects, one per alignment column."""
        from collections import Counter
        if not aln_seqs:
            return []
        seqs_list = list(aln_seqs.values())
        n_cols = len(seqs_list[0])
        counts = []
        for col_idx in range(n_cols):
            c = Counter()
            for seq in seqs_list:
                if col_idx < len(seq):
                    aa = seq[col_idx].upper()
                    if aa not in GAP_CHARS and aa in set(AA_ALPHABET):
                        c[aa] += 1
            counts.append(c)
        return counts

    def _kl(p_counts, q_counts, pc=pseudocount):
        """Asymmetric KL divergence KL(P || Q) in bits."""
        p_total = sum(p_counts.get(a, 0) for a in AA_ALPHABET) + pc * len(AA_ALPHABET)
        q_total = sum(q_counts.get(a, 0) for a in AA_ALPHABET) + pc * len(AA_ALPHABET)
        kl = 0.0
        for a in AA_ALPHABET:
            p = (p_counts.get(a, 0) + pc) / p_total
            q = (q_counts.get(a, 0) + pc) / q_total
            kl += p * math.log2(p / q)
        return float(kl)

    def _jsd(p_counts, q_counts, pc=pseudocount):
        """Symmetric Jensen-Shannon divergence in bits (range [0, 1])."""
        p_total = sum(p_counts.get(a, 0) for a in AA_ALPHABET) + pc * len(AA_ALPHABET)
        q_total = sum(q_counts.get(a, 0) for a in AA_ALPHABET) + pc * len(AA_ALPHABET)
        p_vec = [(p_counts.get(a, 0) + pc) / p_total for a in AA_ALPHABET]
        q_vec = [(q_counts.get(a, 0) + pc) / q_total for a in AA_ALPHABET]
        m_vec = [(p + q) / 2.0 for p, q in zip(p_vec, q_vec)]
        jsd = 0.0
        for p, q, m in zip(p_vec, q_vec, m_vec):
            if m > 0:
                if p > 0:
                    jsd += 0.5 * p * math.log2(p / m)
                if q > 0:
                    jsd += 0.5 * q * math.log2(q / m)
        return float(jsd)

    def _top_aa(counts):
        """Return the most common amino acid at a position, or '-' if all gaps."""
        if not counts:
            return "-"
        return max(counts, key=counts.get)

    # ── Load cluster MSAs and compute per-column counts ───────────────────────
    cluster_counts = {}   # cl_id -> list[Counter]
    cluster_n_seqs = {}   # cl_id -> int

    for cl_id, aln_fp in cluster_aln_paths.items():
        try:
            aln_seqs = read_fasta(aln_fp)
        except Exception as exc:
            print(
                f"[embed/cluster/kl] Could not read {aln_fp}: {exc}",
                file=sys.stderr,
            )
            continue
        n_seqs = len(aln_seqs)
        if n_seqs < min_cluster_size:
            print(
                f"[embed/cluster/kl] Cluster {cl_id} ({n_seqs} seqs) is below "
                f"min_cluster_size={min_cluster_size} — skipping.",
                file=sys.stderr,
            )
            continue
        cluster_counts[cl_id] = _site_counts(aln_seqs)
        cluster_n_seqs[cl_id] = n_seqs

    eligible = sorted(cluster_counts.keys())
    if len(eligible) < 2:
        return pd.DataFrame(), pd.DataFrame()

    # ── Compute pairwise KL / JSD ─────────────────────────────────────────────
    per_position_rows = []

    for cl_a, cl_b in combinations(eligible, 2):
        counts_a = cluster_counts[cl_a]
        counts_b = cluster_counts[cl_b]

        # Alignment lengths must match; take the shorter to be safe
        L = min(len(counts_a), len(counts_b))
        if L == 0:
            continue

        n_a = cluster_n_seqs[cl_a]
        n_b = cluster_n_seqs[cl_b]
        pair_label = f"cluster_{cl_a}:cluster_{cl_b}"

        for pos in range(L):
            ca = counts_a[pos]
            cb = counts_b[pos]
            kl_ab = _kl(ca, cb)
            kl_ba = _kl(cb, ca)
            jsd_val = _jsd(ca, cb)
            per_position_rows.append(
                {
                    "hmm_name": hmm_name,
                    "cluster_A": cl_a,
                    "cluster_B": cl_b,
                    "pair": pair_label,
                    "aln_position": pos + 1,
                    "kl_A_to_B": round(kl_ab, 6),
                    "kl_B_to_A": round(kl_ba, 6),
                    "js_divergence": round(jsd_val, 6),
                    "top_aa_A": _top_aa(ca),
                    "top_aa_B": _top_aa(cb),
                    "n_seqs_A": n_a,
                    "n_seqs_B": n_b,
                }
            )

    if not per_position_rows:
        return pd.DataFrame(), pd.DataFrame()

    per_position_df = pd.DataFrame(per_position_rows)

    # ── Build top-sites summary ───────────────────────────────────────────────
    top_rows = []
    for pair_label, grp in per_position_df.groupby("pair"):
        top = grp.nlargest(top_n_sites, "js_divergence")
        top_rows.append(top)

    top_sites_df = pd.concat(top_rows, ignore_index=True) if top_rows else pd.DataFrame()

    return per_position_df, top_sites_df


def _score_noise_against_hmms(ids, labels, seqs, hmm_dir, hmm_name,
                               summary_dir, threads=4, force=False):
    """Score noise sequences against per-cluster HMMs using ``hmmscan``.

    Writes ``summary/<HMM>.cluster_recruitment.tsv`` with best-HMM hits for
    each noise sequence.  Silently skips if no HMM files are present or
    ``hmmscan`` is not available.
    """
    import tempfile
    import sys

    recruitment_fp = os.path.join(summary_dir, f"{hmm_name}.cluster_recruitment.tsv")
    if os.path.exists(recruitment_fp) and not force:
        return

    # Collect noise protein sequences
    noise_seqs: dict = {}
    for tip, label in zip(ids, labels):
        if label == -1:
            protein = tip.split("|", 1)[1] if "|" in tip else tip
            seq = seqs.get(tip) or seqs.get(protein)
            if seq:
                noise_seqs[protein] = seq

    if not noise_seqs:
        return

    # Collect built HMM files
    hmm_files = sorted(glob.glob(os.path.join(hmm_dir, "cluster_*.hmm")))
    if not hmm_files:
        return

    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            # Write noise sequences to a temp FASTA
            noise_fasta = os.path.join(tmpdir, "noise.faa")
            write_fasta(noise_fasta, noise_seqs)

            # Concatenate all cluster HMMs into a library
            combined_hmm = os.path.join(tmpdir, "cluster_lib.hmm")
            with open(combined_hmm, "wb") as fout:
                for hf in hmm_files:
                    with open(hf, "rb") as fin:
                        fout.write(fin.read())

            # hmmpress the library
            subprocess.run(
                ["hmmpress", "-f", combined_hmm],
                check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
            )

            # Run hmmscan
            tbl_out = os.path.join(tmpdir, "noise_hmmscan.tbl")
            subprocess.run(
                ["hmmscan", "--cpu", str(threads),
                 "--noali", "--tblout", tbl_out, combined_hmm, noise_fasta],
                check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
            )

            # Parse tblout (space-delimited, # = comment)
            rows = []
            with open(tbl_out) as fh:
                for line in fh:
                    if line.startswith("#"):
                        continue
                    parts = line.split()
                    if len(parts) < 9:
                        continue
                    rows.append({
                        "hmm_name": hmm_name,
                        "protein_id": parts[2],
                        "best_cluster_hmm": parts[0],
                        "evalue": parts[4],
                        "bitscore": parts[5],
                    })

            if rows:
                safe_mkdir(summary_dir)
                pd.DataFrame(rows).to_csv(recruitment_fp, sep="\t", index=False)
                print(f"[embed/cluster] Wrote noise recruitment: {recruitment_fp}")

    except Exception as e:
        print(
            f"[embed/cluster] Noise HMM scoring failed for {hmm_name}: {e}",
            file=sys.stderr,
        )


def _validate_emb_cfg(emb_cfg: dict) -> dict:
    """Validate and normalise the embeddings config section.

    Returns a copy of *emb_cfg* with defaults filled in and values normalised.
    Raises ``ValueError`` for invalid option combinations.
    """
    import sys
    cfg = dict(emb_cfg)

    # cluster_on
    cluster_on = str(cfg.get("cluster_on", "PCA")).lower()
    if cluster_on not in _VALID_CLUSTER_ON:
        raise ValueError(
            f"[embed] Invalid 'cluster_on' value '{cfg.get('cluster_on')}'. "
            f"Must be one of: {sorted(_VALID_CLUSTER_ON)}"
        )
    cfg["_cluster_on"] = cluster_on  # normalised lowercase key used internally

    # cluster_method
    cluster_method = str(cfg.get("cluster_method", "hdbscan")).lower()
    if cluster_method not in _VALID_CLUSTER_METHODS:
        raise ValueError(
            f"[embed] Invalid 'cluster_method' value '{cfg.get('cluster_method')}'. "
            f"Must be one of: {sorted(_VALID_CLUSTER_METHODS)}"
        )
    cfg["_cluster_method"] = cluster_method

    # umap_dimensions
    umap_dims = int(cfg.get("umap_dimensions", 2))
    if umap_dims not in (2, 3):
        print(
            f"[embed] Warning: 'umap_dimensions' must be 2 or 3; got {umap_dims}. Defaulting to 2.",
            file=sys.stderr,
        )
        umap_dims = 2
    cfg["_umap_dimensions"] = umap_dims

    # hdbscan_metric – any string is accepted; invalid values produce a runtime
    # error from scikit-learn, which is caught gracefully in _run_hdbscan.
    cfg.setdefault("hdbscan_metric", "euclidean")

    # kNN parameters
    cfg.setdefault("run_knn", True)
    knn_neighbors = int(cfg.get("knn_neighbors", 10))
    if knn_neighbors < 1:
        print(
            f"[embed] Warning: 'knn_neighbors' must be >= 1; got {knn_neighbors}. Defaulting to 10.",
            file=sys.stderr,
        )
        knn_neighbors = 10
    cfg["knn_neighbors"] = knn_neighbors

    # pca_components default
    cfg.setdefault("pca_components", 50)

    # cluster_subworkflow defaults
    sub_cfg = cfg.get("cluster_subworkflow", {})
    if not isinstance(sub_cfg, dict):
        sub_cfg = {}
    sub_cfg.setdefault("enabled", False)
    sub_cfg.setdefault("build_cluster_msas", True)
    sub_cfg.setdefault("seed_membership", "core_only")
    sub_cfg.setdefault("build_cluster_hmms", True)
    sub_cfg.setdefault("classify_noise", True)
    sub_cfg.setdefault("recover_affiliates", True)
    sub_cfg.setdefault("generate_sequence_logos", True)
    sub_cfg.setdefault("logo_format", ["png", "svg"])
    sub_cfg.setdefault("compare_cluster_hmms", False)
    cfg["cluster_subworkflow"] = sub_cfg

    return cfg


def compute_embeddings_for_hmm(hmm_name: str, seqs: dict, emb_cfg: dict, outdir_embeddings: str,
                                force: bool, clades: dict | None, tax_map: dict | None = None,
                                outdir: str | None = None, summary_dir: str | None = None,
                                threads: int = 4):
    """Compute embeddings, PCA, kNN, UMAP (visualization only), clustering, and save plots + TSVs.

    Clustering behaviour is controlled by the following config options:

    * ``cluster_on`` (``"PCA"`` | ``"embeddings"``): whether to cluster the
      PCA-reduced vectors or the raw embedding vectors.
    * ``cluster_method`` (``"hdbscan"`` | ``"leiden"``): clustering algorithm.
    * ``hdbscan_metric`` (default ``"euclidean"``): distance metric for HDBSCAN.
    * ``umap_dimensions`` (``2`` | ``3``): dimensionality of the UMAP projection
      used **only** for visualisation, never for clustering.
    * ``run_knn`` (default ``True``): compute a cosine kNN graph on PCA space.
    * ``knn_neighbors`` (default ``10``): number of neighbors for kNN.
    * ``cluster_subworkflow.enabled`` (default ``False``): when ``True``, run the
      optional cluster-aware subworkflow (per-cluster MSA, logos, HMMs).

    Returns a list of cluster assignment dicts (for clade_assignment.tsv), or
    empty list.
    """
    import sys
    safe_mkdir(outdir_embeddings)

    out_npy = os.path.join(outdir_embeddings, f"{hmm_name}.embeddings.npy")
    out_pca = os.path.join(outdir_embeddings, f"{hmm_name}.pca.tsv")
    out_knn = os.path.join(outdir_embeddings, f"{hmm_name}.knn.tsv")
    out_umap = os.path.join(outdir_embeddings, f"{hmm_name}.umap.tsv")
    out_umap_png = os.path.join(outdir_embeddings, f"{hmm_name}.umap.png")
    out_umap_clust_png = os.path.join(outdir_embeddings, f"{hmm_name}.umap.clustered.png")
    out_meta = os.path.join(outdir_embeddings, f"{hmm_name}.pca.meta.json")
    out_disp = os.path.join(outdir_embeddings, f"{hmm_name}.dispersion.tsv")
    out_vec_tsv = os.path.join(outdir_embeddings, f"{hmm_name}.vectors.tsv")

    if os.path.exists(out_pca) and os.path.exists(out_npy) and os.path.exists(out_umap) and not force:
        return []

    seqs = {k: v.replace(" ", "").replace("\n", "").replace("*", "").replace(".", "") for k, v in seqs.items()}
    if len(seqs) < 3:
        print(f"[embed] Warning: HMM '{hmm_name}' has less than 3 sequences. Skipping embeddings.", file=sys.stderr)
        return []

    # ── Validate config ───────────────────────────────────────────────────
    try:
        emb_cfg = _validate_emb_cfg(emb_cfg)
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        return []

    cluster_on = emb_cfg["_cluster_on"]        # "pca" or "embeddings"
    cluster_method = emb_cfg["_cluster_method"]  # "hdbscan" or "leiden"
    umap_dims = emb_cfg["_umap_dimensions"]      # 2 or 3

    backend = emb_cfg["backend"]
    device = emb_cfg["device"]
    batch_size = int(emb_cfg["batch_size"])
    model_name = emb_cfg["model"]
    repr_layer = emb_cfg.get("repr_layer", None)
    model_dir = emb_cfg.get("model_dir", None)

    try:
        if backend == "esm":
            ids, X = _embed_esm(seqs, model_name=model_name, device=device, batch_size=batch_size, repr_layer=repr_layer, model_dir=model_dir)
        elif backend == "transformers":
            ids, X = _embed_transformers(seqs, model_id_or_path=model_name, device=device, batch_size=batch_size, model_dir=model_dir)
        else:
            print(f"[embed] Error: Unknown backend '{backend}'", file=sys.stderr)
            return []
    except Exception as e:
        print(f"[embed] FAILED {hmm_name}: {e}", file=sys.stderr)
        return []

    X = X.astype(np.float32)
    np.save(out_npy, X)

    # ── PCA ───────────────────────────────────────────────────────────────
    ncomp = int(emb_cfg["pca_components"])
    if ncomp < 2:
        ncomp = 2
    ncomp = min(ncomp, X.shape[1], max(2, X.shape[0] - 1))

    Z, var = _pca_fit_transform(X, n_components=ncomp)

    rows = []
    for tip, coords in zip(ids, Z):
        genome = tip.split("|", 1)[0] if "|" in tip else "Unknown"
        protein = tip.split("|", 1)[1] if "|" in tip else tip
        r = {"hmm": hmm_name, "tip": tip, "genome": genome, "protein": protein}
        for j in range(Z.shape[1]):
            r[f"PC{j+1}"] = float(coords[j])
        rows.append(r)
    pd.DataFrame(rows).to_csv(out_pca, sep="\t", index=False)

    # ── kNN on PCA space ──────────────────────────────────────────────────
    if emb_cfg.get("run_knn", True):
        try:
            from sklearn.neighbors import NearestNeighbors
            knn_neighbors = int(emb_cfg.get("knn_neighbors", 10))
            n_knn = min(knn_neighbors, len(ids) - 1)
            if n_knn < 1:
                import sys as _sys
                print(f"[embed] kNN skipped for {hmm_name}: not enough sequences.", file=_sys.stderr)
            else:
                nbrs = NearestNeighbors(n_neighbors=n_knn, metric="cosine").fit(Z)
                distances, indices = nbrs.kneighbors(Z)
                knn_rows = []
                for tip, dists, idxs in zip(ids, distances, indices):
                    protein = tip.split("|", 1)[1] if "|" in tip else tip
                    for rank, (dist, j) in enumerate(zip(dists, idxs), start=1):
                        neighbor_tip = ids[j]
                        neighbor_protein = neighbor_tip.split("|", 1)[1] if "|" in neighbor_tip else neighbor_tip
                        knn_rows.append({
                            "protein_id": protein,
                            "neighbor_rank": rank,
                            "neighbor_id": neighbor_protein,
                            "distance": float(dist),
                        })
                pd.DataFrame(knn_rows).to_csv(out_knn, sep="\t", index=False)
        except Exception as e:
            import sys as _sys
            print(f"[embed] kNN failed for {hmm_name}: {e}", file=_sys.stderr)

    # ── UMAP (visualization only) ─────────────────────────────────────────
    U = None
    try:
        import umap
        import warnings as _w
        reducer = umap.UMAP(n_components=umap_dims, random_state=42)
        with _w.catch_warnings():
            _w.simplefilter("ignore")
            U = reducer.fit_transform(Z)

        u_rows = []
        for tip, coords in zip(ids, U):
            genome = tip.split("|", 1)[0] if "|" in tip else "Unknown"
            protein = tip.split("|", 1)[1] if "|" in tip else tip
            row = {
                "hmm": hmm_name,
                "tip": tip,
                "genome": genome,
                "protein": protein,
                "UMAP1": float(coords[0]),
                "UMAP2": float(coords[1]),
            }
            if umap_dims >= 3:
                row["UMAP3"] = float(coords[2])
            u_rows.append(row)
        pd.DataFrame(u_rows).to_csv(out_umap, sep="\t", index=False)

        # Save UMAP plot (colored by clades if available)
        _save_umap_plot(U, ids, hmm_name, out_umap_png, clades=clades)
    except ImportError:
        print("[embed] UMAP skip: umap-learn not installed.", file=sys.stderr)
    except Exception as e:
        print(f"[embed] UMAP failed for {hmm_name}: {e}", file=sys.stderr)

    # ── Clustering ────────────────────────────────────────────────────────
    # Select the feature matrix to cluster on (PCA or raw embeddings)
    cluster_assignments = []
    if emb_cfg.get("cluster_embeddings", True):
        X_clust = Z if cluster_on == "pca" else X
        print(f"[embed] Clustering on {'PCA' if cluster_on == 'pca' else 'raw embeddings'} "
              f"using {cluster_method.upper()} for {hmm_name}")

        if cluster_method == "leiden":
            leiden_metric = emb_cfg.get("hdbscan_metric", "cosine")  # reused for kNN graph metric
            labels = _run_leiden(X_clust, metric=leiden_metric)
        else:
            # default: hdbscan
            min_cs = int(emb_cfg.get("hdbscan_min_cluster_size", 5))
            hdb_metric = emb_cfg.get("hdbscan_metric", "euclidean")
            labels = _run_hdbscan(X_clust, min_cluster_size=min_cs, metric=hdb_metric)

        if labels is not None:
            from ..utils.helpers import normalize_genome_id
            unique_clusters = set(labels) - {-1}
            n_clusters = len(unique_clusters)
            noise_count = sum(1 for lb in labels if lb == -1)
            print(f"[embed] {cluster_method.upper()} found {n_clusters} clusters for {hmm_name} "
                  f"({noise_count} noise points)")

            for tip, label in zip(ids, labels):
                genome = tip.split("|", 1)[0] if "|" in tip else "Unknown"
                protein = tip.split("|", 1)[1] if "|" in tip else tip
                taxonomy = "Unknown"
                if tax_map:
                    norm_g = normalize_genome_id(genome)
                    taxonomy = tax_map.get(norm_g, "Unknown")
                cluster_assignments.append({
                    "hmm": hmm_name,
                    "protein": protein,
                    "genome": genome,
                    "cluster_id": int(label),
                    "taxonomy": taxonomy,
                })

            # Save cluster-colored UMAP plot (UMAP used only for visualization)
            if U is not None:
                _save_umap_plot(U, ids, hmm_name, out_umap_clust_png,
                                cluster_labels=labels,
                                title_suffix=f" ({cluster_method.upper()} clusters)")

            # ── Optional cluster subworkflow ───────────────────────────────
            if emb_cfg.get("cluster_subworkflow", {}).get("enabled", False):
                _subworkflow_outdir = outdir if outdir else os.path.dirname(outdir_embeddings)
                _subworkflow_sumdir = summary_dir
                run_cluster_subworkflow(
                    hmm_name=hmm_name,
                    seqs=seqs,
                    ids=ids,
                    Z=Z,
                    labels=labels,
                    emb_cfg=emb_cfg,
                    outdir=_subworkflow_outdir,
                    summary_dir=_subworkflow_sumdir,
                    threads=threads,
                    force=force,
                )

    # ── metadata ──────────────────────────────────────────────────────────
    meta = {
        "hmm": hmm_name,
        "backend": backend,
        "model": model_name,
        "device": device,
        "batch_size": batch_size,
        "repr_layer": repr_layer,
        "pooling": emb_cfg["pooling"],
        "pca_components": int(Z.shape[1]),
        "explained_variance_ratio": var,
        "n_sequences": int(len(ids)),
        "vector_dim": int(X.shape[1]),
        "cluster_on": cluster_on,
        "cluster_method": cluster_method,
        "umap_dimensions": umap_dims,
    }
    write_json(meta, out_meta)

    # optional full vectors TSV
    if emb_cfg.get("write_full_vectors", False):
        vec_rows = []
        for tip, v in zip(ids, X):
            r = {"tip": tip}
            for j, val in enumerate(v):
                r[f"d{j+1}"] = float(val)
            vec_rows.append(r)
        pd.DataFrame(vec_rows).to_csv(out_vec_tsv, sep="\t", index=False)

    # clade dispersion
    if clades:
        Z_df = pd.DataFrame(Z, columns=[f"PC{j+1}" for j in range(Z.shape[1])])
        Z_df["tip"] = ids
        Z_df = Z_df.set_index("tip")

        disp_rows = []
        for cname, tips in clades.items():
            tips_in = [t for t in tips if t in Z_df.index]
            if len(tips_in) < 3:
                continue
            M = Z_df.loc[tips_in].to_numpy()
            centroid = M.mean(axis=0)
            d = np.sqrt(((M - centroid) ** 2).sum(axis=1))
            disp_rows.append({
                "hmm": hmm_name,
                "clade": cname,
                "n": int(len(tips_in)),
                "mean_dist_to_centroid": float(d.mean()),
                "median_dist_to_centroid": float(np.median(d)),
            })
        if disp_rows:
            pd.DataFrame(disp_rows).to_csv(out_disp, sep="\t", index=False)

    return cluster_assignments

def run_embed(cfg, hmm_to_seqs, clades, emb_dir, fasta_dir, hmm_keep,
              force=False, summary_dir=None, tax_map=None):
    print("\n[embed] Computing per-HMM embeddings...")
    emb_cfg = cfg.get("embeddings", {})

    # Resolve model_dir: explicit config value → {outdir}/models
    outdir = cfg.get("output", {}).get("outdir", None)
    if "model_dir" not in emb_cfg or emb_cfg["model_dir"] is None:
        if outdir:
            emb_cfg = dict(emb_cfg)  # shallow copy so we don't mutate the caller's dict
            emb_cfg["model_dir"] = os.path.join(outdir, "models")
            print(f"[embed] Model cache directory: {emb_cfg['model_dir']}")

    threads = int(cfg.get("resources", {}).get("cpu", 4))

    # ensure we have sequence bins
    if not hmm_to_seqs:
        for fp in glob.glob(os.path.join(fasta_dir, "*.faa")):
            hmm = os.path.basename(fp).rsplit(".", 1)[0]
            if hmm_keep is not None and hmm not in hmm_keep:
                continue
            hmm_to_seqs[hmm] = read_fasta(fp)

    all_cluster_assignments = []

    for hmm, seqs in hmm_to_seqs.items():
        if hmm_keep is not None and hmm not in hmm_keep:
            continue
        assignments = compute_embeddings_for_hmm(
            hmm, seqs, emb_cfg, emb_dir, force=force,
            clades=clades, tax_map=tax_map,
            outdir=outdir, summary_dir=summary_dir,
            threads=threads,
        )
        if assignments:
            all_cluster_assignments.extend(assignments)

    # Write combined clade_assignment.tsv
    if all_cluster_assignments and summary_dir:
        from ..utils.helpers import safe_mkdir
        safe_mkdir(summary_dir)
        out_fp = os.path.join(summary_dir, "clade_assignment.tsv")
        pd.DataFrame(all_cluster_assignments).to_csv(out_fp, sep="\t", index=False)
        print(f"[embed] Wrote cluster assignments: {out_fp}  ({len(all_cluster_assignments)} proteins)")


def embed_combined_with_ancestors(cfg, hmm_to_seqs, ancestral_seqs, clades,
                                   emb_dir, fasta_dir, hmm_keep,
                                   force=False, summary_dir=None, tax_map=None):
    """Embed modern + ancestral sequences in a single UMAP/HDBSCAN space.

    Parameters
    ----------
    ancestral_seqs : dict
        {node_id: amino_acid_sequence} from ASR parsing
    """
    print("\n[embed] Computing COMBINED embeddings (modern + ancestral)...")

    emb_cfg = cfg.get("embeddings", {})

    if "model_dir" not in emb_cfg or emb_cfg["model_dir"] is None:
        outdir = cfg.get("output", {}).get("outdir", None)
        if outdir:
            emb_cfg = dict(emb_cfg)
            emb_cfg["model_dir"] = os.path.join(outdir, "models")

    # Combine modern sequences from all HMMs
    combined_seqs = {}
    for hmm, seqs in hmm_to_seqs.items():
        if hmm_keep is not None and hmm not in hmm_keep:
            continue
        for tip, s in seqs.items():
            if tip not in combined_seqs:
                combined_seqs[tip] = s

    # Also load from disk if hmm_to_seqs is empty
    if not combined_seqs:
        for fp in glob.glob(os.path.join(fasta_dir, "*.faa")):
            hmm = os.path.basename(fp).rsplit(".", 1)[0]
            if hmm == "combined_all_hits":
                continue
            if hmm_keep is not None and hmm not in hmm_keep:
                continue
            seqs = read_fasta(fp)
            combined_seqs.update(seqs)

    # Add ancestral sequences with ANC| prefix
    for node_id, seq in ancestral_seqs.items():
        combined_seqs[f"ANC|{node_id}"] = seq

    if len(combined_seqs) < 3:
        import sys
        print("[embed] Combined + ancestral has <3 sequences. Skipping.", file=sys.stderr)
        return

    # Clean sequences
    combined_seqs = {
        k: v.replace(" ", "").replace("\n", "").replace("*", "").replace(".", "")
        for k, v in combined_seqs.items()
    }

    # Run embedding
    backend = emb_cfg.get("backend", "esm")
    device = emb_cfg.get("device", "cuda")
    model_name = emb_cfg.get("model", "esm2_t33_650M_UR50D")
    batch_size = emb_cfg.get("batch_size", 8)
    repr_layer = emb_cfg.get("repr_layer", None)
    model_dir = emb_cfg.get("model_dir", None)

    safe_mkdir(emb_dir)
    prefix = "combined_with_ancestors"
    out_npy = os.path.join(emb_dir, f"{prefix}.embeddings.npy")
    out_umap_png = os.path.join(emb_dir, f"{prefix}.umap.png")
    out_umap_clust_png = os.path.join(emb_dir, f"{prefix}.umap.clustered.png")

    if os.path.exists(out_npy) and not force:
        print(f"[embed] Combined embeddings already exist: {out_npy}")
        return

    print(f"[embed] Embedding {len(combined_seqs)} sequences "
          f"({sum(1 for k in combined_seqs if k.startswith('ANC|'))} ancestral, "
          f"{sum(1 for k in combined_seqs if not k.startswith('ANC|'))} modern)...")

    if backend == "esm":
        ids, X = _embed_esm(combined_seqs, model_name, device, batch_size,
                            repr_layer, model_dir=model_dir)
    else:
        ids, X = _embed_transformers(combined_seqs, model_name, device,
                                      batch_size, model_dir=model_dir)

    np.save(out_npy, X)

    # UMAP
    U = None
    try:
        import umap
        reducer = umap.UMAP(n_components=2, random_state=42, n_neighbors=15)
        U = reducer.fit_transform(X)

        # Save UMAP with ancestor/modern distinction
        _save_combined_umap(U, ids, out_umap_png, title="Combined Modern + Ancestral")

    except ImportError:
        import sys
        print("[embed] umap-learn not installed — skipping UMAP.", file=sys.stderr)
    except Exception as e:
        import sys
        print(f"[embed] UMAP failed for combined: {e}", file=sys.stderr)

    # HDBSCAN
    cluster_assignments = []
    if emb_cfg.get("cluster_embeddings", True):
        min_cs = int(emb_cfg.get("hdbscan_min_cluster_size", 5))
        labels = _run_hdbscan(X, min_cluster_size=min_cs)
        if labels is not None:
            for i, (tip, label) in enumerate(zip(ids, labels)):
                genome = tip.split("|", 1)[0] if "|" in tip else ""
                is_ancestral = tip.startswith("ANC|")
                tax = ""
                if tax_map and genome in tax_map:
                    tax = tax_map[genome]
                cluster_assignments.append({
                    "hmm": "combined",
                    "protein": tip,
                    "genome": genome,
                    "cluster_id": label,
                    "type": "ancestral" if is_ancestral else "modern",
                    "taxonomy": tax,
                })

            if U is not None:
                _save_combined_umap(U, ids, out_umap_clust_png,
                                    cluster_labels=labels,
                                    title="Combined UMAP (HDBSCAN clusters)")

    # Write combined assignments
    if cluster_assignments and summary_dir:
        safe_mkdir(summary_dir)
        out_fp = os.path.join(summary_dir, "clade_assignment.tsv")
        df = pd.DataFrame(cluster_assignments)
        # Append to existing if present
        if os.path.exists(out_fp):
            existing = pd.read_csv(out_fp, sep="\t")
            df = pd.concat([existing, df], ignore_index=True)
        df.to_csv(out_fp, sep="\t", index=False)
        print(f"[embed] Updated clade_assignment.tsv with combined embeddings "
              f"({len(cluster_assignments)} entries)")


def _save_combined_umap(U, ids, out_png, cluster_labels=None, title=""):
    """Save UMAP scatter with distinct markers for modern (circles) vs ancestral (triangles)."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(10, 8))

        is_anc = [i.startswith("ANC|") for i in ids]
        mod_idx = [j for j, a in enumerate(is_anc) if not a]
        anc_idx = [j for j, a in enumerate(is_anc) if a]

        if cluster_labels is not None:
            unique_labels = sorted(set(cluster_labels))
            cmap = plt.cm.get_cmap("tab20", max(len(unique_labels), 1))
            for label in unique_labels:
                color = "lightgrey" if label == -1 else cmap(unique_labels.index(label) % 20)
                lbl = "Noise" if label == -1 else f"Cluster {label}"
                # Modern: circles
                mask_mod = [j for j in mod_idx if cluster_labels[j] == label]
                if mask_mod:
                    ax.scatter(U[mask_mod, 0], U[mask_mod, 1], c=[color],
                              s=20, alpha=0.7, marker="o", edgecolors="none", label=lbl)
                # Ancestral: triangles
                mask_anc = [j for j in anc_idx if cluster_labels[j] == label]
                if mask_anc:
                    ax.scatter(U[mask_anc, 0], U[mask_anc, 1], c=[color],
                              s=40, alpha=0.9, marker="^", edgecolors="black",
                              linewidths=0.5,
                              label=f"{lbl} (ancestral)" if not mask_mod else "")
        else:
            if mod_idx:
                ax.scatter(U[mod_idx, 0], U[mod_idx, 1], c="steelblue",
                          s=20, alpha=0.7, marker="o", edgecolors="none", label="Modern")
            if anc_idx:
                ax.scatter(U[anc_idx, 0], U[anc_idx, 1], c="firebrick",
                          s=40, alpha=0.9, marker="^", edgecolors="black",
                          linewidths=0.5, label="Ancestral")

        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
        ax.set_title(title)
        ax.legend(fontsize=7, loc="best", framealpha=0.7, markerscale=1.5)
        fig.tight_layout()
        fig.savefig(out_png, dpi=200)
        plt.close(fig)
        print(f"[embed] Saved combined UMAP plot: {out_png}")
    except ImportError:
        pass
    except Exception as e:
        import sys
        print(f"[embed] Combined UMAP plot failed: {e}", file=sys.stderr)

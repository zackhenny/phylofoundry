"""
maape.py — MAAPE (Modular Assembly Analysis of Protein Embeddings) step.

Implements the 6-step MAAPE pipeline (https://github.com/Qinlab502/MAAPE) natively,
consuming embedding vectors produced by the embed step and writing evolutionary
network graphs, edge lists, and visualization plots to the maape/ output directory.

Pipeline steps
--------------
1. Load & normalise embeddings (PCA + L2 normalisation).
2. Sliding-window path generation — identify assembly paths between sub-vectors.
3. Co-occurrence matrix → weighted directed edge list.
4. KNN similarity graph construction (FAISS IndexFlatIP).
5. MAAPE network visualisation (directional weights overlaid on KNN graph).
6. Aggregated condensed visualisation (super-graph of clusters).
"""

from __future__ import annotations

import os
import sys
import glob
import pickle
import math
from concurrent.futures import ProcessPoolExecutor
from typing import Any

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Module-level constants
# ---------------------------------------------------------------------------

#: Default number of PCA dimensions used by the MAAPE paper when not otherwise
#: specified.  Applied as the upper bound when ``pca_components`` is ``null``
#: in the config (actual value is ``min(n_sequences, DEFAULT_PCA_COMPONENTS)``).
DEFAULT_PCA_COMPONENTS: int = 110

#: Fixed random seed for reproducible spring-layout placement.
_LAYOUT_SEED: int = 42


# ---------------------------------------------------------------------------
# Step 1 helpers — embedding loading and PCA/L2 normalisation
# ---------------------------------------------------------------------------

def _load_embeddings(hmm: str, emb_dir: str) -> tuple[list[str], np.ndarray]:
    """Load embedding vectors for *hmm* from the embed step outputs.

    Prefers the raw ``.npy`` embedding file produced when
    ``embeddings.write_full_vectors`` is enabled in the embed step.  Falls
    back to the PCA-reduced ``.pca.tsv`` table when the raw file is absent.

    Parameters
    ----------
    hmm:
        HMM / gene-family name.
    emb_dir:
        Path to the ``embeddings/`` directory produced by the embed step.

    Returns
    -------
    tuple[list[str], np.ndarray]
        ``(ids, X)`` where *ids* is the ordered list of sequence identifiers
        and *X* is the ``(n_seqs, n_dims)`` float32 embedding matrix.
    """
    raw_npy = os.path.join(emb_dir, f"{hmm}.embeddings.npy")
    pca_tsv = os.path.join(emb_dir, f"{hmm}.pca.tsv")

    if os.path.exists(raw_npy):
        X = np.load(raw_npy).astype(np.float32)
        # Companion IDs file — written alongside the .npy by the embed step
        ids_file = os.path.join(emb_dir, f"{hmm}.embeddings_ids.txt")
        if os.path.exists(ids_file):
            with open(ids_file) as fh:
                ids = [line.strip() for line in fh if line.strip()]
        else:
            # Fall back: generate synthetic IDs
            ids = [f"seq_{i}" for i in range(len(X))]
        print(f"[maape] Loaded raw embeddings for {hmm}: {X.shape}")
        return ids, X

    if os.path.exists(pca_tsv):
        df = pd.read_csv(pca_tsv, sep="\t")
        # The PCA TSV has a header; the first column is the sequence ID.
        id_col = df.columns[0]
        # Numeric PCA columns (PC1, PC2, … or similar)
        pca_cols = [c for c in df.columns if c != id_col]
        ids = df[id_col].tolist()
        X = df[pca_cols].values.astype(np.float32)
        print(f"[maape] Loaded PCA embeddings for {hmm} from TSV: {X.shape}")
        return ids, X

    print(f"[maape] WARNING: No embedding file found for {hmm} in {emb_dir}. "
          "Skipping.", file=sys.stderr)
    return [], np.empty((0, 0), dtype=np.float32)


def _run_pca_normalize(X: np.ndarray, n_components: int) -> np.ndarray:
    """PCA-reduce *X* to *n_components* dimensions and L2-normalise.

    Parameters
    ----------
    X:
        Input feature matrix ``(n_samples, n_features)``.
    n_components:
        Target dimensionality after PCA.

    Returns
    -------
    np.ndarray
        L2-normalised reduced matrix ``(n_samples, n_components)``.
    """
    from sklearn.decomposition import PCA

    n_samples, n_features = X.shape
    n_comp = min(n_components, n_samples, n_features)

    pca = PCA(n_components=n_comp, random_state=0)
    Z = pca.fit_transform(X).astype(np.float32)

    # L2 normalise each row
    norms = np.linalg.norm(Z, axis=1, keepdims=True)
    norms = np.where(norms == 0, 1.0, norms)
    Z = Z / norms
    return Z


def _l2_normalize(X: np.ndarray) -> np.ndarray:
    """L2-normalise each row of *X* in-place and return the result."""
    norms = np.linalg.norm(X, axis=1, keepdims=True)
    norms = np.where(norms == 0, 1.0, norms)
    return (X / norms).astype(np.float32)


# ---------------------------------------------------------------------------
# Step 2 helper — similarity threshold scaling
# ---------------------------------------------------------------------------

def _compute_similarity_thresholds(
    window_sizes: list[int],
    base_threshold: float,
) -> np.ndarray:
    """Compute per-window-size cosine similarity thresholds.

    MAAPE scales the base threshold (calibrated for ``window_size=5``) to
    larger windows using square-root scaling:

        threshold(w) = base_threshold * sqrt(w / 5)

    Parameters
    ----------
    window_sizes:
        List of sliding-window sizes (e.g. ``[5, 10, 20, 40, 80]``).
    base_threshold:
        Cosine similarity threshold at ``window_size=5``.

    Returns
    -------
    np.ndarray
        Array of thresholds, one per window size.
    """
    thresholds = np.array(
        [base_threshold * math.sqrt(w / 5.0) for w in window_sizes],
        dtype=np.float64,
    )
    return thresholds


# ---------------------------------------------------------------------------
# Step 2 — path generation
# ---------------------------------------------------------------------------

def _generate_paths(
    X: np.ndarray,
    ids: list[str],
    window_sizes: list[int],
    thresholds: np.ndarray,
) -> dict[tuple[int, int], list[tuple[int, int, int]]]:
    """Sliding-window path generation (MAAPE Step 2).

    For each window size ``w``, the embedding of each sequence is segmented
    into overlapping sub-vectors of length ``w``.  All cross-sequence sub-
    vector pairs are compared by cosine similarity; pairs whose similarity
    exceeds the per-window threshold are recorded as assembly *paths*.

    Parameters
    ----------
    X:
        L2-normalised embedding matrix ``(n_seqs, n_dims)``.
    ids:
        Sequence identifiers (length ``n_seqs``).
    window_sizes:
        Sliding window sizes.
    thresholds:
        Cosine similarity thresholds, one per window size.

    Returns
    -------
    dict
        ``{(i, j): [(window_size, pos_i, pos_j), ...]}`` mapping ordered
        ``(i, j)`` index pairs (``i < j``) to the list of assembly paths found
        between sequences *i* and *j*.
    """
    n_seqs, n_dims = X.shape
    paths: dict[tuple[int, int], list[tuple[int, int, int]]] = {}

    for w_idx, w in enumerate(window_sizes):
        thresh = float(thresholds[w_idx])
        if w > n_dims:
            continue

        # Extract all sub-vectors for this window size: shape (n_seqs, n_windows, w)
        n_windows = n_dims - w + 1
        # Build sub-vector matrix: (n_seqs * n_windows, w)
        sub_vecs = np.lib.stride_tricks.sliding_window_view(X, window_shape=w, axis=1)
        # sub_vecs shape: (n_seqs, n_windows, w) → flatten first two dims
        sv = sub_vecs.reshape(-1, w).astype(np.float64)

        # L2-normalise sub-vectors
        sv_norms = np.linalg.norm(sv, axis=1, keepdims=True)
        sv_norms = np.where(sv_norms == 0, 1.0, sv_norms)
        sv_normed = sv / sv_norms  # (n_seqs * n_windows, w)

        # For each pair of sequences compare their sub-vectors.
        # To avoid O(n²) full matrix multiply for large datasets,
        # process in blocks of 64 sequences.
        BLOCK = 64
        for i_start in range(0, n_seqs, BLOCK):
            i_end = min(i_start + BLOCK, n_seqs)
            for j_start in range(i_start, n_seqs, BLOCK):
                j_end = min(j_start + BLOCK, n_seqs)

                for i in range(i_start, i_end):
                    sv_i = sv_normed[i * n_windows:(i + 1) * n_windows]  # (n_windows, w)
                    j_lo = i + 1 if j_start == i_start else j_start
                    for j in range(j_lo, j_end):
                        sv_j = sv_normed[j * n_windows:(j + 1) * n_windows]
                        # Cosine similarities: (n_windows_i, n_windows_j)
                        sims = sv_i @ sv_j.T
                        # Find pairs exceeding threshold
                        pos_pairs = np.argwhere(sims >= thresh)
                        if pos_pairs.size == 0:
                            continue
                        key = (i, j)
                        entry = paths.setdefault(key, [])
                        for pp in pos_pairs:
                            entry.append((w, int(pp[0]), int(pp[1])))

    return paths


# ---------------------------------------------------------------------------
# Step 3 — co-occurrence matrix → directed weighted edge list
# ---------------------------------------------------------------------------

def _calculate_weights_and_edges(
    paths: dict[tuple[int, int], list[tuple[int, int, int]]],
    n_seqs: int,
) -> list[tuple[int, int, float]]:
    """Build the co-occurrence matrix and derive the directed edge list (Step 3).

    Each assembly path ``(w, pos_i, pos_j)`` stored for pair ``(i, j)``
    contributes +1 to ``C[i, j]``.  The directed edge goes from ``i`` to
    ``j`` when ``C[i, j] > C[j, i]`` (sequence *i* contributes more
    sub-vector segments to *j* than vice versa — interpreted as *j* having
    evolved from *i*).  The edge weight is the asymmetry ratio:

        weight = C[source, target] / (C[i,j] + C[j,i])

    Parameters
    ----------
    paths:
        Output of :func:`_generate_paths`.
    n_seqs:
        Total number of sequences.

    Returns
    -------
    list[tuple[int, int, float]]
        ``[(source_idx, target_idx, weight), ...]``.
    """
    # Build co-occurrence matrix
    C = np.zeros((n_seqs, n_seqs), dtype=np.int32)
    for (i, j), path_list in paths.items():
        # Count directional evidence from sub-vector position asymmetry.
        # pos_i < pos_j  → sub-vector of i aligns earlier in j  → i→j direction
        # pos_j < pos_i  → sub-vector of j aligns earlier in i  → j→i direction
        # pos_i == pos_j → symmetric, contributes equally to both directions
        for (w, pi, pj) in path_list:
            if pi < pj:
                C[i, j] += 1
            elif pj < pi:
                C[j, i] += 1
            else:
                C[i, j] += 1
                C[j, i] += 1

    edge_list: list[tuple[int, int, float]] = []
    for i in range(n_seqs):
        for j in range(i + 1, n_seqs):
            total = int(C[i, j]) + int(C[j, i])
            if total == 0:
                continue
            if C[i, j] >= C[j, i]:
                src, tgt = i, j
                w = float(C[i, j]) / total
            else:
                src, tgt = j, i
                w = float(C[j, i]) / total
            edge_list.append((src, tgt, w))

    return edge_list


# ---------------------------------------------------------------------------
# Step 4 — KNN graph construction
# ---------------------------------------------------------------------------

def _build_knn_graph(
    X: np.ndarray,
    ids: list[str],
    k: int,
    threshold: float,
    edge_list: list[tuple[int, int, float]],
) -> "Any":  # networkx.DiGraph
    """Build the KNN similarity graph and overlay directional edge weights (Step 4).

    Uses FAISS ``IndexFlatIP`` (inner-product search on L2-normalised vectors
    equals cosine similarity) to find the *k* nearest neighbours of each
    sequence.  KNN edges whose cosine similarity is below *threshold* are
    discarded.  Directional weights from *edge_list* are then mapped onto the
    surviving KNN edges to produce the final MAAPE directed graph.

    Parameters
    ----------
    X:
        L2-normalised embedding matrix ``(n_seqs, n_dims)``.
    ids:
        Sequence identifiers.
    k:
        Number of nearest neighbours per sequence.
    threshold:
        Minimum cosine similarity for an edge to be retained.
    edge_list:
        Directed edge list ``[(src, tgt, weight)]`` from Step 3.

    Returns
    -------
    networkx.DiGraph
        Directed graph with node attribute ``label`` and edge attributes
        ``cosine_sim``, ``direction_weight``, ``source``, ``target``.
    """
    try:
        import networkx as nx
    except ImportError:
        raise ImportError(
            "[maape] networkx is required for graph construction. "
            "Install it with: pip install networkx"
        )

    try:
        import faiss
        _faiss_available = True
    except ImportError:
        _faiss_available = False
        print(
            "[maape] faiss-cpu not installed — falling back to sklearn NearestNeighbors "
            "for KNN construction (slower for large datasets).",
            file=sys.stderr,
        )

    n_seqs = len(ids)
    X_f32 = X.astype(np.float32)

    # Build a lookup from (i, j) ordered pair → direction_weight
    dir_weight: dict[tuple[int, int], float] = {}
    for (src, tgt, w) in edge_list:
        dir_weight[(src, tgt)] = w

    G = nx.DiGraph()
    for idx, label in enumerate(ids):
        G.add_node(idx, label=label)

    k_eff = min(k + 1, n_seqs)  # +1 because FAISS returns the query itself

    if _faiss_available:
        index = faiss.IndexFlatIP(X_f32.shape[1])
        index.add(X_f32)
        sims, nbr_indices = index.search(X_f32, k_eff)
        # sims[i, 0] == 1.0 is the self-hit; skip it
        for i in range(n_seqs):
            for rank in range(k_eff):
                j = int(nbr_indices[i, rank])
                if j == i:
                    continue
                sim = float(sims[i, rank])
                if sim < threshold:
                    continue
                # Determine direction from edge_list
                if (i, j) in dir_weight:
                    G.add_edge(i, j, cosine_sim=sim,
                               direction_weight=dir_weight[(i, j)])
                elif (j, i) in dir_weight:
                    G.add_edge(j, i, cosine_sim=sim,
                               direction_weight=dir_weight[(j, i)])
                else:
                    # No direction determined — add undirected as i→j
                    G.add_edge(i, j, cosine_sim=sim, direction_weight=0.5)
    else:
        from sklearn.neighbors import NearestNeighbors
        nbrs = NearestNeighbors(n_neighbors=k_eff, metric="cosine", algorithm="auto")
        nbrs.fit(X_f32)
        distances, indices = nbrs.kneighbors(X_f32)
        for i in range(n_seqs):
            for rank in range(k_eff):
                j = int(indices[i, rank])
                if j == i:
                    continue
                sim = 1.0 - float(distances[i, rank])
                if sim < threshold:
                    continue
                if (i, j) in dir_weight:
                    G.add_edge(i, j, cosine_sim=sim,
                               direction_weight=dir_weight[(i, j)])
                elif (j, i) in dir_weight:
                    G.add_edge(j, i, cosine_sim=sim,
                               direction_weight=dir_weight[(j, i)])
                else:
                    G.add_edge(i, j, cosine_sim=sim, direction_weight=0.5)

    return G


# ---------------------------------------------------------------------------
# Step 5 — MAAPE network visualisation
# ---------------------------------------------------------------------------

def _visualize_maape(
    G: "Any",
    ids: list[str],
    out_png: str,
    hmm_name: str,
    color_map: dict[str, str] | None = None,
) -> None:
    """Render and save the MAAPE evolutionary network (Step 5).

    Nodes are coloured by clade/cluster when *color_map* is provided.
    Edge arrows point in the inferred evolutionary direction; edge width
    encodes the direction weight.

    Parameters
    ----------
    G:
        Directed graph from :func:`_build_knn_graph`.
    ids:
        Sequence identifiers aligned with node indices.
    out_png:
        Output PNG path.
    hmm_name:
        Used in the plot title.
    color_map:
        Optional ``{sequence_id: hex_colour}`` dict for node colouring.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.cm as cm
        import networkx as nx
    except ImportError:
        print("[maape] matplotlib/networkx not available — skipping visualisation.",
              file=sys.stderr)
        return

    try:
        fig, ax = plt.subplots(figsize=(12, 10))

        pos = nx.spring_layout(G, seed=_LAYOUT_SEED, k=1.0 / max(1, math.sqrt(len(G))))

        # Node colours
        if color_map:
            node_colors = [color_map.get(ids[n], "#aaaaaa") for n in G.nodes()]
        else:
            node_colors = ["#4a90d9"] * len(G.nodes())

        nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=60,
                               alpha=0.85, ax=ax)

        # Edge widths / colours scaled by direction_weight
        edges = list(G.edges(data=True))
        if edges:
            weights = [d.get("direction_weight", 0.5) for _, _, d in edges]
            edge_widths = [0.5 + 2.5 * w for w in weights]
            edge_colors = [cm.plasma(w) for w in weights]
            nx.draw_networkx_edges(
                G, pos,
                edgelist=[(u, v) for u, v, _ in edges],
                width=edge_widths,
                edge_color=edge_colors,
                arrows=True,
                arrowsize=8,
                arrowstyle="-|>",
                connectionstyle="arc3,rad=0.05",
                ax=ax,
            )

        ax.set_title(
            f"{hmm_name} — MAAPE Evolutionary Network\n"
            f"({len(G.nodes())} sequences, {len(G.edges())} directed edges)",
            fontsize=11,
        )
        ax.axis("off")
        fig.tight_layout()
        os.makedirs(os.path.dirname(out_png), exist_ok=True)
        fig.savefig(out_png, dpi=200)
        plt.close(fig)
        print(f"[maape] Saved MAAPE network plot: {out_png}")
    except Exception as exc:
        print(f"[maape] Visualisation failed for {hmm_name}: {exc}", file=sys.stderr)


# ---------------------------------------------------------------------------
# Step 6 — aggregated condensed visualisation
# ---------------------------------------------------------------------------

def _visualize_aggregated(
    G: "Any",
    ids: list[str],
    labels: list[int] | None,
    out_png: str,
    hmm_name: str,
) -> None:
    """Render a condensed super-graph where nodes are clusters (Step 6).

    Clusters are derived from *labels* (HDBSCAN/Leiden cluster IDs from the
    embed step) or from Louvain/Leiden community detection on *G* when labels
    are unavailable.  Each super-node represents one cluster; super-edges
    summarise the directed flow between clusters.

    Parameters
    ----------
    G:
        MAAPE directed graph from Step 5.
    ids:
        Sequence identifiers.
    labels:
        Optional per-sequence cluster labels (``-1`` = noise).
    out_png:
        Output PNG path.
    hmm_name:
        Used in the plot title.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.cm as cm
        import networkx as nx
    except ImportError:
        print("[maape] matplotlib/networkx not available — skipping aggregated plot.",
              file=sys.stderr)
        return

    try:
        n_seqs = len(ids)

        # Assign cluster membership
        if labels is not None and len(labels) == n_seqs:
            cluster_map = {i: int(labels[i]) for i in range(n_seqs)}
        else:
            # Community detection on G using greedy modularity
            undirected = G.to_undirected()
            try:
                communities = nx.algorithms.community.greedy_modularity_communities(
                    undirected
                )
                cluster_map = {}
                for cid, community in enumerate(communities):
                    for node in community:
                        cluster_map[node] = cid
            except Exception:
                cluster_map = {i: 0 for i in range(n_seqs)}

        # Build super-graph
        S = nx.DiGraph()
        cluster_ids = sorted(set(cluster_map.values()))
        for cid in cluster_ids:
            count = sum(1 for v in cluster_map.values() if v == cid)
            S.add_node(cid, count=count)

        super_edges: dict[tuple[int, int], float] = {}
        for u, v, data in G.edges(data=True):
            cu = cluster_map.get(u, -1)
            cv = cluster_map.get(v, -1)
            if cu == cv or cu == -1 or cv == -1:
                continue
            key = (cu, cv)
            w = data.get("direction_weight", 0.5)
            super_edges[key] = super_edges.get(key, 0.0) + w

        for (cu, cv), total_w in super_edges.items():
            S.add_edge(cu, cv, weight=total_w)

        fig, ax = plt.subplots(figsize=(10, 8))
        pos = nx.spring_layout(S, seed=0)

        node_sizes = [200 + 20 * S.nodes[n].get("count", 1) for n in S.nodes()]
        cmap = plt.cm.get_cmap("tab20", max(len(cluster_ids), 1))
        node_colors = [cmap(i % 20) for i in range(len(cluster_ids))]

        nx.draw_networkx_nodes(S, pos, node_size=node_sizes,
                               node_color=node_colors, alpha=0.9, ax=ax)
        nx.draw_networkx_labels(S, pos, font_size=8, ax=ax)

        s_edges = list(S.edges(data=True))
        if s_edges:
            max_w = max(d.get("weight", 1.0) for _, _, d in s_edges) or 1.0
            e_widths = [0.5 + 3.0 * d.get("weight", 0.5) / max_w
                        for _, _, d in s_edges]
            nx.draw_networkx_edges(
                S, pos,
                edgelist=[(u, v) for u, v, _ in s_edges],
                width=e_widths,
                arrows=True,
                arrowsize=12,
                arrowstyle="-|>",
                connectionstyle="arc3,rad=0.08",
                ax=ax,
            )

        ax.set_title(
            f"{hmm_name} — MAAPE Aggregated Evolutionary Network\n"
            f"({len(cluster_ids)} clusters)",
            fontsize=11,
        )
        ax.axis("off")
        fig.tight_layout()
        os.makedirs(os.path.dirname(out_png), exist_ok=True)
        fig.savefig(out_png, dpi=200)
        plt.close(fig)
        print(f"[maape] Saved aggregated plot: {out_png}")
    except Exception as exc:
        print(f"[maape] Aggregated visualisation failed for {hmm_name}: {exc}",
              file=sys.stderr)


# ---------------------------------------------------------------------------
# Per-HMM orchestrator
# ---------------------------------------------------------------------------

def _load_cluster_labels(hmm: str, emb_dir: str) -> list[int] | None:
    """Load HDBSCAN cluster labels from the embed step TSV outputs.

    Searches for the ``<hmm>.embeddings.tsv`` or ``<hmm>.pca.tsv`` table that
    contains a ``cluster`` column, and returns the labels in the same order as
    the sequence IDs loaded by :func:`_load_embeddings`.

    Returns ``None`` when no cluster labels are found.
    """
    for suffix in ("embeddings.tsv", "pca.tsv"):
        tsv = os.path.join(emb_dir, f"{hmm}.{suffix}")
        if not os.path.exists(tsv):
            continue
        try:
            df = pd.read_csv(tsv, sep="\t")
            if "cluster" in df.columns:
                return df["cluster"].tolist()
        except Exception:
            continue
    return None


def _load_clade_color_map(hmm: str, clade_assign_dir: str | None) -> dict[str, str]:
    """Load clade/cluster assignments and map them to distinct hex colours.

    Returns an empty dict when no assignment file is found.
    """
    if not clade_assign_dir:
        return {}
    tsv = os.path.join(clade_assign_dir, f"{hmm}.tsv")
    if not os.path.exists(tsv):
        return {}
    try:
        import matplotlib.cm as cm
        import matplotlib.colors as mcolors
        df = pd.read_csv(tsv, sep="\t")
        # Expected columns: protein_id (or similar) and clade_name / cluster
        id_col = df.columns[0]
        clade_col = None
        for c in ("clade_name", "cluster", "clade"):
            if c in df.columns:
                clade_col = c
                break
        if clade_col is None:
            return {}
        clades = df[clade_col].unique().tolist()
        cmap = cm.get_cmap("tab20", max(len(clades), 1))
        clade_to_color = {
            clade: mcolors.to_hex(cmap(i % 20))
            for i, clade in enumerate(clades)
        }
        color_map: dict[str, str] = {}
        for _, row in df.iterrows():
            color_map[str(row[id_col])] = clade_to_color.get(row[clade_col], "#aaaaaa")
        return color_map
    except Exception:
        return {}


def run_maape(
    cfg: dict,
    hmm_name: str,
    emb_dir: str,
    maape_dir: str,
    clade_assign_dir: str | None = None,
    force: bool = False,
) -> dict:
    """Orchestrate the full MAAPE pipeline for a single HMM gene family.

    Parameters
    ----------
    cfg:
        Pipeline configuration dict (``cfg["maape"]`` is the MAAPE section).
    hmm_name:
        HMM / gene-family name.
    emb_dir:
        Path to the ``embeddings/`` directory.
    maape_dir:
        Path to the MAAPE output directory.
    clade_assign_dir:
        Optional path to ``clade_assignments/`` for node colouring.
    force:
        Overwrite existing output files.

    Returns
    -------
    dict
        Summary stats: ``{hmm_name, n_nodes, n_edges, mean_edge_weight,
        n_paths}``.
    """
    maape_cfg = cfg.get("maape", {})
    window_sizes = maape_cfg.get("window_sizes", [5, 10, 20, 40, 80])
    knn_k = int(maape_cfg.get("knn_k", 20))
    knn_threshold = float(maape_cfg.get("knn_threshold", 0.5))
    pca_components_cfg = maape_cfg.get("pca_components", None)
    base_threshold = float(maape_cfg.get("similarity_threshold_base", 0.00001))
    reuse_pca = bool(maape_cfg.get("reuse_pca", False))
    generate_aggregated = bool(maape_cfg.get("generate_aggregated", True))
    # When color_scheme is set to a dict/mapping in the config it overrides the
    # auto-detected clade colour map loaded from clade_assignments/.
    cfg_color_scheme = maape_cfg.get("color_scheme", None)

    os.makedirs(maape_dir, exist_ok=True)

    network_pkl = os.path.join(maape_dir, f"{hmm_name}.maape_network.pkl")
    paths_pkl = os.path.join(maape_dir, f"{hmm_name}.paths.pkl")
    edgelist_txt = os.path.join(maape_dir, f"{hmm_name}.edge_list.txt")
    net_png = os.path.join(maape_dir, f"{hmm_name}.maape_network.png")
    agg_png = os.path.join(maape_dir, f"{hmm_name}.maape_aggregated.png")

    # ── Step 1: Load embeddings ──────────────────────────────────────────────
    ids, X_raw = _load_embeddings(hmm_name, emb_dir)
    if len(ids) == 0:
        print(f"[maape] No sequences for {hmm_name} — skipping.", file=sys.stderr)
        return {"hmm_name": hmm_name, "n_nodes": 0, "n_edges": 0,
                "mean_edge_weight": 0.0, "n_paths": 0}

    n_seqs = len(ids)

    if reuse_pca:
        # Use whatever dimensionality the embed step already produced
        X = _l2_normalize(X_raw)
    else:
        # Run MAAPE's own PCA as in the original paper
        n_comp = pca_components_cfg if pca_components_cfg else min(n_seqs, DEFAULT_PCA_COMPONENTS)
        n_comp = min(n_comp, n_seqs, X_raw.shape[1])
        if n_comp < X_raw.shape[1]:
            print(f"[maape] Running PCA({n_comp}) for {hmm_name}...")
            X = _run_pca_normalize(X_raw, n_comp)
        else:
            X = _l2_normalize(X_raw)

    # ── Step 2: Path generation ──────────────────────────────────────────────
    if os.path.exists(paths_pkl) and not force:
        print(f"[maape] Loading cached paths for {hmm_name}...")
        with open(paths_pkl, "rb") as fh:
            paths = pickle.load(fh)
    else:
        print(f"[maape] Generating assembly paths for {hmm_name} "
              f"({n_seqs} seqs, windows={window_sizes})...")
        thresholds = _compute_similarity_thresholds(window_sizes, base_threshold)
        paths = _generate_paths(X, ids, window_sizes, thresholds)
        with open(paths_pkl, "wb") as fh:
            pickle.dump(paths, fh, protocol=4)
        print(f"[maape] Saved paths: {paths_pkl} ({len(paths)} sequence pairs)")

    # ── Step 3: Weights and directed edges ───────────────────────────────────
    edge_list = _calculate_weights_and_edges(paths, n_seqs)

    # Write human-readable edge list (protein IDs).
    # All edge indices are guaranteed to be valid offsets into `ids` because
    # _calculate_weights_and_edges only produces indices in [0, n_seqs).
    if not os.path.exists(edgelist_txt) or force:
        with open(edgelist_txt, "w") as fh:
            fh.write("source\ttarget\tweight\n")
            for (src, tgt, w) in edge_list:
                assert 0 <= src < len(ids) and 0 <= tgt < len(ids), (
                    f"[maape] Edge index out of range for {hmm_name}: "
                    f"src={src}, tgt={tgt}, n_seqs={len(ids)}"
                )
                fh.write(f"{ids[src]}\t{ids[tgt]}\t{w:.6f}\n")
        print(f"[maape] Wrote edge list: {edgelist_txt} ({len(edge_list)} edges)")

    # ── Step 4: KNN graph ────────────────────────────────────────────────────
    if os.path.exists(network_pkl) and not force:
        print(f"[maape] Loading cached network for {hmm_name}...")
        with open(network_pkl, "rb") as fh:
            G = pickle.load(fh)
    else:
        print(f"[maape] Building KNN graph for {hmm_name} (k={knn_k})...")
        G = _build_knn_graph(X, ids, knn_k, knn_threshold, edge_list)
        with open(network_pkl, "wb") as fh:
            pickle.dump(G, fh, protocol=4)
        print(f"[maape] Saved network: {network_pkl}")

    n_edges = len(G.edges())
    mean_weight = (
        float(
            np.mean([d.get("direction_weight", 0.5)
                     for _, _, d in G.edges(data=True)])
        )
        if n_edges > 0 else 0.0
    )

    # ── Step 5: MAAPE network visualisation ──────────────────────────────────
    if not os.path.exists(net_png) or force:
        # Prefer user-supplied color_scheme from config over auto-detected clade map
        if isinstance(cfg_color_scheme, dict):
            color_map: dict[str, str] | None = cfg_color_scheme
        else:
            color_map = _load_clade_color_map(hmm_name, clade_assign_dir) or None
        _visualize_maape(G, ids, net_png, hmm_name, color_map=color_map)

    # ── Step 6: Aggregated visualisation ─────────────────────────────────────
    if generate_aggregated and (not os.path.exists(agg_png) or force):
        cluster_labels = _load_cluster_labels(hmm_name, emb_dir)
        _visualize_aggregated(G, ids, cluster_labels, agg_png, hmm_name)

    print(f"[maape] Complete for {hmm_name}: "
          f"{n_seqs} nodes, {n_edges} edges, mean_weight={mean_weight:.4f}")

    return {
        "hmm_name": hmm_name,
        "n_nodes": n_seqs,
        "n_edges": n_edges,
        "mean_edge_weight": round(mean_weight, 6),
        "n_paths": len(paths),
    }


# ---------------------------------------------------------------------------
# Multi-HMM orchestrator (called from pipeline.py)
# ---------------------------------------------------------------------------

def _worker_maape(args_pack: tuple) -> dict:
    """Worker function for parallel per-HMM MAAPE execution."""
    cfg, hmm_name, emb_dir, maape_dir, clade_assign_dir, force = args_pack
    try:
        return run_maape(cfg, hmm_name, emb_dir, maape_dir, clade_assign_dir, force)
    except Exception as exc:
        import traceback
        print(f"[maape] ERROR for {hmm_name}: {exc}\n{traceback.format_exc()}",
              file=sys.stderr)
        return {"hmm_name": hmm_name, "n_nodes": 0, "n_edges": 0,
                "mean_edge_weight": 0.0, "n_paths": 0, "error": str(exc)}


def run_maape_all(
    cfg: dict,
    emb_dir: str,
    maape_dir: str,
    hmm_list: list[str] | None = None,
    clade_assign_dir: str | None = None,
    force: bool = False,
    n_jobs: int = 1,
) -> list[dict]:
    """Run MAAPE for all HMMs in *emb_dir*, writing results to *maape_dir*.

    Parameters
    ----------
    cfg:
        Pipeline configuration dict.
    emb_dir:
        Path to the ``embeddings/`` directory containing ``.embeddings.npy``
        or ``.pca.tsv`` files.
    maape_dir:
        Path to the MAAPE output directory.
    hmm_list:
        Explicit list of HMM names to process.  When ``None``, the function
        discovers HMMs from the ``.embeddings.npy`` / ``.pca.tsv`` files in
        *emb_dir*.
    clade_assign_dir:
        Optional path to ``clade_assignments/`` for node colouring.
    force:
        Overwrite existing output files.
    n_jobs:
        Number of parallel workers (one per HMM).

    Returns
    -------
    list[dict]
        Per-HMM summary dicts; also written to ``maape/maape_summary.tsv``.
    """
    os.makedirs(maape_dir, exist_ok=True)

    maape_cfg = cfg.get("maape", {})
    per_hmm = bool(maape_cfg.get("per_hmm", True))

    if not per_hmm:
        print("[maape] per_hmm=false — skipping individual HMM processing.")
        return []

    # Discover HMM names from embedding files
    if hmm_list is None:
        npy_files = glob.glob(os.path.join(emb_dir, "*.embeddings.npy"))
        tsv_files = glob.glob(os.path.join(emb_dir, "*.pca.tsv"))
        names_npy = {os.path.basename(f).replace(".embeddings.npy", "")
                     for f in npy_files}
        names_tsv = {os.path.basename(f).replace(".pca.tsv", "")
                     for f in tsv_files}
        hmm_list = sorted(names_npy | names_tsv)

    if not hmm_list:
        print("[maape] No embedding files found in emb_dir — nothing to do.",
              file=sys.stderr)
        return []

    print(f"[maape] Running MAAPE for {len(hmm_list)} HMM(s) "
          f"with {n_jobs} worker(s)...")

    tasks = [
        (cfg, hmm, emb_dir, maape_dir, clade_assign_dir, force)
        for hmm in hmm_list
    ]

    results: list[dict] = []
    if n_jobs > 1 and len(tasks) > 1:
        workers = min(n_jobs, len(tasks))
        with ProcessPoolExecutor(max_workers=workers) as exe:
            for res in exe.map(_worker_maape, tasks):
                results.append(res)
    else:
        for task in tasks:
            results.append(_worker_maape(task))

    # Write summary TSV
    summary_fp = os.path.join(maape_dir, "maape_summary.tsv")
    if results:
        df = pd.DataFrame(results)
        df.to_csv(summary_fp, sep="\t", index=False)
        print(f"[maape] Wrote summary: {summary_fp}")

    return results

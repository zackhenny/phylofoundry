import os
import glob
import warnings
import pandas as pd
import numpy as np
from ..utils.bio import read_fasta
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
    """Generate and save a UMAP scatter plot as PNG."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(8, 6))

        if cluster_labels is not None:
            unique_labels = sorted(set(cluster_labels))
            cmap = plt.cm.get_cmap("tab20", max(len(unique_labels), 1))
            for label in unique_labels:
                mask = [cl == label for cl in cluster_labels]
                pts = U[[i for i, m in enumerate(mask) if m]]
                color = "lightgrey" if label == -1 else cmap(unique_labels.index(label) % 20)
                lbl = "Noise" if label == -1 else f"Cluster {label}"
                ax.scatter(pts[:, 0], pts[:, 1], c=[color], s=20, alpha=0.7, label=lbl, edgecolors="none")
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
                    ax.scatter(pts[:, 0], pts[:, 1], c=[cmap(ci)], s=20, alpha=0.7, label=cn, edgecolors="none")
            # Unassigned
            mask = [a is None for a in assigned]
            pts = U[[i for i, m in enumerate(mask) if m]]
            if len(pts) > 0:
                ax.scatter(pts[:, 0], pts[:, 1], c="lightgrey", s=15, alpha=0.5, label="unassigned", edgecolors="none")
            ax.legend(fontsize=7, loc="best", framealpha=0.7, markerscale=1.5)
        else:
            ax.scatter(U[:, 0], U[:, 1], c="steelblue", s=20, alpha=0.7, edgecolors="none")

        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
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


def _run_hdbscan(X, min_cluster_size=5):
    """Cluster embeddings with HDBSCAN. Returns list of integer labels (-1 = noise)."""
    try:
        from sklearn.cluster import HDBSCAN
        clusterer = HDBSCAN(min_cluster_size=min_cluster_size)
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


def compute_embeddings_for_hmm(hmm_name: str, seqs: dict, emb_cfg: dict, outdir_embeddings: str,
                               force: bool, clades: dict | None, tax_map: dict | None = None):
    """Compute embeddings, PCA, UMAP, HDBSCAN clustering, and save plots + TSVs.
    Returns a list of cluster assignment dicts (for clade_assignment.tsv), or empty list."""
    safe_mkdir(outdir_embeddings)

    out_npy = os.path.join(outdir_embeddings, f"{hmm_name}.embeddings.npy")
    out_pca = os.path.join(outdir_embeddings, f"{hmm_name}.pca.tsv")
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
        import sys
        print(f"[embed] Warning: HMM '{hmm_name}' has less than 3 sequences. Skipping embeddings.", file=sys.stderr)
        return []

    backend = emb_cfg["backend"]
    device = emb_cfg["device"]
    batch_size = int(emb_cfg["batch_size"])
    model_name = emb_cfg["model"]
    repr_layer = emb_cfg.get("repr_layer", None)

    model_dir = emb_cfg.get("model_dir", None)

    try:
        import sys
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

    # PCA
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

    # ── UMAP ──────────────────────────────────────────────────────────────
    U = None
    try:
        import umap
        import warnings as _w
        reducer = umap.UMAP(n_components=2, random_state=42)
        with _w.catch_warnings():
            _w.simplefilter("ignore")
            U = reducer.fit_transform(X)
        
        u_rows = []
        for tip, coords in zip(ids, U):
            genome = tip.split("|", 1)[0] if "|" in tip else "Unknown"
            protein = tip.split("|", 1)[1] if "|" in tip else tip
            u_rows.append({
                "hmm": hmm_name,
                "tip": tip,
                "genome": genome,
                "protein": protein,
                "UMAP1": float(coords[0]),
                "UMAP2": float(coords[1])
            })
        pd.DataFrame(u_rows).to_csv(out_umap, sep="\t", index=False)

        # Save UMAP plot (colored by clades if available)
        _save_umap_plot(U, ids, hmm_name, out_umap_png, clades=clades)
    except ImportError:
        import sys
        print("[embed] UMAP skip: umap-learn not installed.", file=sys.stderr)
    except Exception as e:
        import sys
        print(f"[embed] UMAP failed for {hmm_name}: {e}", file=sys.stderr)

    # ── HDBSCAN clustering ────────────────────────────────────────────────
    cluster_assignments = []
    if emb_cfg.get("cluster_embeddings", True):
        min_cs = int(emb_cfg.get("hdbscan_min_cluster_size", 5))
        labels = _run_hdbscan(X, min_cluster_size=min_cs)
        if labels is not None:
            # Genome ID normalizer for taxonomy lookup
            from ..utils.helpers import normalize_genome_id
            n_clusters = len(set(labels) - {-1})
            print(f"[embed] HDBSCAN found {n_clusters} clusters for {hmm_name} ({sum(1 for l in labels if l == -1)} noise points)")

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

            # Save cluster-colored UMAP plot
            if U is not None:
                _save_umap_plot(U, ids, hmm_name, out_umap_clust_png,
                                cluster_labels=labels, title_suffix=" (HDBSCAN clusters)")

    # metadata
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
    if "model_dir" not in emb_cfg or emb_cfg["model_dir"] is None:
        outdir = cfg.get("output", {}).get("outdir", None)
        if outdir:
            emb_cfg = dict(emb_cfg)  # shallow copy so we don't mutate the caller's dict
            emb_cfg["model_dir"] = os.path.join(outdir, "models")
            print(f"[embed] Model cache directory: {emb_cfg['model_dir']}")

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
            clades=clades, tax_map=tax_map
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

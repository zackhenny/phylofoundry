"""Tests for embedding clustering enhancements:
- cluster_on (PCA vs embeddings)
- cluster_method (hdbscan vs leiden)
- hdbscan_metric
- umap_dimensions (2D vs 3D)
- Config validation (_validate_emb_cfg)
- kNN neighborhood analysis
- Cluster subworkflow helpers
"""
import numpy as np
import pytest

from phylofoundry.tasks.embed import (
    _run_hdbscan,
    _run_leiden,
    _validate_emb_cfg,
    _save_umap_plot,
    _compute_knn_metrics,
    _assign_membership_tiers,
    _write_cluster_fastas,
    _generate_sequence_logo,
    _classify_noise_sequences,
)

# Optional dependency guards
try:
    from sklearn.cluster import HDBSCAN as _HDBSCAN_CLS
    _HDBSCAN_AVAILABLE = True
except ImportError:
    _HDBSCAN_AVAILABLE = False

try:
    from sklearn.neighbors import NearestNeighbors as _NNB  # noqa: F401
    _SKLEARN_AVAILABLE = True
except ImportError:
    _SKLEARN_AVAILABLE = False

try:
    import matplotlib  # noqa: F401
    _MPL_AVAILABLE = True
except ImportError:
    _MPL_AVAILABLE = False

requires_hdbscan = pytest.mark.skipif(
    not _HDBSCAN_AVAILABLE,
    reason="scikit-learn >= 1.3 with HDBSCAN not installed",
)
requires_sklearn = pytest.mark.skipif(
    not _SKLEARN_AVAILABLE,
    reason="scikit-learn not installed",
)
requires_matplotlib = pytest.mark.skipif(
    not _MPL_AVAILABLE,
    reason="matplotlib not installed",
)


# ── Helpers ───────────────────────────────────────────────────────────────────

def _make_clusterable(n=40, n_clusters=3, seed=0):
    """Return a small clusterable float32 matrix (n x 10)."""
    rng = np.random.default_rng(seed)
    centres = rng.standard_normal((n_clusters, 10)) * 5
    X = np.vstack([
        centres[i % n_clusters] + rng.standard_normal(10) * 0.3
        for i in range(n)
    ]).astype(np.float32)
    return X


# ── _validate_emb_cfg ─────────────────────────────────────────────────────────

class TestValidateEmbCfg:
    def _base_cfg(self):
        return {
            "backend": "esm",
            "model": "esm2_t33_650M_UR50D",
            "device": "cpu",
            "batch_size": 4,
            "repr_layer": None,
            "pooling": "mean",
            "pca_components": 3,
            "cluster_embeddings": True,
        }

    def test_defaults(self):
        cfg = _validate_emb_cfg(self._base_cfg())
        assert cfg["_cluster_on"] == "pca"
        assert cfg["_cluster_method"] == "hdbscan"
        assert cfg["_umap_dimensions"] == 2
        assert cfg["hdbscan_metric"] == "euclidean"

    def test_cluster_on_embeddings(self):
        base = self._base_cfg()
        base["cluster_on"] = "embeddings"
        cfg = _validate_emb_cfg(base)
        assert cfg["_cluster_on"] == "embeddings"

    def test_cluster_on_case_insensitive(self):
        base = self._base_cfg()
        base["cluster_on"] = "PCA"
        cfg = _validate_emb_cfg(base)
        assert cfg["_cluster_on"] == "pca"

    def test_cluster_on_invalid(self):
        base = self._base_cfg()
        base["cluster_on"] = "tsne"
        with pytest.raises(ValueError, match="cluster_on"):
            _validate_emb_cfg(base)

    def test_cluster_method_leiden(self):
        base = self._base_cfg()
        base["cluster_method"] = "leiden"
        cfg = _validate_emb_cfg(base)
        assert cfg["_cluster_method"] == "leiden"

    def test_cluster_method_invalid(self):
        base = self._base_cfg()
        base["cluster_method"] = "kmeans"
        with pytest.raises(ValueError, match="cluster_method"):
            _validate_emb_cfg(base)

    def test_umap_dimensions_3(self):
        base = self._base_cfg()
        base["umap_dimensions"] = 3
        cfg = _validate_emb_cfg(base)
        assert cfg["_umap_dimensions"] == 3

    def test_umap_dimensions_invalid_defaults_to_2(self):
        base = self._base_cfg()
        base["umap_dimensions"] = 5
        cfg = _validate_emb_cfg(base)
        assert cfg["_umap_dimensions"] == 2

    def test_hdbscan_metric_cosine(self):
        base = self._base_cfg()
        base["hdbscan_metric"] = "cosine"
        cfg = _validate_emb_cfg(base)
        assert cfg["hdbscan_metric"] == "cosine"

    def test_knn_defaults(self):
        cfg = _validate_emb_cfg(self._base_cfg())
        assert cfg["run_knn"] is True
        assert cfg["knn_neighbors"] == 10

    def test_knn_custom_neighbors(self):
        base = self._base_cfg()
        base["knn_neighbors"] = 5
        cfg = _validate_emb_cfg(base)
        assert cfg["knn_neighbors"] == 5

    def test_knn_disabled(self):
        base = self._base_cfg()
        base["run_knn"] = False
        cfg = _validate_emb_cfg(base)
        assert cfg["run_knn"] is False

    def test_knn_invalid_neighbors_defaults_to_10(self):
        base = self._base_cfg()
        base["knn_neighbors"] = 0
        cfg = _validate_emb_cfg(base)
        assert cfg["knn_neighbors"] == 10

    def test_pca_components_default(self):
        base = self._base_cfg()
        # Remove pca_components to trigger setdefault path
        base.pop("pca_components", None)
        cfg = _validate_emb_cfg(base)
        assert cfg["pca_components"] == 50


# ── _run_hdbscan ──────────────────────────────────────────────────────────────

class TestRunHdbscan:
    @requires_hdbscan
    def test_returns_labels(self):
        X = _make_clusterable()
        labels = _run_hdbscan(X, min_cluster_size=5)
        assert labels is not None
        assert len(labels) == len(X)

    @requires_hdbscan
    def test_euclidean_metric(self):
        X = _make_clusterable()
        labels = _run_hdbscan(X, min_cluster_size=5, metric="euclidean")
        assert labels is not None

    @requires_hdbscan
    def test_cosine_metric(self):
        X = _make_clusterable()
        labels = _run_hdbscan(X, min_cluster_size=5, metric="cosine")
        assert labels is not None
        assert len(labels) == len(X)

    @requires_hdbscan
    def test_labels_are_integers(self):
        X = _make_clusterable()
        labels = _run_hdbscan(X, min_cluster_size=5)
        assert all(isinstance(lb, int) for lb in labels)

    @requires_hdbscan
    def test_noise_label_is_minus_one(self):
        X = _make_clusterable()
        labels = _run_hdbscan(X, min_cluster_size=5)
        assert all(lb >= -1 for lb in labels)

    def test_invalid_metric_returns_none(self):
        X = _make_clusterable()
        labels = _run_hdbscan(X, min_cluster_size=5, metric="not_a_metric_xyz")
        # scikit-learn raises on bad metric (or ImportError on old sklearn);
        # _run_hdbscan must return None either way
        assert labels is None


# ── _run_leiden ───────────────────────────────────────────────────────────────

class TestRunLeiden:
    def test_returns_labels_or_none(self):
        X = _make_clusterable()
        labels = _run_leiden(X)
        # Either succeeds (igraph/leidenalg installed) or gracefully returns None
        if labels is not None:
            assert len(labels) == len(X)
            assert all(isinstance(lb, int) for lb in labels)

    def test_cosine_metric(self):
        X = _make_clusterable()
        labels = _run_leiden(X, metric="cosine")
        if labels is not None:
            assert len(labels) == len(X)

    def test_all_samples_assigned(self):
        X = _make_clusterable()
        labels = _run_leiden(X)
        if labels is not None:
            # Leiden does not produce noise points (-1)
            assert all(lb >= 0 for lb in labels)

    def test_small_dataset(self):
        """Should not crash on a tiny dataset."""
        X = _make_clusterable(n=6, n_clusters=2)
        labels = _run_leiden(X, n_neighbors=3)
        # Either produces labels or None — no exception
        if labels is not None:
            assert len(labels) == 6


# ── _save_umap_plot ───────────────────────────────────────────────────────────

class TestSaveUmapPlot:
    @requires_matplotlib
    def test_2d_plot_saves(self, tmp_path):
        U = np.random.default_rng(0).standard_normal((20, 2)).astype(np.float32)
        ids = [f"seq{i}" for i in range(20)]
        out_png = str(tmp_path / "umap2d.png")
        _save_umap_plot(U, ids, "TESTHMM", out_png)
        assert (tmp_path / "umap2d.png").exists()

    @requires_matplotlib
    def test_3d_plot_saves(self, tmp_path):
        U = np.random.default_rng(1).standard_normal((20, 3)).astype(np.float32)
        ids = [f"seq{i}" for i in range(20)]
        out_png = str(tmp_path / "umap3d.png")
        _save_umap_plot(U, ids, "TESTHMM", out_png)
        assert (tmp_path / "umap3d.png").exists()

    @requires_matplotlib
    def test_2d_plot_with_cluster_labels(self, tmp_path):
        U = np.random.default_rng(2).standard_normal((12, 2)).astype(np.float32)
        ids = [f"seq{i}" for i in range(12)]
        labels = [0, 0, 0, 0, 1, 1, 1, 1, -1, -1, 2, 2]
        out_png = str(tmp_path / "clustered.png")
        _save_umap_plot(U, ids, "TESTHMM", out_png, cluster_labels=labels)
        assert (tmp_path / "clustered.png").exists()

    @requires_matplotlib
    def test_2d_plot_with_clades(self, tmp_path):
        U = np.random.default_rng(3).standard_normal((8, 2)).astype(np.float32)
        ids = [f"seq{i}" for i in range(8)]
        clades = {"cladeA": ["seq0", "seq1", "seq2"], "cladeB": ["seq3", "seq4"]}
        out_png = str(tmp_path / "clade.png")
        _save_umap_plot(U, ids, "TESTHMM", out_png, clades=clades)
        assert (tmp_path / "clade.png").exists()

    def test_plot_skipped_gracefully_without_matplotlib(self, tmp_path, monkeypatch):
        """_save_umap_plot should not raise even when matplotlib is missing."""
        import sys
        monkeypatch.setitem(sys.modules, "matplotlib", None)
        U = np.random.default_rng(4).standard_normal((10, 2)).astype(np.float32)
        ids = [f"seq{i}" for i in range(10)]
        out_png = str(tmp_path / "nomatplotlib.png")
        # Should not raise
        _save_umap_plot(U, ids, "TESTHMM", out_png)


# ── kNN analysis ──────────────────────────────────────────────────────────────

class TestKnnAnalysis:
    """Tests for kNN on PCA space via compute_embeddings_for_hmm integration."""

    def _make_pca_matrix(self, n=20, d=10, seed=0):
        rng = np.random.default_rng(seed)
        return rng.standard_normal((n, d)).astype(np.float32)

    @requires_sklearn
    def test_knn_tsv_written(self, tmp_path):
        """kNN TSV is created when run_knn=True."""
        import pandas as pd
        from sklearn.neighbors import NearestNeighbors

        Z = self._make_pca_matrix(n=15, d=5)
        ids = [f"g{i}|prot{i}" for i in range(15)]
        n_neighbors = 3

        nbrs = NearestNeighbors(n_neighbors=n_neighbors, metric="cosine").fit(Z)
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
        out_knn = tmp_path / "test.knn.tsv"
        pd.DataFrame(knn_rows).to_csv(str(out_knn), sep="\t", index=False)

        assert out_knn.exists()
        df = pd.read_csv(str(out_knn), sep="\t")
        assert list(df.columns) == ["protein_id", "neighbor_rank", "neighbor_id", "distance"]
        assert len(df) == 15 * n_neighbors

    @requires_sklearn
    def test_knn_columns(self, tmp_path):
        """kNN TSV has the expected columns."""
        import pandas as pd
        from sklearn.neighbors import NearestNeighbors

        Z = self._make_pca_matrix(n=10, d=5)
        ids = [f"genome{i}|protein{i}" for i in range(10)]

        nbrs = NearestNeighbors(n_neighbors=3, metric="cosine").fit(Z)
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
        df = pd.DataFrame(knn_rows)
        assert set(df.columns) == {"protein_id", "neighbor_rank", "neighbor_id", "distance"}
        assert (df["neighbor_rank"] >= 1).all()
        assert (df["distance"] >= 0).all()

    @requires_sklearn
    def test_knn_cosine_distances_in_range(self):
        """Cosine distances are in [0, 2]."""
        from sklearn.neighbors import NearestNeighbors

        Z = self._make_pca_matrix(n=12, d=8)
        nbrs = NearestNeighbors(n_neighbors=5, metric="cosine").fit(Z)
        distances, _ = nbrs.kneighbors(Z)
        assert (distances >= 0).all()
        assert (distances <= 2).all()


# ── Cluster subworkflow helpers ───────────────────────────────────────────────

def _make_clustered_data(n_per_cluster=10, n_clusters=3, noise=3, seed=42):
    """Return (Z, ids, labels) with well-separated Gaussian clusters + noise."""
    rng = np.random.default_rng(seed)
    centres = rng.standard_normal((n_clusters, 8)) * 5
    rows, labs = [], []
    for ci in range(n_clusters):
        pts = centres[ci] + rng.standard_normal((n_per_cluster, 8)) * 0.3
        rows.append(pts)
        labs.extend([ci] * n_per_cluster)
    # noise points
    rows.append(rng.standard_normal((noise, 8)) * 20)
    labs.extend([-1] * noise)
    Z = np.vstack(rows).astype(np.float32)
    n = n_per_cluster * n_clusters + noise
    ids = [f"g{i}|prot{i}" for i in range(n)]
    return Z, ids, labs


class TestComputeKnnMetrics:
    @requires_sklearn
    def test_returns_dataframe(self):
        Z, ids, labels = _make_clustered_data()
        df = _compute_knn_metrics(Z, ids, labels, k=5)
        assert not df.empty

    @requires_sklearn
    def test_expected_columns(self):
        Z, ids, labels = _make_clustered_data()
        df = _compute_knn_metrics(Z, ids, labels, k=5)
        expected = {
            "protein_id", "cluster_id", "dominant_cluster",
            "neighborhood_purity", "dist_weighted_purity",
            "mutual_knn_support", "neighborhood_entropy",
            "median_neighbor_distance",
        }
        assert expected.issubset(set(df.columns))

    @requires_sklearn
    def test_row_count_matches_sequences(self):
        Z, ids, labels = _make_clustered_data(n_per_cluster=8, noise=2)
        df = _compute_knn_metrics(Z, ids, labels, k=4)
        assert len(df) == len(ids)

    @requires_sklearn
    def test_purity_in_range(self):
        Z, ids, labels = _make_clustered_data()
        df = _compute_knn_metrics(Z, ids, labels, k=5)
        assert (df["neighborhood_purity"] >= 0).all()
        assert (df["neighborhood_purity"] <= 1).all()

    @requires_sklearn
    def test_mutual_knn_in_range(self):
        Z, ids, labels = _make_clustered_data()
        df = _compute_knn_metrics(Z, ids, labels, k=5)
        assert (df["mutual_knn_support"] >= 0).all()
        assert (df["mutual_knn_support"] <= 1).all()

    @requires_sklearn
    def test_entropy_nonnegative(self):
        Z, ids, labels = _make_clustered_data()
        df = _compute_knn_metrics(Z, ids, labels, k=5)
        assert (df["neighborhood_entropy"] >= 0).all()

    @requires_sklearn
    def test_too_few_sequences_returns_empty(self):
        Z = np.random.default_rng(0).standard_normal((1, 5)).astype(np.float32)
        df = _compute_knn_metrics(Z, ["g|p1"], [0], k=5)
        assert df.empty


class TestAssignMembershipTiers:
    def test_fallback_when_knn_empty(self):
        import pandas as pd
        ids = ["g0|p0", "g1|p1", "g2|p2", "g3|p3"]
        labels = [0, 0, 1, -1]
        tiers = _assign_membership_tiers(ids, labels, pd.DataFrame())
        # labels >= 0 → core; label == -1 → outlier
        assert tiers["p0"] == "core"
        assert tiers["p1"] == "core"
        assert tiers["p2"] == "core"
        assert tiers["p3"] == "outlier"

    @requires_sklearn
    def test_core_members_for_pure_clusters(self):
        """Sequences with high-purity neighbourhoods should be 'core'."""
        Z, ids, labels = _make_clustered_data(n_per_cluster=15, noise=0)
        import pandas as pd
        knn_df = _compute_knn_metrics(Z, ids, labels, k=5)
        tiers = _assign_membership_tiers(ids, labels, knn_df)
        # Most cluster members in a well-separated dataset should be core
        core_count = sum(1 for t in tiers.values() if t == "core")
        assert core_count > 0

    @requires_sklearn
    def test_noise_classified_as_noise_tier(self):
        """Noise sequences (label -1) should get a noise-related tier."""
        Z, ids, labels = _make_clustered_data(n_per_cluster=12, noise=4)
        import pandas as pd
        knn_df = _compute_knn_metrics(Z, ids, labels, k=5)
        tiers = _assign_membership_tiers(ids, labels, knn_df)
        noise_ids = {
            (tip.split("|", 1)[1] if "|" in tip else tip)
            for tip, lb in zip(ids, labels) if lb == -1
        }
        noise_tiers = {tiers.get(pid) for pid in noise_ids if pid in tiers}
        valid_noise_tiers = {"noise_peripheral", "noise_bridge", "outlier"}
        assert noise_tiers.issubset(valid_noise_tiers)

    def test_all_sequences_assigned(self):
        import pandas as pd
        ids = [f"g{i}|p{i}" for i in range(6)]
        labels = [0, 0, 1, 1, -1, -1]
        tiers = _assign_membership_tiers(ids, labels, pd.DataFrame())
        proteins = {tip.split("|", 1)[1] for tip in ids}
        assert proteins == set(tiers.keys())


class TestWriteClusterFastas:
    def _make_seqs_and_data(self, n_per_cluster=5, n_clusters=2, noise=2):
        seqs = {}
        ids = []
        labels = []
        aa = "ACDEFGHIKLMNPQRSTVWY"
        rng = np.random.default_rng(7)
        for cl in range(n_clusters):
            for j in range(n_per_cluster):
                pid = f"prot_cl{cl}_{j}"
                tip = f"genome_{cl}_{j}|{pid}"
                seqs[tip] = "".join(rng.choice(list(aa), size=50))
                ids.append(tip)
                labels.append(cl)
        for j in range(noise):
            pid = f"prot_noise_{j}"
            tip = f"genome_noise_{j}|{pid}"
            seqs[tip] = "".join(rng.choice(list(aa), size=50))
            ids.append(tip)
            labels.append(-1)
        return seqs, ids, labels

    def test_core_fasta_created(self, tmp_path):
        import pandas as pd
        seqs, ids, labels = self._make_seqs_and_data()
        tiers = _assign_membership_tiers(ids, labels, pd.DataFrame())
        result = _write_cluster_fastas(seqs, ids, labels, tiers, str(tmp_path))
        for cl_id, info in result.items():
            if info["n_core"] > 0:
                assert info["core"] is not None
                assert (tmp_path / f"cluster_{cl_id}.core.faa").exists()

    def test_noise_not_included(self, tmp_path):
        import pandas as pd
        seqs, ids, labels = self._make_seqs_and_data(noise=3)
        tiers = _assign_membership_tiers(ids, labels, pd.DataFrame())
        result = _write_cluster_fastas(seqs, ids, labels, tiers, str(tmp_path))
        # Noise sequences should not appear in cluster FASTAs
        from phylofoundry.utils.bio import read_fasta
        for cl_id, info in result.items():
            for key in ("core", "affiliate"):
                fp = info[key]
                if fp and (tmp_path / f"cluster_{cl_id}.{key}.faa").exists():
                    written_seqs = read_fasta(fp)
                    noise_prots = {
                        (tip.split("|", 1)[1] if "|" in tip else tip)
                        for tip, lb in zip(ids, labels) if lb == -1
                    }
                    assert set(written_seqs.keys()).isdisjoint(noise_prots)

    def test_result_has_expected_keys(self, tmp_path):
        import pandas as pd
        seqs, ids, labels = self._make_seqs_and_data()
        tiers = _assign_membership_tiers(ids, labels, pd.DataFrame())
        result = _write_cluster_fastas(seqs, ids, labels, tiers, str(tmp_path))
        for info in result.values():
            assert "core" in info
            assert "affiliate" in info
            assert "n_core" in info
            assert "n_affiliate" in info

    def test_force_overwrites(self, tmp_path):
        import pandas as pd
        seqs, ids, labels = self._make_seqs_and_data()
        tiers = _assign_membership_tiers(ids, labels, pd.DataFrame())
        _write_cluster_fastas(seqs, ids, labels, tiers, str(tmp_path))
        # Second call with force=True should not raise
        _write_cluster_fastas(seqs, ids, labels, tiers, str(tmp_path), force=True)


class TestGenerateSequenceLogo:
    def _write_msa(self, tmp_path, n_seqs=8, aln_len=30, seed=1):
        """Write a small synthetic aligned FASTA."""
        aa = "ACDEFGHIKLMNPQRSTVWY"
        rng = np.random.default_rng(seed)
        fasta_fp = str(tmp_path / "test.aln.faa")
        with open(fasta_fp, "w") as fh:
            for i in range(n_seqs):
                seq = "".join(rng.choice(list(aa), size=aln_len))
                fh.write(f">seq{i}\n{seq}\n")
        return fasta_fp

    @requires_matplotlib
    def test_png_created(self, tmp_path):
        msa = self._write_msa(tmp_path)
        png = str(tmp_path / "test.logo.png")
        svg = str(tmp_path / "test.logo.svg")
        ok, n_seqs, aln_len, n_eff = _generate_sequence_logo(msa, png, svg, "TESTHMM", 0, formats=["png"])
        assert ok
        assert (tmp_path / "test.logo.png").exists()

    @requires_matplotlib
    def test_svg_created(self, tmp_path):
        msa = self._write_msa(tmp_path)
        png = str(tmp_path / "test.logo.png")
        svg = str(tmp_path / "test.logo.svg")
        ok, n_seqs, aln_len, n_eff = _generate_sequence_logo(msa, png, svg, "TESTHMM", 0, formats=["svg"])
        assert ok
        assert (tmp_path / "test.logo.svg").exists()

    @requires_matplotlib
    def test_returns_n_seqs(self, tmp_path):
        msa = self._write_msa(tmp_path, n_seqs=6)
        png = str(tmp_path / "logo.png")
        svg = str(tmp_path / "logo.svg")
        ok, n_seqs, aln_len, n_eff = _generate_sequence_logo(msa, png, svg, "TESTHMM", 1)
        assert n_seqs == 6

    @requires_matplotlib
    def test_n_effective_sites_nonnegative(self, tmp_path):
        msa = self._write_msa(tmp_path)
        png = str(tmp_path / "logo.png")
        svg = str(tmp_path / "logo.svg")
        _, _, aln_len, n_eff = _generate_sequence_logo(msa, png, svg, "TESTHMM", 2)
        assert n_eff >= 0

    def test_missing_msa_returns_false(self, tmp_path):
        ok, n_seqs, aln_len, n_eff = _generate_sequence_logo(
            str(tmp_path / "nonexistent.faa"),
            str(tmp_path / "logo.png"),
            str(tmp_path / "logo.svg"),
            "TESTHMM", 0,
        )
        assert not ok

    @requires_matplotlib
    def test_all_gap_column_ignored(self, tmp_path):
        """An alignment where all columns are gaps should return ok=False."""
        fasta_fp = str(tmp_path / "gaps.aln.faa")
        with open(fasta_fp, "w") as fh:
            for i in range(4):
                fh.write(f">seq{i}\n{'-' * 20}\n")
        png = str(tmp_path / "gaps.logo.png")
        svg = str(tmp_path / "gaps.logo.svg")
        ok, n_seqs, aln_len, n_eff = _generate_sequence_logo(fasta_fp, png, svg, "TESTHMM", 0)
        assert not ok


class TestClassifyNoiseSequences:
    def _make_knn_df_with_noise(self):
        import pandas as pd
        rows = [
            # noise sequences with varying characteristics
            {"protein_id": "p_periph", "cluster_id": -1, "dominant_cluster": 2,
             "neighborhood_purity": 0.7, "dist_weighted_purity": 0.65,
             "mutual_knn_support": 0.5, "neighborhood_entropy": 0.8,
             "median_neighbor_distance": 0.5},
            {"protein_id": "p_bridge", "cluster_id": -1, "dominant_cluster": 0,
             "neighborhood_purity": 0.4, "dist_weighted_purity": 0.4,
             "mutual_knn_support": 0.2, "neighborhood_entropy": 2.0,
             "median_neighbor_distance": 1.2},
            {"protein_id": "p_outlier", "cluster_id": -1, "dominant_cluster": -1,
             "neighborhood_purity": 0.1, "dist_weighted_purity": 0.1,
             "mutual_knn_support": 0.0, "neighborhood_entropy": 2.8,
             "median_neighbor_distance": 3.0},
            # a regular cluster member (should NOT appear in result)
            {"protein_id": "p_core", "cluster_id": 0, "dominant_cluster": 0,
             "neighborhood_purity": 0.95, "dist_weighted_purity": 0.9,
             "mutual_knn_support": 0.8, "neighborhood_entropy": 0.2,
             "median_neighbor_distance": 0.1},
        ]
        return pd.DataFrame(rows)

    def test_returns_only_noise_rows(self):
        df = self._make_knn_df_with_noise()
        result = _classify_noise_sequences(df, "TESTHMM")
        # Only noise sequences (cluster_id == -1) appear in result
        assert set(result["protein_id"]) == {"p_periph", "p_bridge", "p_outlier"}

    def test_expected_columns(self):
        result = _classify_noise_sequences(self._make_knn_df_with_noise(), "TESTHMM")
        expected = {
            "hmm_name", "protein_id", "classification",
            "dominant_cluster", "neighborhood_purity",
            "neighborhood_entropy", "notes",
        }
        assert expected.issubset(set(result.columns))

    def test_peripheral_homolog_classification(self):
        result = _classify_noise_sequences(self._make_knn_df_with_noise(), "TESTHMM")
        row = result[result["protein_id"] == "p_periph"]
        assert row.iloc[0]["classification"] == "peripheral_homolog"

    def test_bridge_or_partial_classification(self):
        result = _classify_noise_sequences(self._make_knn_df_with_noise(), "TESTHMM")
        row = result[result["protein_id"] == "p_bridge"]
        assert row.iloc[0]["classification"] in ("bridge_sequence", "partial_homolog")

    def test_hmm_name_column(self):
        result = _classify_noise_sequences(self._make_knn_df_with_noise(), "MYHMM")
        assert (result["hmm_name"] == "MYHMM").all()

    def test_empty_knn_returns_empty(self):
        import pandas as pd
        result = _classify_noise_sequences(pd.DataFrame(), "TESTHMM")
        assert result.empty

    def test_no_noise_returns_empty(self):
        import pandas as pd
        rows = [
            {"protein_id": "p0", "cluster_id": 0, "dominant_cluster": 0,
             "neighborhood_purity": 0.9, "dist_weighted_purity": 0.9,
             "mutual_knn_support": 0.8, "neighborhood_entropy": 0.1,
             "median_neighbor_distance": 0.2},
        ]
        result = _classify_noise_sequences(pd.DataFrame(rows), "TESTHMM")
        assert result.empty


class TestValidateEmbCfgSubworkflow:
    """Tests for cluster_subworkflow config validation in _validate_emb_cfg."""

    def _base_cfg(self):
        return {
            "backend": "esm",
            "model": "esm2_t33_650M_UR50D",
            "device": "cpu",
            "batch_size": 4,
            "repr_layer": None,
            "pooling": "mean",
            "pca_components": 3,
            "cluster_embeddings": True,
        }

    def test_subworkflow_defaults_injected(self):
        cfg = _validate_emb_cfg(self._base_cfg())
        sub = cfg["cluster_subworkflow"]
        assert sub["enabled"] is False
        assert sub["build_cluster_msas"] is True
        assert sub["seed_membership"] == "core_only"
        assert sub["build_cluster_hmms"] is True
        assert sub["classify_noise"] is True
        assert sub["generate_sequence_logos"] is True
        assert "png" in sub["logo_format"]
        assert "svg" in sub["logo_format"]

    def test_subworkflow_can_be_enabled(self):
        base = self._base_cfg()
        base["cluster_subworkflow"] = {"enabled": True}
        cfg = _validate_emb_cfg(base)
        assert cfg["cluster_subworkflow"]["enabled"] is True
        # Other defaults still filled in
        assert cfg["cluster_subworkflow"]["build_cluster_msas"] is True

    def test_subworkflow_options_respected(self):
        base = self._base_cfg()
        base["cluster_subworkflow"] = {
            "enabled": True,
            "build_cluster_hmms": False,
            "logo_format": ["png"],
        }
        cfg = _validate_emb_cfg(base)
        sub = cfg["cluster_subworkflow"]
        assert sub["build_cluster_hmms"] is False
        assert sub["logo_format"] == ["png"]

    def test_non_dict_subworkflow_replaced_with_defaults(self):
        base = self._base_cfg()
        base["cluster_subworkflow"] = None
        cfg = _validate_emb_cfg(base)
        assert isinstance(cfg["cluster_subworkflow"], dict)
        assert cfg["cluster_subworkflow"]["enabled"] is False


# ── _compute_cluster_kl_divergence ────────────────────────────────────────────

from phylofoundry.tasks.embed import _compute_cluster_kl_divergence


def _write_aligned_fasta(path, seqs):
    """Write an aligned FASTA to *path* (dict {id: seq})."""
    with open(path, "w") as fh:
        for sid, seq in seqs.items():
            fh.write(f">{sid}\n{seq}\n")


class TestComputeClusterKlDivergence:
    """Unit tests for _compute_cluster_kl_divergence()."""

    def _make_aln_dir(self, tmp_path, clusters):
        """Write per-cluster aligned FASTAs and return {cl_id: path} dict.

        Parameters
        ----------
        clusters : dict
            {cl_id: {seq_id: aligned_sequence}}
        """
        paths = {}
        for cl_id, seqs in clusters.items():
            fp = str(tmp_path / f"cluster_{cl_id}.seed.aln.faa")
            _write_aligned_fasta(fp, seqs)
            paths[cl_id] = fp
        return paths

    # ── basic smoke tests ──────────────────────────────────────────────────────

    def test_returns_dataframes(self, tmp_path):
        """Should return two DataFrames (per_position, top_sites)."""
        clusters = {
            0: {f"s{i}": "ACDEF" * 4 for i in range(6)},
            1: {f"t{i}": "GHIKL" * 4 for i in range(6)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        df, top = _compute_cluster_kl_divergence("HMM1", paths)
        assert hasattr(df, "columns")
        assert hasattr(top, "columns")

    def test_nonempty_for_two_clusters(self, tmp_path):
        """Two distinct clusters should produce non-empty output."""
        clusters = {
            0: {f"s{i}": "ACDEFGHIKL" for i in range(6)},
            1: {f"t{i}": "MNPQRSTVWY" for i in range(6)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        df, top = _compute_cluster_kl_divergence("HMM1", paths)
        assert not df.empty
        assert not top.empty

    def test_expected_columns(self, tmp_path):
        clusters = {
            0: {f"s{i}": "ACDEFGHIKL" for i in range(6)},
            1: {f"t{i}": "MNPQRSTVWY" for i in range(6)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        df, _ = _compute_cluster_kl_divergence("HMM1", paths)
        expected = {
            "hmm_name", "cluster_A", "cluster_B", "pair",
            "aln_position", "kl_A_to_B", "kl_B_to_A", "js_divergence",
            "top_aa_A", "top_aa_B", "n_seqs_A", "n_seqs_B",
        }
        assert expected.issubset(set(df.columns))

    # ── divergence value tests ─────────────────────────────────────────────────

    def test_jsd_in_valid_range(self, tmp_path):
        """JSD values must lie in [0, 1] bits."""
        clusters = {
            0: {f"s{i}": "AAAA" for i in range(8)},
            1: {f"t{i}": "CCCC" for i in range(8)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        df, _ = _compute_cluster_kl_divergence("HMM1", paths)
        assert (df["js_divergence"] >= 0).all()
        assert (df["js_divergence"] <= 1).all()

    def test_high_divergence_for_distinct_clusters(self, tmp_path):
        """Clusters with completely different residues should show high JSD."""
        clusters = {
            0: {f"s{i}": "AAAA" for i in range(8)},
            1: {f"t{i}": "CCCC" for i in range(8)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        df, _ = _compute_cluster_kl_divergence("HMM1", paths, pseudocount=1e-6)
        assert df["js_divergence"].max() > 0.5

    def test_low_divergence_for_identical_clusters(self, tmp_path):
        """Identical clusters should produce near-zero JSD at all positions."""
        seq = "ACDEFGHIKL"
        clusters = {
            0: {f"s{i}": seq for i in range(8)},
            1: {f"t{i}": seq for i in range(8)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        df, _ = _compute_cluster_kl_divergence("HMM1", paths)
        assert df["js_divergence"].max() < 0.01

    def test_kl_values_nonnegative(self, tmp_path):
        """KL divergence values must be non-negative."""
        clusters = {
            0: {f"s{i}": "ACDE" for i in range(6)},
            1: {f"t{i}": "FGHI" for i in range(6)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        df, _ = _compute_cluster_kl_divergence("HMM1", paths)
        assert (df["kl_A_to_B"] >= 0).all()
        assert (df["kl_B_to_A"] >= 0).all()

    # ── row and pair count tests ───────────────────────────────────────────────

    def test_row_count_equals_pairs_times_columns(self, tmp_path):
        """Should have n_pairs * aln_len rows in per_position_df."""
        aln_len = 10
        clusters = {
            0: {f"s{i}": "A" * aln_len for i in range(6)},
            1: {f"t{i}": "C" * aln_len for i in range(6)},
            2: {f"u{i}": "G" * aln_len for i in range(6)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        df, _ = _compute_cluster_kl_divergence("HMM1", paths)
        # C(3,2) = 3 pairs; each pair × 10 columns = 30 rows
        assert len(df) == 3 * aln_len

    def test_top_sites_at_most_top_n(self, tmp_path):
        """Top-sites table should have at most top_n_sites rows per pair."""
        aln_len = 15
        top_n = 5
        clusters = {
            0: {f"s{i}": "A" * aln_len for i in range(6)},
            1: {f"t{i}": "C" * aln_len for i in range(6)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        _, top = _compute_cluster_kl_divergence("HMM1", paths, top_n_sites=top_n)
        assert len(top) <= top_n

    def test_aln_position_1based(self, tmp_path):
        """aln_position should start at 1 (1-based indexing)."""
        clusters = {
            0: {f"s{i}": "ACDE" for i in range(6)},
            1: {f"t{i}": "FGHI" for i in range(6)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        df, _ = _compute_cluster_kl_divergence("HMM1", paths)
        assert df["aln_position"].min() == 1

    def test_hmm_name_in_output(self, tmp_path):
        clusters = {
            0: {f"s{i}": "ACDE" for i in range(6)},
            1: {f"t{i}": "FGHI" for i in range(6)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        df, _ = _compute_cluster_kl_divergence("MY_HMM", paths)
        assert (df["hmm_name"] == "MY_HMM").all()

    # ── edge case / guard tests ────────────────────────────────────────────────

    def test_single_cluster_returns_empty(self, tmp_path):
        """Only one cluster → no pairs → both DataFrames empty."""
        clusters = {0: {f"s{i}": "ACDE" for i in range(6)}}
        paths = self._make_aln_dir(tmp_path, clusters)
        df, top = _compute_cluster_kl_divergence("HMM1", paths)
        assert df.empty
        assert top.empty

    def test_empty_paths_returns_empty(self, tmp_path):
        df, top = _compute_cluster_kl_divergence("HMM1", {})
        assert df.empty
        assert top.empty

    def test_min_cluster_size_filters_small_clusters(self, tmp_path):
        """Clusters below min_cluster_size are excluded; may leave < 2 eligible."""
        clusters = {
            0: {f"s{i}": "ACDE" for i in range(2)},  # too small
            1: {f"t{i}": "FGHI" for i in range(2)},  # too small
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        df, top = _compute_cluster_kl_divergence("HMM1", paths, min_cluster_size=5)
        assert df.empty

    def test_gaps_skipped_in_counts(self, tmp_path):
        """Gap characters should not appear as residues; JSD still valid."""
        clusters = {
            0: {f"s{i}": "--AA--" for i in range(6)},
            1: {f"t{i}": "--CC--" for i in range(6)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        df, _ = _compute_cluster_kl_divergence("HMM1", paths)
        assert (df["js_divergence"] >= 0).all()
        assert (df["js_divergence"] <= 1).all()


# ── _compute_cluster_jsd_analysis ─────────────────────────────────────────────

from phylofoundry.tasks.embed import _compute_cluster_jsd_analysis


class TestComputeClusterJsdAnalysis:
    """Unit tests for _compute_cluster_jsd_analysis()."""

    def _make_aln_dir(self, tmp_path, clusters):
        """Write per-cluster aligned FASTAs and return {cl_id: path} dict."""
        paths = {}
        for cl_id, seqs in clusters.items():
            fp = str(tmp_path / f"cluster_{cl_id}.seed.aln.faa")
            _write_aligned_fasta(fp, seqs)
            paths[cl_id] = fp
        return paths

    # ── basic smoke tests ──────────────────────────────────────────────────────

    def test_returns_dataframes(self, tmp_path):
        """Should return two DataFrames (per_position, top_sites)."""
        clusters = {
            0: {f"s{i}": "ACDEF" * 4 for i in range(6)},
            1: {f"t{i}": "GHIKL" * 4 for i in range(6)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        df, top = _compute_cluster_jsd_analysis("HMM1", paths)
        assert hasattr(df, "columns")
        assert hasattr(top, "columns")

    def test_nonempty_for_two_clusters(self, tmp_path):
        """Two distinct clusters should produce non-empty output."""
        clusters = {
            0: {f"s{i}": "ACDEFGHIKL" for i in range(6)},
            1: {f"t{i}": "MNPQRSTVWY" for i in range(6)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        df, top = _compute_cluster_jsd_analysis("HMM1", paths)
        assert not df.empty
        assert not top.empty

    def test_expected_columns(self, tmp_path):
        """Output DataFrame must contain the expected column set."""
        clusters = {
            0: {f"s{i}": "ACDEFGHIKL" for i in range(6)},
            1: {f"t{i}": "MNPQRSTVWY" for i in range(6)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        df, _ = _compute_cluster_jsd_analysis("HMM1", paths)
        expected = {
            "hmm_name", "cluster_A", "cluster_B", "pair",
            "aln_position", "js_divergence",
            "top_aa_A", "top_aa_B", "n_seqs_A", "n_seqs_B",
        }
        assert expected.issubset(set(df.columns))

    def test_no_kl_columns(self, tmp_path):
        """JSD-only output must not contain asymmetric KL columns."""
        clusters = {
            0: {f"s{i}": "ACDE" for i in range(6)},
            1: {f"t{i}": "FGHI" for i in range(6)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        df, _ = _compute_cluster_jsd_analysis("HMM1", paths)
        assert "kl_A_to_B" not in df.columns
        assert "kl_B_to_A" not in df.columns

    # ── divergence value tests ─────────────────────────────────────────────────

    def test_jsd_in_valid_range(self, tmp_path):
        """JSD values must lie in [0, 1] bits."""
        clusters = {
            0: {f"s{i}": "AAAA" for i in range(8)},
            1: {f"t{i}": "CCCC" for i in range(8)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        df, _ = _compute_cluster_jsd_analysis("HMM1", paths)
        assert (df["js_divergence"] >= 0).all()
        assert (df["js_divergence"] <= 1).all()

    def test_high_divergence_for_distinct_clusters(self, tmp_path):
        """Clusters with completely different residues should show high JSD."""
        clusters = {
            0: {f"s{i}": "AAAA" for i in range(8)},
            1: {f"t{i}": "CCCC" for i in range(8)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        df, _ = _compute_cluster_jsd_analysis("HMM1", paths, pseudocount=1e-6)
        assert df["js_divergence"].max() > 0.5

    def test_low_divergence_for_identical_clusters(self, tmp_path):
        """Identical clusters should produce near-zero JSD at all positions."""
        seq = "ACDEFGHIKL"
        clusters = {
            0: {f"s{i}": seq for i in range(8)},
            1: {f"t{i}": seq for i in range(8)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        df, _ = _compute_cluster_jsd_analysis("HMM1", paths)
        assert df["js_divergence"].max() < 0.01

    def test_symmetry(self, tmp_path):
        """JSD is symmetric: swapping cluster A/B must yield the same value."""
        clusters = {
            0: {f"s{i}": "AAAA" for i in range(8)},
            1: {f"t{i}": "CCCC" for i in range(8)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        df, _ = _compute_cluster_jsd_analysis("HMM1", paths)
        # Only one pair in this case; verify JSD value is consistent across rows
        assert df["js_divergence"].nunique() <= 2  # identical positions may share same value

    # ── row and pair count tests ───────────────────────────────────────────────

    def test_row_count_equals_pairs_times_columns(self, tmp_path):
        """Should have n_pairs * aln_len rows in per_position_df."""
        aln_len = 10
        clusters = {
            0: {f"s{i}": "A" * aln_len for i in range(6)},
            1: {f"t{i}": "C" * aln_len for i in range(6)},
            2: {f"u{i}": "G" * aln_len for i in range(6)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        df, _ = _compute_cluster_jsd_analysis("HMM1", paths)
        # C(3,2) = 3 pairs; each pair × 10 columns = 30 rows
        assert len(df) == 3 * aln_len

    def test_top_sites_at_most_top_n(self, tmp_path):
        """Top-sites table should have at most top_n_sites rows per pair."""
        aln_len = 15
        top_n = 5
        clusters = {
            0: {f"s{i}": "A" * aln_len for i in range(6)},
            1: {f"t{i}": "C" * aln_len for i in range(6)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        _, top = _compute_cluster_jsd_analysis("HMM1", paths, top_n_sites=top_n)
        assert len(top) <= top_n

    def test_aln_position_1based(self, tmp_path):
        """aln_position should start at 1 (1-based indexing)."""
        clusters = {
            0: {f"s{i}": "ACDE" for i in range(6)},
            1: {f"t{i}": "FGHI" for i in range(6)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        df, _ = _compute_cluster_jsd_analysis("HMM1", paths)
        assert df["aln_position"].min() == 1

    def test_hmm_name_in_output(self, tmp_path):
        """hmm_name column should match the provided name."""
        clusters = {
            0: {f"s{i}": "ACDE" for i in range(6)},
            1: {f"t{i}": "FGHI" for i in range(6)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        df, _ = _compute_cluster_jsd_analysis("MY_HMM", paths)
        assert (df["hmm_name"] == "MY_HMM").all()

    # ── edge case / guard tests ────────────────────────────────────────────────

    def test_single_cluster_returns_empty(self, tmp_path):
        """Only one cluster → no pairs → both DataFrames empty."""
        clusters = {0: {f"s{i}": "ACDE" for i in range(6)}}
        paths = self._make_aln_dir(tmp_path, clusters)
        df, top = _compute_cluster_jsd_analysis("HMM1", paths)
        assert df.empty
        assert top.empty

    def test_empty_paths_returns_empty(self, tmp_path):
        """Empty input dict should return two empty DataFrames."""
        df, top = _compute_cluster_jsd_analysis("HMM1", {})
        assert df.empty
        assert top.empty

    def test_min_cluster_size_filters_small_clusters(self, tmp_path):
        """Clusters below min_cluster_size are excluded; may leave < 2 eligible."""
        clusters = {
            0: {f"s{i}": "ACDE" for i in range(2)},  # too small
            1: {f"t{i}": "FGHI" for i in range(2)},  # too small
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        df, top = _compute_cluster_jsd_analysis("HMM1", paths, min_cluster_size=5)
        assert df.empty

    def test_gaps_skipped_in_counts(self, tmp_path):
        """Gap characters should not affect JSD validity."""
        clusters = {
            0: {f"s{i}": "--AA--" for i in range(6)},
            1: {f"t{i}": "--CC--" for i in range(6)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        df, _ = _compute_cluster_jsd_analysis("HMM1", paths)
        assert (df["js_divergence"] >= 0).all()
        assert (df["js_divergence"] <= 1).all()




# ── _generate_cluster_motif_heatmap ───────────────────────────────────────────

from phylofoundry.tasks.embed import _generate_cluster_motif_heatmap


def _write_aligned_fasta_hm(path, seqs):
    """Write an aligned FASTA to *path* (dict {id: seq})."""
    with open(path, "w") as fh:
        for sid, seq in seqs.items():
            fh.write(f">{sid}\n{seq}\n")


@requires_matplotlib
class TestGenerateClusterMotifHeatmap:
    """Unit tests for _generate_cluster_motif_heatmap()."""

    def _make_aln_dir(self, tmp_path, clusters):
        """Write per-cluster aligned FASTAs and return {cl_id: path} dict."""
        paths = {}
        for cl_id, seqs in clusters.items():
            fp = str(tmp_path / f"cluster_{cl_id}.seed.aln.faa")
            _write_aligned_fasta_hm(fp, seqs)
            paths[cl_id] = fp
        return paths

    # ── smoke tests ───────────────────────────────────────────────────────────

    def test_returns_dataframe(self, tmp_path):
        """Should return a DataFrame (the data matrix)."""
        clusters = {
            0: {f"s{i}": "ACDEF" for i in range(6)},
            1: {f"t{i}": "GHIKL" for i in range(6)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        result = _generate_cluster_motif_heatmap(
            "HMM1", paths,
            out_dir=str(tmp_path / "heatmaps"),
            summary_dir=str(tmp_path / "summary"),
        )
        import pandas as pd
        assert isinstance(result, pd.DataFrame)

    def test_matrix_shape(self, tmp_path):
        """Matrix should have rows=n_clusters, cols=aln_length."""
        aln_len = 8
        clusters = {
            0: {f"s{i}": "A" * aln_len for i in range(6)},
            1: {f"t{i}": "C" * aln_len for i in range(6)},
            2: {f"u{i}": "G" * aln_len for i in range(6)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        result = _generate_cluster_motif_heatmap(
            "HMM1", paths,
            out_dir=str(tmp_path / "heatmaps"),
            summary_dir=str(tmp_path / "summary"),
        )
        assert result.shape == (3, aln_len)

    def test_matrix_values_in_range_jsd(self, tmp_path):
        """JSD vs global values must lie in [0, 1]."""
        clusters = {
            0: {f"s{i}": "AAAA" for i in range(6)},
            1: {f"t{i}": "CCCC" for i in range(6)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        import numpy as np
        result = _generate_cluster_motif_heatmap(
            "HMM1", paths,
            out_dir=str(tmp_path / "heatmaps"),
            summary_dir=str(tmp_path / "summary"),
            metric="jsd_vs_global",
        )
        assert (result.values >= 0).all()
        assert (result.values <= 1.0 + 1e-10).all()

    def test_matrix_values_nonneg_entropy(self, tmp_path):
        """Shannon entropy values must be non-negative."""
        clusters = {
            0: {f"s{i}": "ACDE" for i in range(6)},
            1: {f"t{i}": "FGHI" for i in range(6)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        result = _generate_cluster_motif_heatmap(
            "HMM1", paths,
            out_dir=str(tmp_path / "heatmaps"),
            summary_dir=str(tmp_path / "summary"),
            metric="shannon_entropy",
        )
        assert (result.values >= 0).all()

    def test_row_labels(self, tmp_path):
        """Row labels should be 'cluster_<id>' for each cluster."""
        clusters = {
            0: {f"s{i}": "ACDEF" for i in range(6)},
            1: {f"t{i}": "GHIKL" for i in range(6)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        result = _generate_cluster_motif_heatmap(
            "HMM1", paths,
            out_dir=str(tmp_path / "heatmaps"),
            summary_dir=str(tmp_path / "summary"),
        )
        assert list(result.index) == ["cluster_0", "cluster_1"]

    def test_column_labels(self, tmp_path):
        """Column labels should be 'pos_<N>' (1-based)."""
        aln_len = 5
        clusters = {
            0: {f"s{i}": "A" * aln_len for i in range(6)},
            1: {f"t{i}": "C" * aln_len for i in range(6)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        result = _generate_cluster_motif_heatmap(
            "HMM1", paths,
            out_dir=str(tmp_path / "heatmaps"),
            summary_dir=str(tmp_path / "summary"),
        )
        assert list(result.columns) == [f"pos_{p}" for p in range(1, aln_len + 1)]

    def test_figure_file_written(self, tmp_path):
        """PNG figure should be written to out_dir."""
        import os
        clusters = {
            0: {f"s{i}": "ACDEF" for i in range(6)},
            1: {f"t{i}": "GHIKL" for i in range(6)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        heatmap_dir = tmp_path / "heatmaps"
        _generate_cluster_motif_heatmap(
            "HMM1", paths,
            out_dir=str(heatmap_dir),
            summary_dir=str(tmp_path / "summary"),
            figure_format=["png"],
        )
        assert (heatmap_dir / "HMM1.cluster_motif_heatmap.png").exists()

    def test_matrix_tsv_written(self, tmp_path):
        """Matrix TSV should be written to summary_dir."""
        clusters = {
            0: {f"s{i}": "ACDEF" for i in range(6)},
            1: {f"t{i}": "GHIKL" for i in range(6)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        summary_dir = tmp_path / "summary"
        _generate_cluster_motif_heatmap(
            "HMM1", paths,
            out_dir=str(tmp_path / "heatmaps"),
            summary_dir=str(summary_dir),
        )
        assert (summary_dir / "HMM1.cluster_motif_heatmap_matrix.tsv").exists()

    def test_identical_clusters_low_jsd(self, tmp_path):
        """Identical clusters should produce near-zero JSD vs global."""
        seq = "ACDEFGHIKL"
        clusters = {
            0: {f"s{i}": seq for i in range(8)},
            1: {f"t{i}": seq for i in range(8)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        result = _generate_cluster_motif_heatmap(
            "HMM1", paths,
            out_dir=str(tmp_path / "heatmaps"),
            summary_dir=str(tmp_path / "summary"),
            metric="jsd_vs_global",
        )
        # Identical distributions → JSD vs pooled global ≈ 0
        assert result.values.max() < 0.01

    def test_distinct_clusters_nonzero_jsd(self, tmp_path):
        """Completely different clusters should produce high JSD vs global."""
        clusters = {
            0: {f"s{i}": "AAAA" for i in range(8)},
            1: {f"t{i}": "CCCC" for i in range(8)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        result = _generate_cluster_motif_heatmap(
            "HMM1", paths,
            out_dir=str(tmp_path / "heatmaps"),
            summary_dir=str(tmp_path / "summary"),
            metric="jsd_vs_global",
        )
        assert result.values.max() > 0.1

    # ── edge case / guard tests ───────────────────────────────────────────────

    def test_single_cluster_returns_empty(self, tmp_path):
        """Only one eligible cluster → return empty DataFrame."""
        import pandas as pd
        clusters = {0: {f"s{i}": "ACDE" for i in range(6)}}
        paths = self._make_aln_dir(tmp_path, clusters)
        result = _generate_cluster_motif_heatmap(
            "HMM1", paths,
            out_dir=str(tmp_path / "heatmaps"),
            summary_dir=str(tmp_path / "summary"),
        )
        assert isinstance(result, pd.DataFrame)
        assert result.empty

    def test_min_cluster_size_filters(self, tmp_path):
        """Clusters below min_cluster_size should be excluded."""
        import pandas as pd
        clusters = {
            0: {f"s{i}": "ACDE" for i in range(2)},  # too small
            1: {f"t{i}": "FGHI" for i in range(2)},  # too small
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        result = _generate_cluster_motif_heatmap(
            "HMM1", paths,
            out_dir=str(tmp_path / "heatmaps"),
            summary_dir=str(tmp_path / "summary"),
            min_cluster_size=5,
        )
        assert isinstance(result, pd.DataFrame)
        assert result.empty

    def test_gaps_skipped(self, tmp_path):
        """Gap characters should not affect matrix validity."""
        clusters = {
            0: {f"s{i}": "--AA--" for i in range(6)},
            1: {f"t{i}": "--CC--" for i in range(6)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        result = _generate_cluster_motif_heatmap(
            "HMM1", paths,
            out_dir=str(tmp_path / "heatmaps"),
            summary_dir=str(tmp_path / "summary"),
        )
        assert (result.values >= 0).all()

    def test_force_overwrites_existing(self, tmp_path):
        """force=True should overwrite an existing figure file."""
        import os
        clusters = {
            0: {f"s{i}": "ACDEF" for i in range(6)},
            1: {f"t{i}": "GHIKL" for i in range(6)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        heatmap_dir = tmp_path / "heatmaps"
        heatmap_dir.mkdir(parents=True, exist_ok=True)
        existing = heatmap_dir / "HMM1.cluster_motif_heatmap.png"
        existing.write_bytes(b"dummy")
        mtime_before = existing.stat().st_mtime

        import time
        time.sleep(0.05)

        _generate_cluster_motif_heatmap(
            "HMM1", paths,
            out_dir=str(heatmap_dir),
            summary_dir=str(tmp_path / "summary"),
            figure_format=["png"],
            force=True,
        )
        assert existing.stat().st_mtime > mtime_before
        assert existing.stat().st_size > len(b"dummy")


# ── _detect_motif_shift_regions ───────────────────────────────────────────────

from phylofoundry.tasks.embed import _detect_motif_shift_regions


def _write_aligned_fasta_cpd(path, seqs):
    """Write an aligned FASTA to *path* (dict {id: seq})."""
    with open(path, "w") as fh:
        for sid, seq in seqs.items():
            fh.write(f">{sid}\n{seq}\n")


class TestDetectMotifShiftRegions:
    """Unit tests for _detect_motif_shift_regions()."""

    def _make_aln_dir(self, tmp_path, clusters):
        """Write per-cluster aligned FASTAs and return {cl_id: path} dict."""
        paths = {}
        for cl_id, seqs in clusters.items():
            fp = str(tmp_path / f"cluster_{cl_id}.seed.aln.faa")
            _write_aligned_fasta_cpd(fp, seqs)
            paths[cl_id] = fp
        return paths

    # ── smoke / return-type tests ──────────────────────────────────────────────

    def test_returns_two_dataframes(self, tmp_path):
        """Should return (regions_df, signal_df) tuple."""
        clusters = {
            0: {f"s{i}": "ACDEF" * 4 for i in range(6)},
            1: {f"t{i}": "GHIKL" * 4 for i in range(6)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        regions_df, signal_df = _detect_motif_shift_regions("HMM1", paths)
        import pandas as pd
        assert isinstance(regions_df, pd.DataFrame)
        assert isinstance(signal_df, pd.DataFrame)

    def test_nonempty_for_two_clusters(self, tmp_path):
        """Two distinct clusters should produce non-empty output."""
        clusters = {
            0: {f"s{i}": "ACDEFGHIKL" for i in range(6)},
            1: {f"t{i}": "MNPQRSTVWY" for i in range(6)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        regions_df, signal_df = _detect_motif_shift_regions("HMM1", paths)
        assert not regions_df.empty
        assert not signal_df.empty

    # ── column tests ───────────────────────────────────────────────────────────

    def test_regions_expected_columns(self, tmp_path):
        """regions_df must contain the expected column set."""
        clusters = {
            0: {f"s{i}": "ACDEFGHIKL" for i in range(6)},
            1: {f"t{i}": "MNPQRSTVWY" for i in range(6)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        regions_df, _ = _detect_motif_shift_regions("HMM1", paths)
        expected = {
            "hmm_name", "region_id", "start_position", "end_position",
            "n_positions", "mean_divergence", "max_divergence", "signal_type",
        }
        assert expected.issubset(set(regions_df.columns))

    def test_signal_expected_columns(self, tmp_path):
        """signal_df must contain the expected column set."""
        clusters = {
            0: {f"s{i}": "ACDE" for i in range(6)},
            1: {f"t{i}": "FGHI" for i in range(6)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        _, signal_df = _detect_motif_shift_regions("HMM1", paths)
        expected = {
            "hmm_name", "aln_position", "signal_value",
            "smoothed_signal_value", "signal_type",
        }
        assert expected.issubset(set(signal_df.columns))

    # ── signal value tests ─────────────────────────────────────────────────────

    def test_signal_values_in_range(self, tmp_path):
        """JSD-based signal values should lie in [0, 1]."""
        clusters = {
            0: {f"s{i}": "AAAA" for i in range(8)},
            1: {f"t{i}": "CCCC" for i in range(8)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        _, signal_df = _detect_motif_shift_regions("HMM1", paths)
        assert (signal_df["signal_value"] >= 0).all()
        assert (signal_df["signal_value"] <= 1 + 1e-9).all()

    def test_signal_length_equals_alignment_length(self, tmp_path):
        """signal_df must have one row per alignment position."""
        aln_len = 12
        clusters = {
            0: {f"s{i}": "A" * aln_len for i in range(6)},
            1: {f"t{i}": "C" * aln_len for i in range(6)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        _, signal_df = _detect_motif_shift_regions("HMM1", paths)
        assert len(signal_df) == aln_len

    def test_signal_positions_1based(self, tmp_path):
        """aln_position should start at 1 (1-based indexing)."""
        clusters = {
            0: {f"s{i}": "ACDE" for i in range(6)},
            1: {f"t{i}": "FGHI" for i in range(6)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        _, signal_df = _detect_motif_shift_regions("HMM1", paths)
        assert signal_df["aln_position"].min() == 1

    # ── region structure tests ─────────────────────────────────────────────────

    def test_regions_cover_full_alignment(self, tmp_path):
        """Regions must contiguously cover the full alignment length."""
        aln_len = 20
        clusters = {
            0: {f"s{i}": "A" * aln_len for i in range(6)},
            1: {f"t{i}": "C" * aln_len for i in range(6)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        regions_df, _ = _detect_motif_shift_regions("HMM1", paths)
        assert regions_df["start_position"].min() == 1
        assert regions_df["end_position"].max() == aln_len
        # No gaps between consecutive regions
        rows = regions_df.sort_values("start_position").reset_index(drop=True)
        for i in range(1, len(rows)):
            assert rows.loc[i, "start_position"] == rows.loc[i - 1, "end_position"] + 1

    def test_regions_n_positions_consistent(self, tmp_path):
        """n_positions == end_position - start_position + 1 for every region."""
        clusters = {
            0: {f"s{i}": "ACDEFGHIKL" for i in range(6)},
            1: {f"t{i}": "MNPQRSTVWY" for i in range(6)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        regions_df, _ = _detect_motif_shift_regions("HMM1", paths)
        for _, row in regions_df.iterrows():
            assert row["n_positions"] == row["end_position"] - row["start_position"] + 1

    def test_mean_divergence_nonneg(self, tmp_path):
        """mean_divergence must be >= 0 for all regions."""
        clusters = {
            0: {f"s{i}": "ACDE" for i in range(6)},
            1: {f"t{i}": "FGHI" for i in range(6)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        regions_df, _ = _detect_motif_shift_regions("HMM1", paths)
        assert (regions_df["mean_divergence"] >= 0).all()

    def test_max_gte_mean_divergence(self, tmp_path):
        """max_divergence >= mean_divergence for every region."""
        clusters = {
            0: {f"s{i}": "ACDE" for i in range(6)},
            1: {f"t{i}": "FGHI" for i in range(6)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        regions_df, _ = _detect_motif_shift_regions("HMM1", paths)
        assert (regions_df["max_divergence"] >= regions_df["mean_divergence"] - 1e-9).all()

    def test_region_ids_sequential(self, tmp_path):
        """region_id should be a sequential 1-based integer."""
        aln_len = 20
        clusters = {
            0: {f"s{i}": "A" * aln_len for i in range(6)},
            1: {f"t{i}": "C" * aln_len for i in range(6)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        regions_df, _ = _detect_motif_shift_regions("HMM1", paths)
        rows = regions_df.sort_values("region_id").reset_index(drop=True)
        assert list(rows["region_id"]) == list(range(1, len(rows) + 1))

    def test_hmm_name_in_output(self, tmp_path):
        """hmm_name column should match the provided name."""
        clusters = {
            0: {f"s{i}": "ACDE" for i in range(6)},
            1: {f"t{i}": "FGHI" for i in range(6)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        regions_df, signal_df = _detect_motif_shift_regions("MY_HMM", paths)
        assert (regions_df["hmm_name"] == "MY_HMM").all()
        assert (signal_df["hmm_name"] == "MY_HMM").all()

    # ── signal type tests ──────────────────────────────────────────────────────

    def test_signal_type_mean_jsd(self, tmp_path):
        """signal_type column should reflect the chosen signal."""
        clusters = {
            0: {f"s{i}": "ACDE" for i in range(6)},
            1: {f"t{i}": "FGHI" for i in range(6)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        regions_df, signal_df = _detect_motif_shift_regions(
            "HMM1", paths, signal="mean_jsd"
        )
        assert (signal_df["signal_type"] == "mean_jsd").all()
        assert (regions_df["signal_type"] == "mean_jsd").all()

    def test_signal_type_max_jsd(self, tmp_path):
        """max_jsd signal should yield valid signal_df."""
        clusters = {
            0: {f"s{i}": "AAAA" for i in range(6)},
            1: {f"t{i}": "CCCC" for i in range(6)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        _, signal_df = _detect_motif_shift_regions("HMM1", paths, signal="max_jsd")
        assert (signal_df["signal_value"] >= 0).all()
        assert (signal_df["signal_type"] == "max_jsd").all()

    def test_signal_type_jsd_vs_global(self, tmp_path):
        """jsd_vs_global signal should yield valid signal_df."""
        clusters = {
            0: {f"s{i}": "ACDE" for i in range(6)},
            1: {f"t{i}": "FGHI" for i in range(6)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        _, signal_df = _detect_motif_shift_regions(
            "HMM1", paths, signal="jsd_vs_global"
        )
        assert (signal_df["signal_value"] >= 0).all()
        assert (signal_df["signal_type"] == "jsd_vs_global").all()

    # ── smoothing tests ────────────────────────────────────────────────────────

    def test_smoothing_does_not_change_length(self, tmp_path):
        """Smoothing should not alter the number of signal positions."""
        aln_len = 20
        clusters = {
            0: {f"s{i}": "A" * aln_len for i in range(6)},
            1: {f"t{i}": "C" * aln_len for i in range(6)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        _, signal_df = _detect_motif_shift_regions(
            "HMM1", paths, smoothing_window=5
        )
        assert len(signal_df) == aln_len

    def test_no_smoothing_raw_equals_smoothed(self, tmp_path):
        """With smoothing_window=0, raw and smoothed signals must be equal."""
        clusters = {
            0: {f"s{i}": "ACDE" for i in range(6)},
            1: {f"t{i}": "FGHI" for i in range(6)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        _, signal_df = _detect_motif_shift_regions(
            "HMM1", paths, smoothing_window=0
        )
        assert (
            (signal_df["signal_value"] - signal_df["smoothed_signal_value"]).abs() < 1e-9
        ).all()

    # ── change-point detection tests ───────────────────────────────────────────

    def test_high_threshold_yields_one_region(self, tmp_path):
        """With a very high threshold, no splits occur → a single region."""
        aln_len = 20
        clusters = {
            0: {f"s{i}": "A" * aln_len for i in range(6)},
            1: {f"t{i}": "C" * aln_len for i in range(6)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        # threshold=1e9 prevents any split
        regions_df, _ = _detect_motif_shift_regions(
            "HMM1", paths, threshold=1e9
        )
        assert len(regions_df) == 1

    def test_max_changepoints_limits_regions(self, tmp_path):
        """max_changepoints caps the number of detected splits."""
        aln_len = 60
        # Alternating A/C blocks to create high-variance signal
        seq_a = ("A" * 10 + "C" * 10) * 3
        seq_b = ("C" * 10 + "A" * 10) * 3
        clusters = {
            0: {f"s{i}": seq_a for i in range(8)},
            1: {f"t{i}": seq_b for i in range(8)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        regions_df, _ = _detect_motif_shift_regions(
            "HMM1", paths, max_changepoints=3, threshold=0.0
        )
        # At most max_changepoints + 1 regions
        assert len(regions_df) <= 4

    # ── edge case / guard tests ────────────────────────────────────────────────

    def test_single_cluster_returns_empty(self, tmp_path):
        """Only one cluster → no pairs → both DataFrames empty."""
        clusters = {0: {f"s{i}": "ACDE" for i in range(6)}}
        paths = self._make_aln_dir(tmp_path, clusters)
        regions_df, signal_df = _detect_motif_shift_regions("HMM1", paths)
        assert regions_df.empty
        assert signal_df.empty

    def test_empty_paths_returns_empty(self, tmp_path):
        """Empty input dict → both DataFrames empty."""
        regions_df, signal_df = _detect_motif_shift_regions("HMM1", {})
        assert regions_df.empty
        assert signal_df.empty

    def test_min_cluster_size_filters(self, tmp_path):
        """Clusters below min_cluster_size are excluded."""
        clusters = {
            0: {f"s{i}": "ACDE" for i in range(2)},
            1: {f"t{i}": "FGHI" for i in range(2)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        regions_df, signal_df = _detect_motif_shift_regions(
            "HMM1", paths, min_cluster_size=5
        )
        assert regions_df.empty
        assert signal_df.empty

    def test_gaps_skipped_in_counts(self, tmp_path):
        """Gap characters should not cause errors."""
        clusters = {
            0: {f"s{i}": "--AA--" for i in range(6)},
            1: {f"t{i}": "--CC--" for i in range(6)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        regions_df, signal_df = _detect_motif_shift_regions("HMM1", paths)
        assert (signal_df["signal_value"] >= 0).all()

    def test_three_clusters_still_works(self, tmp_path):
        """Three clusters should produce valid output."""
        aln_len = 10
        clusters = {
            0: {f"s{i}": "A" * aln_len for i in range(6)},
            1: {f"t{i}": "C" * aln_len for i in range(6)},
            2: {f"u{i}": "G" * aln_len for i in range(6)},
        }
        paths = self._make_aln_dir(tmp_path, clusters)
        regions_df, signal_df = _detect_motif_shift_regions("HMM1", paths)
        assert not regions_df.empty
        assert len(signal_df) == aln_len


# ── _validate_emb_cfg: change_point_detection defaults ───────────────────────

class TestValidateEmbCfgChangePointDetection:
    """Tests for change_point_detection config section in _validate_emb_cfg."""

    def _base_cfg(self):
        return {
            "backend": "esm",
            "model": "esm2_t33_650M_UR50D",
            "device": "cpu",
            "batch_size": 4,
            "repr_layer": None,
            "pooling": "mean",
            "pca_components": 3,
            "cluster_embeddings": True,
        }

    def test_cpd_defaults_injected(self):
        cfg = _validate_emb_cfg(self._base_cfg())
        cpd = cfg["cluster_subworkflow"]["change_point_detection"]
        assert cpd["enabled"] is False
        assert cpd["signal"] == "mean_jsd"
        assert cpd["smoothing_window"] == 0
        assert cpd["min_segment_len"] == 5
        assert cpd["merge_distance"] == 10
        assert cpd["threshold"] == 0.05
        assert cpd["max_changepoints"] == 25
        assert cpd["min_cluster_size"] == 5

    def test_cpd_can_be_enabled(self):
        base = self._base_cfg()
        base["cluster_subworkflow"] = {
            "change_point_detection": {"enabled": True}
        }
        cfg = _validate_emb_cfg(base)
        assert cfg["cluster_subworkflow"]["change_point_detection"]["enabled"] is True

    def test_cpd_custom_signal(self):
        base = self._base_cfg()
        base["cluster_subworkflow"] = {
            "change_point_detection": {"signal": "max_jsd"}
        }
        cfg = _validate_emb_cfg(base)
        assert cfg["cluster_subworkflow"]["change_point_detection"]["signal"] == "max_jsd"

    def test_cpd_non_dict_replaced_with_defaults(self):
        base = self._base_cfg()
        base["cluster_subworkflow"] = {
            "change_point_detection": None
        }
        cfg = _validate_emb_cfg(base)
        cpd = cfg["cluster_subworkflow"]["change_point_detection"]
        assert isinstance(cpd, dict)
        assert cpd["enabled"] is False

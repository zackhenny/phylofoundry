"""Tests for embedding clustering enhancements:
- cluster_on (PCA vs embeddings)
- cluster_method (hdbscan vs leiden)
- hdbscan_metric
- umap_dimensions (2D vs 3D)
- Config validation (_validate_emb_cfg)
"""
import numpy as np
import pytest

from phylofoundry.tasks.embed import (
    _run_hdbscan,
    _run_leiden,
    _validate_emb_cfg,
    _save_umap_plot,
)

# Optional dependency guards
try:
    from sklearn.cluster import HDBSCAN as _HDBSCAN_CLS
    _HDBSCAN_AVAILABLE = True
except ImportError:
    _HDBSCAN_AVAILABLE = False

try:
    import matplotlib  # noqa: F401
    _MPL_AVAILABLE = True
except ImportError:
    _MPL_AVAILABLE = False

requires_hdbscan = pytest.mark.skipif(
    not _HDBSCAN_AVAILABLE,
    reason="scikit-learn >= 1.3 with HDBSCAN not installed",
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


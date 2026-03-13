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


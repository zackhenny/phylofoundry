"""Tests for the aa_composition task module."""

import os
import math
import pytest
import pandas as pd

from phylofoundry.tasks.aa_composition import (
    STANDARD_AAS,
    WEIGHTED_ZC,
    CARBON_NUMBER,
    HYDROPHOBICITY,
    NH2O_QEC,
    THERMOSTABLE_RESIDUES,
    _clean_sequence,
    compute_gene_metrics,
    run_comparative_stats,
    _benjamini_hochberg,
    run_aa_composition,
)


# ---------------------------------------------------------------------------
# _clean_sequence
# ---------------------------------------------------------------------------

class TestCleanSequence:
    def test_removes_gaps(self):
        assert _clean_sequence("A-C.D*E") == list("ACDE")

    def test_removes_nonstandard(self):
        # 'B', 'Z', 'X' are not in AA_ALPHABET
        assert _clean_sequence("ABCXZ") == list("AC")

    def test_removes_initial_met(self):
        result = _clean_sequence("MAACC", remove_initial_met=True)
        assert result[0] == "A"
        assert len(result) == 4

    def test_keeps_internal_met(self):
        result = _clean_sequence("AMCC", remove_initial_met=True)
        assert "M" in result

    def test_no_remove_met(self):
        result = _clean_sequence("MAACC", remove_initial_met=False)
        assert result[0] == "M"
        assert len(result) == 5

    def test_uppercase_conversion(self):
        assert _clean_sequence("acde") == list("ACDE")

    def test_empty_input(self):
        assert _clean_sequence("") == []

    def test_all_nonstandard(self):
        assert _clean_sequence("BXZU") == []


# ---------------------------------------------------------------------------
# compute_gene_metrics
# ---------------------------------------------------------------------------

class TestComputeGeneMetrics:
    def test_basic_composition(self):
        # A simple sequence with known composition
        seq = "AAAC"  # 3 A + 1 C, length = 4 (no Met removal since doesn't start with M)
        row = compute_gene_metrics("g1", seq, "HMM1", compute_pi=False,
                                   remove_initial_met=False)
        assert row["gene_id"] == "g1"
        assert row["hmm"] == "HMM1"
        assert row["len"] == 4
        assert abs(row["aa_A_comp"] - 0.75) < 1e-9
        assert abs(row["aa_C_comp"] - 0.25) < 1e-9
        # All other AAs should be 0.0
        for aa in STANDARD_AAS:
            if aa not in ("A", "C"):
                assert row[f"aa_{aa}_comp"] == 0.0

    def test_zc_single_residue(self):
        # Single A: Zc = WEIGHTED_ZC["A"] / CARBON_NUMBER["A"]
        row = compute_gene_metrics("g1", "A", "HMM1", compute_pi=False,
                                   remove_initial_met=False)
        expected_zc = WEIGHTED_ZC["A"] / CARBON_NUMBER["A"]
        assert abs(row["Zc"] - expected_zc) < 1e-9

    def test_gravy_single_residue(self):
        row = compute_gene_metrics("g1", "I", "HMM1", compute_pi=False,
                                   remove_initial_met=False)
        assert abs(row["gravy"] - HYDROPHOBICITY["I"]) < 1e-9

    def test_nH2O_single_residue(self):
        row = compute_gene_metrics("g1", "W", "HMM1", compute_pi=False,
                                   remove_initial_met=False)
        assert abs(row["nH2O"] - NH2O_QEC["W"]) < 1e-9

    def test_thermo_freq_all_thermo(self):
        # Sequence of all thermostable residues
        seq = "IVYWREL"
        row = compute_gene_metrics("g1", seq, "HMM1", compute_pi=False,
                                   remove_initial_met=False)
        assert abs(row["thermo_freq"] - 1.0) < 1e-9

    def test_thermo_freq_none(self):
        # A and C are not in THERMOSTABLE_RESIDUES
        seq = "AAAACC"
        row = compute_gene_metrics("g1", seq, "HMM1", compute_pi=False,
                                   remove_initial_met=False)
        assert row["thermo_freq"] == 0.0

    def test_empty_sequence_returns_nan(self):
        row = compute_gene_metrics("g1", "", "HMM1", compute_pi=False)
        assert row["len"] == 0
        assert math.isnan(row["Zc"])
        assert math.isnan(row["aa_A_comp"])

    def test_tip_or_clade(self):
        row = compute_gene_metrics("g1", "ACDE", "HMM1", tip_or_clade="cladeA",
                                   compute_pi=False, remove_initial_met=False)
        assert row["tip_or_clade"] == "cladeA"

    def test_sulfur_content_cm(self):
        # C and M are sulfur-containing
        seq = "ACMDEF"  # 1 C + 1 M = 2 S atoms in 6 residues
        row = compute_gene_metrics("g1", seq, "HMM1", compute_pi=False,
                                   remove_initial_met=False)
        assert abs(row["S_content"] - 2 / 6) < 1e-9

    def test_pi_computed_when_requested(self):
        # Requires BioPython — skip if not installed
        pytest.importorskip("Bio")
        seq = "ACDEFGHIKLMNPQRSTVWY"
        row = compute_gene_metrics("g1", seq, "HMM1", compute_pi=True,
                                   remove_initial_met=False)
        assert "pi" in row
        # pI should be a positive float
        assert not math.isnan(row["pi"])
        assert row["pi"] > 0

    def test_pi_nan_when_not_requested(self):
        row = compute_gene_metrics("g1", "ACDE", "HMM1", compute_pi=False,
                                   remove_initial_met=False)
        assert "pi" not in row or math.isnan(row.get("pi", float("nan")))


# ---------------------------------------------------------------------------
# _benjamini_hochberg
# ---------------------------------------------------------------------------

class TestBenjaminiHochberg:
    def test_all_significant(self):
        pvalues = [0.001, 0.002, 0.003]
        q = _benjamini_hochberg(pvalues)
        # All q-values ≤ 1
        assert all(v <= 1.0 for v in q)
        # Monotone: largest p → largest q
        assert q[2] >= q[0]

    def test_empty(self):
        assert _benjamini_hochberg([]) == []

    def test_single_pvalue(self):
        q = _benjamini_hochberg([0.04])
        assert len(q) == 1
        assert abs(q[0] - 0.04) < 1e-9

    def test_identical_pvalues(self):
        pvalues = [0.05, 0.05, 0.05]
        q = _benjamini_hochberg(pvalues)
        assert len(q) == 3


# ---------------------------------------------------------------------------
# run_comparative_stats
# ---------------------------------------------------------------------------

class TestRunComparativeStats:
    def _make_df(self):
        import numpy as np
        rng = np.random.default_rng(42)
        n = 30
        clades = ["A"] * 10 + ["B"] * 10 + ["C"] * 10
        aa_A = rng.normal(0.1, 0.01, n).tolist()
        aa_C = rng.normal(0.05, 0.005, n).tolist()
        return pd.DataFrame({
            "tip_or_clade": clades,
            "aa_A_comp": aa_A,
            "aa_C_comp": aa_C,
        })

    def test_returns_dataframe(self):
        pytest.importorskip("scipy")
        df = self._make_df()
        stats = run_comparative_stats(df, ["aa_A_comp", "aa_C_comp"])
        assert isinstance(stats, pd.DataFrame)
        assert "metric" in stats.columns
        assert "H_stat" in stats.columns
        assert "p_value" in stats.columns
        assert "q_value" in stats.columns
        assert len(stats) == 2

    def test_single_group_returns_empty(self):
        pytest.importorskip("scipy")
        df = pd.DataFrame({
            "tip_or_clade": ["A"] * 5,
            "aa_A_comp": [0.1] * 5,
        })
        stats = run_comparative_stats(df, ["aa_A_comp"])
        assert stats.empty

    def test_q_values_bounded(self):
        pytest.importorskip("scipy")
        df = self._make_df()
        stats = run_comparative_stats(df, ["aa_A_comp", "aa_C_comp"])
        valid_q = stats["q_value"].dropna()
        assert (valid_q >= 0).all() and (valid_q <= 1).all()


# ---------------------------------------------------------------------------
# run_aa_composition (integration)
# ---------------------------------------------------------------------------

class TestRunAaComposition:
    def _write_fasta(self, path, records):
        with open(path, "w") as f:
            for h, s in records.items():
                f.write(f">{h}\n{s}\n")

    def test_basic_run(self, tmp_path):
        """run_aa_composition produces expected output TSV."""
        fasta_dir = tmp_path / "fasta_per_hmm"
        fasta_dir.mkdir()
        summary_dir = tmp_path / "summary"
        summary_dir.mkdir()

        # Write two per-HMM FASTA files
        seqs1 = {
            "gene_A1": "ACDEFGHIKLMNPQRSTVWY",
            "gene_A2": "AAAAAAAAAA",
        }
        seqs2 = {
            "gene_B1": "CCCCCCCCCC",
        }
        self._write_fasta(str(fasta_dir / "HMM1.faa"), seqs1)
        self._write_fasta(str(fasta_dir / "HMM2.faa"), seqs2)

        cfg = {
            "aa_composition": {
                "enabled": True,
                "remove_initial_met": False,
                "compute_pi": False,
                "generate_plots": False,
                "top_n_aas_heatmap": 5,
            }
        }

        df = run_aa_composition(cfg, str(fasta_dir), str(summary_dir),
                                hmm_keep=None, clade_assign_dir=None)

        assert isinstance(df, pd.DataFrame)
        assert len(df) == 3  # 2 + 1 genes
        assert "gene_id" in df.columns
        assert "hmm" in df.columns
        assert "len" in df.columns
        assert "aa_A_comp" in df.columns
        assert "Zc" in df.columns

        # gene_A2 is all-A
        row_a2 = df[df["gene_id"] == "gene_A2"].iloc[0]
        assert abs(row_a2["aa_A_comp"] - 1.0) < 1e-9
        assert row_a2["len"] == 10

        # Output file should exist
        out_tsv = str(summary_dir / "aa_composition" / "aa_comp_per_gene.tsv")
        assert os.path.exists(out_tsv)
        loaded = pd.read_csv(out_tsv, sep="\t")
        assert len(loaded) == 3

    def test_hmm_keep_filters(self, tmp_path):
        fasta_dir = tmp_path / "fasta_per_hmm"
        fasta_dir.mkdir()
        summary_dir = tmp_path / "summary"
        summary_dir.mkdir()

        self._write_fasta(str(fasta_dir / "HMM1.faa"), {"g1": "AAAA"})
        self._write_fasta(str(fasta_dir / "HMM2.faa"), {"g2": "CCCC"})

        cfg = {"aa_composition": {"compute_pi": False, "generate_plots": False,
                                   "remove_initial_met": False}}

        df = run_aa_composition(cfg, str(fasta_dir), str(summary_dir),
                                hmm_keep={"HMM1"})

        assert len(df) == 1
        assert df.iloc[0]["hmm"] == "HMM1"

    def test_empty_fasta_dir_returns_empty_df(self, tmp_path):
        fasta_dir = tmp_path / "empty"
        fasta_dir.mkdir()
        summary_dir = tmp_path / "summary"
        summary_dir.mkdir()

        cfg = {"aa_composition": {"compute_pi": False, "generate_plots": False}}
        df = run_aa_composition(cfg, str(fasta_dir), str(summary_dir))
        assert df.empty

    def test_clade_assignment_loaded(self, tmp_path):
        fasta_dir = tmp_path / "fasta_per_hmm"
        fasta_dir.mkdir()
        summary_dir = tmp_path / "summary"
        summary_dir.mkdir()
        clade_dir = tmp_path / "clade_assignments"
        clade_dir.mkdir()

        self._write_fasta(str(fasta_dir / "HMM1.faa"),
                          {"g1": "AAAA", "g2": "CCCC"})

        # Write per-HMM clade TSV
        clade_tsv = clade_dir / "HMM1.tsv"
        pd.DataFrame({
            "tip": ["g1", "g2"],
            "clade_name": ["cladeX", "cladeY"],
        }).to_csv(str(clade_tsv), sep="\t", index=False)

        cfg = {"aa_composition": {"compute_pi": False, "generate_plots": False,
                                   "remove_initial_met": False}}

        df = run_aa_composition(cfg, str(fasta_dir), str(summary_dir),
                                clade_assign_dir=str(clade_dir))

        assert set(df["tip_or_clade"]) == {"cladeX", "cladeY"}

    def test_column_order(self, tmp_path):
        fasta_dir = tmp_path / "fasta_per_hmm"
        fasta_dir.mkdir()
        summary_dir = tmp_path / "summary"
        summary_dir.mkdir()

        self._write_fasta(str(fasta_dir / "H1.faa"), {"g1": "ACDE"})
        cfg = {"aa_composition": {"compute_pi": False, "generate_plots": False,
                                   "remove_initial_met": False}}
        df = run_aa_composition(cfg, str(fasta_dir), str(summary_dir))

        cols = list(df.columns)
        assert cols[0] == "gene_id"
        assert cols[1] == "hmm"
        assert cols[2] == "len"
        assert cols[3] == "tip_or_clade"
        # AA columns come next
        aa_cols = [c for c in cols if c.startswith("aa_")]
        assert len(aa_cols) == 20

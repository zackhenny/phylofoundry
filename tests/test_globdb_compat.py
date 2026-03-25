"""Tests for GlobDB compatibility options.

Verifies that the pipeline correctly:
- Accepts ``--diamond-db``, ``--combined-faa``, and ``--globdb-taxonomy`` CLI args.
- Applies them to the resolved config.
- Validates that missing files produce helpful errors.
- Skips prep when a prebuilt combined_faa is supplied.
- Uses a prebuilt DIAMOND DB instead of calling makedb.
- Loads proteome sequences from a combined_faa when faa_dir is not set.
"""

from __future__ import annotations

import argparse
import os
from copy import deepcopy
from unittest.mock import MagicMock, patch

import pytest

from phylofoundry.constants import DEFAULT_CONFIG
from phylofoundry.main import _apply_step_args, _build_parser


# ── Helpers ───────────────────────────────────────────────────────────────────


def _parse(argv: list[str]) -> argparse.Namespace:
    """Parse *argv* through the full argument parser (with legacy routing)."""
    from phylofoundry.main import _ALL_SUBCMDS

    if not argv or argv[0] not in _ALL_SUBCMDS:
        argv = ["run"] + argv
    ap = _build_parser()
    return ap.parse_args(argv)


def _valid_cfg(**overrides) -> dict:
    cfg = deepcopy(DEFAULT_CONFIG)
    cfg["inputs"]["faa_dir"] = "/tmp/faa"
    cfg["inputs"]["hmm_input"] = "/tmp/hmm"
    cfg["output"]["outdir"] = "/tmp/out"
    cfg.update(overrides)
    return cfg


# ── CLI argument parsing ───────────────────────────────────────────────────────


class TestCLIArgs:
    """Verify that the new CLI flags are parsed and applied correctly."""

    def test_run_diamond_db_arg(self):
        args = _parse(["run", "--outdir", "/tmp/out", "--diamond_db", "/path/to/db.dmnd"])
        assert args.diamond_db == "/path/to/db.dmnd"

    def test_run_combined_faa_arg(self):
        args = _parse(["run", "--outdir", "/tmp/out", "--combined_faa", "/path/to/combined.faa"])
        assert args.combined_faa == "/path/to/combined.faa"

    def test_run_globdb_taxonomy_arg(self):
        args = _parse(["run", "--outdir", "/tmp/out", "--globdb_taxonomy", "/path/to/tax.tsv"])
        assert args.globdb_taxonomy == "/path/to/tax.tsv"

    def test_hmmer_diamond_db_arg(self):
        args = _parse(["hmmer", "--outdir", "/tmp/out", "--diamond_db", "/path/to/db.dmnd"])
        assert args.diamond_db == "/path/to/db.dmnd"

    def test_hmmer_combined_faa_arg(self):
        args = _parse(["hmmer", "--outdir", "/tmp/out", "--combined_faa", "/path/to/combined.faa"])
        assert args.combined_faa == "/path/to/combined.faa"

    def test_taxonomy_globdb_taxonomy_arg(self):
        args = _parse(["taxonomy", "--outdir", "/tmp/out", "--globdb_taxonomy", "/path/to/tax.tsv"])
        assert args.globdb_taxonomy == "/path/to/tax.tsv"

    def test_prep_combined_faa_arg(self):
        args = _parse(["prep", "--outdir", "/tmp/out", "--combined_faa", "/path/to/combined.faa"])
        assert args.combined_faa == "/path/to/combined.faa"


# ── _apply_step_args propagation ──────────────────────────────────────────────


class TestApplyStepArgs:
    """Verify that CLI args are propagated into the config dict."""

    def test_run_diamond_db_applied_to_cfg(self):
        cfg = _valid_cfg()
        args = _parse(["run", "--outdir", "/tmp/out", "--diamond_db", "/db.dmnd"])
        _apply_step_args(args, cfg)
        assert cfg["inputs"]["diamond_db"] == "/db.dmnd"

    def test_run_combined_faa_applied_to_cfg(self):
        cfg = _valid_cfg()
        args = _parse(["run", "--outdir", "/tmp/out", "--combined_faa", "/combined.faa"])
        _apply_step_args(args, cfg)
        assert cfg["inputs"]["combined_faa"] == "/combined.faa"

    def test_run_globdb_taxonomy_applied_to_cfg(self):
        cfg = _valid_cfg()
        args = _parse(["run", "--outdir", "/tmp/out", "--globdb_taxonomy", "/tax.tsv"])
        _apply_step_args(args, cfg)
        assert cfg["inputs"]["globdb_taxonomy_file"] == "/tax.tsv"

    def test_taxonomy_globdb_taxonomy_applied_to_cfg(self):
        cfg = _valid_cfg()
        args = _parse(["taxonomy", "--outdir", "/tmp/out", "--globdb_taxonomy", "/tax.tsv"])
        _apply_step_args(args, cfg)
        assert cfg["inputs"]["globdb_taxonomy_file"] == "/tax.tsv"

    def test_hmmer_diamond_db_applied_to_cfg(self):
        cfg = _valid_cfg()
        args = _parse(["hmmer", "--outdir", "/tmp/out", "--diamond_db", "/db.dmnd"])
        _apply_step_args(args, cfg)
        assert cfg["inputs"]["diamond_db"] == "/db.dmnd"


# ── DEFAULT_CONFIG includes new fields ────────────────────────────────────────


class TestDefaultConfig:
    def test_diamond_db_in_inputs(self):
        assert "diamond_db" in DEFAULT_CONFIG["inputs"]
        assert DEFAULT_CONFIG["inputs"]["diamond_db"] is None

    def test_combined_faa_in_inputs(self):
        assert "combined_faa" in DEFAULT_CONFIG["inputs"]
        assert DEFAULT_CONFIG["inputs"]["combined_faa"] is None

    def test_globdb_taxonomy_file_in_inputs(self):
        assert "globdb_taxonomy_file" in DEFAULT_CONFIG["inputs"]
        assert DEFAULT_CONFIG["inputs"]["globdb_taxonomy_file"] is None


# ── validate_config with prebuilt files ───────────────────────────────────────


class TestValidateConfigPrebuilt:
    """validate_config should fail fast when prebuilt paths don't exist."""

    def test_missing_combined_faa_raises(self, tmp_path):
        from phylofoundry.config import validate_config

        cfg = _valid_cfg()
        cfg["inputs"]["combined_faa"] = str(tmp_path / "nonexistent.faa")
        with pytest.raises(SystemExit, match="combined_faa"):
            validate_config(cfg)

    def test_missing_diamond_db_raises(self, tmp_path):
        from phylofoundry.config import validate_config

        cfg = _valid_cfg()
        cfg["inputs"]["diamond_db"] = str(tmp_path / "nonexistent.dmnd")
        with pytest.raises(SystemExit, match="diamond_db"):
            validate_config(cfg)

    def test_missing_globdb_taxonomy_raises(self, tmp_path):
        from phylofoundry.config import validate_config

        cfg = _valid_cfg()
        cfg["inputs"]["globdb_taxonomy_file"] = str(tmp_path / "nonexistent.tsv")
        with pytest.raises(SystemExit, match="globdb_taxonomy_file"):
            validate_config(cfg)

    def test_valid_combined_faa_passes(self, tmp_path):
        from phylofoundry.config import validate_config

        faa = tmp_path / "combined.faa"
        faa.write_text(">g1~p1\nMPPP\n")
        cfg = _valid_cfg()
        cfg["inputs"]["combined_faa"] = str(faa)
        # Should not raise
        validate_config(cfg)

    def test_combined_faa_allows_no_faa_dir_in_diamond_mode(self, tmp_path):
        """When combined_faa and diamond_db are provided, faa_dir is not required."""
        from phylofoundry.config import validate_config

        faa = tmp_path / "combined.faa"
        faa.write_text(">g1~p1\nMPPP\n")
        db = tmp_path / "db.dmnd"
        db.write_text("")  # minimal stub
        query = tmp_path / "query.faa"
        query.write_text(">query1\nMPPP\n")

        cfg = _valid_cfg()
        cfg["inputs"]["faa_dir"] = None
        cfg["inputs"]["combined_faa"] = str(faa)
        cfg["inputs"]["diamond_db"] = str(db)
        cfg["inputs"]["diamond_query"] = str(query)
        cfg["diamond"]["enabled"] = True
        # Should not raise
        validate_config(cfg)


# ── Diamond task uses prebuilt DB ─────────────────────────────────────────────


class TestDiamondUsesPrebuiltDB:
    """run_diamond skips makedb when inputs.diamond_db is set."""

    def test_prebuilt_db_skips_makedb(self, tmp_path):
        from phylofoundry.tasks import diamond as diamond_task

        # Write a minimal query FASTA and stub DB file
        query_fp = tmp_path / "query.faa"
        query_fp.write_text(">q1\nMPPP\n")
        db_dmnd = tmp_path / "prebuilt.dmnd"
        db_dmnd.write_text("")  # dummy

        cfg = deepcopy(DEFAULT_CONFIG)
        cfg["inputs"]["diamond_query"] = str(query_fp)
        cfg["inputs"]["diamond_db"] = str(db_dmnd)
        cfg["resources"]["cpu"] = 1

        summary_dir = tmp_path / "summary"
        summary_dir.mkdir()

        with patch.object(diamond_task, "make_diamond_db") as mock_makedb, \
             patch.object(diamond_task, "worker_diamond_search", return_value=None):
            import pandas as pd
            with patch("phylofoundry.tasks.diamond.ProcessPoolExecutor") as mock_pool:
                mock_pool.return_value.__enter__.return_value.map.return_value = []
                diamond_task.run_diamond(
                    cfg, [], "/irrelevant/combined.faa",
                    str(tmp_path), str(summary_dir), force=False
                )

        mock_makedb.assert_not_called()


# ── Pipeline prep skip with prebuilt combined_faa ─────────────────────────────


class TestPipelinePrepSkip:
    """The prep step is skipped when inputs.combined_faa is provided."""

    def test_prep_skipped_when_combined_faa_provided(self, tmp_path):
        """Pipeline logs skip message and marks prep success without running it."""
        pytest.importorskip("Bio", reason="BioPython required for pipeline import")
        from phylofoundry.pipeline import run_pipeline
        from phylofoundry.tasks import prep as prep_task

        combined_faa = tmp_path / "combined.faa"
        combined_faa.write_text(">genome1.faa~prot1\nMPPP\n")

        cfg = deepcopy(DEFAULT_CONFIG)
        cfg["inputs"]["combined_faa"] = str(combined_faa)
        cfg["inputs"]["hmm_input"] = str(tmp_path / "markers")  # won't be read
        cfg["output"]["outdir"] = str(tmp_path / "out")
        cfg["workflow"]["start_at"] = "prep"
        cfg["workflow"]["stop_after"] = "prep"
        cfg["diamond"]["enabled"] = True
        cfg["inputs"]["diamond_query"] = str(tmp_path / "query.faa")
        # Create dummy query file so validator is happy
        (tmp_path / "query.faa").write_text(">q1\nMPPP\n")

        with patch.object(prep_task, "run_prep_diamond_mode") as mock_prep, \
             patch.object(prep_task, "run_prep") as mock_prep_hmm:
            run_pipeline(cfg)

        mock_prep.assert_not_called()
        mock_prep_hmm.assert_not_called()


# ── Load proteomes from combined_faa ─────────────────────────────────────────


class TestLoadProteomesFromCombinedFaa:
    """_load_proteomes_from_combined_faa correctly reconstructs proteome_seqs."""

    def test_basic_parse(self, tmp_path):
        pytest.importorskip("Bio", reason="BioPython required for pipeline import")
        from phylofoundry.pipeline import _load_proteomes_from_combined_faa

        faa = tmp_path / "combined.faa"
        faa.write_text(
            ">genome1.faa~prot1\nMPPP\n"
            ">genome1.faa~prot2\nMAAA\n"
            ">genome2.faa~prot3\nMCCC\n"
        )

        result = _load_proteomes_from_combined_faa(str(faa))

        assert "genome1.faa" in result
        assert "genome2.faa" in result
        assert result["genome1.faa"]["prot1"] == "MPPP"
        assert result["genome1.faa"]["prot2"] == "MAAA"
        assert result["genome2.faa"]["prot3"] == "MCCC"

    def test_no_tilde_fallback(self, tmp_path):
        pytest.importorskip("Bio", reason="BioPython required for pipeline import")
        from phylofoundry.pipeline import _load_proteomes_from_combined_faa

        faa = tmp_path / "combined.faa"
        faa.write_text(">just_a_protein\nMPPP\n")

        result = _load_proteomes_from_combined_faa(str(faa))

        assert "unknown" in result
        assert "just_a_protein" in result["unknown"]

    def test_basic_parse_via_utils(self, tmp_path):
        """Test the underlying logic without importing pipeline (avoids Bio dep)."""
        from phylofoundry.utils.bio import read_fasta

        faa = tmp_path / "combined.faa"
        faa.write_text(
            ">genome1.faa~prot1\nMPPP\n"
            ">genome1.faa~prot2\nMAAA\n"
            ">genome2.faa~prot3\nMCCC\n"
        )

        # Replicate the same logic as _load_proteomes_from_combined_faa
        flat = read_fasta(str(faa))
        proteome_seqs: dict = {}
        for header, seq in flat.items():
            if "~" in header:
                genome, protein = header.split("~", 1)
            else:
                genome = "unknown"
                protein = header
            proteome_seqs.setdefault(genome, {})[protein] = seq

        assert "genome1.faa" in proteome_seqs
        assert proteome_seqs["genome1.faa"]["prot1"] == "MPPP"
        assert proteome_seqs["genome2.faa"]["prot3"] == "MCCC"

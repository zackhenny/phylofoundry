"""Tests for the input provenance and artifact-linking system (--input-run)."""

from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path

import pytest


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_fake_run(tmp_path: Path, run_id: str = "test-run-001") -> Path:
    """Create a minimal fake PhyloFoundry output directory with checkpoint DB."""
    outdir = tmp_path / "fake_run"
    outdir.mkdir()

    # Create the logs/ directory and checkpoint DB
    logs = outdir / "logs"
    logs.mkdir()

    import sqlite3
    db_path = logs / "checkpoint.db"
    conn = sqlite3.connect(str(db_path))
    conn.executescript(
        """
        CREATE TABLE IF NOT EXISTS runs (
            run_id TEXT PRIMARY KEY,
            start_time TEXT,
            end_time TEXT,
            outdir TEXT,
            config_path TEXT,
            ndjson_path TEXT,
            sqlite_path TEXT
        );
        CREATE TABLE IF NOT EXISTS steps (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            run_id TEXT,
            step_name TEXT,
            sequence INTEGER,
            state TEXT,
            input_fingerprint TEXT,
            params TEXT,
            started_at TEXT,
            completed_at TEXT,
            outputs TEXT,
            error TEXT,
            trace TEXT
        );
        CREATE TABLE IF NOT EXISTS artifacts (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            run_id TEXT,
            step_name TEXT,
            path TEXT,
            sha256 TEXT,
            size_bytes INTEGER,
            created_at TEXT
        );
        """
    )
    conn.execute(
        "INSERT INTO runs (run_id, start_time, outdir) VALUES (?, ?, ?)",
        (run_id, "2024-01-15T10:00:00Z", str(outdir)),
    )
    conn.execute(
        "INSERT INTO steps (run_id, step_name, state) VALUES (?, ?, ?)",
        (run_id, "phylo", "success"),
    )
    conn.commit()
    conn.close()

    # Write a minimal config.yaml with _run_meta
    config_yaml = f"""
_run_meta:
  run_id: "{run_id}"
  start_time: "2024-01-15T10:00:00Z"
  ndjson_path: "{logs / 'pipeline.ndjson'}"
  sqlite_path: "{db_path}"
  outdir: "{outdir}"
output:
  outdir: "{outdir}"
"""
    (outdir / "config.yaml").write_text(config_yaml)

    # Create some fake artifact directories that hyphy needs
    (outdir / "trees_iqtree").mkdir()
    (outdir / "codon_alignments").mkdir()
    (outdir / "summary").mkdir()

    return outdir


# ---------------------------------------------------------------------------
# InputRunResolver tests
# ---------------------------------------------------------------------------


class TestInputRunResolver:
    def test_valid_directory(self, tmp_path):
        from phylofoundry.provenance import InputRunResolver

        outdir = _make_fake_run(tmp_path)
        resolver = InputRunResolver(str(outdir))
        assert resolver.directory == str(outdir)

    def test_invalid_directory_raises(self, tmp_path):
        from phylofoundry.provenance import InputRunError, InputRunResolver

        with pytest.raises(InputRunError, match="does not exist"):
            InputRunResolver(str(tmp_path / "nonexistent"))

    def test_run_id_extracted(self, tmp_path):
        from phylofoundry.provenance import InputRunResolver

        outdir = _make_fake_run(tmp_path, run_id="my-run-123")
        resolver = InputRunResolver(str(outdir))
        assert resolver.run_id == "my-run-123"

    def test_provenance_record_structure(self, tmp_path):
        from phylofoundry.provenance import InputRunResolver

        outdir = _make_fake_run(tmp_path, run_id="prov-run")
        resolver = InputRunResolver(str(outdir))
        record = resolver.provenance_record()
        assert record["input_run_dir"] == str(outdir)
        assert record["input_run_id"] == "prov-run"
        assert "input_run_config" in record

    def test_validate_missing_artifacts_raises(self, tmp_path):
        from phylofoundry.provenance import InputRunError, InputRunResolver

        # Empty directory — no artifacts present
        empty = tmp_path / "empty"
        empty.mkdir()
        (empty / "config.yaml").write_text(
            '_run_meta:\n  run_id: "x"\noutput:\n  outdir: "x"\n'
        )
        resolver = InputRunResolver(str(empty))
        with pytest.raises(InputRunError, match="missing required artifacts"):
            resolver.validate_for_step("hyphy")

    def test_validate_passes_when_artifacts_present(self, tmp_path):
        from phylofoundry.provenance import InputRunResolver

        outdir = _make_fake_run(tmp_path)
        resolver = InputRunResolver(str(outdir))
        # Should not raise — trees_iqtree and codon_alignments exist
        resolver.validate_for_step("hyphy")

    def test_artifact_hashes_file_artifacts(self, tmp_path):
        from phylofoundry.provenance import InputRunResolver

        outdir = _make_fake_run(tmp_path)
        # Create a real file artifact (best_hits_tsv)
        summary = outdir / "summary"
        best_hits = summary / "best_hits.competitive.tsv"
        best_hits.write_text("hit1\thit2\n")

        resolver = InputRunResolver(str(outdir))
        hashes = resolver.artifact_hashes("taxonomy_integrate")
        # best_hits_tsv is a file — should have a hash
        assert "best_hits_tsv" in hashes
        assert len(hashes["best_hits_tsv"]) == 64  # SHA-256 hex

    def test_artifact_hashes_skips_directories(self, tmp_path):
        from phylofoundry.provenance import InputRunResolver

        outdir = _make_fake_run(tmp_path)
        resolver = InputRunResolver(str(outdir))
        hashes = resolver.artifact_hashes("hyphy")
        # trees_dir and codon_dir are directories — should not be hashed
        assert "trees_dir" not in hashes
        assert "codon_dir" not in hashes


# ---------------------------------------------------------------------------
# apply_input_run tests
# ---------------------------------------------------------------------------


class TestApplyInputRun:
    def test_no_op_when_not_set(self, tmp_path):
        from phylofoundry.provenance import apply_input_run

        cfg = {"_checkpoint": {}}
        apply_input_run(cfg)
        assert "_input_run" not in cfg

    def test_records_provenance(self, tmp_path):
        from phylofoundry.provenance import apply_input_run

        outdir = _make_fake_run(tmp_path, run_id="ir-run-42")
        cfg = {"_checkpoint": {"input_run": str(outdir)}, "workflow": {}}
        apply_input_run(cfg)
        assert cfg["_input_run"]["input_run_id"] == "ir-run-42"
        assert cfg["_input_paths"]["input_run_dir"] == str(outdir)

    def test_raises_on_missing_artifacts(self, tmp_path):
        from phylofoundry.provenance import apply_input_run, InputRunError

        empty = tmp_path / "empty"
        empty.mkdir()
        (empty / "config.yaml").write_text(
            '_run_meta:\n  run_id: "y"\noutput:\n  outdir: "y"\n'
        )
        cfg = {
            "_checkpoint": {"input_run": str(empty)},
            "workflow": {"start_at": "hyphy"},
        }
        with pytest.raises(InputRunError, match="missing required artifacts"):
            apply_input_run(cfg, step_name="hyphy")


# ---------------------------------------------------------------------------
# list_runs tests
# ---------------------------------------------------------------------------


class TestListRuns:
    def test_empty_directory(self, tmp_path, capsys):
        from phylofoundry.provenance import list_runs

        list_runs(str(tmp_path))
        out = capsys.readouterr().out
        assert "No PhyloFoundry runs found" in out

    def test_finds_run(self, tmp_path, capsys):
        from phylofoundry.provenance import list_runs

        _make_fake_run(tmp_path, run_id="found-run-99")
        list_runs(str(tmp_path))
        out = capsys.readouterr().out
        assert "found-run-99" in out

    def test_json_output(self, tmp_path, capsys):
        from phylofoundry.provenance import list_runs

        _make_fake_run(tmp_path, run_id="json-run")
        list_runs(str(tmp_path), as_json=True)
        out = capsys.readouterr().out
        data = json.loads(out)
        assert isinstance(data, list)
        run_ids = [r.get("run_id") for r in data]
        assert "json-run" in run_ids

    def test_nonexistent_directory(self, tmp_path, capsys):
        from phylofoundry.provenance import list_runs

        list_runs(str(tmp_path / "does_not_exist"))
        out = capsys.readouterr().out
        assert "No PhyloFoundry runs found" in out


# ---------------------------------------------------------------------------
# CLI: list-runs subcommand
# ---------------------------------------------------------------------------


class TestListRunsCLI:
    def test_list_runs_help(self):
        import subprocess
        import sys

        result = subprocess.run(
            [sys.executable, "-m", "phylofoundry.main", "list-runs", "--help"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "--json" in result.stdout

    def test_input_run_in_help(self):
        import subprocess
        import sys

        result = subprocess.run(
            [sys.executable, "-m", "phylofoundry.main", "hyphy", "--help"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "--input-run" in result.stdout


# ---------------------------------------------------------------------------
# Checkpoint provenance recording
# ---------------------------------------------------------------------------


class TestCheckpointProvenanceRecording:
    def test_input_run_in_ndjson(self, tmp_path):
        """start_run should include input_run provenance in the NDJSON event."""
        import json as _json
        from phylofoundry.checkpoint import Checkpointer, new_run_id

        outdir = tmp_path / "new_run"
        outdir.mkdir()
        prior = _make_fake_run(tmp_path, run_id="prior-run")

        cfg = {
            "_input_run": {
                "input_run_dir": str(prior),
                "input_run_id": "prior-run",
                "input_run_config": str(prior / "config.yaml"),
            },
            "output": {"outdir": str(outdir)},
        }

        ckpt = Checkpointer(str(outdir))
        run_id = new_run_id()
        ckpt.start_run(run_id, cfg)

        ndjson_path = outdir / "logs" / "pipeline.ndjson"
        with open(ndjson_path) as fh:
            events = [_json.loads(line) for line in fh if line.strip()]

        run_start = next(e for e in events if e.get("event") == "run_start")
        assert "input_run" in run_start
        assert run_start["input_run"]["input_run_id"] == "prior-run"

    def test_input_run_in_config_yaml(self, tmp_path):
        """start_run should include input_run provenance in the config snapshot."""
        from ruamel.yaml import YAML
        from phylofoundry.checkpoint import Checkpointer, new_run_id

        outdir = tmp_path / "new_run2"
        outdir.mkdir()
        prior = _make_fake_run(tmp_path, run_id="prior-run-2")

        cfg = {
            "_input_run": {
                "input_run_dir": str(prior),
                "input_run_id": "prior-run-2",
                "input_run_config": None,
            },
            "output": {"outdir": str(outdir)},
        }

        ckpt = Checkpointer(str(outdir))
        run_id = new_run_id()
        ckpt.start_run(run_id, cfg)

        yml = YAML()
        with open(outdir / "config.yaml") as fh:
            saved = yml.load(fh)

        assert saved["_run_meta"]["input_run"]["input_run_id"] == "prior-run-2"

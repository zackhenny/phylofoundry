"""Tests for the NDJSON+SQLite checkpoint and resume system."""

from __future__ import annotations

import json
import os
import sqlite3
import tempfile
from copy import deepcopy
from pathlib import Path

import pytest


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def tmp_outdir(tmp_path):
    d = tmp_path / "outdir"
    d.mkdir()
    return str(d)


@pytest.fixture
def checkpointer(tmp_outdir):
    from phylofoundry.checkpoint import Checkpointer

    return Checkpointer(tmp_outdir)


@pytest.fixture
def run_id():
    from phylofoundry.checkpoint import new_run_id

    return new_run_id()


# ---------------------------------------------------------------------------
# compute_fingerprint
# ---------------------------------------------------------------------------


class TestComputeFingerprint:
    def test_deterministic(self, tmp_path):
        from phylofoundry.checkpoint import compute_fingerprint

        f = tmp_path / "a.txt"
        f.write_bytes(b"hello")
        fp1 = compute_fingerprint({"k": 1}, [str(f)])
        fp2 = compute_fingerprint({"k": 1}, [str(f)])
        assert fp1 == fp2

    def test_different_params_differ(self):
        from phylofoundry.checkpoint import compute_fingerprint

        fp1 = compute_fingerprint({"k": 1}, [])
        fp2 = compute_fingerprint({"k": 2}, [])
        assert fp1 != fp2

    def test_different_file_content_differs(self, tmp_path):
        from phylofoundry.checkpoint import compute_fingerprint

        f = tmp_path / "b.txt"
        f.write_bytes(b"version1")
        fp1 = compute_fingerprint({}, [str(f)])
        f.write_bytes(b"version2")
        fp2 = compute_fingerprint({}, [str(f)])
        assert fp1 != fp2

    def test_missing_file_handled(self):
        from phylofoundry.checkpoint import compute_fingerprint

        # Should not raise; missing files are silently skipped
        fp = compute_fingerprint({}, ["/nonexistent/path/file.txt"])
        assert isinstance(fp, str)
        assert len(fp) == 40

    def test_returns_40_char_hex(self):
        from phylofoundry.checkpoint import compute_fingerprint

        fp = compute_fingerprint({"x": "y"}, [])
        assert len(fp) == 40
        assert all(c in "0123456789abcdef" for c in fp)


# ---------------------------------------------------------------------------
# Checkpointer — run lifecycle
# ---------------------------------------------------------------------------


class TestCheckpointerRunLifecycle:
    def test_start_run_writes_ndjson(self, checkpointer, run_id, tmp_outdir):
        checkpointer.start_run(run_id, {"inputs": {}, "_test": True})
        assert os.path.isfile(checkpointer.ndjson_path)
        lines = Path(checkpointer.ndjson_path).read_text().splitlines()
        records = [json.loads(ln) for ln in lines if ln]
        assert any(r.get("event") == "run_start" and r.get("run_id") == run_id for r in records)

    def test_start_run_writes_sqlite(self, checkpointer, run_id):
        checkpointer.start_run(run_id, {})
        conn = sqlite3.connect(checkpointer.db_path)
        row = conn.execute("SELECT run_id FROM runs WHERE run_id=?", (run_id,)).fetchone()
        conn.close()
        assert row is not None
        assert row[0] == run_id

    def test_start_run_writes_config_snapshot(self, checkpointer, run_id, tmp_outdir):
        cfg = {"inputs": {"faa_dir": "/tmp/faa"}}
        checkpointer.start_run(run_id, cfg)
        snapshot_path = os.path.join(tmp_outdir, "config.yaml")
        assert os.path.isfile(snapshot_path)

    def test_config_snapshot_contains_run_meta(self, checkpointer, run_id, tmp_outdir):
        checkpointer.start_run(run_id, {})
        snapshot_path = os.path.join(tmp_outdir, "config.yaml")
        from ruamel.yaml import YAML
        yml = YAML()
        with open(snapshot_path) as fh:
            data = yml.load(fh)
        assert "_run_meta" in data
        assert data["_run_meta"]["run_id"] == run_id

    def test_end_run_writes_ndjson(self, checkpointer, run_id):
        checkpointer.start_run(run_id, {})
        checkpointer.end_run(run_id, success=True)
        lines = Path(checkpointer.ndjson_path).read_text().splitlines()
        records = [json.loads(ln) for ln in lines if ln]
        assert any(r.get("event") == "run_end" for r in records)

    def test_end_run_sets_end_time(self, checkpointer, run_id):
        checkpointer.start_run(run_id, {})
        checkpointer.end_run(run_id, success=True)
        conn = sqlite3.connect(checkpointer.db_path)
        row = conn.execute("SELECT end_time FROM runs WHERE run_id=?", (run_id,)).fetchone()
        conn.close()
        assert row[0] is not None


# ---------------------------------------------------------------------------
# Checkpointer — step lifecycle
# ---------------------------------------------------------------------------


class TestCheckpointerStepLifecycle:
    def _setup(self, checkpointer, run_id):
        checkpointer.start_run(run_id, {})

    def test_step_pending_to_success(self, checkpointer, run_id):
        self._setup(checkpointer, run_id)
        checkpointer.record_step_pending(run_id, "prep", params={"k": "v"})
        checkpointer.record_step_running(run_id, "prep")
        checkpointer.record_step_success(run_id, "prep", input_fingerprint="abc123")

        last = checkpointer.latest_step(run_id, "prep")
        assert last is not None
        assert last["state"] == "success"
        assert last["input_fingerprint"] == "abc123"

    def test_step_failure_records_error(self, checkpointer, run_id):
        self._setup(checkpointer, run_id)
        checkpointer.record_step_pending(run_id, "hmmer")
        checkpointer.record_step_running(run_id, "hmmer")
        checkpointer.record_step_failure(run_id, "hmmer", error="something broke", trace="traceback here")

        last = checkpointer.latest_step(run_id, "hmmer")
        assert last["state"] == "failed"
        assert last["error"] == "something broke"
        assert last["trace"] == "traceback here"

    def test_step_skipped(self, checkpointer, run_id):
        self._setup(checkpointer, run_id)
        checkpointer.record_step_skipped(run_id, "phylo", reason="fingerprint match")
        lines = Path(checkpointer.ndjson_path).read_text().splitlines()
        records = [json.loads(ln) for ln in lines if ln]
        skip_records = [r for r in records if r.get("event") == "step_skipped"]
        assert any(r["step"] == "phylo" for r in skip_records)

    def test_latest_step_returns_none_when_absent(self, checkpointer, run_id):
        self._setup(checkpointer, run_id)
        result = checkpointer.latest_step(run_id, "nonexistent_step")
        assert result is None

    def test_all_step_states(self, checkpointer, run_id):
        self._setup(checkpointer, run_id)
        checkpointer.record_step_pending(run_id, "prep")
        checkpointer.record_step_running(run_id, "prep")
        checkpointer.record_step_success(run_id, "prep")
        checkpointer.record_step_pending(run_id, "hmmer")
        checkpointer.record_step_running(run_id, "hmmer")
        checkpointer.record_step_failure(run_id, "hmmer", error="oops")

        states = checkpointer.all_step_states(run_id)
        assert states["prep"] == "success"
        assert states["hmmer"] == "failed"

    def test_ndjson_contains_all_transitions(self, checkpointer, run_id):
        self._setup(checkpointer, run_id)
        checkpointer.record_step_pending(run_id, "extract")
        checkpointer.record_step_running(run_id, "extract")
        checkpointer.record_step_success(run_id, "extract")

        lines = Path(checkpointer.ndjson_path).read_text().splitlines()
        events = [json.loads(ln)["event"] for ln in lines if ln]
        assert "step_pending" in events
        assert "step_running" in events
        assert "step_success" in events


# ---------------------------------------------------------------------------
# Checkpointer — artifact recording
# ---------------------------------------------------------------------------


class TestCheckpointerArtifacts:
    def test_artifacts_recorded_on_success(self, checkpointer, run_id, tmp_path):
        checkpointer.start_run(run_id, {})
        # Create a real file to hash
        f = tmp_path / "result.tsv"
        f.write_text("col1\tcol2\nA\tB\n")

        checkpointer.record_step_pending(run_id, "extract")
        checkpointer.record_step_running(run_id, "extract")
        checkpointer.record_step_success(run_id, "extract", outputs=[str(f)])

        conn = sqlite3.connect(checkpointer.db_path)
        rows = conn.execute(
            "SELECT path, sha256 FROM artifacts WHERE run_id=? AND step_name=?",
            (run_id, "extract"),
        ).fetchall()
        conn.close()
        assert len(rows) == 1
        assert rows[0][0] == str(f)
        assert len(rows[0][1]) == 64  # full SHA-256 hex

    def test_missing_artifact_file_skipped(self, checkpointer, run_id):
        checkpointer.start_run(run_id, {})
        checkpointer.record_step_pending(run_id, "phylo")
        checkpointer.record_step_running(run_id, "phylo")
        # Passing a non-existent path should not raise
        checkpointer.record_step_success(run_id, "phylo", outputs=["/no/such/file.treefile"])
        conn = sqlite3.connect(checkpointer.db_path)
        count = conn.execute(
            "SELECT COUNT(*) FROM artifacts WHERE run_id=?", (run_id,)
        ).fetchone()[0]
        conn.close()
        # No artifact should be recorded for missing file
        assert count == 0


# ---------------------------------------------------------------------------
# load_checkpoint_meta / locate_saved_config
# ---------------------------------------------------------------------------


class TestLoadCheckpointMeta:
    def test_locates_yaml_config(self, tmp_outdir):
        from phylofoundry.checkpoint import locate_saved_config

        yaml_path = os.path.join(tmp_outdir, "config.yaml")
        Path(yaml_path).write_text("_run_meta:\n  run_id: run-xyz\n")
        result = locate_saved_config(tmp_outdir)
        assert result == yaml_path

    def test_falls_back_to_json(self, tmp_outdir):
        from phylofoundry.checkpoint import locate_saved_config

        json_path = os.path.join(tmp_outdir, "config.json")
        Path(json_path).write_text(json.dumps({"_run_meta": {"run_id": "run-abc"}}))
        result = locate_saved_config(tmp_outdir)
        assert result == json_path

    def test_yaml_preferred_over_json(self, tmp_outdir):
        from phylofoundry.checkpoint import locate_saved_config

        yaml_path = os.path.join(tmp_outdir, "config.yaml")
        json_path = os.path.join(tmp_outdir, "config.json")
        Path(yaml_path).write_text("_run_meta:\n  run_id: run-yaml\n")
        Path(json_path).write_text(json.dumps({"_run_meta": {"run_id": "run-json"}}))
        result = locate_saved_config(tmp_outdir)
        assert result == yaml_path

    def test_returns_none_when_absent(self, tmp_outdir):
        from phylofoundry.checkpoint import locate_saved_config

        result = locate_saved_config(tmp_outdir)
        assert result is None

    def test_load_checkpoint_meta_from_directory(self, tmp_outdir):
        from phylofoundry.checkpoint import Checkpointer, load_checkpoint_meta

        ckptr = Checkpointer(tmp_outdir)
        run_id = "test-run-123"
        ckptr.start_run(run_id, {})

        meta = load_checkpoint_meta(tmp_outdir)
        assert meta is not None
        assert meta.run_id == run_id
        assert meta.outdir == tmp_outdir

    def test_load_checkpoint_meta_from_yaml_path(self, tmp_outdir):
        from phylofoundry.checkpoint import Checkpointer, load_checkpoint_meta

        ckptr = Checkpointer(tmp_outdir)
        run_id = "direct-path-test"
        ckptr.start_run(run_id, {})

        yaml_path = os.path.join(tmp_outdir, "config.yaml")
        meta = load_checkpoint_meta(yaml_path)
        assert meta is not None
        assert meta.run_id == run_id

    def test_load_returns_none_for_missing_path(self):
        from phylofoundry.checkpoint import load_checkpoint_meta

        result = load_checkpoint_meta("/nonexistent/dir/that/does/not/exist")
        assert result is None


# ---------------------------------------------------------------------------
# build_resume_plan
# ---------------------------------------------------------------------------


class TestBuildResumePlan:
    def test_skip_on_matching_fingerprint(self, checkpointer, run_id):
        checkpointer.start_run(run_id, {})
        checkpointer.record_step_pending(run_id, "prep", input_fingerprint="fp-prep")
        checkpointer.record_step_running(run_id, "prep")
        checkpointer.record_step_success(run_id, "prep", input_fingerprint="fp-prep")

        from phylofoundry.checkpoint import build_resume_plan

        plan = build_resume_plan(checkpointer, run_id, ["prep"], {"prep": "fp-prep"})
        assert plan["prep"] == "skip"

    def test_run_on_mismatched_fingerprint(self, checkpointer, run_id):
        checkpointer.start_run(run_id, {})
        checkpointer.record_step_pending(run_id, "prep", input_fingerprint="fp-old")
        checkpointer.record_step_running(run_id, "prep")
        checkpointer.record_step_success(run_id, "prep", input_fingerprint="fp-old")

        from phylofoundry.checkpoint import build_resume_plan

        plan = build_resume_plan(checkpointer, run_id, ["prep"], {"prep": "fp-new"})
        assert plan["prep"] == "run"

    def test_run_when_no_prior_record(self, checkpointer, run_id):
        checkpointer.start_run(run_id, {})

        from phylofoundry.checkpoint import build_resume_plan

        plan = build_resume_plan(checkpointer, run_id, ["prep"], {"prep": "fp-x"})
        assert plan["prep"] == "run"

    def test_resume_on_interrupted_running(self, checkpointer, run_id):
        checkpointer.start_run(run_id, {})
        checkpointer.record_step_pending(run_id, "hmmer")
        checkpointer.record_step_running(run_id, "hmmer")
        # Never recorded as success/failure — simulating an interrupted run

        from phylofoundry.checkpoint import build_resume_plan

        plan = build_resume_plan(checkpointer, run_id, ["hmmer"], {"hmmer": "fp-h"})
        assert plan["hmmer"] == "resume"

    def test_force_always_runs(self, checkpointer, run_id):
        checkpointer.start_run(run_id, {})
        checkpointer.record_step_pending(run_id, "prep", input_fingerprint="fp-match")
        checkpointer.record_step_running(run_id, "prep")
        checkpointer.record_step_success(run_id, "prep", input_fingerprint="fp-match")

        from phylofoundry.checkpoint import build_resume_plan

        plan = build_resume_plan(
            checkpointer, run_id, ["prep"], {"prep": "fp-match"}, force=True
        )
        assert plan["prep"] == "run"

    def test_failed_step_reruns(self, checkpointer, run_id):
        checkpointer.start_run(run_id, {})
        checkpointer.record_step_pending(run_id, "extract")
        checkpointer.record_step_running(run_id, "extract")
        checkpointer.record_step_failure(run_id, "extract", error="boom")

        from phylofoundry.checkpoint import build_resume_plan

        plan = build_resume_plan(checkpointer, run_id, ["extract"], {"extract": "fp-e"})
        assert plan["extract"] == "run"


# ---------------------------------------------------------------------------
# NDJSON atomicity — _atomic_append_ndjson
# ---------------------------------------------------------------------------


class TestAtomicAppend:
    def test_multiple_appends_accumulate(self, tmp_path):
        from phylofoundry.checkpoint import _atomic_append_ndjson

        p = str(tmp_path / "test.ndjson")
        for i in range(5):
            _atomic_append_ndjson(p, {"seq": i})
        lines = Path(p).read_text().splitlines()
        assert len(lines) == 5
        for i, ln in enumerate(lines):
            assert json.loads(ln)["seq"] == i

    def test_idempotent_on_existing_file(self, tmp_path):
        from phylofoundry.checkpoint import _atomic_append_ndjson

        p = str(tmp_path / "existing.ndjson")
        Path(p).write_text(json.dumps({"existing": True}) + "\n")
        _atomic_append_ndjson(p, {"new": True})
        lines = Path(p).read_text().splitlines()
        assert len(lines) == 2
        assert json.loads(lines[0])["existing"] is True
        assert json.loads(lines[1])["new"] is True


# ---------------------------------------------------------------------------
# CLI flag presence
# ---------------------------------------------------------------------------


class TestCLIResumeFlags:
    """Verify that --resume, --resume-from, --no-resume are registered."""

    def test_resume_flag_parsed(self):
        from phylofoundry.main import _build_parser

        ap = _build_parser()
        args = ap.parse_args(["run", "--resume"])
        assert args.resume is True

    def test_no_resume_flag_parsed(self):
        from phylofoundry.main import _build_parser

        ap = _build_parser()
        args = ap.parse_args(["run", "--no-resume"])
        assert args.no_resume is True

    def test_resume_from_flag_parsed(self):
        from phylofoundry.main import _build_parser

        ap = _build_parser()
        args = ap.parse_args(["run", "--resume-from", "/tmp/previous_run"])
        assert args.resume_from == "/tmp/previous_run"

    def test_force_still_present(self):
        from phylofoundry.main import _build_parser

        ap = _build_parser()
        args = ap.parse_args(["run", "--force"])
        assert args.force is True

    def test_resume_flags_on_embed_subcommand(self):
        from phylofoundry.main import _build_parser

        ap = _build_parser()
        args = ap.parse_args(["embed", "--resume"])
        assert args.resume is True

    def test_defaults_are_false(self):
        from phylofoundry.main import _build_parser

        ap = _build_parser()
        args = ap.parse_args(["run"])
        assert args.resume is False
        assert args.no_resume is False
        assert args.resume_from is None


# ---------------------------------------------------------------------------
# SubtaskProgressTracker
# ---------------------------------------------------------------------------


class TestSubtaskProgressTracker:
    """Tests for per-subtask NDJSON progress tracking within pipeline steps."""

    def test_fresh_tracker_nothing_done(self, tmp_path):
        from phylofoundry.checkpoint import SubtaskProgressTracker

        log = str(tmp_path / ".progress.ndjson")
        tracker = SubtaskProgressTracker(log, resume=True)
        assert not tracker.is_done("hmm_A")
        assert tracker.completed_count() == 0

    def test_record_and_query_success(self, tmp_path):
        from phylofoundry.checkpoint import SubtaskProgressTracker

        log = str(tmp_path / ".progress.ndjson")
        tracker = SubtaskProgressTracker(log, resume=False)
        tracker.record_running("hmm_A")
        assert not tracker.is_done("hmm_A")
        tracker.record_success("hmm_A")
        assert tracker.is_done("hmm_A")
        assert tracker.completed_count() == 1

    def test_record_failure_not_marked_done(self, tmp_path):
        from phylofoundry.checkpoint import SubtaskProgressTracker

        log = str(tmp_path / ".progress.ndjson")
        tracker = SubtaskProgressTracker(log, resume=False)
        tracker.record_running("hmm_B")
        tracker.record_failure("hmm_B", error="OOM")
        assert not tracker.is_done("hmm_B")

    def test_record_skipped_marks_done(self, tmp_path):
        from phylofoundry.checkpoint import SubtaskProgressTracker

        log = str(tmp_path / ".progress.ndjson")
        tracker = SubtaskProgressTracker(log, resume=False)
        tracker.record_skipped("hmm_C", reason="output exists")
        assert tracker.is_done("hmm_C")

    def test_persist_and_reload_on_resume(self, tmp_path):
        from phylofoundry.checkpoint import SubtaskProgressTracker

        log = str(tmp_path / ".progress.ndjson")

        # First run: complete hmm_A and hmm_B; hmm_C interrupted
        t1 = SubtaskProgressTracker(log, resume=False)
        t1.record_running("hmm_A")
        t1.record_success("hmm_A")
        t1.record_running("hmm_B")
        t1.record_success("hmm_B")
        t1.record_running("hmm_C")
        # hmm_C interrupted — no success record

        # Second run: resume=True should see A and B done, C not
        t2 = SubtaskProgressTracker(log, resume=True)
        assert t2.is_done("hmm_A")
        assert t2.is_done("hmm_B")
        assert not t2.is_done("hmm_C")
        assert t2.completed_count() == 2

    def test_no_resume_ignores_existing_log(self, tmp_path):
        from phylofoundry.checkpoint import SubtaskProgressTracker

        log = str(tmp_path / ".progress.ndjson")

        # Pre-populate with a successful subtask
        t1 = SubtaskProgressTracker(log, resume=False)
        t1.record_success("hmm_X")

        # resume=False → should NOT see hmm_X as done
        t2 = SubtaskProgressTracker(log, resume=False)
        assert not t2.is_done("hmm_X")

    def test_cleanup_removes_log(self, tmp_path):
        from phylofoundry.checkpoint import SubtaskProgressTracker

        log = str(tmp_path / ".progress.ndjson")
        tracker = SubtaskProgressTracker(log, resume=False)
        tracker.record_success("hmm_A")
        assert os.path.exists(log)
        tracker.cleanup()
        assert not os.path.exists(log)

    def test_cleanup_idempotent_no_file(self, tmp_path):
        from phylofoundry.checkpoint import SubtaskProgressTracker

        log = str(tmp_path / ".progress.ndjson")
        tracker = SubtaskProgressTracker(log, resume=False)
        # No records written — file may not exist
        tracker.cleanup()  # should not raise

    def test_ndjson_format_correct(self, tmp_path):
        from phylofoundry.checkpoint import SubtaskProgressTracker

        log = str(tmp_path / ".progress.ndjson")
        tracker = SubtaskProgressTracker(log, resume=False)
        tracker.record_running("hmm_A")
        tracker.record_success("hmm_A", n_sequences=42)

        lines = Path(log).read_text().splitlines()
        assert len(lines) == 2
        rec_run = json.loads(lines[0])
        rec_ok = json.loads(lines[1])
        assert rec_run["subtask"] == "hmm_A"
        assert rec_run["state"] == "running"
        assert rec_ok["subtask"] == "hmm_A"
        assert rec_ok["state"] == "success"
        assert rec_ok["n_sequences"] == 42
        assert "timestamp" in rec_ok

    def test_multiple_subtasks(self, tmp_path):
        from phylofoundry.checkpoint import SubtaskProgressTracker

        log = str(tmp_path / ".progress.ndjson")
        tracker = SubtaskProgressTracker(log, resume=False)
        hmms = [f"hmm_{i}" for i in range(10)]
        for h in hmms:
            tracker.record_running(h)
            tracker.record_success(h)

        assert tracker.completed_count() == 10
        for h in hmms:
            assert tracker.is_done(h)

    def test_corrupted_log_line_tolerated(self, tmp_path):
        from phylofoundry.checkpoint import SubtaskProgressTracker

        log = tmp_path / ".progress.ndjson"
        # Write one valid line and one corrupted line
        log.write_text(
            '{"subtask":"hmm_A","state":"success","timestamp":"2024-01-01T00:00:00Z"}\n'
            'NOT VALID JSON !!!\n'
            '{"subtask":"hmm_B","state":"success","timestamp":"2024-01-01T00:00:00Z"}\n'
        )
        tracker = SubtaskProgressTracker(str(log), resume=True)
        assert tracker.is_done("hmm_A")
        assert tracker.is_done("hmm_B")
        assert tracker.completed_count() == 2

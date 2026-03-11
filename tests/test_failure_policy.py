"""Tests for failure_policy.py and pipeline failure handling integration."""

import os
import csv
import tempfile
import pytest

from phylofoundry.failure_policy import (
    apply_failure_policy,
    build_reverse_dependency_map,
    get_impacted_steps,
)
from phylofoundry.execution_schema import PlannedStep, StepState
from phylofoundry.logging_utils import (
    append_pipeline_log,
    initialize_step_status,
    pipeline_log_path,
    step_status_path,
    update_step_status,
)
from phylofoundry.execution_schema import ExecutionPlan


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def _make_plan_steps(names_states: dict) -> list:
    """Build a list of PlannedStep from {name: StepState}."""
    return [PlannedStep(name=n, state=s) for n, s in names_states.items()]


# ---------------------------------------------------------------------------
# build_reverse_dependency_map
# ---------------------------------------------------------------------------

class TestBuildReverseDependencyMap:
    def test_known_steps_returns_map(self):
        """Registry-backed steps produce a populated reverse map."""
        rev = build_reverse_dependency_map(["prep", "hmmer", "extract"])
        # hmmer depends on prep → prep should have hmmer as downstream
        assert "hmmer" in rev["prep"]
        # extract depends on hmmer → hmmer should have extract as downstream
        assert "extract" in rev["hmmer"]
        # extract has no downstream in this subset
        assert rev["extract"] == []

    def test_unknown_step_skipped(self):
        """Steps not in REGISTRY produce empty downstream lists (no KeyError)."""
        rev = build_reverse_dependency_map(["nonexistent_step"])
        assert rev == {"nonexistent_step": []}

    def test_all_keys_present(self):
        steps = ["prep", "hmmer", "extract", "phylo"]
        rev = build_reverse_dependency_map(steps)
        for name in steps:
            assert name in rev


# ---------------------------------------------------------------------------
# get_impacted_steps
# ---------------------------------------------------------------------------

class TestGetImpactedSteps:
    def test_direct_downstream(self):
        rev = build_reverse_dependency_map(["prep", "hmmer"])
        impacted = get_impacted_steps("prep", rev)
        assert "hmmer" in impacted
        assert "prep" not in impacted  # failed step not in impacted set

    def test_transitive_downstream(self):
        # prep → hmmer → extract
        rev = build_reverse_dependency_map(["prep", "hmmer", "extract"])
        impacted = get_impacted_steps("prep", rev)
        assert "hmmer" in impacted
        assert "extract" in impacted

    def test_unrelated_not_impacted(self):
        # synteny depends on hmmer, not on phylo
        rev = build_reverse_dependency_map(
            ["prep", "hmmer", "extract", "phylo", "synteny"]
        )
        impacted = get_impacted_steps("phylo", rev)
        # synteny is not downstream of phylo
        assert "synteny" not in impacted

    def test_no_downstream(self):
        rev = {"step_a": [], "step_b": []}
        assert get_impacted_steps("step_a", rev) == set()

    def test_cycle_guard(self):
        """Cycles in the reverse map should not cause infinite loops."""
        rev = {"a": ["b"], "b": ["a"]}
        impacted = get_impacted_steps("a", rev)
        assert "b" in impacted


# ---------------------------------------------------------------------------
# apply_failure_policy
# ---------------------------------------------------------------------------

class TestApplyFailurePolicy:
    def test_pending_downstream_blocked(self):
        steps = _make_plan_steps({
            "prep": StepState.FAILED,
            "hmmer": StepState.PENDING,
            "extract": StepState.PENDING,
        })
        rev = build_reverse_dependency_map(["prep", "hmmer", "extract"])
        apply_failure_policy(steps, "prep", rev)
        states = {s.name: s.state for s in steps}
        assert states["hmmer"] == StepState.BLOCKED
        assert states["extract"] == StepState.BLOCKED

    def test_disabled_step_not_blocked(self):
        """Disabled steps must stay disabled, not transition to blocked."""
        steps = _make_plan_steps({
            "prep": StepState.FAILED,
            "hmmer": StepState.DISABLED,
        })
        rev = build_reverse_dependency_map(["prep", "hmmer"])
        apply_failure_policy(steps, "prep", rev)
        states = {s.name: s.state for s in steps}
        assert states["hmmer"] == StepState.DISABLED

    def test_success_step_not_overwritten(self):
        """Already-succeeded steps must not be changed to blocked."""
        steps = _make_plan_steps({
            "prep": StepState.FAILED,
            "hmmer": StepState.SUCCESS,
        })
        rev = build_reverse_dependency_map(["prep", "hmmer"])
        apply_failure_policy(steps, "prep", rev)
        states = {s.name: s.state for s in steps}
        assert states["hmmer"] == StepState.SUCCESS

    def test_unrelated_step_unaffected(self):
        """Steps not in the downstream graph of the failed step are unchanged."""
        steps = _make_plan_steps({
            "prep": StepState.FAILED,
            "hmmer": StepState.PENDING,
            "synteny": StepState.PENDING,  # synteny depends on hmmer, not prep
        })
        # Minimal reverse map: only prep → hmmer (no synteny edge)
        rev = {"prep": ["hmmer"], "hmmer": [], "synteny": []}
        apply_failure_policy(steps, "prep", rev)
        states = {s.name: s.state for s in steps}
        assert states["hmmer"] == StepState.BLOCKED
        assert states["synteny"] == StepState.PENDING  # unrelated


# ---------------------------------------------------------------------------
# append_pipeline_log
# ---------------------------------------------------------------------------

class TestAppendPipelineLog:
    def test_creates_log_file(self, tmp_path):
        logs_dir = str(tmp_path / "logs")
        os.makedirs(logs_dir, exist_ok=True)
        append_pipeline_log(logs_dir, "test message")
        log_path = pipeline_log_path(logs_dir)
        assert os.path.exists(log_path)
        with open(log_path) as fh:
            content = fh.read()
        assert "test message" in content

    def test_appends_multiple_messages(self, tmp_path):
        logs_dir = str(tmp_path / "logs")
        os.makedirs(logs_dir, exist_ok=True)
        append_pipeline_log(logs_dir, "first")
        append_pipeline_log(logs_dir, "second")
        with open(pipeline_log_path(logs_dir)) as fh:
            content = fh.read()
        assert "first" in content
        assert "second" in content

    def test_timestamp_prefix_present(self, tmp_path):
        logs_dir = str(tmp_path / "logs")
        os.makedirs(logs_dir, exist_ok=True)
        append_pipeline_log(logs_dir, "msg")
        with open(pipeline_log_path(logs_dir)) as fh:
            line = fh.readline()
        # Timestamps are wrapped in [...], e.g. [2026-03-11T22:59:47Z]
        assert line.startswith("[")


# ---------------------------------------------------------------------------
# Integration: step_status.tsv reflects failure and blocked states
# ---------------------------------------------------------------------------

class TestStepStatusIntegration:
    def _read_status(self, path):
        rows = {}
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                rows[row["step"]] = row["state"]
        return rows

    def test_status_updated_to_failed_then_blocked(self, tmp_path):
        """Simulates what pipeline does: mark failed, then block downstream."""
        logs_dir = str(tmp_path / "logs")
        os.makedirs(logs_dir, exist_ok=True)

        plan_steps = _make_plan_steps({
            "prep": StepState.PENDING,
            "hmmer": StepState.PENDING,
            "extract": StepState.PENDING,
        })
        plan = ExecutionPlan(
            steps=plan_steps,
            logs_dir=logs_dir,
            enabled_step_names=["prep", "hmmer", "extract"],
            disabled_step_names=[],
        )
        status_path = step_status_path(logs_dir)
        initialize_step_status(plan, status_path)

        # Simulate prep running then failing
        update_step_status(status_path, "prep", "running")
        update_step_status(status_path, "prep", StepState.FAILED.value)

        # Build reverse map and apply failure policy
        rev = build_reverse_dependency_map(["prep", "hmmer", "extract"])
        apply_failure_policy(plan_steps, "prep", rev)

        # Update TSV for blocked steps
        for s in plan_steps:
            if s.state == StepState.BLOCKED:
                update_step_status(status_path, s.name, StepState.BLOCKED.value)

        states = self._read_status(status_path)
        assert states["prep"] == "failed"
        assert states["hmmer"] == "blocked"
        assert states["extract"] == "blocked"

    def test_unrelated_step_remains_pending(self, tmp_path):
        """Steps not downstream of the failing step stay pending."""
        logs_dir = str(tmp_path / "logs")
        os.makedirs(logs_dir, exist_ok=True)

        plan_steps = _make_plan_steps({
            "prep": StepState.PENDING,
            "hmmer": StepState.PENDING,
            "synteny": StepState.PENDING,
        })
        plan = ExecutionPlan(
            steps=plan_steps,
            logs_dir=logs_dir,
            enabled_step_names=["prep", "hmmer", "synteny"],
            disabled_step_names=[],
        )
        status_path = step_status_path(logs_dir)
        initialize_step_status(plan, status_path)

        update_step_status(status_path, "prep", StepState.FAILED.value)
        # Minimal rev map: prep blocks hmmer only (synteny depends on hmmer per registry
        # but we isolate the test using a custom map)
        rev = {"prep": ["hmmer"], "hmmer": [], "synteny": []}
        apply_failure_policy(plan_steps, "prep", rev)
        for s in plan_steps:
            if s.state == StepState.BLOCKED:
                update_step_status(status_path, s.name, StepState.BLOCKED.value)

        states = self._read_status(status_path)
        assert states["prep"] == "failed"
        assert states["hmmer"] == "blocked"
        assert states["synteny"] == "pending"  # unrelated branch unaffected

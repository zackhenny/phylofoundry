"""Logging utilities for PhyloFoundry.

Provides helpers to create the ``logs/`` directory and write structured
outputs such as ``execution_plan.json`` and ``step_status.tsv``.
"""

from __future__ import annotations

import csv
import json
import os
from typing import TYPE_CHECKING, List

if TYPE_CHECKING:
    from .execution_schema import ExecutionPlan

# ---------------------------------------------------------------------------
# Path helpers
# ---------------------------------------------------------------------------


def ensure_logs_dir(outdir: str) -> str:
    """Create ``outdir/logs/`` if it does not exist and return its path."""
    logs_dir = os.path.join(outdir, "logs")
    os.makedirs(logs_dir, exist_ok=True)
    return logs_dir


def pipeline_log_path(logs_dir: str) -> str:
    """Return the path for the main pipeline log file."""
    return os.path.join(logs_dir, "pipeline.log")


def step_log_path(logs_dir: str, step_name: str) -> str:
    """Return the path for an individual step log file."""
    return os.path.join(logs_dir, f"{step_name}.log")


def execution_plan_path(logs_dir: str) -> str:
    """Return the path for the execution plan JSON file."""
    return os.path.join(logs_dir, "execution_plan.json")


def step_status_path(logs_dir: str) -> str:
    """Return the path for the step status TSV file."""
    return os.path.join(logs_dir, "step_status.tsv")


# ---------------------------------------------------------------------------
# Write helpers
# ---------------------------------------------------------------------------


def write_execution_plan(plan: "ExecutionPlan", path: str) -> None:
    """Serialize *plan* to a JSON file at *path*.

    The parent directory is created if it does not exist.
    """
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        json.dump(plan.as_dict(), fh, indent=2)


def initialize_step_status(plan: "ExecutionPlan", path: str) -> None:
    """Write an initial ``step_status.tsv`` from the execution plan.

    Columns: ``step``, ``state``, ``log_path``.

    States reflect the initial plan — enabled steps are ``pending`` and
    disabled steps are ``disabled``.

    The parent directory is created if it does not exist.
    """
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["step", "state", "log_path"])
        for step in plan.steps:
            writer.writerow([step.name, step.state.value, step.log_path or ""])


def append_pipeline_log(logs_dir: str, message: str) -> None:
    """Append a timestamped *message* to ``logs/pipeline.log``.

    The file is created if it does not exist.  Each message is written on its
    own line, prefixed with an ISO-8601 timestamp.
    """
    import datetime

    path = pipeline_log_path(logs_dir)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    timestamp = datetime.datetime.now(datetime.timezone.utc).strftime(
        "%Y-%m-%dT%H:%M:%SZ"
    )
    with open(path, "a") as fh:
        fh.write(f"[{timestamp}] {message}\n")


def append_step_log(logs_dir: str, step_name: str, message: str) -> None:
    """Append a timestamped *message* to ``logs/<step>.log``.

    The file is created if it does not exist.  Each message is written on its
    own line, prefixed with an ISO-8601 timestamp.
    """
    import datetime

    path = step_log_path(logs_dir, step_name)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    timestamp = datetime.datetime.now(datetime.timezone.utc).strftime(
        "%Y-%m-%dT%H:%M:%SZ"
    )
    with open(path, "a") as fh:
        fh.write(f"[{timestamp}] {message}\n")


def update_step_status(path: str, step_name: str, new_state: str) -> None:
    """Update the state for a single step in an existing ``step_status.tsv``.

    Reads the entire file, updates the matching row in-place, and rewrites
    the file. No-ops silently if the file does not exist.
    """
    if not os.path.exists(path):
        return
    rows: List[List[str]] = []
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        rows = list(reader)
    for i, row in enumerate(rows):
        if i == 0:
            continue  # skip header
        if row and row[0] == step_name:
            row[1] = new_state
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerows(rows)

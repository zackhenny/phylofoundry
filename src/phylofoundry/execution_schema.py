"""Execution plan schema for PhyloFoundry.

Defines dataclasses for describing an execution plan, per-step states,
and execution results. All types are designed to be serialization-friendly
(plain dicts/strings) for writing to JSON and TSV outputs.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import List, Optional


class StepState(str, Enum):
    """The state of a planned or completed pipeline step."""

    PENDING = "pending"
    """Scheduled to run but not yet started."""
    RUNNING = "running"
    """Currently executing."""
    SUCCESS = "success"
    """Completed successfully."""
    FAILED = "failed"
    """Completed with an error."""
    BLOCKED = "blocked"
    """Blocked because an upstream dependency failed."""
    DISABLED = "disabled"
    """Not enabled in this run (config flag or out-of-range)."""
    SKIPPED = "skipped"
    """Skipped because outputs already exist and force=False."""


@dataclass
class PlannedStep:
    """A single step as it appears in an execution plan."""

    name: str
    state: StepState
    dependencies: List[str] = field(default_factory=list)
    required_tools: List[str] = field(default_factory=list)
    log_path: Optional[str] = None
    notes: str = ""


@dataclass
class ExecutionPlan:
    """Complete plan for a pipeline run.

    Built by :func:`~phylofoundry.execution_planner.build_execution_plan` and
    written to ``logs/execution_plan.json`` before execution begins.
    """

    steps: List[PlannedStep] = field(default_factory=list)
    logs_dir: str = ""
    enabled_step_names: List[str] = field(default_factory=list)
    disabled_step_names: List[str] = field(default_factory=list)
    required_tools: List[str] = field(default_factory=list)

    def as_dict(self) -> dict:
        """Serialize to a plain dict suitable for JSON output."""
        return {
            "logs_dir": self.logs_dir,
            "enabled_steps": self.enabled_step_names,
            "disabled_steps": self.disabled_step_names,
            "required_tools": self.required_tools,
            "steps": [
                {
                    "name": s.name,
                    "state": s.state.value,
                    "dependencies": s.dependencies,
                    "required_tools": s.required_tools,
                    "log_path": s.log_path,
                    "notes": s.notes,
                }
                for s in self.steps
            ],
        }


@dataclass
class StepResult:
    """Result of executing a single pipeline step."""

    name: str
    state: StepState
    message: str = ""


@dataclass
class ExecutionResult:
    """Aggregate result of a complete pipeline run."""

    step_results: List[StepResult] = field(default_factory=list)

    def overall_success(self) -> bool:
        """Return ``True`` if no step finished in a FAILED state."""
        failed = [r for r in self.step_results if r.state == StepState.FAILED]
        return len(failed) == 0

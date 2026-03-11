"""Execution planner for PhyloFoundry.

Builds an :class:`~phylofoundry.execution_schema.ExecutionPlan` from a
resolved configuration dict and an optional ordered list of workflow steps.

The plan captures:

- which steps are enabled vs disabled for this run
- required (non-optional) external tools across enabled steps
- per-step log file paths
- step dependency edges from the registry
"""

from __future__ import annotations

import os
from typing import Dict, List, Optional

from .constants import STEPS
from .execution_schema import ExecutionPlan, PlannedStep, StepState
from .logging_utils import step_log_path
from .workflow_registry import REGISTRY

# Default execution order — mirrors STEPS in constants.py
DEFAULT_WORKFLOW_ORDER: List[str] = list(STEPS)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _get_config_value(cfg: dict, key_path: str):
    """Retrieve a value from a nested config dict by dot-separated key path.

    Returns ``None`` if any intermediate key is missing.
    """
    parts = key_path.split(".")
    node = cfg
    for part in parts:
        if not isinstance(node, dict):
            return None
        node = node.get(part)
    return node


def _is_step_enabled(
    step_name: str,
    cfg: dict,
    start_at: Optional[str],
    stop_after: Optional[str],
) -> bool:
    """Return ``True`` if *step_name* should be included in the execution plan.

    A step is excluded when:

    - it falls outside the ``[start_at, stop_after]`` range, OR
    - it is optional and its ``enabled_config_key`` resolves to a falsy value.
    """
    # Range check — mirrors pipeline.step_in_range to avoid a circular import.
    all_steps = DEFAULT_WORKFLOW_ORDER
    try:
        i = all_steps.index(step_name)
        i0 = all_steps.index(start_at) if start_at else 0
        i1 = all_steps.index(stop_after) if stop_after else len(all_steps) - 1
    except ValueError:
        # Step not in the known order; include it for forward compatibility.
        return True
    if not (i0 <= i <= i1):
        return False

    step_def = REGISTRY.get(step_name)
    if step_def is None:
        # Unknown step — assume enabled for forward compatibility.
        return True

    if step_def.optional and step_def.enabled_config_key:
        enabled = _get_config_value(cfg, step_def.enabled_config_key)
        if not enabled:
            return False

    return True


def _build_downstream_map(order: List[str]) -> Dict[str, List[str]]:
    """Return a map from each step name to the list of its direct downstream steps.

    Only steps present in *order* are included as downstream targets.
    """
    downstream: Dict[str, List[str]] = {name: [] for name in order}
    for name in order:
        step_def = REGISTRY.get(name)
        if step_def is None:
            continue
        for dep in step_def.dependencies:
            if dep in downstream:
                downstream[dep].append(name)
    return downstream


def _collect_required_tools(enabled_steps: List[str]) -> List[str]:
    """Collect the deduplicated list of non-optional tools needed by *enabled_steps*."""
    tools: List[str] = []
    seen: set = set()
    for name in enabled_steps:
        step_def = REGISTRY.get(name)
        if step_def is None:
            continue
        for req in step_def.tool_requirements:
            if not req.optional and req.name not in seen:
                tools.append(req.name)
                seen.add(req.name)
    return tools


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def build_execution_plan(
    cfg: dict,
    workflow_order: Optional[List[str]] = None,
) -> ExecutionPlan:
    """Build an :class:`~phylofoundry.execution_schema.ExecutionPlan` for *cfg*.

    Parameters
    ----------
    cfg:
        Resolved pipeline configuration dict (output of
        :func:`~phylofoundry.config.resolve_config`).
    workflow_order:
        Ordered list of step names to consider.  Defaults to
        :data:`DEFAULT_WORKFLOW_ORDER`.

    Returns
    -------
    ExecutionPlan
        A fully populated plan with enabled/disabled steps, required tools,
        and per-step log paths.
    """
    if workflow_order is None:
        workflow_order = DEFAULT_WORKFLOW_ORDER

    outdir = cfg.get("output", {}).get("outdir") or ""
    logs_dir = os.path.join(outdir, "logs") if outdir else "logs"

    start_at = cfg.get("workflow", {}).get("start_at")
    stop_after = cfg.get("workflow", {}).get("stop_after")

    enabled: List[str] = []
    disabled: List[str] = []
    for name in workflow_order:
        if _is_step_enabled(name, cfg, start_at, stop_after):
            enabled.append(name)
        else:
            disabled.append(name)

    required_tools = _collect_required_tools(enabled)

    planned_steps: List[PlannedStep] = []
    for name in workflow_order:
        step_def = REGISTRY.get(name)
        deps = step_def.dependencies if step_def else []
        tools = (
            [r.name for r in step_def.tool_requirements if not r.optional]
            if step_def
            else []
        )
        state = StepState.PENDING if name in enabled else StepState.DISABLED
        planned_steps.append(
            PlannedStep(
                name=name,
                state=state,
                dependencies=deps,
                required_tools=tools,
                log_path=step_log_path(logs_dir, name),
            )
        )

    return ExecutionPlan(
        steps=planned_steps,
        logs_dir=logs_dir,
        enabled_step_names=enabled,
        disabled_step_names=disabled,
        required_tools=required_tools,
    )

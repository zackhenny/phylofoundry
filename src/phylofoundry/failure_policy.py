"""Failure policy scaffolding for PhyloFoundry.

Provides helpers for understanding the downstream impact of step failures
and for propagating blocked-step states through the execution plan.

This module is intentionally minimal in Phase 1 of the architecture plan.
It defines the data structures and helper functions needed to support
dependency-aware failure continuation in Phase 4.
"""

from __future__ import annotations

from typing import Dict, List, Set

from .workflow_registry import REGISTRY


def build_reverse_dependency_map(step_names: List[str]) -> Dict[str, List[str]]:
    """Return a map from each step to its direct downstream dependents.

    Parameters
    ----------
    step_names:
        Ordered list of step names to include.

    Returns
    -------
    dict
        ``{step_name: [directly_downstream_step, ...]}``
    """
    reverse: Dict[str, List[str]] = {name: [] for name in step_names}
    for name in step_names:
        step_def = REGISTRY.get(name)
        if step_def is None:
            continue
        for dep in step_def.dependencies:
            if dep in reverse:
                reverse[dep].append(name)
    return reverse


def get_impacted_steps(
    failed_step: str,
    reverse_map: Dict[str, List[str]],
) -> Set[str]:
    """Return the full set of steps transitively blocked by *failed_step*.

    Uses a breadth-first traversal of the reverse dependency graph.

    Parameters
    ----------
    failed_step:
        The name of the step that failed.
    reverse_map:
        Output of :func:`build_reverse_dependency_map`.

    Returns
    -------
    set
        Step names that should be marked ``BLOCKED`` (excludes *failed_step*
        itself).
    """
    blocked: Set[str] = set()
    queue = list(reverse_map.get(failed_step, []))
    while queue:
        step = queue.pop(0)
        if step in blocked:
            continue
        blocked.add(step)
        queue.extend(reverse_map.get(step, []))
    return blocked


def apply_failure_policy(
    plan_steps: list,
    failed_step: str,
    reverse_map: Dict[str, List[str]],
) -> None:
    """Mark steps in *plan_steps* as ``BLOCKED`` if they are downstream of *failed_step*.

    Only steps currently in ``PENDING`` state are transitioned to ``BLOCKED``;
    disabled, already-failed, or already-succeeded steps are left unchanged.

    This mutates the ``state`` attribute of each affected
    :class:`~phylofoundry.execution_schema.PlannedStep` in place.

    Parameters
    ----------
    plan_steps:
        List of :class:`~phylofoundry.execution_schema.PlannedStep` objects.
    failed_step:
        Name of the step that failed.
    reverse_map:
        Output of :func:`build_reverse_dependency_map`.
    """
    from .execution_schema import StepState

    impacted = get_impacted_steps(failed_step, reverse_map)
    for step in plan_steps:
        if step.name in impacted and step.state == StepState.PENDING:
            step.state = StepState.BLOCKED

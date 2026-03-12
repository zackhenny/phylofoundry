"""Preflight, planning, and doctor utilities for PhyloFoundry.

Provides human-readable CLI output for four inspection modes:

- :func:`print_list_steps`   -- list every registered workflow step
- :func:`print_plan`         -- show which steps are enabled for a given config
- :func:`validate_config_only` -- validate config without running the pipeline
- :func:`run_doctor`         -- diagnose tool / environment issues
"""

from __future__ import annotations

import shutil
import sys
from typing import Optional

from .constants import STEPS
from .execution_planner import build_execution_plan
from .workflow_registry import REGISTRY

# ── ANSI colour helpers (disabled on non-TTY) ─────────────────────────────────

_USE_COLOUR = sys.stdout.isatty()

_GREEN  = "\033[32m" if _USE_COLOUR else ""
_YELLOW = "\033[33m" if _USE_COLOUR else ""
_RED    = "\033[31m" if _USE_COLOUR else ""
_BOLD   = "\033[1m"  if _USE_COLOUR else ""
_RESET  = "\033[0m"  if _USE_COLOUR else ""


def _ok(msg: str) -> str:
    return f"{_GREEN}✔{_RESET}  {msg}"


def _warn(msg: str) -> str:
    return f"{_YELLOW}⚠{_RESET}  {msg}"


def _fail(msg: str) -> str:
    return f"{_RED}✘{_RESET}  {msg}"


def _section(title: str) -> None:
    print(f"\n{_BOLD}{title}{_RESET}")
    print("─" * (len(title) + 2))


# ---------------------------------------------------------------------------
# Public helpers
# ---------------------------------------------------------------------------


def print_list_steps() -> None:
    """Print all registered workflow steps with descriptions and metadata."""
    _section("Known workflow steps")
    for name in STEPS:
        step = REGISTRY.get(name)
        if step is None:
            print(f"  {_YELLOW}{name}{_RESET}  (not registered)")
            continue

        optional_tag = f"  {_YELLOW}[optional]{_RESET}" if step.optional else ""
        print(f"  {_BOLD}{name}{_RESET}{optional_tag}")
        print(f"      {step.description}")

        if step.dependencies:
            print(f"      depends-on : {', '.join(step.dependencies)}")

        required_tools = [t for t in step.tool_requirements if not t.optional]
        optional_tools = [t for t in step.tool_requirements if t.optional]
        if required_tools:
            print(f"      tools      : {', '.join(t.name for t in required_tools)}")
        if optional_tools:
            print(f"      opt-tools  : {', '.join(t.name for t in optional_tools)}")

        if step.enabled_config_key:
            print(f"      enabled-key: {step.enabled_config_key}")

        print()


def print_plan(cfg: dict) -> None:
    """Print the execution plan that would result from *cfg*.

    Builds the plan via :func:`~phylofoundry.execution_planner.build_execution_plan`
    and renders it as a human-readable table showing each step's enabled/disabled
    status, required tools, and dependencies.
    """
    plan = build_execution_plan(cfg)

    _section("Execution plan")

    print(f"  Enabled  : {len(plan.enabled_step_names)} step(s)")
    print(f"  Disabled : {len(plan.disabled_step_names)} step(s)")
    if plan.required_tools:
        print(f"  Required tools: {', '.join(plan.required_tools)}")
    print()

    for step in plan.steps:
        from .execution_schema import StepState
        if step.state == StepState.PENDING:
            marker = _ok("ENABLED ")
        else:
            marker = _warn("DISABLED")

        deps_str = f"  ← {', '.join(step.dependencies)}" if step.dependencies else ""
        tools_str = (
            f"  [tools: {', '.join(step.required_tools)}]"
            if step.required_tools
            else ""
        )
        print(f"  {marker}  {_BOLD}{step.name:<22}{_RESET}{deps_str}{tools_str}")

    print()
    if plan.disabled_step_names:
        print(
            f"  {_YELLOW}Tip:{_RESET} disabled optional steps can be enabled via their config flag "
            "(e.g. embeddings.enabled=true)."
        )


def validate_config_only(cfg: dict) -> bool:
    """Validate *cfg* and print the result.

    Returns
    -------
    bool
        ``True`` if the config is valid, ``False`` otherwise.
    """
    _section("Config validation")

    errors: list[str] = []

    faa_arg = cfg.get("inputs", {}).get("faa_dir")
    hmm_arg = cfg.get("inputs", {}).get("hmm_input")
    outdir  = cfg.get("output", {}).get("outdir")

    if not faa_arg:
        errors.append("inputs.faa_dir is not set")
    if not hmm_arg:
        errors.append("inputs.hmm_input is not set")
    if not outdir:
        errors.append("output.outdir is not set")

    start_at   = cfg.get("workflow", {}).get("start_at")
    stop_after = cfg.get("workflow", {}).get("stop_after")
    if start_at and start_at not in STEPS:
        errors.append(f"workflow.start_at '{start_at}' is not a known step")
    if stop_after and stop_after not in STEPS:
        errors.append(f"workflow.stop_after '{stop_after}' is not a known step")
    if start_at and stop_after:
        if STEPS.index(start_at) > STEPS.index(stop_after):
            errors.append(
                f"workflow.start_at '{start_at}' comes after workflow.stop_after '{stop_after}'"
            )

    if errors:
        for err in errors:
            print(f"  {_fail(err)}")
        print()
        return False

    print(f"  {_ok('Config is valid.')}")
    print()
    return True


def run_doctor(cfg: Optional[dict] = None) -> bool:
    """Diagnose tool and environment issues.

    Checks every tool required by any registered step (and every optional
    tool) and reports which are available / missing.  If *cfg* is provided,
    also runs a config validation.

    Returns
    -------
    bool
        ``True`` if no required tools are missing (optional warnings are
        non-fatal), ``False`` otherwise.
    """
    _section("Environment / tool doctor")

    all_ok = True

    # ── Collect all tools from the registry ──────────────────────────────
    required_tools: list[tuple[str, str]] = []   # (exe, step_name)
    optional_tools: list[tuple[str, str]] = []

    for name in STEPS:
        step = REGISTRY.get(name)
        if step is None:
            continue
        for req in step.tool_requirements:
            entry = (req.name, name)
            if req.optional:
                optional_tools.append(entry)
            else:
                required_tools.append(entry)

    # Deduplicate while preserving first-seen order
    seen: set[str] = set()
    unique_required: list[tuple[str, str]] = []
    for exe, sname in required_tools:
        if exe not in seen:
            seen.add(exe)
            unique_required.append((exe, sname))

    seen_opt: set[str] = set()
    unique_optional: list[tuple[str, str]] = []
    for exe, sname in optional_tools:
        if exe not in seen_opt and exe not in seen:
            seen_opt.add(exe)
            unique_optional.append((exe, sname))

    print("  Required tools")
    for exe, sname in unique_required:
        if shutil.which(exe):
            print(f"    {_ok(exe)}")
        else:
            print(f"    {_fail(exe)}  (first needed by step: {sname})")
            all_ok = False

    if unique_optional:
        print()
        print("  Optional tools")
        for exe, sname in unique_optional:
            if shutil.which(exe):
                print(f"    {_ok(exe)}")
            else:
                print(f"    {_warn(exe + ' not found')}  (optional; step: {sname})")

    # ── Python packages ───────────────────────────────────────────────────
    _python_package_checks()

    # ── Config validation (if cfg supplied) ───────────────────────────────
    if cfg is not None:
        print()
        config_ok = validate_config_only(cfg)
        if not config_ok:
            all_ok = False

    print()
    if all_ok:
        print(f"  {_ok('All required tools found. Environment looks good.')}")
    else:
        print(
            f"  {_fail('One or more required tools are missing or the config is invalid.')}"
            "\n  Install the missing tools and re-run."
        )

    return all_ok


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _python_package_checks() -> None:
    """Check optional Python package dependencies and print results."""
    packages = [
        ("torch",    "PyTorch (needed for embed / score_motifs / discover_motifs)"),
        ("skbio",    "scikit-bio (needed for conservation_metrics)"),
        ("pandas",   "pandas (core dependency)"),
        ("Bio",      "Biopython (core dependency)"),
    ]

    print()
    print("  Python packages")
    for mod, label in packages:
        try:
            __import__(mod)
            print(f"    {_ok(label)}")
        except ImportError:
            print(f"    {_warn(label + ' — not installed')}")

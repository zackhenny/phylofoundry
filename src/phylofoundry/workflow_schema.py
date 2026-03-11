"""Workflow metadata schema for PhyloFoundry.

Defines dataclasses for describing workflow steps, their artifacts,
tool requirements, and the overall registry.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional


class CheckpointStrategy(str, Enum):
    """How a step determines whether it can skip execution."""

    NONE = "none"
    OUTPUT_EXISTS = "output_exists"
    FORCE_FLAG = "force_flag"


class ArtifactKind(str, Enum):
    """Classification of artifacts produced or consumed by a step."""

    TABLE = "table"
    FASTA = "fasta"
    ALIGNMENT = "alignment"
    TREE = "tree"
    EMBEDDING = "embedding"
    PLOT = "plot"
    JSON = "json"
    DIRECTORY = "directory"
    LOG = "log"


@dataclass
class ArtifactSpec:
    """Describes an artifact produced or consumed by a step."""

    name: str
    kind: ArtifactKind
    path_template: str
    """Path relative to outdir; may contain ``{step}`` or similar placeholders."""
    description: str = ""
    optional: bool = False


@dataclass
class ToolRequirement:
    """An external tool required by a step."""

    name: str
    """Executable name, e.g. ``"hmmscan"``."""
    config_key: Optional[str] = None
    """Dot-separated config path that overrides the default binary name."""
    optional: bool = False
    """If True, the step can degrade gracefully without this tool."""
    description: str = ""


@dataclass
class StepDefinition:
    """Full contract for a single workflow step."""

    name: str
    description: str
    inputs: List[str] = field(default_factory=list)
    """Names of upstream step outputs consumed by this step."""
    outputs: List[ArtifactSpec] = field(default_factory=list)
    tool_requirements: List[ToolRequirement] = field(default_factory=list)
    dependencies: List[str] = field(default_factory=list)
    """Names of upstream steps that must complete before this step runs."""
    optional: bool = False
    """If True, this step has an ``enabled`` config flag and may be skipped."""
    enabled_config_key: Optional[str] = None
    """Dot-separated config path for the boolean enabled flag, e.g. ``"embeddings.enabled"``."""
    checkpoint_strategy: CheckpointStrategy = CheckpointStrategy.OUTPUT_EXISTS
    notes: str = ""
    """Freeform developer notes (architecture intentions, future changes, etc.)."""


@dataclass
class WorkflowRegistry:
    """Registry of all known workflow steps.

    Steps are stored in a dict keyed by name, and the dependency graph is
    validated at construction time via :meth:`validate`.
    """

    steps: Dict[str, StepDefinition] = field(default_factory=dict)
    """Mapping of step name to :class:`StepDefinition`."""

    def register(self, step: StepDefinition) -> None:
        """Register a step definition."""
        self.steps[step.name] = step

    def validate(self) -> None:
        """Validate that all declared dependencies refer to registered steps.

        Raises
        ------
        ValueError
            If any dependency references an unknown step name.
        """
        errors = []
        for name, step in self.steps.items():
            for dep in step.dependencies:
                if dep not in self.steps:
                    errors.append(
                        f"Step '{name}' declares dependency on unknown step '{dep}'"
                    )
        if errors:
            raise ValueError(
                "Workflow registry validation failed:\n" + "\n".join(errors)
            )

    def get(self, name: str) -> Optional[StepDefinition]:
        """Return a step definition by name, or ``None`` if not found."""
        return self.steps.get(name)

    def ordered_names(self, order: List[str]) -> List[str]:
        """Filter an ordered list to only registered step names."""
        return [n for n in order if n in self.steps]

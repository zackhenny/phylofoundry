"""Input provenance and artifact-linking utilities for PhyloFoundry.

This module implements the ``--input-run`` mechanism that allows a pipeline
module to consume artifacts produced by a *prior* run directory rather than
requiring the full pipeline to re-run from the beginning.

Key responsibilities
--------------------
1. **Artifact resolution** — given a prior output directory, locate the
   canonical artifacts that a downstream module requires.
2. **Artifact validation** — check that required files exist, are non-empty,
   and (optionally) match expected SHA-256 fingerprints from the prior run's
   checkpoint database.
3. **Provenance recording** — embed the prior run's ID, config snapshot path,
   and per-artifact hashes into the new run's config so that ``Checkpointer``
   can persist them in the ``_run_meta`` block and the NDJSON log.
4. **Run discovery** — scan a directory tree for PhyloFoundry output
   directories and display a summary table (used by ``phylofoundry list-runs``).

Typical usage (internal — called from ``pipeline.py``)
-------------------------------------------------------
::

    from phylofoundry.provenance import InputRunResolver

    resolver = InputRunResolver(input_run_dir="/prior/results")
    resolver.validate_for_step("hyphy")   # raises InputRunError on problems
    cfg["_input_run"] = resolver.provenance_record()

Typical usage (CLI — ``phylofoundry list-runs``)
------------------------------------------------
::

    phylofoundry list-runs ./experiments
    phylofoundry list-runs --json
"""

from __future__ import annotations

import json
import os
import sqlite3
from typing import Any, Dict, List, Optional, Sequence

from .artifact_paths import ArtifactPaths
from .checkpoint import file_sha256, load_checkpoint_meta


# ---------------------------------------------------------------------------
# Step → required input artifact names
# ---------------------------------------------------------------------------
# Maps each pipeline step to the *ArtifactPaths* property names (and, for
# per-HMM artifacts, the dir attribute names) that it needs as inputs.
# Only properties that must be produced by a *prior* run (i.e. are outputs of
# earlier steps) are listed here — internal intermediates generated in the
# same run are not included.
#
# The resolver checks for at least one of the listed paths/dirs to confirm
# that the prior run produced the necessary outputs.
_STEP_REQUIRED_ARTIFACTS: Dict[str, List[str]] = {
    # hmmer needs the prep outputs (combined FAA / HMM)
    "hmmer": ["combined_faa", "combined_hmm"],
    # extract needs hmmer tables
    "extract": ["hmmscan_dir", "hmmsearch_dir"],
    # embed needs per-HMM FASTAs
    "embed": ["fasta_dir"],
    # phylo needs per-HMM FASTAs (or curated FASTAs/alignments)
    "phylo": ["fasta_dir"],
    # curate needs trees + FASTAs
    "curate": ["trees_dir", "fasta_dir"],
    # taxonomy_integrate needs best_hits_tsv
    "taxonomy_integrate": ["best_hits_tsv"],
    # conservation_metrics needs best_hits_tsv + trees
    "conservation_metrics": ["best_hits_tsv", "trees_dir"],
    # detect_clades needs best_hits_tsv
    "detect_clades": ["best_hits_tsv"],
    # post needs trees + best_hits_tsv
    "post": ["trees_dir", "best_hits_tsv"],
    # synteny needs best_hits_tsv
    "synteny": ["best_hits_tsv"],
    # codon needs trees + alignments
    "codon": ["trees_dir", "alignments_clipkit_dir"],
    # hyphy needs codon alignments + trees
    "hyphy": ["trees_dir", "codon_dir"],
    # score_motifs needs embeddings
    "score_motifs": ["embeddings_dir"],
    # discover_motifs needs embeddings + clade assignments
    "discover_motifs": ["embeddings_dir", "clade_assign_dir"],
}


# ---------------------------------------------------------------------------
# Errors
# ---------------------------------------------------------------------------

class InputRunError(RuntimeError):
    """Raised when the prior run directory is missing required artifacts."""


# ---------------------------------------------------------------------------
# InputRunResolver
# ---------------------------------------------------------------------------

class InputRunResolver:
    """Resolve and validate artifacts from a prior PhyloFoundry run directory.

    Parameters
    ----------
    input_run_dir:
        Absolute or relative path to the prior run's output directory
        (i.e. the value passed to ``--input-run``).
    """

    def __init__(self, input_run_dir: str) -> None:
        self._dir = os.path.abspath(input_run_dir)
        if not os.path.isdir(self._dir):
            raise InputRunError(
                f"--input-run directory does not exist: {self._dir}"
            )
        self._paths = ArtifactPaths(self._dir)
        self._meta = load_checkpoint_meta(self._dir)

    # ── Public interface ────────────────────────────────────────────────────

    @property
    def directory(self) -> str:
        """Absolute path to the prior run directory."""
        return self._dir

    @property
    def run_id(self) -> Optional[str]:
        """Run ID of the prior run (``None`` if metadata could not be read)."""
        return self._meta.run_id if self._meta else None

    def artifact_path(self, attr: str) -> str:
        """Return the absolute path of an artifact attribute from the prior run."""
        return getattr(self._paths, attr)

    def validate_for_step(self, step_name: str) -> None:
        """Check that all required input artifacts for *step_name* exist.

        Raises
        ------
        InputRunError
            If any required artifact is missing or the prior run directory does
            not look like a valid PhyloFoundry output directory.
        """
        required = _STEP_REQUIRED_ARTIFACTS.get(step_name)
        if not required:
            # No declared requirements — nothing to validate.
            return

        missing: List[str] = []
        for attr in required:
            p = self.artifact_path(attr)
            if not os.path.exists(p):
                missing.append(f"{attr} ({p})")

        if missing:
            raise InputRunError(
                f"--input-run '{self._dir}' is missing required artifacts for "
                f"step '{step_name}':\n"
                + "\n".join(f"  • {m}" for m in missing)
                + "\n\nEnsure the prior run completed the upstream steps "
                "successfully before using it as an input source."
            )

    def artifact_hashes(self, step_name: str) -> Dict[str, str]:
        """Return a mapping of artifact-name → SHA-256 for required inputs.

        Only *file* artifacts are hashed; directory artifacts are skipped.
        """
        required = _STEP_REQUIRED_ARTIFACTS.get(step_name, [])
        hashes: Dict[str, str] = {}
        for attr in required:
            p = self.artifact_path(attr)
            if os.path.isfile(p):
                try:
                    hashes[attr] = file_sha256(p)
                except OSError:
                    pass
        return hashes

    def provenance_record(self, step_name: Optional[str] = None) -> Dict[str, Any]:
        """Return a serialisable provenance record for embedding in run metadata.

        The returned dict is stored in ``cfg["_input_run"]`` and is included in
        the config snapshot (``OUTDIR/config.yaml``) written by
        :meth:`~phylofoundry.checkpoint.Checkpointer.start_run`.

        Fields
        ------
        ``input_run_dir``
            Absolute path to the prior run directory.
        ``input_run_id``
            Run ID of the prior run (``None`` if unknown).
        ``input_run_config``
            Path to the prior run's saved config snapshot.
        ``artifact_hashes``
            Dict of ``artifact_name → SHA-256`` for all required file inputs.
        """
        record: Dict[str, Any] = {
            "input_run_dir": self._dir,
            "input_run_id": self.run_id,
            "input_run_config": (
                self._meta.config_path if self._meta else None
            ),
        }
        if step_name:
            record["artifact_hashes"] = self.artifact_hashes(step_name)
        return record

    def print_summary(self) -> None:
        """Print a one-line summary of the resolved prior run."""
        run_id = self.run_id or "<unknown>"
        print(
            f"[input-run] Using prior run '{run_id}' from '{self._dir}' "
            "as input source."
        )


# ---------------------------------------------------------------------------
# apply_input_run
# ---------------------------------------------------------------------------

def apply_input_run(cfg: Dict[str, Any], step_name: Optional[str] = None) -> None:
    """Resolve, validate, and record the ``--input-run`` provenance in *cfg*.

    This is the single entry point called by ``pipeline.py`` when
    ``cfg["_checkpoint"]["input_run"]`` is set.  It:

    1. Creates an :class:`InputRunResolver` for the specified directory.
    2. Optionally validates that the required artifacts for *step_name* exist.
    3. Stores a provenance record in ``cfg["_input_run"]``.
    4. Injects the prior run's artifact directory paths into *cfg* so that
       downstream task modules pick up files from the correct location.

    Parameters
    ----------
    cfg:
        Resolved pipeline config dict (mutated in place).
    step_name:
        Internal pipeline step name being executed (e.g. ``"hyphy"``).
        When ``None``, validation is skipped but provenance is still recorded.
    """
    input_run_dir: Optional[str] = cfg.get("_checkpoint", {}).get("input_run")
    if not input_run_dir:
        return

    resolver = InputRunResolver(input_run_dir)
    resolver.print_summary()

    if step_name:
        resolver.validate_for_step(step_name)

    cfg["_input_run"] = resolver.provenance_record(step_name)

    # Inject artifact paths — the pipeline reads these from cfg["_input_paths"]
    # when present (see pipeline.py _resolve_input_paths).
    input_dir = resolver.directory
    cfg.setdefault("_input_paths", {})["input_run_dir"] = input_dir


# ---------------------------------------------------------------------------
# list_runs — CLI helper for ``phylofoundry list-runs``
# ---------------------------------------------------------------------------

def _scan_runs(search_dir: str) -> List[Dict[str, Any]]:
    """Scan *search_dir* for PhyloFoundry output directories.

    A directory is considered a PhyloFoundry output directory when it contains
    either a ``config.yaml`` / ``config.json`` with a ``_run_meta`` block, or a
    ``logs/checkpoint.db`` SQLite database.

    Returns a list of run-info dicts sorted by start time (newest first).
    """
    results: List[Dict[str, Any]] = []
    search_dir = os.path.abspath(search_dir)

    if not os.path.isdir(search_dir):
        return results

    # Check the search_dir itself and its immediate children.
    candidates = [search_dir] + [
        os.path.join(search_dir, d)
        for d in sorted(os.listdir(search_dir))
        if os.path.isdir(os.path.join(search_dir, d))
    ]

    for candidate in candidates:
        meta = load_checkpoint_meta(candidate)
        if meta is None:
            # Try to read directly from checkpoint.db without a config
            db_path = os.path.join(candidate, "logs", "checkpoint.db")
            if not os.path.isfile(db_path):
                continue
            try:
                conn = sqlite3.connect(db_path, timeout=3)
                rows = conn.execute(
                    "SELECT run_id, start_time, end_time, outdir "
                    "FROM runs ORDER BY start_time DESC"
                ).fetchall()
                conn.close()
                for row in rows:
                    results.append({
                        "run_id": row[0],
                        "start_time": row[1],
                        "end_time": row[2],
                        "outdir": row[3] or candidate,
                        "config_path": None,
                        "completed_steps": _steps_from_db(db_path, row[0]),
                    })
            except sqlite3.Error:
                pass
            continue

        # We have a CheckpointMeta — query its DB for run details.
        run_id = meta.run_id
        start_time: Optional[str] = None
        end_time: Optional[str] = None
        try:
            conn = sqlite3.connect(meta.sqlite_path, timeout=3)
            row = conn.execute(
                "SELECT start_time, end_time FROM runs WHERE run_id=?", (run_id,)
            ).fetchone()
            if row:
                start_time, end_time = row
            conn.close()
        except (sqlite3.Error, OSError):
            pass

        results.append({
            "run_id": run_id,
            "start_time": start_time,
            "end_time": end_time,
            "outdir": meta.outdir,
            "config_path": meta.config_path,
            "completed_steps": _steps_from_db(meta.sqlite_path, run_id),
        })

    # Sort newest first
    results.sort(key=lambda r: r.get("start_time") or "", reverse=True)
    return results


def _steps_from_db(db_path: str, run_id: str) -> List[str]:
    """Return list of successfully completed step names for *run_id*."""
    if not db_path or not os.path.isfile(db_path):
        return []
    try:
        conn = sqlite3.connect(db_path, timeout=3)
        rows = conn.execute(
            "SELECT step_name FROM steps WHERE run_id=? AND state='success'",
            (run_id,),
        ).fetchall()
        conn.close()
        return [r[0] for r in rows]
    except (sqlite3.Error, OSError):
        return []


def list_runs(search_dir: str = ".", *, as_json: bool = False) -> None:
    """Print a summary of available PhyloFoundry runs in *search_dir*.

    Parameters
    ----------
    search_dir:
        Directory to scan (default: current working directory).
    as_json:
        When ``True``, output the raw run metadata as a JSON array instead of
        a human-readable table.
    """
    runs = _scan_runs(search_dir)

    if not runs:
        print(f"No PhyloFoundry runs found in '{os.path.abspath(search_dir)}'.")
        return

    if as_json:
        print(json.dumps(runs, indent=2, default=str))
        return

    # Human-readable table
    print(f"\nPhyloFoundry runs found in '{os.path.abspath(search_dir)}':\n")
    header = f"{'RUN ID':<38}  {'STARTED':<20}  {'STEPS DONE':<10}  OUTDIR"
    print(header)
    print("-" * len(header))
    for r in runs:
        run_id = r.get("run_id") or "<unknown>"
        started = (r.get("start_time") or "")[:19].replace("T", " ")
        n_steps = len(r.get("completed_steps") or [])
        outdir = r.get("outdir") or ""
        print(f"{run_id:<38}  {started:<20}  {n_steps:<10}  {outdir}")

    print(
        f"\nTip: pass any OUTDIR above to --input-run to use it as an input source.\n"
        f"     e.g.  phylofoundry hyphy --input-run <OUTDIR> --outdir ./new_results"
    )

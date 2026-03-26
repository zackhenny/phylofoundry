"""Checkpoint and logging utilities for PhyloFoundry.

Implements:
- Append-only NDJSON log per pipeline run (atomic writes, fsync per entry).
- SQLite database for fast lookup of runs, steps, and artifacts.
- Content fingerprinting: each step hashes its parameters + input artifacts +
  pipeline version so that re-runs can detect unchanged steps and skip them.
- Config snapshot writing: saves a copy of the resolved config to OUTDIR at
  run start so ``--resume`` can reload it automatically.
- Atomic writes using a tmp file + ``os.replace``, WAL mode for SQLite.

Typical usage
-------------
Run start::

    checkpointer = Checkpointer(outdir)
    checkpointer.start_run(run_id, cfg)

Per-step::

    checkpointer.record_step_pending(run_id, "hmmer", params={...})
    checkpointer.record_step_running(run_id, "hmmer")
    checkpointer.record_step_success(run_id, "hmmer", outputs=[...])
    # or on failure:
    checkpointer.record_step_failure(run_id, "hmmer", error="...", trace="...")

Resume::

    meta = load_checkpoint_meta(outdir)
    checkpointer = Checkpointer(outdir)
    step_info = checkpointer.latest_step(run_id, "hmmer")
    fp = compute_fingerprint(params, input_files)
    if step_info and step_info["state"] == "success" and step_info["input_fingerprint"] == fp:
        # skip step
        ...
"""

from __future__ import annotations

import hashlib
import json
import os
import sqlite3
import tempfile
import traceback
import uuid
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional

# ---------------------------------------------------------------------------
# Version sentinel — bump when step logic changes and old checkpoints should
# be treated as stale.
# ---------------------------------------------------------------------------
CHECKPOINT_VERSION = "1"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _now_iso() -> str:
    """Return the current UTC time as an ISO-8601 string."""
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _atomic_append_ndjson(path: str, record: dict) -> None:
    """Append *record* as a JSON line to *path* atomically.

    Uses a temporary file in the same directory + ``os.replace`` to avoid
    partial writes on crash, then fsyncs the containing directory to ensure
    durability.
    """
    line = json.dumps(record, separators=(",", ":"), default=str) + "\n"
    dir_ = os.path.dirname(os.path.abspath(path))
    os.makedirs(dir_, exist_ok=True)
    # Append-only: read existing content, write new file atomically.
    existing = b""
    if os.path.exists(path):
        with open(path, "rb") as fh:
            existing = fh.read()
    fd, tmp = tempfile.mkstemp(dir=dir_)
    try:
        with os.fdopen(fd, "wb") as fh:
            fh.write(existing)
            fh.write(line.encode())
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp, path)
    except Exception:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise
    # fsync the directory so the rename is visible
    try:
        dfd = os.open(dir_, os.O_RDONLY)
        os.fsync(dfd)
        os.close(dfd)
    except OSError:
        pass


# ---------------------------------------------------------------------------
# Fingerprinting
# ---------------------------------------------------------------------------


def compute_fingerprint(params: Dict[str, Any], input_files: List[str]) -> str:
    """Compute a SHA-256 fingerprint for a step.

    The fingerprint is computed over:
    - The serialised *params* dict (sorted keys for stability).
    - The SHA-256 hash of each *input_files* path (if the file exists).
    - The :data:`CHECKPOINT_VERSION` sentinel.

    Parameters
    ----------
    params:
        A dict of parameters that uniquely describe this step's configuration.
    input_files:
        A list of file paths whose content should be included in the hash.
        Missing files are skipped (they will cause a mismatch on the next run).

    Returns
    -------
    str
        40-character hex prefix of the full SHA-256 digest.  This provides
        ~160 bits of collision resistance, which is sufficient for step
        deduplication.  Use :func:`file_sha256` for full-length artifact hashes.
    """
    h = hashlib.sha256()
    h.update(CHECKPOINT_VERSION.encode())
    h.update(json.dumps(params, sort_keys=True, default=str).encode())
    for fp in sorted(input_files):
        if os.path.isfile(fp):
            fh = hashlib.sha256()
            with open(fp, "rb") as f:
                for chunk in iter(lambda: f.read(65536), b""):
                    fh.update(chunk)
            h.update(fp.encode())
            h.update(fh.digest())
    return h.hexdigest()[:40]


def file_sha256(path: str) -> str:
    """Return the full SHA-256 hex digest of *path*."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


# ---------------------------------------------------------------------------
# Per-subtask progress tracker
# ---------------------------------------------------------------------------


class SubtaskProgressTracker:
    """Track fine-grained per-subtask progress within a pipeline step.

    Long-running steps (e.g. ``embed``, ``phylo``) iterate over many
    independent sub-tasks (one per HMM profile, genome, etc.).  When a step
    is interrupted mid-way—due to OOM, SIGKILL, scheduler wall-time, or any
    other reason—this tracker lets the step resume from exactly where it
    stopped on the next ``--resume`` run.

    Progress is persisted in a hidden NDJSON file (e.g.
    ``.embed_progress.ndjson``) inside the step's output directory.  Each
    line is a JSON record with at minimum ``{"subtask": "<id>",
    "state": "running"|"success"|"failed", "timestamp": "..."}``.

    **Typical usage**::

        tracker = SubtaskProgressTracker(
            os.path.join(emb_dir, ".embed_progress.ndjson"),
            resume=resume,
        )
        for hmm in hmms:
            if tracker.is_done(hmm):
                print(f"[embed] RESUME: skipping {hmm}")
                continue
            tracker.record_running(hmm)
            try:
                process(hmm)
                tracker.record_success(hmm)
            except Exception as exc:
                tracker.record_failure(hmm, error=str(exc))
                raise
        tracker.cleanup()  # remove progress file on full success

    Parameters
    ----------
    log_path:
        Path to the NDJSON progress file.  The file is created automatically
        on first write; its containing directory must already exist.
    resume:
        When ``True`` (default), load any existing progress from *log_path*
        so that already-completed sub-tasks can be skipped.  When ``False``
        (e.g. ``--force`` mode) the tracker starts fresh and ignores any
        existing file.
    """

    def __init__(self, log_path: str, resume: bool = True) -> None:
        self.log_path = log_path
        self._done: set = set()
        if resume:
            self._load()

    # ── Loading ──────────────────────────────────────────────────────────────

    def _load(self) -> None:
        """Populate *_done* from any existing progress log."""
        if not os.path.exists(self.log_path):
            return
        try:
            with open(self.log_path) as fh:
                for line in fh:
                    line = line.strip()
                    if not line:
                        continue
                    try:
                        rec = json.loads(line)
                        if rec.get("state") == "success":
                            self._done.add(rec["subtask"])
                    except (json.JSONDecodeError, KeyError):
                        pass
        except OSError:
            pass

    # ── Queries ──────────────────────────────────────────────────────────────

    def is_done(self, subtask_id: str) -> bool:
        """Return ``True`` if *subtask_id* was previously recorded as successful."""
        return subtask_id in self._done

    def completed_count(self) -> int:
        """Return the number of successfully completed sub-tasks."""
        return len(self._done)

    # ── Recording ────────────────────────────────────────────────────────────

    def record_running(self, subtask_id: str, **extra: Any) -> None:
        """Append a ``running`` record for *subtask_id*."""
        self._append({"subtask": subtask_id, "state": "running", **extra})

    def record_success(self, subtask_id: str, **extra: Any) -> None:
        """Append a ``success`` record and mark *subtask_id* as done."""
        self._done.add(subtask_id)
        self._append({"subtask": subtask_id, "state": "success", **extra})

    def record_failure(self, subtask_id: str, error: str = "", **extra: Any) -> None:
        """Append a ``failed`` record for *subtask_id*."""
        self._append({"subtask": subtask_id, "state": "failed", "error": error, **extra})

    def record_skipped(self, subtask_id: str, reason: str = "", **extra: Any) -> None:
        """Append a ``skipped`` record (e.g. output already exists) and mark done."""
        self._done.add(subtask_id)
        self._append({"subtask": subtask_id, "state": "skipped", "reason": reason, **extra})

    # ── Cleanup ──────────────────────────────────────────────────────────────

    def cleanup(self) -> None:
        """Remove the progress log file after the full step completes successfully.

        Idempotent — does nothing when the file does not exist.
        """
        try:
            if os.path.exists(self.log_path):
                os.remove(self.log_path)
        except OSError:
            pass

    # ── Internal ─────────────────────────────────────────────────────────────

    def _append(self, record: dict) -> None:
        record["timestamp"] = _now_iso()
        _atomic_append_ndjson(self.log_path, record)


# ---------------------------------------------------------------------------
# SQLite schema
# ---------------------------------------------------------------------------

_SCHEMA = """
PRAGMA journal_mode=WAL;

CREATE TABLE IF NOT EXISTS runs (
    run_id      TEXT PRIMARY KEY,
    start_time  TEXT NOT NULL,
    end_time    TEXT,
    outdir      TEXT,
    config_path TEXT,
    ndjson_path TEXT,
    sqlite_path TEXT
);

CREATE TABLE IF NOT EXISTS steps (
    id                INTEGER PRIMARY KEY AUTOINCREMENT,
    run_id            TEXT NOT NULL,
    step_name         TEXT NOT NULL,
    sequence          INTEGER DEFAULT 0,
    state             TEXT NOT NULL,
    input_fingerprint TEXT,
    params            TEXT,
    outputs           TEXT,
    error             TEXT,
    trace             TEXT,
    started_at        TEXT,
    completed_at      TEXT,
    FOREIGN KEY (run_id) REFERENCES runs(run_id)
);

CREATE INDEX IF NOT EXISTS idx_steps_run_step ON steps (run_id, step_name);

CREATE TABLE IF NOT EXISTS artifacts (
    id          INTEGER PRIMARY KEY AUTOINCREMENT,
    run_id      TEXT NOT NULL,
    step_name   TEXT NOT NULL,
    path        TEXT NOT NULL,
    sha256      TEXT,
    size_bytes  INTEGER,
    created_at  TEXT,
    FOREIGN KEY (run_id) REFERENCES runs(run_id)
);
"""


def _open_db(db_path: str) -> sqlite3.Connection:
    """Open (or create) the SQLite checkpoint database at *db_path*.

    ``check_same_thread=False`` is set so the connection can be used from
    the main thread after being created in a constructor.  PhyloFoundry
    pipeline steps run sequentially; concurrent multi-threaded access to
    the same connection is not expected.  For parallel execution, callers
    should open separate connections.
    """
    os.makedirs(os.path.dirname(os.path.abspath(db_path)), exist_ok=True)
    conn = sqlite3.connect(db_path, timeout=30, check_same_thread=False)
    conn.row_factory = sqlite3.Row
    conn.executescript(_SCHEMA)
    conn.commit()
    return conn


# ---------------------------------------------------------------------------
# Checkpoint metadata
# ---------------------------------------------------------------------------


@dataclass
class CheckpointMeta:
    """Paths and identifiers for an existing checkpoint."""

    run_id: str
    outdir: str
    ndjson_path: str
    sqlite_path: str
    config_path: str
    extra: Dict[str, Any] = field(default_factory=dict)


def locate_saved_config(outdir: str) -> Optional[str]:
    """Return the path to a saved config inside *outdir*, or ``None``.

    Prefers ``config.yaml`` (canonical); falls back to ``config.json`` for
    legacy runs.
    """
    for name in ("config.yaml", "config.json"):
        p = os.path.join(outdir, name)
        if os.path.isfile(p):
            return p
    return None


def load_checkpoint_meta(outdir_or_path: str) -> Optional[CheckpointMeta]:
    """Load checkpoint metadata from *outdir_or_path*.

    Accepts either:
    - A path to a saved ``config.yaml`` / ``config.json``, or
    - A path to an output directory that contains one.

    Returns ``None`` when no checkpoint can be found.
    """
    if os.path.isdir(outdir_or_path):
        cfg_path = locate_saved_config(outdir_or_path)
        outdir = outdir_or_path
    elif os.path.isfile(outdir_or_path):
        cfg_path = outdir_or_path
        outdir = os.path.dirname(outdir_or_path)
    else:
        return None

    if cfg_path is None:
        return None

    meta: Dict[str, Any] = {}
    try:
        if cfg_path.endswith(".yaml") or cfg_path.endswith(".yml"):
            from ruamel.yaml import YAML
            yml = YAML()
            with open(cfg_path) as fh:
                loaded = yml.load(fh) or {}
            meta = loaded.get("_run_meta", {})
        else:
            with open(cfg_path) as fh:
                loaded = json.load(fh)
            meta = loaded.get("_run_meta", {})
    except Exception:
        meta = {}

    run_id = meta.get("run_id")
    if not run_id:
        # No embedded metadata — look for the latest checkpoint DB
        db_path = os.path.join(outdir, "logs", "checkpoint.db")
        ndjson_path = os.path.join(outdir, "logs", "pipeline.ndjson")
        if not os.path.isfile(db_path):
            return None
        # Try to retrieve latest run_id from DB
        try:
            conn = sqlite3.connect(db_path, timeout=5)
            row = conn.execute(
                "SELECT run_id FROM runs ORDER BY start_time DESC LIMIT 1"
            ).fetchone()
            conn.close()
            if row:
                run_id = row[0]
            else:
                return None
        except Exception:
            return None

    ndjson_path = meta.get(
        "ndjson_path", os.path.join(outdir, "logs", "pipeline.ndjson")
    )
    db_path = meta.get(
        "sqlite_path", os.path.join(outdir, "logs", "checkpoint.db")
    )
    return CheckpointMeta(
        run_id=run_id,
        outdir=outdir,
        ndjson_path=ndjson_path,
        sqlite_path=db_path,
        config_path=cfg_path,
        extra=meta,
    )


# ---------------------------------------------------------------------------
# Checkpointer
# ---------------------------------------------------------------------------


class Checkpointer:
    """High-level checkpoint manager for a pipeline run.

    Parameters
    ----------
    outdir:
        Pipeline output directory.  NDJSON and SQLite files are stored in
        ``outdir/logs/``.
    """

    def __init__(self, outdir: str) -> None:
        self._outdir = outdir
        self._logs_dir = os.path.join(outdir, "logs")
        os.makedirs(self._logs_dir, exist_ok=True)
        self._ndjson_path = os.path.join(self._logs_dir, "pipeline.ndjson")
        self._db_path = os.path.join(self._logs_dir, "checkpoint.db")
        self._conn = _open_db(self._db_path)

    # ── Properties ──────────────────────────────────────────────────────────

    @property
    def ndjson_path(self) -> str:
        return self._ndjson_path

    @property
    def db_path(self) -> str:
        return self._db_path

    # ── Run lifecycle ────────────────────────────────────────────────────────

    def start_run(self, run_id: str, cfg: Dict[str, Any]) -> None:
        """Record the start of a new pipeline run and write a config snapshot.

        The config is saved to ``OUTDIR/config.yaml`` with an embedded
        ``_run_meta`` block containing the run ID and checkpoint paths.
        """
        now = _now_iso()
        # Embed run metadata into the config snapshot
        meta: Dict[str, Any] = {
            "run_id": run_id,
            "start_time": now,
            "ndjson_path": self._ndjson_path,
            "sqlite_path": self._db_path,
            "outdir": self._outdir,
        }
        # Include input provenance when --input-run was supplied.
        input_run_record = cfg.get("_input_run")
        if input_run_record:
            meta["input_run"] = input_run_record
        self._write_config_snapshot(cfg, meta)

        # Register run in SQLite
        self._conn.execute(
            "INSERT OR REPLACE INTO runs "
            "(run_id, start_time, outdir, config_path, ndjson_path, sqlite_path) "
            "VALUES (?, ?, ?, ?, ?, ?)",
            (
                run_id,
                now,
                self._outdir,
                os.path.join(self._outdir, "config.yaml"),
                self._ndjson_path,
                self._db_path,
            ),
        )
        self._conn.commit()

        # Write NDJSON entry — include input_run provenance for auditability.
        ndjson_event: Dict[str, Any] = {
            "event": "run_start",
            "run_id": run_id,
            "timestamp": now,
            "outdir": self._outdir,
        }
        if input_run_record:
            ndjson_event["input_run"] = input_run_record
        self._append(ndjson_event)

    def end_run(self, run_id: str, success: bool) -> None:
        """Record the end of a pipeline run."""
        now = _now_iso()
        self._conn.execute(
            "UPDATE runs SET end_time=? WHERE run_id=?", (now, run_id)
        )
        self._conn.commit()
        self._append({
            "event": "run_end",
            "run_id": run_id,
            "timestamp": now,
            "success": success,
        })

    # ── Step lifecycle ───────────────────────────────────────────────────────

    def record_step_pending(
        self,
        run_id: str,
        step_name: str,
        params: Optional[Dict[str, Any]] = None,
        input_fingerprint: Optional[str] = None,
    ) -> None:
        """Record that *step_name* is about to be scheduled."""
        now = _now_iso()
        seq = self._next_sequence(run_id, step_name)
        self._conn.execute(
            "INSERT INTO steps "
            "(run_id, step_name, sequence, state, input_fingerprint, params, started_at) "
            "VALUES (?, ?, ?, 'pending', ?, ?, ?)",
            (
                run_id,
                step_name,
                seq,
                input_fingerprint,
                json.dumps(params or {}, default=str),
                now,
            ),
        )
        self._conn.commit()
        self._append({
            "event": "step_pending",
            "run_id": run_id,
            "step": step_name,
            "sequence": seq,
            "input_fingerprint": input_fingerprint,
            "params": params or {},
            "timestamp": now,
        })

    def record_step_running(self, run_id: str, step_name: str) -> None:
        """Record that *step_name* has started executing."""
        now = _now_iso()
        self._conn.execute(
            "UPDATE steps SET state='running', started_at=? "
            "WHERE run_id=? AND step_name=? AND id=("
            "  SELECT id FROM steps WHERE run_id=? AND step_name=? "
            "  ORDER BY id DESC LIMIT 1"
            ")",
            (now, run_id, step_name, run_id, step_name),
        )
        self._conn.commit()
        self._append({
            "event": "step_running",
            "run_id": run_id,
            "step": step_name,
            "timestamp": now,
        })

    def record_step_success(
        self,
        run_id: str,
        step_name: str,
        outputs: Optional[List[str]] = None,
        input_fingerprint: Optional[str] = None,
    ) -> None:
        """Record that *step_name* completed successfully."""
        now = _now_iso()
        outputs_json = json.dumps(outputs or [], default=str)
        self._conn.execute(
            "UPDATE steps SET state='success', completed_at=?, outputs=?, "
            "input_fingerprint=COALESCE(?, input_fingerprint) "
            "WHERE run_id=? AND step_name=? AND id=("
            "  SELECT id FROM steps WHERE run_id=? AND step_name=? "
            "  ORDER BY id DESC LIMIT 1"
            ")",
            (
                now,
                outputs_json,
                input_fingerprint,
                run_id,
                step_name,
                run_id,
                step_name,
            ),
        )
        self._conn.commit()
        self._append({
            "event": "step_success",
            "run_id": run_id,
            "step": step_name,
            "outputs": outputs or [],
            "input_fingerprint": input_fingerprint,
            "timestamp": now,
        })
        if outputs:
            self._record_artifacts(run_id, step_name, outputs, now)

    def record_step_failure(
        self,
        run_id: str,
        step_name: str,
        error: str = "",
        trace: str = "",
    ) -> None:
        """Record that *step_name* failed."""
        now = _now_iso()
        self._conn.execute(
            "UPDATE steps SET state='failed', completed_at=?, error=?, trace=? "
            "WHERE run_id=? AND step_name=? AND id=("
            "  SELECT id FROM steps WHERE run_id=? AND step_name=? "
            "  ORDER BY id DESC LIMIT 1"
            ")",
            (now, error, trace, run_id, step_name, run_id, step_name),
        )
        self._conn.commit()
        self._append({
            "event": "step_failed",
            "run_id": run_id,
            "step": step_name,
            "error": error,
            "trace": trace,
            "timestamp": now,
        })

    def record_step_skipped(self, run_id: str, step_name: str, reason: str = "") -> None:
        """Record that *step_name* was skipped (fingerprint match)."""
        now = _now_iso()
        self._append({
            "event": "step_skipped",
            "run_id": run_id,
            "step": step_name,
            "reason": reason,
            "timestamp": now,
        })

    # ── Resume queries ───────────────────────────────────────────────────────

    def latest_step(self, run_id: str, step_name: str) -> Optional[Dict[str, Any]]:
        """Return the most recent checkpoint record for *step_name* in *run_id*.

        Returns ``None`` when no record exists.
        """
        row = self._conn.execute(
            "SELECT state, input_fingerprint, params, outputs, error, trace, "
            "started_at, completed_at "
            "FROM steps WHERE run_id=? AND step_name=? ORDER BY id DESC LIMIT 1",
            (run_id, step_name),
        ).fetchone()
        if row is None:
            return None
        return dict(row)

    def all_step_states(self, run_id: str) -> Dict[str, str]:
        """Return a mapping of step_name → state for all steps in *run_id*."""
        rows = self._conn.execute(
            "SELECT step_name, state FROM steps WHERE run_id=? "
            "GROUP BY step_name HAVING id=MAX(id)",
            (run_id,),
        ).fetchall()
        return {r["step_name"]: r["state"] for r in rows}

    def latest_run(self) -> Optional[Dict[str, Any]]:
        """Return the most recently started run's metadata, or ``None``."""
        row = self._conn.execute(
            "SELECT * FROM runs ORDER BY start_time DESC LIMIT 1"
        ).fetchone()
        return dict(row) if row else None

    # ── Internal helpers ─────────────────────────────────────────────────────

    def _append(self, record: dict) -> None:
        _atomic_append_ndjson(self._ndjson_path, record)

    def _next_sequence(self, run_id: str, step_name: str) -> int:
        row = self._conn.execute(
            "SELECT COUNT(*) FROM steps WHERE run_id=? AND step_name=?",
            (run_id, step_name),
        ).fetchone()
        return (row[0] if row else 0) + 1

    def _record_artifacts(
        self, run_id: str, step_name: str, paths: List[str], ts: str
    ) -> None:
        for p in paths:
            if not os.path.isfile(p):
                continue
            try:
                sha = file_sha256(p)
                size = os.path.getsize(p)
            except OSError:
                sha, size = "", 0
            self._conn.execute(
                "INSERT INTO artifacts (run_id, step_name, path, sha256, size_bytes, created_at) "
                "VALUES (?, ?, ?, ?, ?, ?)",
                (run_id, step_name, p, sha, size, ts),
            )
        self._conn.commit()

    def _write_config_snapshot(self, cfg: Dict[str, Any], meta: Dict[str, Any]) -> None:
        """Write *cfg* + *meta* to ``OUTDIR/config.yaml`` (canonical) atomically."""
        from copy import deepcopy

        snapshot = deepcopy(cfg)
        snapshot["_run_meta"] = meta

        yaml_path = os.path.join(self._outdir, "config.yaml")
        os.makedirs(self._outdir, exist_ok=True)

        dir_ = os.path.dirname(os.path.abspath(yaml_path))
        fd, tmp = tempfile.mkstemp(dir=dir_, suffix=".yaml")
        try:
            from ruamel.yaml import YAML
            yml = YAML()
            yml.default_flow_style = False
            with os.fdopen(fd, "w") as fh:
                yml.dump(snapshot, fh)
            os.replace(tmp, yaml_path)
        except Exception:
            try:
                os.unlink(tmp)
            except OSError:
                pass
            raise


# ---------------------------------------------------------------------------
# Resume helpers
# ---------------------------------------------------------------------------


def new_run_id() -> str:
    """Generate a new unique run ID."""
    return str(uuid.uuid4())


def build_resume_plan(
    checkpointer: "Checkpointer",
    run_id: str,
    steps: List[str],
    step_fingerprints: Dict[str, str],
    force: bool = False,
) -> Dict[str, str]:
    """Determine the action for each step when resuming.

    Parameters
    ----------
    checkpointer:
        An open :class:`Checkpointer` connected to the previous run.
    run_id:
        The run ID of the previous run to resume.
    steps:
        Ordered list of step names in the pipeline.
    step_fingerprints:
        Mapping of step_name → fingerprint for the *current* pipeline config.
    force:
        When ``True`` every step is marked for re-run regardless of checkpoint.

    Returns
    -------
    dict
        Mapping of step_name → action string:
        ``"skip"`` | ``"run"`` | ``"resume"``
    """
    plan: Dict[str, str] = {}
    for step in steps:
        if force:
            plan[step] = "run"
            continue
        last = checkpointer.latest_step(run_id, step)
        if last is None:
            plan[step] = "run"
        elif last["state"] == "success":
            current_fp = step_fingerprints.get(step)
            if current_fp and last["input_fingerprint"] == current_fp:
                plan[step] = "skip"
            else:
                plan[step] = "run"
        elif last["state"] in ("running", "pending"):
            # Interrupted — attempt resume; implementation falls back to re-run
            plan[step] = "resume"
        else:
            # failed or unknown
            plan[step] = "run"
    return plan


def print_resume_summary(plan: Dict[str, str]) -> None:
    """Print a human-readable summary of the resume plan to stdout."""
    skip = [s for s, a in plan.items() if a == "skip"]
    run = [s for s, a in plan.items() if a == "run"]
    resume = [s for s, a in plan.items() if a == "resume"]

    print("\n[checkpoint] Resume plan:")
    if skip:
        print(f"  SKIP  ({len(skip)}): {', '.join(skip)}")
    if resume:
        print(f"  RESUME({len(resume)}): {', '.join(resume)}")
    if run:
        print(f"  RUN   ({len(run)}): {', '.join(run)}")
    print()

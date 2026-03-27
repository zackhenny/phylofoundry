import os
import glob
import traceback as _traceback
from .artifact_paths import ArtifactPaths
from .constants import STEPS
from .utils.helpers import safe_mkdir, write_json
from .utils.bio import read_fasta
from .tasks import prep, hmmer, extract, embed, phylo, curate, post, synteny, codon, hyphy, asr
from .tasks import taxonomy_integrate, conservation_metrics, detect_clades, aa_composition
from .execution_planner import build_execution_plan
from .execution_schema import StepState
from .failure_policy import apply_failure_policy, build_reverse_dependency_map
from .logging_utils import (
    append_pipeline_log,
    append_step_log,
    ensure_logs_dir,
    execution_plan_path,
    initialize_step_status,
    step_status_path,
    update_step_status,
    write_execution_plan,
)
from .checkpoint import (
    Checkpointer,
    CheckpointMeta,
    build_resume_plan,
    compute_fingerprint,
    load_checkpoint_meta,
    new_run_id,
    print_resume_summary,
)


def step_in_range(step, start_at, stop_after):
    """Check whether `step` falls within [start_at, stop_after]."""
    i = STEPS.index(step)
    i0 = STEPS.index(start_at) if start_at else 0
    i1 = STEPS.index(stop_after) if stop_after else len(STEPS) - 1
    return i0 <= i <= i1


def load_manifest(hmm_manifest_fp):
    if not hmm_manifest_fp:
        return None
    keep = set()
    with open(hmm_manifest_fp) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            keep.add(s)
    return keep


def _load_proteomes_lazy(genomes, faa_dir):
    """Load all proteomes into a dict keyed by genome filename."""
    proteome_seqs = {}
    for g in genomes:
        proteome_seqs[g] = read_fasta(os.path.join(faa_dir, g))
    return proteome_seqs


def _load_proteomes_from_combined_faa(combined_faa_path):
    """Parse a combined proteomes FASTA to reconstruct per-genome protein dicts.

    The combined FASTA produced by the prep step (or provided prebuilt) uses
    headers of the form ``>genome~protein_id``.  This function rebuilds the
    nested ``{genome: {protein_id: sequence}}`` mapping consumed by the
    ``extract`` step.

    Parameters
    ----------
    combined_faa_path : str
        Path to the combined proteomes FASTA file.

    Returns
    -------
    dict[str, dict[str, str]]
        Mapping ``genome → {protein_id → sequence}``.
    """
    flat = read_fasta(combined_faa_path)
    proteome_seqs: dict = {}
    for header, seq in flat.items():
        if "~" in header:
            genome, protein = header.split("~", 1)
        else:
            # Fall back: treat entire header as protein under a synthetic genome key
            genome = "unknown"
            protein = header
        proteome_seqs.setdefault(genome, {})[protein] = seq
    return proteome_seqs


def run_pipeline(cfg):
    # ── Setup ──────────────────────────────────────────────────────────────
    faa_arg = cfg["inputs"]["faa_dir"]
    hmm_arg = cfg["inputs"].get("hmm_input")
    outdir = cfg["output"]["outdir"]

    # Prebuilt GlobDB-compatible inputs (all optional).
    prebuilt_combined_faa = cfg["inputs"].get("combined_faa")
    prebuilt_diamond_db = cfg["inputs"].get("diamond_db")

    # Optional working directory for large intermediate files.
    # When set, combined_proteomes.faa, combined.hmm, DIAMOND databases, and
    # embedding vectors are written here instead of outdir.
    workdir_cfg = cfg["output"].get("workdir")
    work_base = os.path.abspath(workdir_cfg) if workdir_cfg else outdir

    start_at = cfg["workflow"]["start_at"]
    stop_after = cfg["workflow"]["stop_after"]
    force = bool(cfg["workflow"]["force"])

    # ── Detect search mode ────────────────────────────────────────────────
    use_diamond = cfg.get("diamond", {}).get("enabled", False)
    if use_diamond:
        diamond_query = cfg["inputs"].get("diamond_query")
        if not diamond_query:
            raise SystemExit(
                "diamond.enabled=true but inputs.diamond_query is not set. "
                "Provide a path to a FASTA file or directory of FASTA files."
            )

    phy_cfg = cfg["phylo"]
    emb_cfg = cfg["embeddings"]
    post_cfg = cfg["post"]
    tax_int_cfg = cfg.get("taxonomy_integrate", {})
    cons_met_cfg = cfg.get("conservation_metrics", {})
    det_cla_cfg = cfg.get("detect_clades", {})
    aa_comp_cfg = cfg.get("aa_composition", {})
    synteny_cfg = cfg["synteny"]
    codon_cfg = cfg["codon"]
    hyphy_cfg = cfg["hyphy"]
    motif_cfg = cfg.get("motifs", {})
    discover_cfg = cfg.get("discover", {})

    safe_mkdir(outdir)
    if work_base != outdir:
        safe_mkdir(work_base)
        print(f"[pipeline] Working directory: {work_base}")

    # ── Checkpoint / resume setup ─────────────────────────────────────────
    ckpt_cfg = cfg.get("_checkpoint", {})
    _resume = ckpt_cfg.get("resume", False)
    _resume_from = ckpt_cfg.get("resume_from", None)
    _no_resume = ckpt_cfg.get("no_resume", False)

    checkpointer = Checkpointer(outdir)
    resume_meta: "CheckpointMeta | None" = None
    resume_plan: "dict[str, str]" = {}

    if not _no_resume and not force:
        if _resume_from:
            resume_meta = load_checkpoint_meta(_resume_from)
            if resume_meta is None:
                print(
                    f"[checkpoint] WARNING: --resume-from '{_resume_from}' could not "
                    f"be loaded. Starting a fresh run."
                )
        elif _resume:
            resume_meta = load_checkpoint_meta(outdir)
            if resume_meta is None:
                print(
                    f"[checkpoint] WARNING: No checkpoint found in '{outdir}'. "
                    f"Starting a fresh run."
                )

    if resume_meta is not None:
        print(f"[checkpoint] Resuming run '{resume_meta.run_id}' from '{resume_meta.outdir}'")
        run_id = resume_meta.run_id
        # Open the previous checkpointer if the DB is in a different directory
        if os.path.abspath(resume_meta.outdir) != os.path.abspath(outdir):
            prev_checkpointer = Checkpointer(resume_meta.outdir)
        else:
            prev_checkpointer = checkpointer

        # Build per-step fingerprints from current config
        step_fps: dict = {}
        for s in STEPS:
            params_for_fp = {
                "step": s,
                "workflow": {k: v for k, v in cfg["workflow"].items()
                             if k not in ("force",)},
            }
            step_fps[s] = compute_fingerprint(params_for_fp, [])

        resume_plan = build_resume_plan(
            prev_checkpointer, run_id, STEPS, step_fps, force=force
        )
        print_resume_summary(resume_plan)
    else:
        run_id = new_run_id()

    # Start (or resume) the run and write a config snapshot to OUTDIR
    checkpointer.start_run(run_id, cfg)
    print(f"[checkpoint] Run ID: {run_id}")

    # ── Input provenance (--input-run) ────────────────────────────────────
    # When --input-run is set, validate that the prior run produced the
    # artifacts needed by the active step(s) and record full provenance.
    _input_run_dir: "str | None" = ckpt_cfg.get("input_run")
    if _input_run_dir:
        from .provenance import apply_input_run, InputRunError
        _active_step = cfg["workflow"].get("start_at")
        try:
            apply_input_run(cfg, step_name=_active_step)
        except InputRunError as exc:
            raise SystemExit(f"\n[input-run] ERROR: {exc}") from exc

    # Helper: should this step be skipped due to resume?
    def _is_resume_skip(step_name: str) -> bool:
        return resume_plan.get(step_name) == "skip"

    # ── Execution planning and logging scaffold ────────────────────────────
    # Build the plan before any work begins so that the logs/ directory and
    # status files are always present even if an early step fails.
    exec_plan = build_execution_plan(cfg)
    logs_dir = ensure_logs_dir(outdir)
    write_execution_plan(exec_plan, execution_plan_path(logs_dir))
    status_path = step_status_path(logs_dir)
    initialize_step_status(exec_plan, status_path)

    # ── Failure policy setup ───────────────────────────────────────────────
    # Build the reverse dependency map once so that any step failure can
    # propagate blocked states to downstream dependent steps without stopping
    # unrelated branches.
    all_step_names = [s.name for s in exec_plan.steps]
    reverse_map = build_reverse_dependency_map(all_step_names)

    def _is_blocked(step_name: str) -> bool:
        """Return True if step_name has been marked BLOCKED in the plan."""
        for s in exec_plan.steps:
            if s.name == step_name:
                return s.state == StepState.BLOCKED
        return False

    def _log_step_event(step_name: str, event: str) -> None:
        """Write step-level and pipeline-level lifecycle messages."""
        append_pipeline_log(logs_dir, f"{event}: {step_name}")
        append_step_log(logs_dir, step_name, event)

    def _on_step_failure(step_name: str, exc: Exception) -> None:
        """Record failure and propagate BLOCKED state to downstream steps."""
        update_step_status(status_path, step_name, StepState.FAILED.value)
        # Also mark the step as failed in the plan so _is_blocked stays correct.
        for s in exec_plan.steps:
            if s.name == step_name:
                s.state = StepState.FAILED
                break
        msg_fail = f"STEP FAILED: {step_name} — {exc}"
        print(f"\n[pipeline] {msg_fail}")
        append_pipeline_log(logs_dir, msg_fail)
        append_step_log(logs_dir, step_name, f"FAILED: {exc}")
        checkpointer.record_step_failure(
            run_id, step_name, error=str(exc),
            trace=_traceback.format_exc()
        )
        # Record which steps were already blocked before this failure so we
        # only log newly blocked steps.
        previously_blocked = {s.name for s in exec_plan.steps if s.state == StepState.BLOCKED}
        # Propagate blocked state to enabled downstream dependents.
        apply_failure_policy(exec_plan.steps, step_name, reverse_map)
        for s in exec_plan.steps:
            if s.state == StepState.BLOCKED and s.name not in previously_blocked:
                update_step_status(status_path, s.name, StepState.BLOCKED.value)
                blocked_msg = f"BLOCKED: {s.name} (upstream failure: {step_name})"
                print(f"[pipeline] {blocked_msg}")
                append_pipeline_log(logs_dir, blocked_msg)

    # ── Output structure ───────────────────────────────────────────────────
    paths = ArtifactPaths(outdir)

    # When --input-run is specified, read upstream artifacts from the prior
    # run directory instead of the current outdir.
    _input_run_dir_abs = cfg.get("_input_paths", {}).get("input_run_dir")
    input_paths = ArtifactPaths(_input_run_dir_abs) if _input_run_dir_abs else paths

    hmmscan_dir = input_paths.hmmscan_dir
    hmmsearch_dir = input_paths.hmmsearch_dir
    fasta_dir = input_paths.fasta_dir
    raw_fasta_dir = input_paths.fasta_dir
    aln_dir = input_paths.alignments_hmm_dir
    clipkit_dir = input_paths.alignments_clipkit_dir
    tree_dir = input_paths.trees_dir
    summary_dir = input_paths.summary_dir
    motifs_dir = paths.motifs_dir
    discover_dir = paths.discover_dir
    post_dir = input_paths.conservation_metrics_dir
    codon_dir = input_paths.codon_dir
    hyphy_dir = paths.hyphy_dir  # always write hyphy output to new outdir
    # Embedding vectors can be large; place them in work_base when a separate
    # working directory is configured.
    emb_dir = os.path.join(work_base, "embeddings")
    clade_assign_dir = input_paths.clade_assign_dir
    curated_dir = input_paths.curated_dir

    for d in paths.all_output_dirs():
        safe_mkdir(d)
    # Ensure emb_dir exists even when it lives outside outdir
    safe_mkdir(emb_dir)

    write_json(cfg, paths.resolved_config_json)

    hmm_keep = load_manifest(cfg["workflow"]["hmm_manifest"])

    # ── Input discovery ────────────────────────────────────────────────────
    # faa_dir is required unless a prebuilt combined_faa (or diamond_db in
    # DIAMOND mode) has been supplied — in those cases genome FASTAs are
    # not needed for the prep / DB-build steps.
    faa_dir = None
    genomes: list[str] = []
    if faa_arg:
        faa_abs = os.path.abspath(faa_arg)
        if os.path.isfile(faa_abs):
            if not faa_abs.endswith(".faa"):
                raise SystemExit(f"inputs.faa_dir points to a file but not .faa: {faa_abs}")
            faa_dir = os.path.dirname(faa_abs) or "."
            genomes = [os.path.basename(faa_abs)]
        else:
            if not os.path.isdir(faa_abs):
                raise SystemExit(f"inputs.faa_dir must be a directory or a single .faa file: {faa_abs}")
            faa_dir = faa_abs
            genomes = sorted(f for f in os.listdir(faa_dir) if f.endswith(".faa"))
        if not genomes:
            raise SystemExit("No .faa inputs found.")
    elif not prebuilt_combined_faa and not (use_diamond and prebuilt_diamond_db):
        raise SystemExit(
            "inputs.faa_dir is not set and no prebuilt combined_faa or diamond_db "
            "was provided. Supply --faa_dir, --combined_faa, or --diamond_db."
        )

    # HMM input discovery — skipped when DIAMOND mode is active
    hmm_dir = None
    hmm_files = []
    hmm_input_mode = None
    if not use_diamond:
        if not hmm_arg:
            raise SystemExit("inputs.hmm_input is not set.")
        hmm_abs = os.path.abspath(hmm_arg)
        if os.path.isfile(hmm_abs):
            if not hmm_abs.endswith(".hmm"):
                raise SystemExit(f"inputs.hmm_input points to a file but not .hmm: {hmm_abs}")
            hmm_dir = os.path.dirname(hmm_abs) or "."
            hmm_files = [os.path.basename(hmm_abs)]
            hmm_input_mode = "file"
        else:
            if not os.path.isdir(hmm_abs):
                raise SystemExit(f"inputs.hmm_input must be a directory or .hmm file: {hmm_abs}")
            hmm_dir = hmm_abs
            hmm_files = sorted(f for f in os.listdir(hmm_dir) if f.endswith(".hmm"))
            hmm_input_mode = "dir"
        if not hmm_files:
            raise SystemExit("No .hmm files found.")

    # Large intermediate files go to work_base (may equal outdir when no
    # separate working directory is configured).
    # Use prebuilt combined_faa if provided; otherwise the pipeline builds it.
    if prebuilt_combined_faa:
        combined_faa = os.path.abspath(prebuilt_combined_faa)
        print(f"[pipeline] Using prebuilt combined proteomes FASTA: {combined_faa}")
    else:
        combined_faa = os.path.join(work_base, "combined_proteomes.faa")
    combined_hmm = os.path.join(work_base, "combined.hmm")

    # ── STEP: prep ─────────────────────────────────────────────────────────
    if step_in_range("prep", start_at, stop_after):
        if _is_blocked("prep"):
            print("[pipeline] Skipping blocked step: prep")
        elif _is_resume_skip("prep"):
            print("[pipeline] SKIP: prep (checkpoint match)")
            update_step_status(status_path, "prep", "success")
            checkpointer.record_step_skipped(run_id, "prep", "fingerprint match")
        elif prebuilt_combined_faa:
            print(
                "[pipeline] Skipping prep build: using prebuilt combined_faa "
                f"({combined_faa})"
            )
            update_step_status(status_path, "prep", "success")
            _log_step_event("prep", "SKIPPED (prebuilt combined_faa provided)")
            checkpointer.record_step_skipped(run_id, "prep", "prebuilt combined_faa provided")
        else:
            checkpointer.record_step_pending(run_id, "prep")
            update_step_status(status_path, "prep", "running")
            _log_step_event("prep", "START")
            checkpointer.record_step_running(run_id, "prep")
            try:
                if use_diamond:
                    prep.run_prep_diamond_mode(cfg, genomes, faa_dir, combined_faa, force)
                else:
                    prep.run_prep(cfg, genomes, faa_dir, hmm_input_mode, hmm_dir,
                                  hmm_files, combined_faa, combined_hmm, force)
                update_step_status(status_path, "prep", "success")
                _log_step_event("prep", "SUCCESS")
                checkpointer.record_step_success(run_id, "prep")
            except Exception as exc:
                _on_step_failure("prep", exc)
    if stop_after == "prep":
        checkpointer.end_run(run_id, success=True)
        return

    # ── STEP: hmmer (or diamond) ───────────────────────────────────────────
    scan_df = None
    search_df = None
    best_df = None
    if step_in_range("hmmer", start_at, stop_after):
        if _is_blocked("hmmer"):
            print("[pipeline] Skipping blocked step: hmmer")
        elif _is_resume_skip("hmmer"):
            print("[pipeline] SKIP: hmmer (checkpoint match)")
            update_step_status(status_path, "hmmer", "success")
            checkpointer.record_step_skipped(run_id, "hmmer", "fingerprint match")
        else:
            checkpointer.record_step_pending(run_id, "hmmer")
            update_step_status(status_path, "hmmer", "running")
            _log_step_event("hmmer", "START")
            checkpointer.record_step_running(run_id, "hmmer")
            try:
                if use_diamond:
                    from .tasks import diamond as diamond_task
                    scan_df, search_df, best_df = diamond_task.run_diamond(
                        cfg, genomes, combined_faa, work_base, summary_dir, force
                    )
                else:
                    scan_df, search_df, best_df = hmmer.run_hmmer(
                        cfg, genomes, faa_dir, hmm_files, hmm_dir, combined_hmm,
                        combined_faa, outdir, summary_dir, hmmscan_dir, hmmsearch_dir,
                        hmm_keep, force,
                        resume=(resume_plan.get("hmmer") == "resume"),
                    )
                # Optional: clean up combined FAA to reclaim disk space
                prep_cfg = cfg.get("prep", {})
                if prep_cfg.get("cleanup_combined_faa", False) and os.path.exists(combined_faa):
                    os.remove(combined_faa)
                    print("[pipeline] Removed combined_proteomes.faa (prep.cleanup_combined_faa=true).")
                update_step_status(status_path, "hmmer", "success")
                _log_step_event("hmmer", "SUCCESS")
                checkpointer.record_step_success(run_id, "hmmer")
            except Exception as exc:
                _on_step_failure("hmmer", exc)
    if stop_after == "hmmer":
        checkpointer.end_run(run_id, success=True)
        return

    # ── Name mapping for hmmalign ──────────────────────────────────────────
    # In DIAMOND mode there are no HMM profiles, so this mapping is empty.
    # The phylo step will use MAFFT alignment in that case.
    name_to_hmm_path = {}
    if not use_diamond:
        for hf in hmm_files:
            fp = os.path.join(hmm_dir, hf)
            try:
                with open(fp) as f:
                    for line in f:
                        if line.startswith("NAME"):
                            name_to_hmm_path[line.split()[1]] = fp
                            break
            except OSError:
                continue

    # ── Helper: load hit DataFrames from disk if we skipped hmmer ──────────
    def _ensure_hit_dfs():
        nonlocal scan_df, search_df
        if scan_df is not None and search_df is not None:
            return
        import pandas as pd
        if use_diamond:
            hits_diamond_tsv = os.path.join(summary_dir, "diamond_hits.filtered.tsv")
            search_df = pd.read_csv(hits_diamond_tsv, sep="\t") if os.path.exists(hits_diamond_tsv) else pd.DataFrame()
            scan_df = pd.DataFrame()
            return
        hits_scan_tsv = os.path.join(summary_dir, "hmmscan_hits.filtered.tsv")
        hits_search_tsv = os.path.join(summary_dir, "hmmsearch_hits.filtered.tsv")
        scan_df = pd.read_csv(hits_scan_tsv, sep="\t") if os.path.exists(hits_scan_tsv) else pd.DataFrame()
        search_df = pd.read_csv(hits_search_tsv, sep="\t") if os.path.exists(hits_search_tsv) else pd.DataFrame()

    # ── STEP: extract ──────────────────────────────────────────────────────
    hmm_to_seqs = {}
    if step_in_range("extract", start_at, stop_after):
        if _is_blocked("extract"):
            print("[pipeline] Skipping blocked step: extract")
        elif _is_resume_skip("extract"):
            print("[pipeline] SKIP: extract (checkpoint match)")
            update_step_status(status_path, "extract", "success")
            checkpointer.record_step_skipped(run_id, "extract", "fingerprint match")
        else:
            checkpointer.record_step_pending(run_id, "extract")
            update_step_status(status_path, "extract", "running")
            _log_step_event("extract", "START")
            checkpointer.record_step_running(run_id, "extract")
            try:
                _ensure_hit_dfs()
                # When per-genome FAA files are available, load proteins from
                # them directly.  When only a prebuilt combined_faa was
                # supplied (faa_dir is None), parse it instead.
                if faa_dir and genomes:
                    proteome_seqs = _load_proteomes_lazy(genomes, faa_dir)
                else:
                    print(
                        "[pipeline] faa_dir not set; loading proteome sequences "
                        f"from combined_faa: {combined_faa}"
                    )
                    proteome_seqs = _load_proteomes_from_combined_faa(combined_faa)
                hmm_to_seqs = extract.run_extract(
                    cfg, scan_df, search_df, fasta_dir, hmm_keep, proteome_seqs, force
                )
                del proteome_seqs  # free memory after extraction
                update_step_status(status_path, "extract", "success")
                _log_step_event("extract", "SUCCESS")
                checkpointer.record_step_success(run_id, "extract")
            except Exception as exc:
                _on_step_failure("extract", exc)
    if stop_after == "extract":
        checkpointer.end_run(run_id, success=True)
        return

    # ── Taxonomy map (used by embed + post) ──────────────────────────────
    tax_map = {}
    if (cfg["inputs"].get("gtdb_dir") or cfg["inputs"].get("taxonomy_file")
            or cfg["inputs"].get("globdb_taxonomy_file")):
        tax_map = post._load_taxonomy(
            cfg["inputs"].get("gtdb_dir"),
            cfg["inputs"].get("taxonomy_file"),
            cfg["inputs"].get("globdb_taxonomy_file"),
        )
        if tax_map:
            print(f"[pipeline] Loaded taxonomy for {len(tax_map)} genomes.")

    # ── STEP: embed ────────────────────────────────────────────────────────
    if step_in_range("embed", start_at, stop_after) and emb_cfg.get("enabled", False):
        if _is_blocked("embed"):
            print("[pipeline] Skipping blocked step: embed")
        elif _is_resume_skip("embed"):
            print("[pipeline] SKIP: embed (checkpoint match)")
            update_step_status(status_path, "embed", "success")
            checkpointer.record_step_skipped(run_id, "embed", "fingerprint match")
        else:
            checkpointer.record_step_pending(run_id, "embed")
            update_step_status(status_path, "embed", "running")
            _log_step_event("embed", "START")
            checkpointer.record_step_running(run_id, "embed")
            try:
                clades = None
                if post_cfg.get("clades_tsv", None):
                    try:
                        clades = post.load_clades_tsv(post_cfg["clades_tsv"])
                    except Exception:
                        clades = None
                embed.run_embed(cfg, hmm_to_seqs, clades, emb_dir, fasta_dir, hmm_keep,
                                force, summary_dir=summary_dir, tax_map=tax_map,
                                resume=(resume_plan.get("embed") == "resume"))
                update_step_status(status_path, "embed", "success")
                _log_step_event("embed", "SUCCESS")
                checkpointer.record_step_success(run_id, "embed")
            except Exception as exc:
                _on_step_failure("embed", exc)
    if stop_after == "embed":
        checkpointer.end_run(run_id, success=True)
        return

    # ── STEP: phylo ────────────────────────────────────────────────────────
    # In DIAMOND mode there are no HMM profiles for hmmalign, so MAFFT must be
    # used for alignment.  Auto-enable it here so users don't need to set it.
    if use_diamond:
        phy_cfg = dict(phy_cfg)  # shallow copy to avoid mutating shared cfg
        phy_cfg["mafft"] = True
        cfg = dict(cfg)
        cfg["phylo"] = phy_cfg

    if step_in_range("phylo", start_at, stop_after):
        if _is_blocked("phylo"):
            print("[pipeline] Skipping blocked step: phylo")
        elif _is_resume_skip("phylo"):
            print("[pipeline] SKIP: phylo (checkpoint match)")
            update_step_status(status_path, "phylo", "success")
            checkpointer.record_step_skipped(run_id, "phylo", "fingerprint match")
        else:
            checkpointer.record_step_pending(run_id, "phylo")
            update_step_status(status_path, "phylo", "running")
            _log_step_event("phylo", "START")
            checkpointer.record_step_running(run_id, "phylo")
            try:
                phylo.run_phylo(cfg, hmm_to_seqs, fasta_dir, aln_dir, clipkit_dir,
                                tree_dir, name_to_hmm_path, hmm_keep, force,
                                resume=(resume_plan.get("phylo") == "resume"))

                # ── ASR parsing (after phylo, if ASR was enabled) ─────────────────
                if not phy_cfg.get("no_asr", False):
                    all_ancestral = asr.run_asr_parse(cfg, tree_dir, fasta_dir,
                                                       hmm_keep, force)

                    # Combined embedding with ancestors if available
                    if all_ancestral and emb_cfg.get("enabled", False):
                        # Merge all ancestral seqs across HMMs
                        merged_ancestral = {}
                        for hmm_name, anc_seqs in all_ancestral.items():
                            for node_id, seq in anc_seqs.items():
                                merged_ancestral[f"{hmm_name}_{node_id}"] = seq

                        if merged_ancestral:
                            clades = None
                            if post_cfg.get("clades_tsv", None):
                                try:
                                    clades = post.load_clades_tsv(post_cfg["clades_tsv"])
                                except Exception:
                                    clades = None
                            embed.embed_combined_with_ancestors(
                                cfg, hmm_to_seqs, merged_ancestral, clades,
                                emb_dir, fasta_dir, hmm_keep,
                                force=force, summary_dir=summary_dir, tax_map=tax_map
                            )

                update_step_status(status_path, "phylo", "success")
                _log_step_event("phylo", "SUCCESS")
                checkpointer.record_step_success(run_id, "phylo")
            except Exception as exc:
                _on_step_failure("phylo", exc)

    if stop_after == "phylo":
        checkpointer.end_run(run_id, success=True)
        return

    # ── STEP: curate ───────────────────────────────────────────────────────
    if step_in_range("curate", start_at, stop_after):
        if _is_blocked("curate"):
            print("[pipeline] Skipping blocked step: curate")
        elif _is_resume_skip("curate"):
            print("[pipeline] SKIP: curate (checkpoint match)")
            update_step_status(status_path, "curate", "success")
            checkpointer.record_step_skipped(run_id, "curate", "fingerprint match")
        else:
            checkpointer.record_step_pending(run_id, "curate")
            update_step_status(status_path, "curate", "running")
            _log_step_event("curate", "START")
            checkpointer.record_step_running(run_id, "curate")
            try:
                curate.run_curate(cfg, tree_dir, fasta_dir, clipkit_dir, emb_dir,
                                  summary_dir, hmm_keep, force,
                                  curated_dir=curated_dir)
                update_step_status(status_path, "curate", "success")
                _log_step_event("curate", "SUCCESS")
                checkpointer.record_step_success(run_id, "curate")
            except Exception as exc:
                _on_step_failure("curate", exc)
    if stop_after == "curate":
        checkpointer.end_run(run_id, success=True)
        return

    # ── Prefer curated overlay outputs for downstream steps ────────────────
    # If the curate step produced overlay artifacts, redirect downstream
    # directory pointers so that post, codon, hyphy, etc. consume the curated
    # (and potentially pruned) versions rather than the raw pipeline outputs.
    _cur_tree_dir    = os.path.join(curated_dir, "trees")
    _cur_fasta_dir   = os.path.join(curated_dir, "fasta_per_hmm")
    _cur_clipkit_dir = os.path.join(curated_dir, "alignments_clipkit")
    if os.path.isdir(_cur_tree_dir) and glob.glob(os.path.join(_cur_tree_dir, "*.treefile")):
        tree_dir    = _cur_tree_dir
        fasta_dir   = _cur_fasta_dir
        clipkit_dir = _cur_clipkit_dir

    # ── STEP: taxonomy_integrate ───────────────────────────────────────────
    if step_in_range("taxonomy_integrate", start_at, stop_after) and tax_int_cfg.get("enabled", False):
        if _is_blocked("taxonomy_integrate"):
            print("[pipeline] Skipping blocked step: taxonomy_integrate")
        elif _is_resume_skip("taxonomy_integrate"):
            print("[pipeline] SKIP: taxonomy_integrate (checkpoint match)")
            update_step_status(status_path, "taxonomy_integrate", "success")
            checkpointer.record_step_skipped(run_id, "taxonomy_integrate", "fingerprint match")
        else:
            checkpointer.record_step_pending(run_id, "taxonomy_integrate")
            update_step_status(status_path, "taxonomy_integrate", "running")
            _log_step_event("taxonomy_integrate", "START")
            checkpointer.record_step_running(run_id, "taxonomy_integrate")
            try:
                taxonomy_integrate.run_taxonomy_integrate(cfg, summary_dir)
                update_step_status(status_path, "taxonomy_integrate", "success")
                _log_step_event("taxonomy_integrate", "SUCCESS")
                checkpointer.record_step_success(run_id, "taxonomy_integrate")
            except Exception as exc:
                _on_step_failure("taxonomy_integrate", exc)
    if stop_after == "taxonomy_integrate":
        checkpointer.end_run(run_id, success=True)
        return

    # ── STEP: conservation_metrics ─────────────────────────────────────────
    if step_in_range("conservation_metrics", start_at, stop_after) and cons_met_cfg.get("enabled", False):
        if _is_blocked("conservation_metrics"):
            print("[pipeline] Skipping blocked step: conservation_metrics")
        elif _is_resume_skip("conservation_metrics"):
            print("[pipeline] SKIP: conservation_metrics (checkpoint match)")
            update_step_status(status_path, "conservation_metrics", "success")
            checkpointer.record_step_skipped(run_id, "conservation_metrics", "fingerprint match")
        else:
            checkpointer.record_step_pending(run_id, "conservation_metrics")
            update_step_status(status_path, "conservation_metrics", "running")
            _log_step_event("conservation_metrics", "START")
            checkpointer.record_step_running(run_id, "conservation_metrics")
            try:
                conservation_metrics.run_conservation_metrics(
                    cfg, tree_dir, clipkit_dir, aln_dir, post_dir, hmm_keep
                )
                update_step_status(status_path, "conservation_metrics", "success")
                _log_step_event("conservation_metrics", "SUCCESS")
                checkpointer.record_step_success(run_id, "conservation_metrics")
            except Exception as exc:
                _on_step_failure("conservation_metrics", exc)
    if stop_after == "conservation_metrics":
        checkpointer.end_run(run_id, success=True)
        return

    # ── STEP: detect_clades ────────────────────────────────────────────────
    if step_in_range("detect_clades", start_at, stop_after) and det_cla_cfg.get("enabled", False):
        if _is_blocked("detect_clades"):
            print("[pipeline] Skipping blocked step: detect_clades")
        elif _is_resume_skip("detect_clades"):
            print("[pipeline] SKIP: detect_clades (checkpoint match)")
            update_step_status(status_path, "detect_clades", "success")
            checkpointer.record_step_skipped(run_id, "detect_clades", "fingerprint match")
        else:
            checkpointer.record_step_pending(run_id, "detect_clades")
            update_step_status(status_path, "detect_clades", "running")
            _log_step_event("detect_clades", "START")
            checkpointer.record_step_running(run_id, "detect_clades")
            try:
                detect_clades.run_detect_clades(
                    cfg, tree_dir, emb_dir, summary_dir, hmm_keep, clade_assign_dir
                )
                update_step_status(status_path, "detect_clades", "success")
                _log_step_event("detect_clades", "SUCCESS")
                checkpointer.record_step_success(run_id, "detect_clades")
            except Exception as exc:
                _on_step_failure("detect_clades", exc)
    if stop_after == "detect_clades":
        checkpointer.end_run(run_id, success=True)
        return

    # ── STEP: aa_composition ───────────────────────────────────────────────
    if step_in_range("aa_composition", start_at, stop_after) and aa_comp_cfg.get("enabled", False):
        if _is_blocked("aa_composition"):
            print("[pipeline] Skipping blocked step: aa_composition")
        elif _is_resume_skip("aa_composition"):
            print("[pipeline] SKIP: aa_composition (checkpoint match)")
            update_step_status(status_path, "aa_composition", "success")
            checkpointer.record_step_skipped(run_id, "aa_composition", "fingerprint match")
        else:
            checkpointer.record_step_pending(run_id, "aa_composition")
            update_step_status(status_path, "aa_composition", "running")
            _log_step_event("aa_composition", "START")
            checkpointer.record_step_running(run_id, "aa_composition")
            try:
                # AA composition is computed on raw post-filter hits, not
                # curate overlay outputs.
                aa_composition.run_aa_composition(
                    cfg, raw_fasta_dir, summary_dir,
                    hmm_keep=hmm_keep,
                    clade_assign_dir=clade_assign_dir,
                    force=force,
                )
                update_step_status(status_path, "aa_composition", "success")
                _log_step_event("aa_composition", "SUCCESS")
                checkpointer.record_step_success(run_id, "aa_composition")
            except Exception as exc:
                _on_step_failure("aa_composition", exc)
    if stop_after == "aa_composition":
        checkpointer.end_run(run_id, success=True)
        return

    # ── STEP: post ─────────────────────────────────────────────────────────
    if step_in_range("post", start_at, stop_after) and post_cfg.get("enabled", False):
        if _is_blocked("post"):
            print("[pipeline] Skipping blocked step: post")
        elif _is_resume_skip("post"):
            print("[pipeline] SKIP: post (checkpoint match)")
            update_step_status(status_path, "post", "success")
            checkpointer.record_step_skipped(run_id, "post", "fingerprint match")
        else:
            checkpointer.record_step_pending(run_id, "post")
            update_step_status(status_path, "post", "running")
            _log_step_event("post", "START")
            checkpointer.record_step_running(run_id, "post")
            try:
                post.run_post(cfg, tree_dir, clipkit_dir, aln_dir, post_dir, summary_dir, hmm_keep, force,
                              clade_assign_dir=clade_assign_dir)
                update_step_status(status_path, "post", "success")
                _log_step_event("post", "SUCCESS")
                checkpointer.record_step_success(run_id, "post")
            except Exception as exc:
                _on_step_failure("post", exc)
    if stop_after == "post":
        checkpointer.end_run(run_id, success=True)
        return

    # ── STEP: synteny ──────────────────────────────────────────────────────
    synteny_dir = os.path.join(outdir, "synteny")
    safe_mkdir(synteny_dir)

    if step_in_range("synteny", start_at, stop_after) and synteny_cfg.get("enabled", False):
        if _is_blocked("synteny"):
            print("[pipeline] Skipping blocked step: synteny")
        elif _is_resume_skip("synteny"):
            print("[pipeline] SKIP: synteny (checkpoint match)")
            update_step_status(status_path, "synteny", "success")
            checkpointer.record_step_skipped(run_id, "synteny", "fingerprint match")
        else:
            checkpointer.record_step_pending(run_id, "synteny")
            update_step_status(status_path, "synteny", "running")
            _log_step_event("synteny", "START")
            checkpointer.record_step_running(run_id, "synteny")
            try:
                _ensure_hit_dfs()
                synteny.run_synteny(cfg, synteny_dir, tree_dir, scan_df, search_df, hmm_keep, force,
                                    clade_assign_dir=clade_assign_dir)
                update_step_status(status_path, "synteny", "success")
                _log_step_event("synteny", "SUCCESS")
                checkpointer.record_step_success(run_id, "synteny")
            except Exception as exc:
                _on_step_failure("synteny", exc)
    if stop_after == "synteny":
        checkpointer.end_run(run_id, success=True)
        return

    # ── STEP: codon ────────────────────────────────────────────────────────
    if step_in_range("codon", start_at, stop_after) and codon_cfg.get("enabled", False):
        if _is_blocked("codon"):
            print("[pipeline] Skipping blocked step: codon")
        elif _is_resume_skip("codon"):
            print("[pipeline] SKIP: codon (checkpoint match)")
            update_step_status(status_path, "codon", "success")
            checkpointer.record_step_skipped(run_id, "codon", "fingerprint match")
        else:
            checkpointer.record_step_pending(run_id, "codon")
            update_step_status(status_path, "codon", "running")
            _log_step_event("codon", "START")
            checkpointer.record_step_running(run_id, "codon")
            try:
                codon.run_codon(cfg, tree_dir, clipkit_dir, aln_dir, codon_dir, hmm_keep, force)
                update_step_status(status_path, "codon", "success")
                _log_step_event("codon", "SUCCESS")
                checkpointer.record_step_success(run_id, "codon")
            except Exception as exc:
                _on_step_failure("codon", exc)
    if stop_after == "codon":
        checkpointer.end_run(run_id, success=True)
        return

    # ── STEP: hyphy ────────────────────────────────────────────────────────
    if step_in_range("hyphy", start_at, stop_after) and hyphy_cfg.get("enabled", False):
        if _is_blocked("hyphy"):
            print("[pipeline] Skipping blocked step: hyphy")
        elif _is_resume_skip("hyphy"):
            print("[pipeline] SKIP: hyphy (checkpoint match)")
            update_step_status(status_path, "hyphy", "success")
            checkpointer.record_step_skipped(run_id, "hyphy", "fingerprint match")
        else:
            checkpointer.record_step_pending(run_id, "hyphy")
            update_step_status(status_path, "hyphy", "running")
            _log_step_event("hyphy", "START")
            checkpointer.record_step_running(run_id, "hyphy")
            try:
                hyphy.run_hyphy(cfg, codon_dir, tree_dir, hyphy_dir, hmm_keep, force,
                                clade_assign_dir=clade_assign_dir)
                update_step_status(status_path, "hyphy", "success")
                _log_step_event("hyphy", "SUCCESS")
                checkpointer.record_step_success(run_id, "hyphy")
            except Exception as exc:
                _on_step_failure("hyphy", exc)
    if stop_after == "hyphy":
        checkpointer.end_run(run_id, success=True)
        return

    # ── STEP: score_motifs ─────────────────────────────────────────────────
    if step_in_range("score_motifs", start_at, stop_after) and motif_cfg.get("enabled", False):
        if _is_blocked("score_motifs"):
            print("[pipeline] Skipping blocked step: score_motifs")
        elif _is_resume_skip("score_motifs"):
            print("[pipeline] SKIP: score_motifs (checkpoint match)")
            update_step_status(status_path, "score_motifs", "success")
            checkpointer.record_step_skipped(run_id, "score_motifs", "fingerprint match")
        else:
            checkpointer.record_step_pending(run_id, "score_motifs")
            update_step_status(status_path, "score_motifs", "running")
            _log_step_event("score_motifs", "START")
            checkpointer.record_step_running(run_id, "score_motifs")
            try:
                from .tasks import score_motifs
                score_motifs.score_motifs(cfg, fasta_dir, summary_dir, motifs_dir, hmm_keep, force)
                update_step_status(status_path, "score_motifs", "success")
                _log_step_event("score_motifs", "SUCCESS")
                checkpointer.record_step_success(run_id, "score_motifs")
            except Exception as exc:
                _on_step_failure("score_motifs", exc)
    if stop_after == "score_motifs":
        checkpointer.end_run(run_id, success=True)
        return

    # ── STEP: discover_motifs ──────────────────────────────────────────────
    if step_in_range("discover_motifs", start_at, stop_after) and discover_cfg.get("enabled", False):
        if _is_blocked("discover_motifs"):
            print("[pipeline] Skipping blocked step: discover_motifs")
        elif _is_resume_skip("discover_motifs"):
            print("[pipeline] SKIP: discover_motifs (checkpoint match)")
            update_step_status(status_path, "discover_motifs", "success")
            checkpointer.record_step_skipped(run_id, "discover_motifs", "fingerprint match")
        else:
            checkpointer.record_step_pending(run_id, "discover_motifs")
            update_step_status(status_path, "discover_motifs", "running")
            _log_step_event("discover_motifs", "START")
            checkpointer.record_step_running(run_id, "discover_motifs")
            try:
                from .tasks import discover_motifs
                discover_motifs.discover_motifs(cfg, fasta_dir, summary_dir, discover_dir, hmm_keep, force,
                                         clade_assign_dir=clade_assign_dir)
                update_step_status(status_path, "discover_motifs", "success")
                _log_step_event("discover_motifs", "SUCCESS")
                checkpointer.record_step_success(run_id, "discover_motifs")
            except Exception as exc:
                _on_step_failure("discover_motifs", exc)

    checkpointer.end_run(run_id, success=True)
    print("\nPipeline complete.")

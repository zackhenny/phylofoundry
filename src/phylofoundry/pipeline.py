import os
import glob
from .constants import STEPS
from .utils.helpers import safe_mkdir, write_json
from .utils.bio import read_fasta
from .tasks import prep, hmmer, extract, embed, phylo, curate, post, synteny, codon, hyphy, asr
from .tasks import taxonomy_integrate, conservation_metrics, detect_clades
from .execution_planner import build_execution_plan
from .execution_schema import StepState
from .failure_policy import apply_failure_policy, build_reverse_dependency_map
from .logging_utils import (
    append_pipeline_log,
    ensure_logs_dir,
    execution_plan_path,
    initialize_step_status,
    step_status_path,
    update_step_status,
    write_execution_plan,
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


def run_pipeline(cfg):
    # ── Setup ──────────────────────────────────────────────────────────────
    faa_arg = cfg["inputs"]["faa_dir"]
    hmm_arg = cfg["inputs"]["hmm_input"]
    outdir = cfg["output"]["outdir"]

    start_at = cfg["workflow"]["start_at"]
    stop_after = cfg["workflow"]["stop_after"]
    force = bool(cfg["workflow"]["force"])

    phy_cfg = cfg["phylo"]
    emb_cfg = cfg["embeddings"]
    post_cfg = cfg["post"]
    tax_int_cfg = cfg.get("taxonomy_integrate", {})
    cons_met_cfg = cfg.get("conservation_metrics", {})
    det_cla_cfg = cfg.get("detect_clades", {})
    synteny_cfg = cfg["synteny"]
    codon_cfg = cfg["codon"]
    hyphy_cfg = cfg["hyphy"]
    motif_cfg = cfg.get("motifs", {})
    discover_cfg = cfg.get("discover", {})

    safe_mkdir(outdir)

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
    hmmscan_dir = os.path.join(outdir, "hmmscan_tbl")
    hmmsearch_dir = os.path.join(outdir, "hmmsearch_tbl")
    fasta_dir = os.path.join(outdir, "fasta_per_hmm")
    aln_dir = os.path.join(outdir, "alignments_hmm")
    clipkit_dir = os.path.join(outdir, "alignments_clipkit")
    tree_dir = os.path.join(outdir, "trees_iqtree")
    summary_dir = os.path.join(outdir, "summary")
    post_dir = os.path.join(summary_dir, "post_scikitbio")
    codon_dir = os.path.join(outdir, "codon_alignments")
    hyphy_dir = os.path.join(summary_dir, "hyphy")
    emb_dir = os.path.join(outdir, "embeddings")
    clade_assign_dir = os.path.join(outdir, "clade_assignments")
    curated_dir = os.path.join(outdir, "curated")

    for d in [hmmscan_dir, hmmsearch_dir, fasta_dir, aln_dir,
              clipkit_dir, tree_dir, summary_dir, post_dir,
              codon_dir, hyphy_dir, emb_dir, clade_assign_dir, curated_dir]:
        safe_mkdir(d)

    write_json(cfg, os.path.join(summary_dir, "resolved_config.json"))

    hmm_keep = load_manifest(cfg["workflow"]["hmm_manifest"])

    # ── Input discovery ────────────────────────────────────────────────────
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

    combined_faa = os.path.join(outdir, "combined_proteomes.faa")
    combined_hmm = os.path.join(outdir, "combined.hmm")

    # ── STEP: prep ─────────────────────────────────────────────────────────
    if step_in_range("prep", start_at, stop_after):
        if _is_blocked("prep"):
            print("[pipeline] Skipping blocked step: prep")
        else:
            update_step_status(status_path, "prep", "running")
            append_pipeline_log(logs_dir, "START: prep")
            try:
                prep.run_prep(cfg, genomes, faa_dir, hmm_input_mode, hmm_dir,
                              hmm_files, combined_faa, combined_hmm, force)
                update_step_status(status_path, "prep", "success")
                append_pipeline_log(logs_dir, "SUCCESS: prep")
            except Exception as exc:
                _on_step_failure("prep", exc)
    if stop_after == "prep":
        return

    # ── STEP: hmmer ────────────────────────────────────────────────────────
    scan_df = None
    search_df = None
    best_df = None
    if step_in_range("hmmer", start_at, stop_after):
        if _is_blocked("hmmer"):
            print("[pipeline] Skipping blocked step: hmmer")
        else:
            update_step_status(status_path, "hmmer", "running")
            append_pipeline_log(logs_dir, "START: hmmer")
            try:
                scan_df, search_df, best_df = hmmer.run_hmmer(
                    cfg, genomes, faa_dir, hmm_files, hmm_dir, combined_hmm,
                    combined_faa, outdir, summary_dir, hmmscan_dir, hmmsearch_dir,
                    hmm_keep, force
                )
                # Optional: clean up combined FAA to reclaim disk space
                prep_cfg = cfg.get("prep", {})
                if prep_cfg.get("cleanup_combined_faa", False) and os.path.exists(combined_faa):
                    os.remove(combined_faa)
                    print("[pipeline] Removed combined_proteomes.faa (prep.cleanup_combined_faa=true).")
                update_step_status(status_path, "hmmer", "success")
                append_pipeline_log(logs_dir, "SUCCESS: hmmer")
            except Exception as exc:
                _on_step_failure("hmmer", exc)
    if stop_after == "hmmer":
        return

    # ── Name mapping for hmmalign ──────────────────────────────────────────
    name_to_hmm_path = {}
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
        hits_scan_tsv = os.path.join(summary_dir, "hmmscan_hits.filtered.tsv")
        hits_search_tsv = os.path.join(summary_dir, "hmmsearch_hits.filtered.tsv")
        scan_df = pd.read_csv(hits_scan_tsv, sep="\t") if os.path.exists(hits_scan_tsv) else pd.DataFrame()
        search_df = pd.read_csv(hits_search_tsv, sep="\t") if os.path.exists(hits_search_tsv) else pd.DataFrame()

    # ── STEP: extract ──────────────────────────────────────────────────────
    hmm_to_seqs = {}
    if step_in_range("extract", start_at, stop_after):
        if _is_blocked("extract"):
            print("[pipeline] Skipping blocked step: extract")
        else:
            update_step_status(status_path, "extract", "running")
            append_pipeline_log(logs_dir, "START: extract")
            try:
                _ensure_hit_dfs()
                proteome_seqs = _load_proteomes_lazy(genomes, faa_dir)
                hmm_to_seqs = extract.run_extract(
                    cfg, scan_df, search_df, fasta_dir, hmm_keep, proteome_seqs, force
                )
                del proteome_seqs  # free memory after extraction
                update_step_status(status_path, "extract", "success")
                append_pipeline_log(logs_dir, "SUCCESS: extract")
            except Exception as exc:
                _on_step_failure("extract", exc)
    if stop_after == "extract":
        return

    # ── Taxonomy map (used by embed + post) ──────────────────────────────
    tax_map = {}
    if cfg["inputs"].get("gtdb_dir") or cfg["inputs"].get("taxonomy_file"):
        tax_map = post._load_taxonomy(
            cfg["inputs"].get("gtdb_dir"),
            cfg["inputs"].get("taxonomy_file"),
        )
        if tax_map:
            print(f"[pipeline] Loaded taxonomy for {len(tax_map)} genomes.")

    # ── STEP: embed ────────────────────────────────────────────────────────
    if step_in_range("embed", start_at, stop_after) and emb_cfg.get("enabled", False):
        if _is_blocked("embed"):
            print("[pipeline] Skipping blocked step: embed")
        else:
            update_step_status(status_path, "embed", "running")
            append_pipeline_log(logs_dir, "START: embed")
            try:
                clades = None
                if post_cfg.get("clades_tsv", None):
                    try:
                        clades = post.load_clades_tsv(post_cfg["clades_tsv"])
                    except Exception:
                        clades = None
                embed.run_embed(cfg, hmm_to_seqs, clades, emb_dir, fasta_dir, hmm_keep,
                                force, summary_dir=summary_dir, tax_map=tax_map)
                update_step_status(status_path, "embed", "success")
                append_pipeline_log(logs_dir, "SUCCESS: embed")
            except Exception as exc:
                _on_step_failure("embed", exc)
    if stop_after == "embed":
        return

    # ── STEP: phylo ────────────────────────────────────────────────────────
    if step_in_range("phylo", start_at, stop_after):
        if _is_blocked("phylo"):
            print("[pipeline] Skipping blocked step: phylo")
        else:
            update_step_status(status_path, "phylo", "running")
            append_pipeline_log(logs_dir, "START: phylo")
            try:
                phylo.run_phylo(cfg, hmm_to_seqs, fasta_dir, aln_dir, clipkit_dir,
                                tree_dir, name_to_hmm_path, hmm_keep, force)

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
                append_pipeline_log(logs_dir, "SUCCESS: phylo")
            except Exception as exc:
                _on_step_failure("phylo", exc)

    if stop_after == "phylo":
        return

    # ── STEP: curate ───────────────────────────────────────────────────────
    if step_in_range("curate", start_at, stop_after):
        if _is_blocked("curate"):
            print("[pipeline] Skipping blocked step: curate")
        else:
            update_step_status(status_path, "curate", "running")
            append_pipeline_log(logs_dir, "START: curate")
            try:
                curate.run_curate(cfg, tree_dir, fasta_dir, clipkit_dir, emb_dir,
                                  summary_dir, hmm_keep, force,
                                  curated_dir=curated_dir)
                update_step_status(status_path, "curate", "success")
                append_pipeline_log(logs_dir, "SUCCESS: curate")
            except Exception as exc:
                _on_step_failure("curate", exc)
    if stop_after == "curate":
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
        else:
            update_step_status(status_path, "taxonomy_integrate", "running")
            append_pipeline_log(logs_dir, "START: taxonomy_integrate")
            try:
                taxonomy_integrate.run_taxonomy_integrate(cfg, summary_dir)
                update_step_status(status_path, "taxonomy_integrate", "success")
                append_pipeline_log(logs_dir, "SUCCESS: taxonomy_integrate")
            except Exception as exc:
                _on_step_failure("taxonomy_integrate", exc)
    if stop_after == "taxonomy_integrate":
        return

    # ── STEP: conservation_metrics ─────────────────────────────────────────
    if step_in_range("conservation_metrics", start_at, stop_after) and cons_met_cfg.get("enabled", False):
        if _is_blocked("conservation_metrics"):
            print("[pipeline] Skipping blocked step: conservation_metrics")
        else:
            update_step_status(status_path, "conservation_metrics", "running")
            append_pipeline_log(logs_dir, "START: conservation_metrics")
            try:
                conservation_metrics.run_conservation_metrics(
                    cfg, tree_dir, clipkit_dir, aln_dir, post_dir, hmm_keep
                )
                update_step_status(status_path, "conservation_metrics", "success")
                append_pipeline_log(logs_dir, "SUCCESS: conservation_metrics")
            except Exception as exc:
                _on_step_failure("conservation_metrics", exc)
    if stop_after == "conservation_metrics":
        return

    # ── STEP: detect_clades ────────────────────────────────────────────────
    if step_in_range("detect_clades", start_at, stop_after) and det_cla_cfg.get("enabled", False):
        if _is_blocked("detect_clades"):
            print("[pipeline] Skipping blocked step: detect_clades")
        else:
            update_step_status(status_path, "detect_clades", "running")
            append_pipeline_log(logs_dir, "START: detect_clades")
            try:
                detect_clades.run_detect_clades(
                    cfg, tree_dir, emb_dir, summary_dir, hmm_keep, clade_assign_dir
                )
                update_step_status(status_path, "detect_clades", "success")
                append_pipeline_log(logs_dir, "SUCCESS: detect_clades")
            except Exception as exc:
                _on_step_failure("detect_clades", exc)
    if stop_after == "detect_clades":
        return

    # ── STEP: post ─────────────────────────────────────────────────────────
    if step_in_range("post", start_at, stop_after) and post_cfg.get("enabled", False):
        if _is_blocked("post"):
            print("[pipeline] Skipping blocked step: post")
        else:
            update_step_status(status_path, "post", "running")
            append_pipeline_log(logs_dir, "START: post")
            try:
                post.run_post(cfg, tree_dir, clipkit_dir, aln_dir, post_dir, summary_dir, hmm_keep, force,
                              clade_assign_dir=clade_assign_dir)
                update_step_status(status_path, "post", "success")
                append_pipeline_log(logs_dir, "SUCCESS: post")
            except Exception as exc:
                _on_step_failure("post", exc)
    if stop_after == "post":
        return

    # ── STEP: synteny ──────────────────────────────────────────────────────
    synteny_dir = os.path.join(outdir, "synteny")
    safe_mkdir(synteny_dir)

    if step_in_range("synteny", start_at, stop_after) and synteny_cfg.get("enabled", False):
        if _is_blocked("synteny"):
            print("[pipeline] Skipping blocked step: synteny")
        else:
            update_step_status(status_path, "synteny", "running")
            append_pipeline_log(logs_dir, "START: synteny")
            try:
                _ensure_hit_dfs()
                synteny.run_synteny(cfg, synteny_dir, tree_dir, scan_df, search_df, hmm_keep, force,
                                    clade_assign_dir=clade_assign_dir)
                update_step_status(status_path, "synteny", "success")
                append_pipeline_log(logs_dir, "SUCCESS: synteny")
            except Exception as exc:
                _on_step_failure("synteny", exc)
    if stop_after == "synteny":
        return

    # ── STEP: codon ────────────────────────────────────────────────────────
    if step_in_range("codon", start_at, stop_after) and codon_cfg.get("enabled", False):
        if _is_blocked("codon"):
            print("[pipeline] Skipping blocked step: codon")
        else:
            update_step_status(status_path, "codon", "running")
            append_pipeline_log(logs_dir, "START: codon")
            try:
                codon.run_codon(cfg, tree_dir, clipkit_dir, aln_dir, codon_dir, hmm_keep, force)
                update_step_status(status_path, "codon", "success")
                append_pipeline_log(logs_dir, "SUCCESS: codon")
            except Exception as exc:
                _on_step_failure("codon", exc)
    if stop_after == "codon":
        return

    # ── STEP: hyphy ────────────────────────────────────────────────────────
    if step_in_range("hyphy", start_at, stop_after) and hyphy_cfg.get("enabled", False):
        if _is_blocked("hyphy"):
            print("[pipeline] Skipping blocked step: hyphy")
        else:
            update_step_status(status_path, "hyphy", "running")
            append_pipeline_log(logs_dir, "START: hyphy")
            try:
                hyphy.run_hyphy(cfg, codon_dir, tree_dir, hyphy_dir, hmm_keep, force,
                                clade_assign_dir=clade_assign_dir)
                update_step_status(status_path, "hyphy", "success")
                append_pipeline_log(logs_dir, "SUCCESS: hyphy")
            except Exception as exc:
                _on_step_failure("hyphy", exc)
    if stop_after == "hyphy":
        return

    # ── STEP: score_motifs ─────────────────────────────────────────────────
    if step_in_range("score_motifs", start_at, stop_after) and motif_cfg.get("enabled", False):
        if _is_blocked("score_motifs"):
            print("[pipeline] Skipping blocked step: score_motifs")
        else:
            update_step_status(status_path, "score_motifs", "running")
            append_pipeline_log(logs_dir, "START: score_motifs")
            try:
                from .tasks import score_motifs
                score_motifs.score_motifs(cfg, fasta_dir, summary_dir, hmm_keep, force)
                update_step_status(status_path, "score_motifs", "success")
                append_pipeline_log(logs_dir, "SUCCESS: score_motifs")
            except Exception as exc:
                _on_step_failure("score_motifs", exc)
    if stop_after == "score_motifs":
        return

    # ── STEP: discover_motifs ──────────────────────────────────────────────
    if step_in_range("discover_motifs", start_at, stop_after) and discover_cfg.get("enabled", False):
        if _is_blocked("discover_motifs"):
            print("[pipeline] Skipping blocked step: discover_motifs")
        else:
            update_step_status(status_path, "discover_motifs", "running")
            append_pipeline_log(logs_dir, "START: discover_motifs")
            try:
                from .tasks import discover_motifs
                discover_motifs.discover_motifs(cfg, fasta_dir, summary_dir, hmm_keep, force,
                                         clade_assign_dir=clade_assign_dir)
                update_step_status(status_path, "discover_motifs", "success")
                append_pipeline_log(logs_dir, "SUCCESS: discover_motifs")
            except Exception as exc:
                _on_step_failure("discover_motifs", exc)

    print("\nPipeline complete.")

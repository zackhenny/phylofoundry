import argparse
import sys
from .config import resolve_config, validate_config, STEPS


def main():
    ap = argparse.ArgumentParser(
        description="PhyloFoundry: Competitive HMM pipeline with JSON config + embeddings"
    )
    ap.add_argument("--config", default=None, help="JSON config file")
    ap.add_argument("--dump_default_config", action="store_true",
                    help="Print default config JSON and exit")

    # CLI overrides
    ap.add_argument("--faa_dir", default=None, help="Override inputs.faa_dir")
    ap.add_argument("--hmm_dir", default=None,
                    help="Override inputs.hmm_input (dir or single .hmm)")
    ap.add_argument("--diamond_query", default=None,
                    help="Override inputs.diamond_query (FASTA file/dir for DIAMOND blastp)")
    ap.add_argument("--diamond_mode", action="store_true",
                    help="Enable DIAMOND search mode (use protein FASTA queries instead of HMMs)")
    ap.add_argument("--outdir", default=None, help="Override output.outdir")
    ap.add_argument("--cpu", type=int, default=None, help="Override resources.cpu")
    ap.add_argument("--start_at", choices=STEPS, default=None,
                    help="Override workflow.start_at")
    ap.add_argument("--stop_after", choices=STEPS, default=None,
                    help="Override workflow.stop_after")
    ap.add_argument("--force", action="store_true",
                    help="Override workflow.force=True")

    # Combined tree flag
    ap.add_argument("--combined", action="store_true",
                    help="Enable combined tree from all HMMs (phylo.combined_tree)")

    # Motif scoring flags
    ap.add_argument("--motifs", default=None,
                    help="Comma-separated motif list for attention scoring "
                         "(e.g., HPEVY,HPEVF)")

    # Preflight / inspection modes
    ap.add_argument("--list-steps", action="store_true",
                    help="List all known workflow steps and exit")
    ap.add_argument("--plan", action="store_true",
                    help="Show the execution plan for the given config and exit "
                         "(no pipeline steps are run)")
    ap.add_argument("--validate-config", action="store_true",
                    help="Validate the config without running the pipeline and exit")
    ap.add_argument("--doctor", action="store_true",
                    help="Check tool availability and environment health, then exit")

    args = ap.parse_args()

    # ── Modes that do not require a full config ────────────────────────────
    if args.list_steps:
        from .preflight import print_list_steps
        print_list_steps()
        sys.exit(0)

    cfg = resolve_config(args)
    if cfg is None:  # dump_default_config was handled
        sys.exit(0)

    # Apply new CLI overrides to config
    if args.combined:
        cfg["phylo"]["combined_tree"] = True

    if args.motifs:
        motif_list = [m.strip() for m in args.motifs.split(",") if m.strip()]
        if motif_list:
            cfg.setdefault("motifs", {})["enabled"] = True
            cfg["motifs"]["motif_list"] = motif_list

    if args.diamond_query is not None:
        cfg["inputs"]["diamond_query"] = args.diamond_query

    if args.diamond_mode:
        cfg.setdefault("diamond", {})["enabled"] = True

    # ── Modes that require a resolved config but skip execution ───────────
    if args.plan:
        from .preflight import print_plan
        print_plan(cfg)
        sys.exit(0)

    if args.validate_config:
        from .preflight import validate_config_only
        ok = validate_config_only(cfg)
        sys.exit(0 if ok else 1)

    if args.doctor:
        from .preflight import run_doctor
        ok = run_doctor(cfg)
        sys.exit(0 if ok else 1)

    # ── Normal pipeline execution ─────────────────────────────────────────
    validate_config(cfg)

    # Check dependencies
    from .utils.helpers import check_dependencies

    deps = ["hmmscan", "hmmsearch"]

    if cfg.get("diamond", {}).get("enabled", False):
        # In diamond mode, hmmscan/hmmsearch are not needed; replace with diamond
        deps = ["diamond"]

    if cfg["phylo"]:
        deps.append(cfg["phylo"].get("iqtree_bin", "iqtree"))
        deps.append("mafft")
        deps.append("clipkit")

    if cfg["codon"].get("enabled", False):
        deps.append(cfg["codon"].get("pal2nal_cmd", "pal2nal.pl"))

    if cfg["hyphy"].get("enabled", False):
        deps.append(cfg["hyphy"].get("hyphy_bin", "hyphy"))

    if cfg["synteny"].get("enabled", False):
        sim_method = cfg["synteny"].get("similarity", {}).get("method", "diamond")
        if sim_method == "diamond":
            deps.append(cfg["synteny"].get("similarity", {}).get(
                "diamond_bin", "diamond"))
        elif sim_method == "mmseqs":
            deps.append(cfg["synteny"].get("similarity", {}).get(
                "mmseqs_bin", "mmseqs"))

    check_dependencies(deps)

    from .pipeline import run_pipeline
    run_pipeline(cfg)


if __name__ == "__main__":
    main()

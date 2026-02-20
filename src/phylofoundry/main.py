import argparse
import sys
from .config import resolve_config, validate_config, STEPS
from .pipeline import run_pipeline


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

    # Motif discovery flags
    ap.add_argument("--standard-clade", type=int, default=None,
                    help="HDBSCAN cluster ID for the standard reference clade")
    ap.add_argument("--novel-clade", type=int, default=None,
                    help="HDBSCAN cluster ID for the novel comparison clade")

    args = ap.parse_args()

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

    if args.standard_clade is not None:
        cfg.setdefault("discover", {})["standard_clade"] = args.standard_clade
    if args.novel_clade is not None:
        cfg.setdefault("discover", {})["novel_clade"] = args.novel_clade
    if args.standard_clade is not None and args.novel_clade is not None:
        cfg.setdefault("discover", {})["enabled"] = True

    validate_config(cfg)

    # Check dependencies
    from .utils.helpers import check_dependencies

    deps = ["hmmscan", "hmmsearch"]

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

    run_pipeline(cfg)


if __name__ == "__main__":
    main()

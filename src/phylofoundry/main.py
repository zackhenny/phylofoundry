import argparse
import os
import sys

from .config import (
    STEPS,
    config_explain,
    config_template,
    load_config_file,
    resolve_config,
    validate_config,
)
from .methods import run_ha_sites
from .pipeline import run_pipeline


def _build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description="PhyloFoundry pipeline")
    ap.add_argument("--version", action="store_true", help="Print version and exit")

    sub = ap.add_subparsers(dest="subcommand")
    cfgp = sub.add_parser("config", help="Config utilities")
    cfgsub = cfgp.add_subparsers(dest="config_cmd")
    t = cfgsub.add_parser("template", help="Print YAML template")
    t.add_argument("--mode", choices=["minimal", "full"], default="minimal")
    e = cfgsub.add_parser("explain", help="Explain configuration sections")
    e.add_argument("path", nargs="?", default=None)
    v = cfgsub.add_parser("validate", help="Validate config file")
    v.add_argument("file")

    hap = sub.add_parser("ha", help="Run standalone HA-site analysis")
    hap.add_argument("--config", required=True, help="JSON or YAML config file")
    scope = hap.add_mutually_exclusive_group(required=True)
    scope.add_argument("--hmm", help="Single HMM ID to process")
    scope.add_argument("--all", action="store_true", help="Process all HMMs")

    ap.add_argument("--config", default=None, help="JSON or YAML config file")
    ap.add_argument("--dump_default_config", action="store_true", help="Print default config JSON and exit")
    ap.add_argument("--faa_dir", default=None, help="Override inputs.faa_dir")
    ap.add_argument("--hmm_dir", default=None, help="Override inputs.hmm_input")
    ap.add_argument("--outdir", default=None, help="Override output.outdir")
    ap.add_argument("--cpu", type=int, default=None, help="Override resources.cpu")
    ap.add_argument("--start_at", choices=STEPS, default=None)
    ap.add_argument("--stop_after", choices=STEPS, default=None)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--combined", action="store_true")
    ap.add_argument("--motifs", default=None)
    ap.add_argument("--set", dest="set_values", action="append", default=[], help="Override config key(s): section.key=value")
    ap.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    return ap


def main():
    ap = _build_parser()
    args = ap.parse_args()

    if args.version:
        from . import __version__

        print(__version__)
        return

    if args.subcommand == "config":
        if args.config_cmd == "template":
            print(config_template(args.mode))
            return
        if args.config_cmd == "explain":
            print(config_explain(args.path))
            return
        if args.config_cmd == "validate":
            cfg = load_config_file(args.file)
            validate_config(cfg)
            print("Config valid")
            return
        ap.error("Missing config subcommand")

    if args.subcommand == "ha":
        cfg = validate_config(load_config_file(args.config))
        outdir = cfg["output"]["outdir"]
        fasta_dir = os.path.join(outdir, "fasta_per_hmm")
        clipkit_dir = os.path.join(outdir, "alignments_clipkit")
        run_ha_sites(
            cfg,
            fasta_dir,
            os.path.join(outdir, "embeddings"),
            os.path.join(outdir, "summary"),
            os.path.join(outdir, "qc"),
            hmm_keep=None if args.all else {args.hmm},
            alignment_dir=clipkit_dir if os.path.exists(clipkit_dir) else None,
        )
        return

    cfg = resolve_config(args)
    if cfg is None:
        sys.exit(0)

    if args.combined:
        cfg["phylo"]["combined_tree"] = True

    if args.motifs:
        motif_list = [m.strip() for m in args.motifs.split(",") if m.strip()]
        if motif_list:
            cfg.setdefault("motifs", {})["enabled"] = True
            cfg["motifs"]["motif_list"] = motif_list

    cfg = validate_config(cfg)
    cfg.setdefault("logging", {})["level"] = args.log_level

    from .utils.helpers import check_dependencies

    deps = ["hmmscan", "hmmsearch"]
    check_dependencies(deps)
    run_pipeline(cfg)


if __name__ == "__main__":
    main()

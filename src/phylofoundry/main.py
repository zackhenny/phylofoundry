import argparse
import sys
from .config import resolve_config, validate_config, STEPS
from .pipeline import run_pipeline

def main():
    ap = argparse.ArgumentParser(description="PhyloFoundry: Competitive HMM pipeline with JSON config + embeddings")
    ap.add_argument("--config", default=None, help="JSON config file")
    ap.add_argument("--dump_default_config", action="store_true", help="Print default config JSON and exit")

    # CLI overrides
    ap.add_argument("--faa_dir", default=None, help="Override inputs.faa_dir")
    ap.add_argument("--hmm_dir", default=None, help="Override inputs.hmm_input (dir or single .hmm)")
    ap.add_argument("--outdir", default=None, help="Override output.outdir")
    ap.add_argument("--cpu", type=int, default=None, help="Override resources.cpu")
    ap.add_argument("--start_at", choices=STEPS, default=None, help="Override workflow.start_at")
    ap.add_argument("--stop_after", choices=STEPS, default=None, help="Override workflow.stop_after")
    ap.add_argument("--force", action="store_true", help="Override workflow.force=True")

    # Parse common args
    args = ap.parse_args()
    
    cfg = resolve_config(args)
    if cfg is None: # dump_default_config was handled
        sys.exit(0)

    validate_config(cfg)
    
    # Check dependencies
    from .utils.helpers import check_dependencies
    check_dependencies(["hmmscan", "hmmsearch", "mafft", "clipkit", "iqtree", "pal2nal.pl", "hyphy"])

    run_pipeline(cfg)

if __name__ == "__main__":
    main()

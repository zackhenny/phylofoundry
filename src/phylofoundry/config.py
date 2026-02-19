import json
import os
import argparse
import shutil
from copy import deepcopy
from pathlib import Path
from .constants import DEFAULT_CONFIG, STEPS
from .utils.helpers import load_json_config, write_json

def deep_update(base: dict, updates: dict) -> dict:
    """Recursively update a dictionary."""
    for k, v in updates.items():
        if isinstance(v, dict) and isinstance(base.get(k), dict):
            deep_update(base[k], v)
        else:
            base[k] = v
    return base

def resolve_config(args: argparse.Namespace) -> dict:
    """Combine default config, JSON config, and CLI overrides."""
    
    if hasattr(args, 'dump_default_config') and args.dump_default_config:
        print(json.dumps(DEFAULT_CONFIG, indent=2, sort_keys=True))
        return None

    cfg = deepcopy(DEFAULT_CONFIG)
    if args.config:
        deep_update(cfg, load_json_config(args.config))

    # CLI overrides
    if args.faa_dir is not None:
        cfg["inputs"]["faa_dir"] = args.faa_dir
    if args.hmm_dir is not None:
        cfg["inputs"]["hmm_input"] = args.hmm_dir
    if args.outdir is not None:
        cfg["output"]["outdir"] = args.outdir
    if args.cpu is not None:
        cfg["resources"]["cpu"] = int(args.cpu)
    
    # Auto-detect SLURM_CPUS_PER_TASK if cpu not explicitly set in CLI (though CLI default is None, so check if it's default from config)
    # Actually, better logic: if args.cpu IS set, use it. If NOT set, check SLURM. If SLURM not set, use config default.
    if args.cpu is None:
        slurm_cpus = os.environ.get("SLURM_CPUS_PER_TASK")
        if slurm_cpus:
            try:
                cfg["resources"]["cpu"] = int(slurm_cpus)
                print(f"Auto-detected SLURM_CPUS_PER_TASK: {slurm_cpus}")
            except ValueError:
                pass

    if args.start_at is not None:
        cfg["workflow"]["start_at"] = args.start_at
    if args.stop_after is not None:
        cfg["workflow"]["stop_after"] = args.stop_after
    if args.force:
        cfg["workflow"]["force"] = True

    # Auto-detect IQ-TREE binary if default "iqtree" is not found but v2/v3 are
    # Only if user hasn't overridden it in config file (we check if it's still default)
    # Note: merge logic might have overwritten it. If it's still "iqtree", we check.
    current_bin = cfg["phylo"].get("iqtree_bin", "iqtree")
    if current_bin == "iqtree" and not shutil.which("iqtree"):
        for cand in ["iqtree2", "iqtree3"]:
            if shutil.which(cand):
                cfg["phylo"]["iqtree_bin"] = cand
                print(f"Auto-detected IQ-TREE binary: {cand}")
                break

    return cfg

def validate_config(cfg: dict):
    """Validate required configuration fields."""
    faa_arg = cfg["inputs"]["faa_dir"]
    hmm_arg = cfg["inputs"]["hmm_input"]
    outdir = cfg["output"]["outdir"]
    
    if not faa_arg or not hmm_arg or not outdir:
        raise SystemExit("Config must specify inputs.faa_dir, inputs.hmm_input, output.outdir (or pass via CLI).")

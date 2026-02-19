import os
import subprocess
import glob
import re
import json
import shutil

def safe_mkdir(p):
    os.makedirs(p, exist_ok=True)

def run_cmd(cmd, quiet=False, shell=False):
    if not quiet:
        if shell and isinstance(cmd, str):
            print("Running:", cmd)
        else:
            print("Running:", " ".join(cmd) if isinstance(cmd, list) else cmd)
            
    stdout_target = subprocess.DEVNULL if quiet else None
    stderr_target = subprocess.DEVNULL if quiet else None
    subprocess.run(cmd, check=True, stdout=stdout_target, stderr=stderr_target, shell=shell)

def normalize_genome_id(x: str) -> str:
    if x is None:
        return x
    x = str(x)
    x = re.sub(r"\.gz$", "", x, flags=re.IGNORECASE)
    while True:
        new = re.sub(r"\.(faa|fna|fa|fasta|ffn|cds)$", "", x, flags=re.IGNORECASE)
        if new == x:
            break
        x = new
    return x

def find_cds_fasta_for_genome(cds_dir, genome):
    if not cds_dir:
        return None
    genome_norm = normalize_genome_id(genome)
    patterns = [
        f"{genome}*.fna", f"{genome}*.ffn", f"{genome}*.cds*.fna",
        f"{genome_norm}*.fna", f"{genome_norm}*.ffn", f"{genome_norm}*.cds*.fna",
        f"*{genome}*.fna", f"*{genome_norm}*.fna",
    ]
    for pat in patterns:
        hits = glob.glob(os.path.join(cds_dir, pat))
        if hits:
            return sorted(hits)[0]
    return None

def check_dependencies(executables: list):
    """Check if executables are in PATH."""
    missing = []
    for exe in executables:
        if not shutil.which(exe):
            missing.append(exe)
    if missing:
        print(f"WARNING: The following executables are not found in PATH: {', '.join(missing)}")

def load_json_config(path: str) -> dict:
    """Load JSON config from file."""
    with open(path) as f:
        return json.load(f)

def write_json(obj: dict, out_fp: str):
    """Write dictionary to JSON file."""
    os.makedirs(os.path.dirname(out_fp), exist_ok=True)
    with open(out_fp, "w") as f:
        json.dump(obj, f, indent=2, sort_keys=True)

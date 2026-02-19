import os
import glob
from ..utils.helpers import run_cmd

def run_hyphy_wrapper(hyphy_bin, test_name, codon_aln_fp, tree_fp, out_json):
    cmd = [hyphy_bin, test_name, "--alignment", codon_aln_fp, "--tree", tree_fp, "--output", out_json]
    try:
        run_cmd(cmd, quiet=True, shell=False)
        return True
    except Exception as e:
        import sys
        print(f"[hyphy] test {test_name} failed: {e}", file=sys.stderr)
        return False

def run_hyphy(cfg, codon_dir, tree_dir, hyphy_dir, hmm_keep, force=False):
    print("\n[hyphy] Running HyPhy tests (generic wrapper)...")
    
    hyphy_cfg = cfg.get("hyphy", {})
    tests = [x.strip() for x in str(hyphy_cfg.get("hyphy_tests", "")).split(",") if x.strip()]

    hmm_names = sorted([os.path.basename(x).split(".")[0] for x in glob.glob(os.path.join(codon_dir, "*.codon.fasta"))])
    if hmm_keep is not None:
        hmm_names = [h for h in hmm_names if h in hmm_keep]

    for hmm in hmm_names:
        codon_aln_fp = os.path.join(codon_dir, f"{hmm}.codon.fasta")
        tree_fp = os.path.join(tree_dir, f"{hmm}.treefile")
        if not os.path.exists(tree_fp):
            continue

        for test in tests:
            out_json = os.path.join(hyphy_dir, f"{hmm}.{test}.json")
            if os.path.exists(out_json) and not force:
                continue
            _ = run_hyphy_wrapper(hyphy_cfg.get("hyphy_bin", "hyphy"), test, codon_aln_fp, tree_fp, out_json)

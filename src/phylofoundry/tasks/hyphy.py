import os
import sys
import glob
from ..utils.helpers import run_cmd

def run_hyphy_wrapper(hyphy_bin, test_name, codon_aln_fp, tree_fp, out_json, extra_args=None):
    """Run a HyPhy test.

    Note: RELAX requires the Newick tree to have branches labeled with
    ``{test}`` and ``{reference}`` markers.  Call sites should verify the
    tree is labeled before invoking this function for RELAX.
    """
    builtins = {"meme", "absrel", "relax", "fel", "busted", "slac", "fade", "fubar"}
    exe_test = test_name.lower() if test_name.lower() in builtins else test_name
    cmd = [hyphy_bin, exe_test, "--alignment", codon_aln_fp, "--tree", tree_fp, "--output", out_json]
    if extra_args:
        cmd.extend(extra_args)
    try:
        run_cmd(cmd, quiet=True, shell=False)
        return True
    except Exception as e:
        print(f"[hyphy] test {test_name} failed: {e}", file=sys.stderr)
        return False

def run_hyphy(cfg, codon_dir, tree_dir, hyphy_dir, hmm_keep, force=False):
    print("\n[hyphy] Running HyPhy tests (generic wrapper)...")
    
    hyphy_cfg = cfg.get("hyphy", {})
    tests = [x.strip() for x in str(hyphy_cfg.get("hyphy_tests", "")).split(",") if x.strip()]
    hyphy_args = hyphy_cfg.get("hyphy_args", {}) # e.g. {"MEME": ["--branches", "All"]}

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

            # RELAX requires branch labels ({test} / {reference}) in the Newick tree.
            # Without them the analysis will fail or hang waiting for interactive input.
            if test.upper() == "RELAX":
                with open(tree_fp) as _tf:
                    tree_content_lower = _tf.read().lower()
                if "{test}" not in tree_content_lower or "{reference}" not in tree_content_lower:
                    print(
                        f"[hyphy] Skipping RELAX for {hmm}: tree has no branch labels. "
                        "RELAX requires branches labeled with {test} and {reference} in "
                        "the Newick string. Please provide a labeled tree.",
                        file=sys.stderr,
                    )
                    continue
            
            # Fetch any test-specific arguments from the config
            extra_args = hyphy_args.get(test, [])
            if isinstance(extra_args, str):
                import shlex
                extra_args = shlex.split(extra_args)
                
            _ = run_hyphy_wrapper(hyphy_cfg.get("hyphy_bin", "hyphy"), test, codon_aln_fp, tree_fp, out_json, extra_args=extra_args)

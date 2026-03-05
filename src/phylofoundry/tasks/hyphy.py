import os
import sys
import glob
import copy
import subprocess
import pandas as pd


def _import_post_helpers():
    from .post import load_clades_tsv, tree_load, load_per_hmm_clades

    return load_clades_tsv, tree_load, load_per_hmm_clades


def _iter_nodes(root):
    if hasattr(root, "preorder"):
        return list(root.preorder(include_self=True))
    stack = [root]
    out = []
    while stack:
        node = stack.pop()
        out.append(node)
        for child in reversed(getattr(node, "children", []) or []):
            stack.append(child)
    return out


def _subtree_nodes(node):
    return _iter_nodes(node)


def _append_branch_label(node, label):
    if node is None:
        return
    marker = f"{{{label}}}"
    current = "" if node.name is None else str(node.name)
    if marker in current:
        return
    node.name = f"{current}{marker}" if current else marker


def _clone_tree(tree):
    if hasattr(tree, "copy"):
        try:
            return tree.copy()
        except TypeError:
            return tree.copy(deep=True)
    return copy.deepcopy(tree)


def _tree_tip_names(tree):
    return {n.name for n in tree.tips() if n.name is not None}


def _build_labeled_tree_for_clade(tree, tips_present, test_name, label_mode="crown", relax_label_reference=True):
    """Return a labeled copy of tree for a given clade and HyPhy test."""
    test_upper = str(test_name).upper()
    labeled = _clone_tree(tree)
    mrca = labeled.lca(tips_present)
    if mrca is None:
        return None

    subtree = _subtree_nodes(mrca)
    fg_nodes = [n for n in subtree if not n.is_tip()]

    # crown: label branches inside subtree only (exclude incoming branch to MRCA)
    # stem: include crown + incoming branch to MRCA (represented by MRCA node label)
    if str(label_mode).lower() != "stem":
        fg_nodes = [n for n in fg_nodes if n is not mrca]
        if not fg_nodes and not mrca.is_tip():
            # Ensure a non-empty foreground set for small crown clades.
            fg_nodes = [mrca]
    elif mrca.parent is None:
        # root has no incoming branch to label
        fg_nodes = [n for n in fg_nodes if n is not mrca]

    if test_upper in {"ABSREL", "BUSTED"}:
        for node in fg_nodes:
            _append_branch_label(node, "FG")
    elif test_upper == "RELAX":
        fg_set = set(fg_nodes)
        for node in fg_nodes:
            _append_branch_label(node, "test")
        if relax_label_reference:
            for node in _iter_nodes(labeled):
                if node.is_tip() or node in fg_set:
                    continue
                _append_branch_label(node, "reference")

    return labeled


def _serialize_tree(tree, out_fp):
    os.makedirs(os.path.dirname(out_fp), exist_ok=True)
    tree.write(out_fp, format="newick")


def _normalize_args(extra_args):
    if extra_args is None:
        return []
    if isinstance(extra_args, str):
        import shlex

        return shlex.split(extra_args)
    return list(extra_args)


def _set_or_replace_arg(args, flag, value):
    args = list(args)
    if flag in args:
        idx = args.index(flag)
        if idx + 1 < len(args):
            args[idx + 1] = value
            return args
    args.extend([flag, value])
    return args


def _ensure_clade_test_args(test_name, base_args, relax_label_reference=True):
    args = _normalize_args(base_args)
    test_upper = str(test_name).upper()

    if test_upper in {"ABSREL", "BUSTED"}:
        args = _set_or_replace_arg(args, "--branches", "FG")
    elif test_upper == "RELAX":
        args = _set_or_replace_arg(args, "--test", "test")
        if relax_label_reference:
            args = _set_or_replace_arg(args, "--reference", "reference")

    return args


def _tree_has_relax_labels(tree_fp):
    with open(tree_fp) as tf:
        tree_content_lower = tf.read().lower()
    return "{test}" in tree_content_lower and "{reference}" in tree_content_lower


def run_hyphy_wrapper(hyphy_bin, test_name, codon_aln_fp, tree_fp, out_json, extra_args=None, log_fp=None):
    """Run a HyPhy test and optionally write a combined stdout/stderr log."""
    builtins = {"meme", "absrel", "relax", "fel", "busted", "slac", "fade", "fubar"}
    exe_test = test_name.lower() if test_name.lower() in builtins else test_name
    cmd = [hyphy_bin, exe_test, "--alignment", codon_aln_fp, "--tree", tree_fp, "--output", out_json]
    if extra_args:
        cmd.extend(extra_args)

    os.makedirs(os.path.dirname(out_json), exist_ok=True)
    if log_fp:
        os.makedirs(os.path.dirname(log_fp), exist_ok=True)

    try:
        proc = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        if log_fp:
            with open(log_fp, "w") as lf:
                lf.write(proc.stdout or "")
        return True, cmd
    except subprocess.CalledProcessError as e:
        if log_fp:
            with open(log_fp, "w") as lf:
                lf.write(e.stdout or "")
                if e.stderr:
                    lf.write("\n")
                    lf.write(e.stderr)
        print(f"[hyphy] test {test_name} failed: {e}", file=sys.stderr)
        return False, cmd


def run_hyphy(cfg, codon_dir, tree_dir, hyphy_dir, hmm_keep, force=False, clade_assign_dir=None):
    print("\n[hyphy] Running HyPhy tests (generic wrapper)...")

    hyphy_cfg = cfg.get("hyphy", {})
    tests = [x.strip() for x in str(hyphy_cfg.get("hyphy_tests", "")).split(",") if x.strip()]
    hyphy_args = hyphy_cfg.get("hyphy_args", {})

    use_detected_clades = bool(hyphy_cfg.get("use_detected_clades", True))
    min_clade_size = int(hyphy_cfg.get("min_clade_size", 4))
    label_mode = str(hyphy_cfg.get("label_mode", "crown")).lower()
    relax_label_reference = bool(hyphy_cfg.get("relax_label_reference", True))

    summary_dir = os.path.dirname(hyphy_dir)
    results_dir = os.path.dirname(summary_dir)
    labeled_tree_root = os.path.join(results_dir, "trees_labeled")
    detected_clades_fp = os.path.join(summary_dir, "detected_clades.tsv")

    clades = None
    tree_load = None
    load_per_hmm = None
    global_clades = None
    if use_detected_clades:
        load_clades_tsv, tree_load, load_per_hmm = _import_post_helpers()
        if os.path.exists(detected_clades_fp):
            global_clades = load_clades_tsv(detected_clades_fp)
            print(f"[hyphy] Loaded {len(global_clades)} detected clades from {detected_clades_fp}")
        # Mark that we should attempt per-HMM loading
        clades = global_clades

    hmm_names = sorted([os.path.basename(x).replace(".codon.fasta", "") for x in glob.glob(os.path.join(codon_dir, "*.codon.fasta"))])
    if hmm_keep is not None:
        hmm_names = [h for h in hmm_names if h in hmm_keep]

    qc_rows = []

    for hmm in hmm_names:
        codon_aln_fp = os.path.join(codon_dir, f"{hmm}.codon.fasta")
        tree_fp = os.path.join(tree_dir, f"{hmm}.treefile")
        if not os.path.exists(tree_fp):
            continue

        if not clades:
            for test in tests:
                out_json = os.path.join(hyphy_dir, f"{hmm}.{test}.json")
                if os.path.exists(out_json) and not force:
                    continue

                if test.upper() == "RELAX" and not _tree_has_relax_labels(tree_fp):
                    print(
                        f"[hyphy] Skipping RELAX for {hmm}: tree has no branch labels. "
                        "RELAX requires branches labeled with {test} and {reference} in "
                        "the Newick string. Please provide a labeled tree.",
                        file=sys.stderr,
                    )
                    continue

                extra_args = _normalize_args(hyphy_args.get(test, []))
                run_hyphy_wrapper(
                    hyphy_cfg.get("hyphy_bin", "hyphy"),
                    test,
                    codon_aln_fp,
                    tree_fp,
                    out_json,
                    extra_args=extra_args,
                    log_fp=None,
                )
            continue

        tree = tree_load(tree_fp)
        tip_names = _tree_tip_names(tree)

        # Load per-HMM clades if available, otherwise fall back to global with prefix filtering
        hmm_clades = None
        if use_detected_clades and load_per_hmm and clade_assign_dir:
            # Check for combined mode
            combined_tree = cfg.get("phylo", {}).get("combined_tree", False)
            if combined_tree:
                hmm_clades = load_per_hmm(clade_assign_dir, "combined")
            else:
                hmm_clades = load_per_hmm(clade_assign_dir, hmm)

        if hmm_clades is None and clades:
            # Fall back to global clades with prefix filtering (backward compat)
            hmm_clades = {cn: tips for cn, tips in clades.items()
                          if "|" not in cn or cn.split("|", 1)[0] == hmm}

        if not hmm_clades:
            continue
        for clade_name in sorted(hmm_clades):
            clade_tips = hmm_clades[clade_name]
            tips_present = [t for t in clade_tips if t in tip_names]
            if len(tips_present) < min_clade_size:
                for test in tests:
                    qc_rows.append(
                        {
                            "hmm": hmm,
                            "test": test,
                            "clade_name": clade_name,
                            "n_tips": len(tips_present),
                            "labeled_tree_path": "",
                            "out_json_path": "",
                            "status": "skipped_min_clade_size",
                        }
                    )
                continue

            for test in tests:
                test_upper = test.upper()
                # MEME remains whole-tree unless user explicitly supplies non-All branches.
                if test_upper == "MEME":
                    extra_args = _normalize_args(hyphy_args.get(test, []))
                    if "--branches" not in extra_args:
                        extra_args = _set_or_replace_arg(extra_args, "--branches", "All")
                    out_json = os.path.join(hyphy_dir, f"{hmm}.{test}.json")
                    if os.path.exists(out_json) and not force:
                        continue
                    ok, _ = run_hyphy_wrapper(
                        hyphy_cfg.get("hyphy_bin", "hyphy"),
                        test,
                        codon_aln_fp,
                        tree_fp,
                        out_json,
                        extra_args=extra_args,
                        log_fp=None,
                    )
                    qc_rows.append(
                        {
                            "hmm": hmm,
                            "test": test,
                            "clade_name": "ALL",
                            "n_tips": len(tip_names),
                            "labeled_tree_path": tree_fp,
                            "out_json_path": out_json,
                            "status": "ok" if ok else "failed",
                        }
                    )
                    continue

                if test_upper not in {"ABSREL", "BUSTED", "RELAX"}:
                    continue

                labeled_tree = _build_labeled_tree_for_clade(
                    tree,
                    tips_present,
                    test_name=test,
                    label_mode=label_mode,
                    relax_label_reference=relax_label_reference,
                )
                if labeled_tree is None:
                    qc_rows.append(
                        {
                            "hmm": hmm,
                            "test": test,
                            "clade_name": clade_name,
                            "n_tips": len(tips_present),
                            "labeled_tree_path": "",
                            "out_json_path": "",
                            "status": "skipped_no_mrca",
                        }
                    )
                    continue

                labeled_tree_fp = os.path.join(labeled_tree_root, hmm, f"{clade_name}.{test_upper}.nwk")
                _serialize_tree(labeled_tree, labeled_tree_fp)

                out_json = os.path.join(hyphy_dir, hmm, test_upper, f"{clade_name}.json")
                out_log = os.path.join(hyphy_dir, hmm, test_upper, f"{clade_name}.log")
                if os.path.exists(out_json) and not force:
                    qc_rows.append(
                        {
                            "hmm": hmm,
                            "test": test_upper,
                            "clade_name": clade_name,
                            "n_tips": len(tips_present),
                            "labeled_tree_path": labeled_tree_fp,
                            "out_json_path": out_json,
                            "status": "skipped_exists",
                        }
                    )
                    continue

                extra_args = _ensure_clade_test_args(
                    test,
                    hyphy_args.get(test, []),
                    relax_label_reference=relax_label_reference,
                )
                ok, _ = run_hyphy_wrapper(
                    hyphy_cfg.get("hyphy_bin", "hyphy"),
                    test,
                    codon_aln_fp,
                    labeled_tree_fp,
                    out_json,
                    extra_args=extra_args,
                    log_fp=out_log,
                )
                qc_rows.append(
                    {
                        "hmm": hmm,
                        "test": test_upper,
                        "clade_name": clade_name,
                        "n_tips": len(tips_present),
                        "labeled_tree_path": labeled_tree_fp,
                        "out_json_path": out_json,
                        "status": "ok" if ok else "failed",
                    }
                )

    if qc_rows:
        pd.DataFrame(qc_rows).to_csv(os.path.join(hyphy_dir, "clade_runs.tsv"), sep="\t", index=False)

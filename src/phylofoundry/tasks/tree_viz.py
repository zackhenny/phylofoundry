"""Tree visualization step using ggtree (R).

Responsibilities
----------------
- For every IQ-TREE ``.treefile`` in *tree_dir*, call the bundled
  ``scripts/plot_iqtree_ggtree.R`` script via ``Rscript`` to produce
  annotated phylogenetic tree plots.
- Optionally overlay clade designations (from ``clade_assignments/``) and
  taxonomy annotations (from ``summary/genome_taxonomy.tsv``).
- Produce per-HMM PNG/PDF/SVG outputs in *viz_out_dir*.
- Failures are non-fatal: a warning is printed and the pipeline continues.

This module is invoked by the pipeline orchestrator after the ``phylo``,
``detect_clades``, and ``taxonomy_integrate`` steps have been run (or
skipped).  It is also automatically triggered when the ``synteny`` step
is enabled, because synteny needs the per-HMM tree panels.
"""

import glob
import os
import shutil
import subprocess
import sys

from ..utils.helpers import safe_mkdir


def _find_rscript() -> str | None:
    """Return the path to ``Rscript``, or *None* if it is not found."""
    return shutil.which("Rscript")


def _r_script_path() -> str:
    """Return the absolute path to ``plot_iqtree_ggtree.R``."""
    # The scripts/ directory lives at the repository root, two levels above
    # this file (src/phylofoundry/tasks/tree_viz.py →
    #   src/phylofoundry/tasks/ → src/phylofoundry/ → src/ → repo root).
    tasks_dir = os.path.dirname(__file__)
    repo_root = os.path.abspath(
        os.path.join(tasks_dir, "..", "..", "..", "..")
    )
    return os.path.join(repo_root, "scripts", "plot_iqtree_ggtree.R")


def run_tree_viz(
    cfg,
    tree_dir: str,
    summary_dir: str,
    clade_assign_dir: str,
    viz_out_dir: str,
    hmm_keep,
    force: bool = False,
) -> None:
    """Generate annotated ggtree plots for every IQ-TREE treefile.

    Parameters
    ----------
    cfg:
        Resolved pipeline configuration dict.
    tree_dir:
        Directory containing ``*.treefile`` outputs from the phylo step.
    summary_dir:
        Directory containing ``genome_taxonomy.tsv`` (taxonomy_integrate
        output) and ``best_hits.competitive.tsv``.
    clade_assign_dir:
        Directory containing per-HMM ``<hmm>.clades.tsv`` files from the
        detect_clades step.
    viz_out_dir:
        Root output directory; per-HMM plots are written directly here.
    hmm_keep:
        Optional set of HMM names to restrict processing to, or ``None``
        to process all HMMs.
    force:
        When *True*, regenerate plots even if they already exist.
    """
    print("\n[tree_viz] Generating ggtree-based tree visualizations...")

    rscript_bin = _find_rscript()
    if rscript_bin is None:
        print(
            "[tree_viz] WARNING: Rscript not found in PATH — "
            "skipping ggtree visualization.",
            file=sys.stderr,
        )
        return

    r_script = _r_script_path()
    if not os.path.exists(r_script):
        print(
            f"[tree_viz] WARNING: R script not found at {r_script} — "
            "skipping ggtree visualization.",
            file=sys.stderr,
        )
        return

    safe_mkdir(viz_out_dir)

    viz_cfg = cfg.get("phylo", {}).get("tree_viz", {})
    formats = viz_cfg.get("formats", ["png"])
    if isinstance(formats, list):
        formats_str = ",".join(formats)
    else:
        formats_str = str(formats)

    bootstrap_min = float(viz_cfg.get("bootstrap_min", 80))
    show_tip_labels = str(viz_cfg.get("show_tip_labels", "auto"))
    width = float(viz_cfg.get("width", 10))
    height = float(viz_cfg.get("height") or 0)
    color_palette = str(viz_cfg.get("color_palette", "Set3"))
    tax_level = str(viz_cfg.get("tax_level", "genus"))

    # Locate taxonomy annotation file (prefer genome_taxonomy.tsv)
    tax_tsv = os.path.join(summary_dir, "genome_taxonomy.tsv")
    if not os.path.exists(tax_tsv):
        tax_tsv = os.path.join(summary_dir, "best_hits.with_taxonomy.tsv")
    if not os.path.exists(tax_tsv):
        tax_tsv = None

    # Collect treefile list
    tree_files = sorted(glob.glob(os.path.join(tree_dir, "*.treefile")))
    if hmm_keep is not None:
        tree_files = [
            fp for fp in tree_files
            if os.path.basename(fp).replace(".treefile", "") in hmm_keep
        ]

    if not tree_files:
        print("[tree_viz] No treefiles found; nothing to visualize.")
        return

    n_done = 0
    n_failed = 0

    for tree_fp in tree_files:
        hmm_name = os.path.basename(tree_fp).replace(".treefile", "")

        # Check if output already exists (use first format as sentinel)
        first_fmt = formats[0] if isinstance(formats, list) else formats
        out_check = os.path.join(viz_out_dir, f"{hmm_name}.tree.{first_fmt}")
        if os.path.exists(out_check) and not force:
            continue

        # Per-HMM clade assignment file
        clades_tsv = os.path.join(clade_assign_dir, f"{hmm_name}.clades.tsv")
        if not os.path.exists(clades_tsv):
            clades_tsv = None

        cmd = [
            rscript_bin, "--vanilla", r_script,
            "--treefile", tree_fp,
            "--outdir", viz_out_dir,
            "--hmm_name", hmm_name,
            "--format", formats_str,
            "--width", str(width),
            "--height", str(height),
            "--bootstrap_min", str(bootstrap_min),
            "--show_tip_labels", show_tip_labels,
            "--color_palette", color_palette,
            "--tax_level", tax_level,
        ]
        if clades_tsv:
            cmd += ["--clades_tsv", clades_tsv]
        if tax_tsv:
            cmd += ["--taxonomy_tsv", tax_tsv]

        try:
            result = subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            if result.returncode != 0:
                print(
                    f"[tree_viz] WARNING: ggtree failed for {hmm_name}:\n"
                    f"{result.stderr[-2000:]}",
                    file=sys.stderr,
                )
                n_failed += 1
            else:
                n_done += 1
        except Exception as exc:
            print(
                f"[tree_viz] WARNING: subprocess error for {hmm_name}: {exc}",
                file=sys.stderr,
            )
            n_failed += 1

    print(
        f"[tree_viz] Done. Plots: {n_done} succeeded, {n_failed} failed. "
        f"Output: {viz_out_dir}"
    )

"""Taxonomy integration step.

Responsibilities:
- Load taxonomy from ``gtdb_dir`` or ``taxonomy_file``.
- Write ``summary/genome_taxonomy.tsv``.
- Write ``summary/best_hits.with_taxonomy.tsv``.

This module contains logic extracted from the monolithic ``post`` step so that
taxonomy integration can run independently of conservation metrics and clade
detection.
"""

import os

import pandas as pd

from .post import _load_taxonomy
from ..utils.helpers import normalize_genome_id


def run_taxonomy_integrate(cfg, summary_dir):
    """Load taxonomy and write summary TSV files.

    Parameters
    ----------
    cfg:
        Resolved pipeline configuration dict.
    summary_dir:
        Path to the ``summary/`` output directory.

    Returns
    -------
    dict
        Mapping of normalised genome ID → lineage string.  Empty dict when no
        taxonomy source is configured or when no entries are loaded.
    """
    print("\n[taxonomy_integrate] Loading and integrating taxonomy...")

    gtdb_dir = cfg["inputs"].get("gtdb_dir")
    tax_file = cfg["inputs"].get("taxonomy_file")
    globdb_tax_file = cfg["inputs"].get("globdb_taxonomy_file")

    if not (gtdb_dir or tax_file or globdb_tax_file):
        print("[taxonomy_integrate] No taxonomy source configured; skipping.")
        return {}

    tax_map = _load_taxonomy(gtdb_dir, tax_file, globdb_tax_file)
    if not tax_map:
        print("[taxonomy_integrate] No taxonomy data loaded.")
        return {}

    print(f"[taxonomy_integrate] Loaded taxonomy for {len(tax_map)} genomes.")

    # 1. Write genome → taxonomy map
    tax_rows = [{"genome": k, "taxonomy": v} for k, v in tax_map.items()]
    genome_tax_fp = os.path.join(summary_dir, "genome_taxonomy.tsv")
    pd.DataFrame(tax_rows).to_csv(genome_tax_fp, sep="\t", index=False)
    print(f"[taxonomy_integrate] Wrote {genome_tax_fp}")

    # 2. Annotate best_hits with taxonomy
    best_hits_fp = os.path.join(summary_dir, "best_hits.competitive.tsv")
    if os.path.exists(best_hits_fp):
        df = pd.read_csv(best_hits_fp, sep="\t")
        if "genome" in df.columns:
            def _lookup_taxonomy(g_filename):
                return tax_map.get(normalize_genome_id(str(g_filename)), "Unknown")

            df["taxonomy"] = df["genome"].apply(_lookup_taxonomy)
            out_fp = os.path.join(summary_dir, "best_hits.with_taxonomy.tsv")
            df.to_csv(out_fp, sep="\t", index=False)
            print(f"[taxonomy_integrate] Wrote {out_fp}")

    return tax_map

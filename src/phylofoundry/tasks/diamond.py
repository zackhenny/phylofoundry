"""DIAMOND BLAST search task for PhyloFoundry.

Provides an alternative to the HMMER-based homology search. When
``diamond.enabled`` is ``True`` in the config, the pipeline replaces the
``hmmer`` step with DIAMOND blastp searches against the combined proteome
database built by the ``prep`` step.

The output DataFrames share the same schema expected by ``extract.py`` and all
downstream steps:
    ``[genome, protein, hmm, bitscore, evalue, coverage]``

The ``hmm`` column is repurposed as the query group name (basename of the query
FASTA file without extension) so that ``extract``, ``phylo``, and every other
downstream step works without modification.
"""

import os
import subprocess
import pandas as pd
from concurrent.futures import ProcessPoolExecutor

from ..utils.helpers import safe_mkdir

# Column names for DIAMOND --outfmt 6 output
_DIAMOND_COLS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen",
]
_DIAMOND_OUTFMT_FIELDS = (
    "qseqid sseqid pident length mismatch gapopen "
    "qstart qend sstart send evalue bitscore qlen slen"
)


def make_diamond_db(combined_faa, db_path, force=False):
    """Build a DIAMOND database from *combined_faa*.

    Parameters
    ----------
    combined_faa : str
        Path to the combined proteomes FASTA (``combined_proteomes.faa``).
    db_path : str
        Destination path prefix for the DIAMOND DB (without ``.dmnd``).
    force : bool
        When ``True``, rebuild even if the ``.dmnd`` file already exists.
    """
    dmnd_file = db_path + ".dmnd"
    if os.path.exists(dmnd_file) and not force:
        print("[diamond] DB already exists, skipping makedb.")
        return
    cmd = ["diamond", "makedb", "--in", combined_faa, "--db", db_path]
    result = subprocess.run(cmd, capture_output=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"diamond makedb failed (exit {result.returncode}):\n"
            + result.stderr.decode(errors="replace")
        )


def load_query_fastas(query_input):
    """Resolve *query_input* to a ``{query_name: fasta_path}`` mapping.

    Handles three cases:

    * A single ``.faa``/``.fasta``/``.fa``/``.fas`` file — one query named
      after the file basename (without extension).
    * A directory of such files — one query per file.
    * Anything else raises ``ValueError``.

    Parameters
    ----------
    query_input : str
        Path to a FASTA file or a directory of FASTA files.

    Returns
    -------
    dict[str, str]
        Mapping from query name → absolute path to FASTA file.
    """
    _EXTENSIONS = (".faa", ".fasta", ".fa", ".fas")
    queries: dict[str, str] = {}

    if os.path.isfile(query_input):
        name = os.path.splitext(os.path.basename(query_input))[0]
        queries[name] = os.path.abspath(query_input)
    elif os.path.isdir(query_input):
        for fname in sorted(os.listdir(query_input)):
            if fname.endswith(_EXTENSIONS):
                name = os.path.splitext(fname)[0]
                queries[name] = os.path.abspath(os.path.join(query_input, fname))
        if not queries:
            raise ValueError(
                f"No FASTA files (.faa/.fasta/.fa/.fas) found in directory: {query_input}"
            )
    else:
        raise ValueError(
            f"inputs.diamond_query must be a FASTA file or directory: {query_input}"
        )
    return queries


def parse_diamond_tab(tab_file):
    """Parse a DIAMOND ``--outfmt 6`` tabular output file.

    Parameters
    ----------
    tab_file : str
        Path to the DIAMOND tabular output file.

    Returns
    -------
    pandas.DataFrame
        DataFrame with columns matching ``_DIAMOND_COLS``.  Returns an empty
        DataFrame when the file does not exist or is empty.
    """
    if not os.path.exists(tab_file) or os.path.getsize(tab_file) == 0:
        return pd.DataFrame(columns=_DIAMOND_COLS)
    try:
        return pd.read_csv(tab_file, sep="\t", header=None, names=_DIAMOND_COLS)
    except Exception:
        return pd.DataFrame(columns=_DIAMOND_COLS)


def worker_diamond_search(args_pack):
    """Run one DIAMOND blastp search for a single query FASTA.

    This function is designed for use with ``ProcessPoolExecutor.map()``.

    Parameters
    ----------
    args_pack : tuple
        ``(query_name, query_fasta, db_path, out_tsv, diamond_cfg)``

    Returns
    -------
    pandas.DataFrame or None
        DataFrame with columns
        ``[genome, protein, hmm, bitscore, evalue, coverage, pident]``
        or ``None`` when the search yields no hits or fails.
    """
    query_name, query_fasta, db_path, out_tsv, diamond_cfg = args_pack

    sensitivity = diamond_cfg.get("sensitivity", "sensitive")

    cmd = [
        "diamond", "blastp",
        "--db", db_path,
        "--query", query_fasta,
        "--out", out_tsv,
        "--outfmt", "6",
    ] + _DIAMOND_OUTFMT_FIELDS.split() + [
        "--evalue", str(diamond_cfg.get("max_evalue", 1e-5)),
        "--max-target-seqs", str(diamond_cfg.get("max_target_seqs", 500)),
        f"--{sensitivity}",
        "--threads", "1",
        "-b", str(diamond_cfg.get("block_size", 2.0)),
        "-c", str(diamond_cfg.get("index_chunks", 4)),
    ]

    try:
        result = subprocess.run(cmd, capture_output=True)
        if result.returncode != 0:
            return None

        df = parse_diamond_tab(out_tsv)
        if df.empty:
            return None

        # Compute coverage as alignment length / query length
        df["coverage"] = df["length"] / df["qlen"]

        # Decompose genome~proteinID from sseqid (format produced by prep step)
        df["genome"] = df["sseqid"].apply(
            lambda x: x.split("~")[0] if "~" in x else "Unknown"
        )
        df["protein"] = df["sseqid"].apply(
            lambda x: x.split("~", 1)[1] if "~" in x else x
        )

        # Reuse the "hmm" column name for downstream compatibility
        df["hmm"] = query_name

        return df[["genome", "protein", "hmm", "bitscore", "evalue", "coverage", "pident"]]
    except Exception:
        return None


def apply_diamond_filtering(df, diamond_cfg):
    """Filter DIAMOND hits by identity, coverage, and e-value.

    Parameters
    ----------
    df : pandas.DataFrame
        Hit DataFrame with at least columns
        ``[pident, coverage, evalue]``.
    diamond_cfg : dict
        ``diamond`` section of the PhyloFoundry config.

    Returns
    -------
    pandas.DataFrame
        Filtered copy of *df*.
    """
    if df.empty:
        return df

    df = df.copy()

    min_identity = float(diamond_cfg.get("min_identity", 30.0))
    min_coverage = float(diamond_cfg.get("min_coverage", 0.5))
    max_evalue = float(diamond_cfg.get("max_evalue", 1e-5))

    if "pident" in df.columns:
        df = df[df["pident"] >= min_identity]
    if "coverage" in df.columns:
        df = df[df["coverage"] >= min_coverage]
    if "evalue" in df.columns:
        df = df[df["evalue"] <= max_evalue]

    return df.copy()


def run_diamond(cfg, genomes, combined_faa, outdir, summary_dir, force=False):
    """Main orchestrator: build DB, run all queries in parallel, filter hits.

    Steps
    -----
    1. Build DIAMOND DB from *combined_faa*, or use a prebuilt DB when
       ``cfg["inputs"]["diamond_db"]`` is set.
    2. Load query FASTAs from ``cfg["inputs"]["diamond_query"]``.
    3. Run parallel blastp searches (one per query, ``--threads 1`` each).
    4. Load cached results for queries not re-run.
    5. Concatenate all hits, apply filtering, save to ``summary/``.
    6. Compute best hits using ``hmmer.best_hits()``.
    7. Return ``(scan_df, search_df, best_df)`` where ``scan_df`` is always
       empty (DIAMOND has no scan/search distinction).

    Parameters
    ----------
    cfg : dict
        Resolved PhyloFoundry configuration dictionary.  When
        ``cfg["inputs"]["diamond_db"]`` is set, that prebuilt ``.dmnd`` file
        is used and ``make_diamond_db`` is skipped.
    genomes : list[str]
        List of genome filenames (not used directly, kept for API parity).
    combined_faa : str
        Path to the combined proteomes FASTA (used only when no prebuilt DB is
        configured).
    outdir : str
        Pipeline output directory.
    summary_dir : str
        Path to the ``summary/`` subdirectory.
    force : bool
        When ``True``, re-run all searches even if outputs already exist.

    Returns
    -------
    tuple[pandas.DataFrame, pandas.DataFrame, pandas.DataFrame]
        ``(scan_df, search_df, best_df)`` where ``scan_df`` is empty,
        ``search_df`` contains all filtered hits, and ``best_df`` contains
        the single best hit per (genome, protein) pair.
    """
    print("\n[diamond] Running DIAMOND blastp searches (with caching)...")

    diamond_cfg = cfg.get("diamond", {})
    cpu = int(cfg["resources"]["cpu"])

    diamond_dir = os.path.join(outdir, "diamond")
    safe_mkdir(diamond_dir)

    # Use a prebuilt DIAMOND DB if provided; otherwise build one from combined_faa.
    prebuilt_db = cfg["inputs"].get("diamond_db")
    if prebuilt_db:
        from ..utils.helpers import resolve_dmnd_path
        db_path = resolve_dmnd_path(prebuilt_db)
        print(f"[diamond] Using prebuilt DIAMOND database: {prebuilt_db}")
    else:
        db_path = os.path.join(diamond_dir, "combined_proteomes")
        make_diamond_db(combined_faa, db_path, force)

    # Resolve query FASTAs
    query_input = cfg["inputs"]["diamond_query"]
    queries = load_query_fastas(query_input)

    if not queries:
        raise RuntimeError("No query FASTA files found for DIAMOND search.")

    print(f"[diamond] Found {len(queries)} query file(s): {', '.join(sorted(queries))}")

    # Identify which queries need to be (re-)run
    run_tasks = []
    cached_names = []
    for qname, qpath in queries.items():
        out_tsv = os.path.join(diamond_dir, f"{qname}.diamond.tsv")
        if os.path.exists(out_tsv) and not force:
            cached_names.append(qname)
        else:
            run_tasks.append((qname, qpath, db_path, out_tsv, diamond_cfg))

    if cached_names:
        print(f"[diamond] Using cached results for: {', '.join(sorted(cached_names))}")

    all_hits = []

    # Run searches in parallel
    if run_tasks:
        with ProcessPoolExecutor(max_workers=cpu) as exe:
            for res in exe.map(worker_diamond_search, run_tasks):
                if res is not None:
                    all_hits.append(res)

    # Load cached results
    for qname in cached_names:
        out_tsv = os.path.join(diamond_dir, f"{qname}.diamond.tsv")
        df = parse_diamond_tab(out_tsv)
        if not df.empty:
            df["coverage"] = df["length"] / df["qlen"]
            df["genome"] = df["sseqid"].apply(
                lambda x: x.split("~")[0] if "~" in x else "Unknown"
            )
            df["protein"] = df["sseqid"].apply(
                lambda x: x.split("~", 1)[1] if "~" in x else x
            )
            df["hmm"] = qname
            all_hits.append(df[["genome", "protein", "hmm", "bitscore", "evalue", "coverage", "pident"]])

    if all_hits:
        search_df = pd.concat(all_hits, ignore_index=True)
    else:
        search_df = pd.DataFrame(
            columns=["genome", "protein", "hmm", "bitscore", "evalue", "coverage", "pident"]
        )

    # Apply filtering
    hits_tsv = os.path.join(summary_dir, "diamond_hits.filtered.tsv")
    if os.path.exists(hits_tsv) and not force and not run_tasks:
        search_df = pd.read_csv(hits_tsv, sep="\t")
    else:
        search_df = apply_diamond_filtering(search_df, diamond_cfg)
        if not search_df.empty:
            search_df.to_csv(hits_tsv, sep="\t", index=False)

    print(
        f"[diamond] {len(search_df)} hits after filtering "
        f"({search_df['hmm'].nunique() if not search_df.empty else 0} queries with hits)."
    )

    # Best hits (reuse hmmer.best_hits — same schema)
    from .hmmer import best_hits
    best_hits_tsv = os.path.join(summary_dir, "best_hits.competitive.tsv")
    if os.path.exists(best_hits_tsv) and not force and not run_tasks:
        best_df = pd.read_csv(best_hits_tsv, sep="\t")
    else:
        best_df = best_hits(search_df) if not search_df.empty else pd.DataFrame()
        if not best_df.empty:
            best_df["source"] = "diamond"
            best_df.to_csv(best_hits_tsv, sep="\t", index=False)

    # scan_df is always empty in DIAMOND mode (no hmmscan equivalent)
    return pd.DataFrame(), search_df, best_df

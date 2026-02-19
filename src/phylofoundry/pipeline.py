import os
import glob
from .constants import STEPS
from .utils.helpers import safe_mkdir, write_json
from .utils.bio import read_fasta
from .tasks import prep, hmmer, extract, embed, phylo, post, codon, hyphy


def step_in_range(step, start_at, stop_after):
    """Check whether `step` falls within [start_at, stop_after]."""
    i = STEPS.index(step)
    i0 = STEPS.index(start_at) if start_at else 0
    i1 = STEPS.index(stop_after) if stop_after else len(STEPS) - 1
    return i0 <= i <= i1


def load_manifest(hmm_manifest_fp):
    if not hmm_manifest_fp:
        return None
    keep = set()
    with open(hmm_manifest_fp) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            keep.add(s)
    return keep


def _load_proteomes_lazy(genomes, faa_dir):
    """Generator-style loader: loads one genome at a time into a shared dict.
    This avoids loading ALL proteomes at once for very large datasets.
    However, for extract/embed we do need all in memory. So we return a dict."""
    proteome_seqs = {}
    for g in genomes:
        proteome_seqs[g] = read_fasta(os.path.join(faa_dir, g))
    return proteome_seqs


def run_pipeline(cfg):
    # ── Setup ──────────────────────────────────────────────────────────────
    faa_arg = cfg["inputs"]["faa_dir"]
    hmm_arg = cfg["inputs"]["hmm_input"]
    outdir = cfg["output"]["outdir"]

    start_at = cfg["workflow"]["start_at"]
    stop_after = cfg["workflow"]["stop_after"]
    force = bool(cfg["workflow"]["force"])

    phy_cfg = cfg["phylo"]
    emb_cfg = cfg["embeddings"]
    post_cfg = cfg["post"]
    codon_cfg = cfg["codon"]
    hyphy_cfg = cfg["hyphy"]

    safe_mkdir(outdir)

    # ── Output structure ───────────────────────────────────────────────────
    hmmscan_dir = os.path.join(outdir, "hmmscan_tbl")
    hmmsearch_dir = os.path.join(outdir, "hmmsearch_tbl")
    fasta_dir = os.path.join(outdir, "fasta_per_hmm")
    aln_dir = os.path.join(outdir, "alignments_hmm")
    clipkit_dir = os.path.join(outdir, "alignments_clipkit")
    tree_dir = os.path.join(outdir, "trees_iqtree")
    summary_dir = os.path.join(outdir, "summary")
    post_dir = os.path.join(summary_dir, "post_scikitbio")
    codon_dir = os.path.join(outdir, "codon_alignments")
    hyphy_dir = os.path.join(summary_dir, "hyphy")
    emb_dir = os.path.join(outdir, "embeddings")

    for d in [hmmscan_dir, hmmsearch_dir, fasta_dir, aln_dir,
              clipkit_dir, tree_dir, summary_dir, post_dir,
              codon_dir, hyphy_dir, emb_dir]:
        safe_mkdir(d)

    write_json(cfg, os.path.join(summary_dir, "resolved_config.json"))

    hmm_keep = load_manifest(cfg["workflow"]["hmm_manifest"])

    # ── Input discovery ────────────────────────────────────────────────────
    faa_abs = os.path.abspath(faa_arg)
    if os.path.isfile(faa_abs):
        if not faa_abs.endswith(".faa"):
            raise SystemExit(f"inputs.faa_dir points to a file but not .faa: {faa_abs}")
        faa_dir = os.path.dirname(faa_abs) or "."
        genomes = [os.path.basename(faa_abs)]
    else:
        if not os.path.isdir(faa_abs):
            raise SystemExit(f"inputs.faa_dir must be a directory or a single .faa file: {faa_abs}")
        faa_dir = faa_abs
        genomes = sorted(f for f in os.listdir(faa_dir) if f.endswith(".faa"))
    if not genomes:
        raise SystemExit("No .faa inputs found.")

    hmm_abs = os.path.abspath(hmm_arg)
    if os.path.isfile(hmm_abs):
        if not hmm_abs.endswith(".hmm"):
            raise SystemExit(f"inputs.hmm_input points to a file but not .hmm: {hmm_abs}")
        hmm_dir = os.path.dirname(hmm_abs) or "."
        hmm_files = [os.path.basename(hmm_abs)]
        hmm_input_mode = "file"
    else:
        if not os.path.isdir(hmm_abs):
            raise SystemExit(f"inputs.hmm_input must be a directory or .hmm file: {hmm_abs}")
        hmm_dir = hmm_abs
        hmm_files = sorted(f for f in os.listdir(hmm_dir) if f.endswith(".hmm"))
        hmm_input_mode = "dir"
    if not hmm_files:
        raise SystemExit("No .hmm files found.")

    combined_faa = os.path.join(outdir, "combined_proteomes.faa")
    combined_hmm = os.path.join(outdir, "combined.hmm")

    # ── STEP: prep ─────────────────────────────────────────────────────────
    if step_in_range("prep", start_at, stop_after):
        prep.run_prep(cfg, genomes, faa_dir, hmm_input_mode, hmm_dir,
                      hmm_files, combined_faa, combined_hmm, force)
    if stop_after == "prep":
        return

    # ── STEP: hmmer ────────────────────────────────────────────────────────
    scan_df = None
    search_df = None
    best_df = None
    if step_in_range("hmmer", start_at, stop_after):
        scan_df, search_df, best_df = hmmer.run_hmmer(
            cfg, genomes, faa_dir, hmm_files, hmm_dir, combined_hmm,
            combined_faa, outdir, summary_dir, hmmscan_dir, hmmsearch_dir,
            hmm_keep, force
        )
    if stop_after == "hmmer":
        return

    # ── Name mapping for hmmalign ──────────────────────────────────────────
    name_to_hmm_path = {}
    for hf in hmm_files:
        fp = os.path.join(hmm_dir, hf)
        try:
            with open(fp) as f:
                for line in f:
                    if line.startswith("NAME"):
                        name_to_hmm_path[line.split()[1]] = fp
                        break
        except OSError:
            continue

    # ── Helper: load hit DataFrames from disk if we skipped hmmer ──────────
    def _ensure_hit_dfs():
        nonlocal scan_df, search_df
        if scan_df is not None and search_df is not None:
            return
        import pandas as pd
        hits_scan_tsv = os.path.join(summary_dir, "hmmscan_hits.filtered.tsv")
        hits_search_tsv = os.path.join(summary_dir, "hmmsearch_hits.filtered.tsv")
        scan_df = pd.read_csv(hits_scan_tsv, sep="\t") if os.path.exists(hits_scan_tsv) else pd.DataFrame()
        search_df = pd.read_csv(hits_search_tsv, sep="\t") if os.path.exists(hits_search_tsv) else pd.DataFrame()

    # ── STEP: extract ──────────────────────────────────────────────────────
    # Only load proteomes into memory when actually needed
    hmm_to_seqs = {}
    if step_in_range("extract", start_at, stop_after):
        _ensure_hit_dfs()
        proteome_seqs = _load_proteomes_lazy(genomes, faa_dir)
        hmm_to_seqs = extract.run_extract(
            cfg, scan_df, search_df, fasta_dir, hmm_keep, proteome_seqs, force
        )
        del proteome_seqs  # free memory after extraction
    if stop_after == "extract":
        return

    # ── STEP: embed ────────────────────────────────────────────────────────
    if step_in_range("embed", start_at, stop_after) and emb_cfg.get("enabled", False):
        clades = None
        if post_cfg.get("clades_tsv", None):
            try:
                clades = post.load_clades_tsv(post_cfg["clades_tsv"])
            except Exception:
                clades = None
        embed.run_embed(cfg, hmm_to_seqs, clades, emb_dir, fasta_dir, hmm_keep, force)
    if stop_after == "embed":
        return

    # ── STEP: phylo ────────────────────────────────────────────────────────
    if step_in_range("phylo", start_at, stop_after):
        phylo.run_phylo(cfg, hmm_to_seqs, fasta_dir, aln_dir, clipkit_dir,
                        tree_dir, name_to_hmm_path, hmm_keep, force)
    if stop_after == "phylo":
        return

    # ── STEP: post ─────────────────────────────────────────────────────────
    if step_in_range("post", start_at, stop_after) and post_cfg.get("enabled", False):
        post.run_post(cfg, tree_dir, clipkit_dir, aln_dir, post_dir, hmm_keep, force)
    if stop_after == "post":
        return

    # ── STEP: codon ────────────────────────────────────────────────────────
    if step_in_range("codon", start_at, stop_after) and codon_cfg.get("enabled", False):
        codon.run_codon(cfg, tree_dir, clipkit_dir, aln_dir, codon_dir, hmm_keep, force)
    if stop_after == "codon":
        return

    # ── STEP: hyphy ────────────────────────────────────────────────────────
    if step_in_range("hyphy", start_at, stop_after) and hyphy_cfg.get("enabled", False):
        hyphy.run_hyphy(cfg, codon_dir, tree_dir, hyphy_dir, hmm_keep, force)

    print("\nPipeline complete.")

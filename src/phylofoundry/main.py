import argparse
import sys
from .config import resolve_config, validate_config, STEPS


# ── Helpers ────────────────────────────────────────────────────────────────────

def _print_default_config_yaml() -> None:
    """Print the bundled annotated default config YAML to stdout.

    Searches for config/config.yaml starting from the repository root
    (resolved relative to this file's location) then falls back to a plain
    YAML dump of DEFAULT_CONFIG when the template cannot be located.
    """
    import io
    from copy import deepcopy
    from pathlib import Path
    from ruamel.yaml import YAML
    from .constants import DEFAULT_CONFIG

    # The source tree layout is: src/phylofoundry/main.py → repo root is 3 levels up.
    # The installed-package layout may differ, so we try both.
    _this = Path(__file__).resolve()
    _candidates = [
        _this.parents[3] / "config" / "config.yaml",  # editable / source install
        _this.parents[2] / "config" / "config.yaml",  # flat installed layout
    ]
    for _template in _candidates:
        if _template.exists():
            print(_template.read_text())
            return
    # Fallback: dump DEFAULT_CONFIG as plain YAML (no comments)
    yml = YAML()
    buf = io.StringIO()
    yml.dump(deepcopy(DEFAULT_CONFIG), buf)
    print(buf.getvalue())

# ── Constants ──────────────────────────────────────────────────────────────────

# Maps CLI subcommand names to internal pipeline step names.
# Hyphenated names are preferred for multi-word steps.
_STEP_SUBCMD_MAP: dict[str, str] = {
    "prep": "prep",
    "hmmer": "hmmer",
    "extract": "extract",
    "embed": "embed",
    "maape": "maape",
    "phylo": "phylo",
    "curate": "curate",
    "taxonomy": "taxonomy_integrate",
    "conservation": "conservation_metrics",
    "detect-clades": "detect_clades",
    "aa-composition": "aa_composition",
    "post": "post",
    "tree-viz": "tree_viz",
    "synteny": "synteny",
    "codon": "codon",
    "hyphy": "hyphy",
    "score-motifs": "score_motifs",
    "discover-motifs": "discover_motifs",
}

# Maps internal step names to their (config_section, key) that enables them.
_STEP_ENABLE_MAP: dict[str, tuple[str, str]] = {
    "embed": ("embeddings", "enabled"),
    "maape": ("maape", "enabled"),
    "curate": ("curate", "enabled"),
    "taxonomy_integrate": ("taxonomy_integrate", "enabled"),
    "conservation_metrics": ("conservation_metrics", "enabled"),
    "detect_clades": ("detect_clades", "enabled"),
    "aa_composition": ("aa_composition", "enabled"),
    "post": ("post", "enabled"),
    "tree_viz": ("phylo", "tree_viz.enabled"),  # enabled via phylo.tree_viz.enabled
    "synteny": ("synteny", "enabled"),
    "codon": ("codon", "enabled"),
    "hyphy": ("hyphy", "enabled"),
    "score_motifs": ("motifs", "enabled"),
    "discover_motifs": ("discover", "enabled"),
}

# All recognised top-level subcommand tokens (used for legacy routing).
_ALL_SUBCMDS: frozenset[str] = frozenset(
    {"run", "list-steps", "list-runs", "plan", "validate", "doctor", "dump-config"}
    | set(_STEP_SUBCMD_MAP.keys())
)


# ── Argument-parser helpers ────────────────────────────────────────────────────

def _add_common_args(parser: argparse.ArgumentParser) -> None:
    """Add shared arguments used by most subparsers."""
    parser.add_argument("--config", default=None,
                        help="Path to a YAML config file (config/config.yaml) or "
                             "legacy JSON config file.  Values in the file are merged "
                             "with the built-in defaults; CLI flags take the highest "
                             "precedence.  Run 'phylofoundry dump-config' to see the "
                             "full annotated default config.")
    parser.add_argument("--faa_dir", default=None,
                        help="Directory of per-genome protein FASTA (.faa) files, or "
                             "path to a single .faa file.  Overrides inputs.faa_dir in "
                             "the config.  Required for the prep step; not needed when "
                             "running a later step that reads existing pipeline outputs.")
    parser.add_argument("--hmm_dir", default=None,
                        help="Directory of HMM profile files (.hmm), or path to a "
                             "single .hmm file.  Overrides inputs.hmm_input in the "
                             "config.  Required for prep/hmmer/extract; may be omitted "
                             "for standalone embed (--fasta_dir) or standalone phylo "
                             "with --mafft.")
    parser.add_argument("--outdir", default=None,
                        help="Output directory for all pipeline artifacts (trees, "
                             "alignments, embeddings, logs, etc.).  Overrides "
                             "output.outdir in the config.  When running a module "
                             "independently this should point to an existing pipeline "
                             "output directory so the module can find its inputs.")
    parser.add_argument("--work_dir", default=None, dest="work_dir",
                        metavar="WORK_DIR",
                        help="Optional scratch directory for large *intermediate* files "
                             "(combined_proteomes.faa, combined.hmm, DIAMOND databases). "
                             "Keeps bulky intermediates off the main output filesystem. "
                             "Overrides output.workdir in the config.  When not set, "
                             "intermediate files are written to --outdir.")
    parser.add_argument("--cpu", type=int, default=None,
                        help="Number of CPU threads to use.  Overrides resources.cpu in "
                             "the config.  Defaults to 8 (or SLURM_CPUS_PER_TASK if "
                             "that environment variable is set).")
    parser.add_argument("--force", action="store_true",
                        help="Re-run all steps unconditionally, ignoring existing "
                             "checkpoints and output files.  Overrides "
                             "workflow.force=True.")
    # ── Resume / checkpoint flags ──────────────────────────────────────────
    parser.add_argument("--resume", action="store_true",
                        help="Resume a previous run stored in --outdir: reads the "
                             "saved config snapshot, reconnects to the checkpoint "
                             "database, and skips steps whose content fingerprint "
                             "matches the prior run.")
    parser.add_argument("--resume-from", dest="resume_from", default=None,
                        metavar="RUN_ID_OR_PATH",
                        help="Explicitly specify the run ID or directory of a prior "
                             "run to resume from.  Overrides --resume (which auto-"
                             "detects from --outdir).")
    parser.add_argument("--no-resume", dest="no_resume", action="store_true",
                        help="Disable automatic resume detection.  The pipeline starts "
                             "fresh even if checkpoint data exists in --outdir.")
    # ── Input provenance flag ──────────────────────────────────────────────
    parser.add_argument("--input-run", dest="input_run", default=None,
                        metavar="PRIOR_OUTDIR",
                        help="Read input artifacts from a *different* pipeline output "
                             "directory (PRIOR_OUTDIR) instead of --outdir.  Use this "
                             "to chain separate pipeline stages: run the first stage "
                             "with its own --outdir, then pass that directory as "
                             "--input-run to the next stage.  Example: "
                             "'phylofoundry phylo --outdir ./phylo_out' followed by "
                             "'phylofoundry hyphy --input-run ./phylo_out "
                             "--outdir ./hyphy_out'.  The prior run's config snapshot "
                             "and artifact checksums are recorded in the new run's "
                             "provenance.")


def _build_parser() -> argparse.ArgumentParser:
    """Return the fully configured argument parser with all subcommands."""
    ap = argparse.ArgumentParser(
        prog="phylofoundry",
        description=(
            "PhyloFoundry: Modular bioinformatic pipeline for competitive HMM\n"
            "analysis, phylogenetics, and protein language model embeddings.\n\n"
            "Run the full pipeline:\n"
            "  phylofoundry run --faa_dir ./proteomes --hmm_dir ./markers "
            "--outdir ./results\n\n"
            "Run a single module:\n"
            "  phylofoundry embed --outdir ./results --cpu 8\n"
            "  phylofoundry phylo --outdir ./results --cpu 16\n\n"
            "Inspect and validate:\n"
            "  phylofoundry list-steps\n"
            "  phylofoundry plan  --config config/config.yaml\n"
            "  phylofoundry validate --config config/config.yaml\n"
            "  phylofoundry doctor"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    sub = ap.add_subparsers(dest="subcommand", metavar="<module>")

    # ── run ────────────────────────────────────────────────────────────────
    p_run = sub.add_parser(
        "run",
        help="Run the full pipeline (or a range of steps with --start_at/--stop_after)",
        description="Run the full PhyloFoundry pipeline from start to finish.",
    )
    _add_common_args(p_run)
    p_run.add_argument("--dump_default_config", action="store_true",
                       help="Print default config JSON and exit")
    p_run.add_argument("--start_at", choices=STEPS, default=None,
                       help="Override workflow.start_at")
    p_run.add_argument("--stop_after", choices=STEPS, default=None,
                       help="Override workflow.stop_after")
    p_run.add_argument("--diamond_query", default=None,
                       help="Override inputs.diamond_query "
                            "(FASTA file/dir for DIAMOND blastp)")
    p_run.add_argument("--diamond_mode", action="store_true",
                       help="Enable DIAMOND search mode "
                            "(use protein FASTA queries instead of HMMs)")
    p_run.add_argument("--diamond_db", default=None,
                       help="Path to a prebuilt DIAMOND .dmnd database "
                            "(overrides inputs.diamond_db; skips makedb step)")
    p_run.add_argument("--combined_faa", default=None,
                       help="Path to a prebuilt combined proteomes FASTA "
                            "(overrides inputs.combined_faa; skips prep build)")
    p_run.add_argument("--globdb_taxonomy", default=None,
                       help="Path to a GlobDB-style headerless taxonomy TSV "
                            "(col1=genome_id, col2=GTDB taxonomy; overrides "
                            "inputs.globdb_taxonomy_file)")
    p_run.add_argument("--combined", action="store_true",
                       help="Enable combined tree from all HMMs "
                            "(phylo.combined_tree)")
    p_run.add_argument("--motifs", default=None,
                       help="Comma-separated motif list for attention scoring "
                            "(e.g., HPEVY,HPEVF)")
    # Legacy inspection flags preserved on 'run' for backward compatibility
    p_run.add_argument("--list-steps", action="store_true",
                       help="List all known workflow steps and exit")
    p_run.add_argument("--plan", action="store_true",
                       help="Show the execution plan for the given config and exit "
                            "(no pipeline steps are run)")
    p_run.add_argument("--validate-config", action="store_true",
                       help="Validate the config without running the pipeline and exit")
    p_run.add_argument("--doctor", action="store_true",
                       help="Check tool availability and environment health, then exit")

    # ── prep ───────────────────────────────────────────────────────────────
    p_prep = sub.add_parser(
        "prep",
        help="Combine FAA and HMM inputs into pipeline inputs",
        description=(
            "Combines all per-genome protein FASTA (.faa) files into a single\n"
            "combined_proteomes.faa and all HMM profile files into a single\n"
            "combined.hmm database, then runs hmmpress to build the binary\n"
            "indexes (.h3f/.h3i/.h3m/.h3p) needed by hmmscan/hmmsearch.\n\n"
            "FASTA HEADER FORMAT\n"
            "-------------------\n"
            "  Each sequence in combined_proteomes.faa is renamed:\n"
            "    >{genome_filename}~{original_protein_id}\n"
            "  e.g. >GCF_000001405.faa~WP_000123456.1\n"
            "  Trailing stop-codon characters ('*') are stripped; 'X' and\n"
            "  other ambiguous residues are kept.\n\n"
            "TYPICAL USE\n"
            "-----------\n"
            "  phylofoundry prep --faa_dir ./proteomes --hmm_dir ./markers \\\n"
            "                    --outdir ./results\n\n"
            "  --faa_dir  directory containing one .faa file per genome\n"
            "  --hmm_dir  directory containing one or more .hmm profile files\n\n"
            "  Outputs written to --outdir (or --work_dir if set):\n"
            "    combined_proteomes.faa   all proteomes concatenated\n"
            "    combined.hmm             all HMM profiles concatenated\n"
            "    combined.hmm.h3{f,i,m,p} hmmpress binary indices\n\n"
            "REQUIREMENTS\n"
            "------------\n"
            "  hmmpress  (HMMER ≥3.3 — part of the phylofoundry conda env)\n\n"
            "USING A CONFIG FILE\n"
            "-------------------\n"
            "  phylofoundry prep --config config/config.yaml --outdir ./results\n\n"
            "  Relevant config keys:\n"
            "    inputs.faa_dir             directory of per-genome .faa files\n"
            "    inputs.hmm_input           directory or single .hmm profile file\n"
            "    output.outdir              output/results directory\n"
            "    output.workdir             scratch dir for large intermediates\n"
            "    prep.cleanup_combined_faa  delete combined_proteomes.faa after\n"
            "                               downstream steps finish (saves disk)\n\n"
            "SKIPPING PREP\n"
            "-------------\n"
            "  Pass --combined_faa to supply a pre-built combined FASTA and skip\n"
            "  this step entirely:\n"
            "    phylofoundry hmmer --combined_faa ./combined.faa --outdir ./results"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    _add_common_args(p_prep)
    p_prep.add_argument("--combined_faa", default=None,
                        help="Path to a pre-built combined proteomes FASTA.  When "
                             "supplied, the prep build step is skipped and this file "
                             "is used directly.  Overrides inputs.combined_faa in the "
                             "config.")
    p_hmmer = sub.add_parser(
        "hmmer",
        help="Run HMM scan/search (or DIAMOND blastp) against proteomes",
        description=(
            "Runs hmmscan (per-genome) and/or hmmsearch (per-HMM) using\n"
            "PyHMMER to identify HMM hits across all input proteomes, then\n"
            "applies competitive filtering to produce a single best-hit per\n"
            "(genome, protein) pair.  Alternatively runs DIAMOND blastp\n"
            "when --diamond_mode is given.\n\n"
            "COMPETITIVE FILTERING\n"
            "---------------------\n"
            "  For each (genome, protein) pair, the HMM with the highest\n"
            "  bitscore 'wins'.  A hit is included only when:\n"
            "    bitscore  ≥ filtering.global_min_score  (default: 25.0)\n"
            "    coverage  ≥ filtering.min_coverage       (default: 0.5)\n"
            "    e-value   ≤ filtering.max_evalue          (if use_evalue=true)\n"
            "  Override with --min_score and --min_coverage, or set\n"
            "  filtering.disable_bitscore_filter/disable_coverage_filter=true.\n\n"
            "TYPICAL USE (after prep)\n"
            "------------------------\n"
            "  phylofoundry hmmer --outdir ./results\n\n"
            "  Expects in --outdir from a previous prep run:\n"
            "    combined_proteomes.faa\n"
            "    combined.hmm\n\n"
            "  Outputs written to --outdir:\n"
            "    hmmscan_tbl/                    per-genome domain hit tables\n"
            "    hmmsearch_tbl/                  per-HMM domain hit tables\n"
            "    summary/best_hits.tsv           all filtered hits\n"
            "    summary/best_hits.competitive.tsv  one best hit per protein\n"
            "      columns: genome, protein, hmm, bitscore, evalue, coverage\n\n"
            "SKIPPING PREP (pre-built inputs)\n"
            "--------------------------------\n"
            "  phylofoundry hmmer --combined_faa ./combined.faa \\\n"
            "                     --hmm_dir ./markers --outdir ./results\n\n"
            "DIAMOND MODE (alternative to HMMER)\n"
            "------------------------------------\n"
            "  phylofoundry hmmer --diamond_mode \\\n"
            "                     --faa_dir ./proteomes \\\n"
            "                     --diamond_query ./queries.faa \\\n"
            "                     --outdir ./results\n\n"
            "  Use --diamond_db to skip the makedb step with a prebuilt .dmnd file.\n"
            "  DIAMOND sensitivity is set via diamond.sensitivity in the config\n"
            "  (fast | mid-sensitive | sensitive | more-sensitive | very-sensitive\n"
            "  | ultra-sensitive; default: sensitive).\n\n"
            "REQUIREMENTS\n"
            "------------\n"
            "  HMMER mode:   PyHMMER Python package (no external binary needed)\n"
            "  DIAMOND mode: diamond binary ≥2.0 on PATH\n\n"
            "USING A CONFIG FILE\n"
            "-------------------\n"
            "  phylofoundry hmmer --config config/config.yaml --outdir ./results\n\n"
            "  Relevant config keys:\n"
            "    inputs.faa_dir             directory of per-genome .faa files\n"
            "    inputs.hmm_input           HMM profiles (dir or single file)\n"
            "    inputs.combined_faa        prebuilt combined proteomes FASTA\n"
            "    inputs.diamond_query       query FASTA for DIAMOND mode\n"
            "    inputs.diamond_db          prebuilt DIAMOND .dmnd database\n"
            "    hmmer.run_scan             enable hmmscan (default: true)\n"
            "    hmmer.run_search           enable hmmsearch (default: true)\n"
            "    filtering.global_min_score minimum bitscore to keep a hit\n"
            "    filtering.min_coverage     minimum alignment coverage\n"
            "    filtering.max_evalue       maximum e-value (used if use_evalue)\n"
            "    filtering.use_evalue       filter by e-value instead of bitscore\n"
            "    filtering.disable_bitscore_filter  skip bitscore threshold\n"
            "    filtering.disable_coverage_filter  skip coverage threshold\n"
            "    diamond.*                  DIAMOND-specific settings\n"
            "    output.outdir, resources.cpu"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    _add_common_args(p_hmmer)
    p_hmmer.add_argument("--diamond_query", default=None,
                         help="FASTA file or directory of FASTA files to use as "
                              "DIAMOND blastp queries.  Only required in --diamond_mode. "
                              "Overrides inputs.diamond_query in the config.")
    p_hmmer.add_argument("--diamond_mode", action="store_true",
                         help="Run DIAMOND blastp instead of HMMER.  Requires "
                              "--diamond_query (protein FASTA) and --faa_dir (target "
                              "proteomes).  Use --diamond_db to skip the makedb step.")
    p_hmmer.add_argument("--diamond_db", default=None,
                         help="Path to a pre-built DIAMOND .dmnd database file. "
                              "Skips the makedb step.  Overrides inputs.diamond_db.")
    p_hmmer.add_argument("--combined_faa", default=None,
                         help="Path to a pre-built combined proteomes FASTA.  Skips "
                              "the prep build step.  Overrides inputs.combined_faa.")
    p_hmmer.add_argument("--min_score", type=float, default=None,
                         help="Minimum bitscore threshold for retaining a hit "
                              "(default: 25.0).  Hits below this value are discarded "
                              "before competitive assignment.  Overrides "
                              "filtering.global_min_score in the config.")
    p_hmmer.add_argument("--min_coverage", type=float, default=None,
                         help="Minimum alignment coverage fraction [0-1] required to "
                              "keep a hit (default: 0.5 = 50%% of the query length "
                              "aligned).  Overrides filtering.min_coverage in the "
                              "config.")

    # ── extract ────────────────────────────────────────────────────────────
    p_extract = sub.add_parser(
        "extract",
        help="Extract sequences for HMM hits into per-HMM FASTA files",
        description=(
            "Reads the competitive best-hits table produced by the hmmer step\n"
            "and extracts the matching protein sequences from combined_proteomes.faa\n"
            "into individual per-HMM FASTA files under fasta_per_hmm/.  Each file\n"
            "contains all proteins whose single best-HMM assignment is that family.\n\n"
            "SEQUENCE HEADER FORMAT\n"
            "----------------------\n"
            "  Each sequence in the per-HMM FASTA is named:\n"
            "    >{genome_filename}|{original_protein_id}\n"
            "  e.g. >GCF_000001405.faa|WP_000123456.1\n"
            "  This label format is used throughout downstream steps as the\n"
            "  'tip label' in phylogenetic trees.\n\n"
            "TYPICAL USE (after hmmer)\n"
            "-------------------------\n"
            "  phylofoundry extract --outdir ./results\n\n"
            "  Expects in --outdir from a previous hmmer run:\n"
            "    combined_proteomes.faa\n"
            "    summary/best_hits.competitive.tsv  (competitive best-hit table)\n\n"
            "  Outputs written to --outdir:\n"
            "    fasta_per_hmm/<HMM_NAME>.faa   one FASTA file per HMM family\n\n"
            "FILTERING\n"
            "---------\n"
            "  The same bitscore/coverage thresholds used in the hmmer step apply.\n"
            "  To keep ALL hits (not just the best per protein), set\n"
            "  phylo.keep_all_hits=true in the config.\n\n"
            "USING A CONFIG FILE\n"
            "-------------------\n"
            "  phylofoundry extract --config config/config.yaml --outdir ./results\n\n"
            "  Relevant config keys:\n"
            "    inputs.faa_dir              directory of per-genome .faa files\n"
            "    output.outdir               output directory\n"
            "    filtering.global_min_score  minimum bitscore (default: 25.0)\n"
            "    filtering.min_coverage      minimum alignment coverage (default: 0.5)\n"
            "    phylo.keep_all_hits         keep all hits, not just competitive best"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    _add_common_args(p_extract)
    p_extract.add_argument("--min_score", type=float, default=None,
                           help="Minimum bitscore threshold for retaining a hit "
                                "(default: 25.0).  Overrides filtering.global_min_score "
                                "in the config.")
    p_extract.add_argument("--min_coverage", type=float, default=None,
                           help="Minimum alignment coverage fraction [0–1] for "
                                "retaining a hit (default: 0.5).  Overrides "
                                "filtering.min_coverage in the config.")

    # ── embed ──────────────────────────────────────────────────────────────
    p_embed = sub.add_parser(
        "embed",
        help="Generate protein language model embeddings (ESM-2 or HuggingFace)",
        description=(
            "Generates per-sequence mean-pooled embeddings using ESM-2 (default)\n"
            "or any HuggingFace AutoModel-compatible protein language model.  Each\n"
            "sequence is represented as a fixed-length vector (e.g. 1280 dims for\n"
            "esm2_t33_650M_UR50D) formed by averaging the token representations of\n"
            "a specified transformer layer.  Optionally runs PCA dimensionality\n"
            "reduction, UMAP projection, and HDBSCAN clustering.\n\n"
            "ESM-2 MODEL OPTIONS\n"
            "-------------------\n"
            "  esm2_t6_8M_UR50D       6 layers,   8M params,  320-dim  (fast/small)\n"
            "  esm2_t12_35M_UR50D    12 layers,  35M params,  480-dim\n"
            "  esm2_t30_150M_UR50D   30 layers, 150M params,  640-dim\n"
            "  esm2_t33_650M_UR50D   33 layers, 650M params, 1280-dim  (default)\n"
            "  esm2_t36_3B_UR50D     36 layers,   3B params, 2560-dim  (high quality)\n"
            "  esm2_t48_15B_UR50D    48 layers,  15B params, 5120-dim  (very large)\n\n"
            "TYPICAL USE (after extract in a pipeline run)\n"
            "---------------------------------------------\n"
            "  phylofoundry embed --outdir ./results --cpu 8\n\n"
            "  Expects per-HMM FASTA files in --outdir/fasta_per_hmm/.\n"
            "  Outputs written to --outdir/embeddings/:\n"
            "    <HMM>.embeddings.npy     raw embedding matrix (N × D)\n"
            "    <HMM>.pca.tsv            PCA coordinates + metadata\n"
            "    <HMM>.umap.tsv           UMAP 2-D coordinates\n"
            "    <HMM>.umap.png           UMAP scatter plot coloured by cluster\n"
            "    <HMM>.knn.tsv            k-nearest-neighbour distance table\n"
            "    <HMM>.dispersion.tsv     per-HDBSCAN-cluster dispersion metrics\n\n"
            "STANDALONE USE (skipping prep/hmmer/extract)\n"
            "--------------------------------------------\n"
            "  phylofoundry embed --fasta_dir ./my_fastas --outdir ./embed_out \\\n"
            "                     --model esm2_t12_35M_UR50D --device cuda\n\n"
            "  --fasta_dir should contain one .faa file per gene family.\n"
            "  Embeddings land in --outdir/embeddings/.\n\n"
            "CHAINING FROM A PRIOR PIPELINE RUN\n"
            "-----------------------------------\n"
            "  phylofoundry embed --input-run ./prior_pipeline_out \\\n"
            "                     --outdir ./embed_out --cpu 16\n\n"
            "  Reads fasta_per_hmm/ from the prior run directory.\n\n"
            "REQUIREMENTS\n"
            "------------\n"
            "  backend=esm:          'esm' Python package (fair-esm, pip install esm)\n"
            "  backend=transformers: transformers, torch Python packages\n"
            "  UMAP/clustering:      umap-learn, hdbscan Python packages\n\n"
            "USING A CONFIG FILE\n"
            "-------------------\n"
            "  phylofoundry embed --config config/config.yaml --outdir ./results\n\n"
            "  Relevant config keys:\n"
            "    embeddings.enabled             must be true (auto-set standalone)\n"
            "    embeddings.backend             'esm' (default) or 'transformers'\n"
            "    embeddings.model               model name/path\n"
            "    embeddings.device              'cuda' or 'cpu'\n"
            "    embeddings.batch_size          sequences per inference batch\n"
            "    embeddings.repr_layer          layer index to extract (null = last)\n"
            "    embeddings.pca_components      number of PCA dimensions (default: 3)\n"
            "    embeddings.cluster_embeddings  run HDBSCAN clustering (default: true)\n"
            "    embeddings.hdbscan_min_cluster_size  min cluster size (default: 5)\n"
            "    embeddings.write_full_vectors  write full-dim TSV in addition to .npy\n"
            "    output.outdir, resources.cpu"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    _add_common_args(p_embed)
    p_embed.add_argument("--model", default=None,
                         help="ESM-2 or HuggingFace model name or path.  Examples: "
                              "esm2_t6_8M_UR50D, esm2_t12_35M_UR50D, "
                              "esm2_t33_650M_UR50D (default), "
                              "facebook/esm2_t6_8M_UR50D.  "
                              "Overrides embeddings.model in the config.")
    p_embed.add_argument("--device", default=None,
                         help="Compute device for inference: 'cuda' (GPU, recommended) "
                              "or 'cpu'.  Defaults to 'cuda' when a GPU is available. "
                              "Overrides embeddings.device in the config.")
    p_embed.add_argument("--batch_size", type=int, default=None,
                         help="Number of sequences per inference batch.  Larger values "
                              "are faster on GPU but use more memory.  Typical range: "
                              "4–32.  Overrides embeddings.batch_size in the config.")
    p_embed.add_argument("--backend", default=None,
                         choices=["esm", "transformers"],
                         help="Embedding backend.  'esm' uses the fair-esm library "
                              "(default, requires 'esm' Python package).  'transformers' "
                              "uses HuggingFace Transformers and supports any "
                              "AutoModel-compatible protein LM.  Overrides "
                              "embeddings.backend in the config.")
    p_embed.add_argument("--repr_layer", type=int, default=None,
                         help="ESM-2 transformer layer whose token representations are "
                              "mean-pooled to form the embedding vector.  Defaults to "
                              "the last layer (e.g. layer 33 for esm2_t33_650M_UR50D). "
                              "Earlier layers capture more local sequence features; "
                              "later layers capture more global structure.  Overrides "
                              "embeddings.repr_layer in the config.")
    p_embed.add_argument("--fasta_dir", default=None,
                         help="Directory containing per-HMM .faa files to embed "
                              "directly, bypassing the prep/hmmer/extract pipeline "
                              "stages.  Each file should be named <HMM_NAME>.faa.  "
                              "Use this when you already have per-family FASTAs and "
                              "do not have (or need) a full prior pipeline run. "
                              "Overrides inputs.fasta_dir in the config.")

    # ── maape ──────────────────────────────────────────────────────────────
    p_maape = sub.add_parser(
        "maape",
        help="Run MAAPE evolutionary network analysis on protein embeddings",
        description=(
            "Constructs a weighted, directed KNN similarity network from protein\n"
            "language model embeddings using the MAAPE algorithm (Modular Assembly\n"
            "Analysis of Protein Embeddings).  Reveals evolutionary relationships\n"
            "and directional flow between protein subfamilies in embedding space.\n\n"
            "ALGORITHM OVERVIEW (6 steps)\n"
            "-----------------------------\n"
            "  1. Load & normalise embeddings (PCA reduction + L2 normalisation)\n"
            "  2. Sliding-window sub-vector path generation across PCA dimensions\n"
            "  3. Co-occurrence matrix → weighted directed edge list\n"
            "  4. KNN similarity graph construction (FAISS IndexFlatIP)\n"
            "  5. MAAPE network visualisation (directional weights on KNN graph)\n"
            "  6. Aggregated condensed visualisation (cluster super-graph)\n\n"
            "TYPICAL USE (after embed in a pipeline run)\n"
            "-------------------------------------------\n"
            "  phylofoundry maape --outdir ./results --cpu 8\n\n"
            "  Expects embedding files in --outdir/embeddings/.\n"
            "  Outputs written to --outdir/maape/:\n"
            "    <HMM>.maape_network.pkl      serialised NetworkX directed graph\n"
            "    <HMM>.edge_list.txt          weighted directed edge list (TSV)\n"
            "    <HMM>.maape_network.png      per-sequence evolutionary network plot\n"
            "    <HMM>.maape_aggregated.png   condensed cluster-level plot (Step 6)\n"
            "    <HMM>.paths.pkl              intermediate assembly paths\n"
            "    maape_summary.tsv            per-HMM network statistics summary\n\n"
            "COLOUR CODING\n"
            "-------------\n"
            "  Nodes are coloured by HDBSCAN cluster if clade_assignments/ exists\n"
            "  (from detect-clades step), otherwise by cluster index.  Override\n"
            "  colours via maape.color_scheme in the config.\n\n"
            "REQUIREMENTS\n"
            "------------\n"
            "  faiss-cpu (or faiss-gpu), networkx, matplotlib Python packages\n\n"
            "USING A CONFIG FILE\n"
            "-------------------\n"
            "  phylofoundry maape --config config/config.yaml --outdir ./results\n\n"
            "  Relevant config keys:\n"
            "    maape.enabled                    must be true (auto-set standalone)\n"
            "    maape.window_sizes               window sizes for sub-vectors\n"
            "                                     (default: [5, 10, 20, 40, 80])\n"
            "    maape.knn_k                      neighbours per node in KNN graph\n"
            "                                     (default: 20)\n"
            "    maape.knn_threshold              min cosine similarity for an edge\n"
            "                                     (default: 0.5)\n"
            "    maape.pca_components             PCA dimensions before KNN\n"
            "                                     (null = auto, up to 110)\n"
            "    maape.similarity_threshold_base  base cosine threshold at window=5\n"
            "                                     (default: 0.00001)\n"
            "    maape.generate_aggregated        run Step 6 super-graph (default: true)\n"
            "    maape.per_hmm                    process each HMM separately (default: true)\n"
            "    maape.reuse_pca                  reuse PCA from embed step (default: false)\n"
            "    maape.color_scheme               colour map: null (auto) or dict\n"
            "    output.outdir, resources.cpu"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    _add_common_args(p_maape)
    p_maape.add_argument("--knn_k", type=int, default=None,
                         help="Number of nearest neighbours for KNN graph "
                              "construction (default: 20).  Overrides maape.knn_k.")
    p_maape.add_argument("--knn_threshold", type=float, default=None,
                         help="Minimum cosine similarity for KNN edges to be "
                              "retained (default: 0.5).  Overrides maape.knn_threshold.")
    p_maape.add_argument("--window_sizes", default=None,
                         help="Comma-separated list of sliding window sizes for "
                              "sub-vector path generation (e.g. '5,10,20,40,80').  "
                              "Overrides maape.window_sizes.")
    p_maape.add_argument("--no_aggregated", action="store_true",
                         help="Skip the Step 6 aggregated condensed visualisation.  "
                              "Overrides maape.generate_aggregated=false.")

    # ── phylo ──────────────────────────────────────────────────────────────
    p_phylo = sub.add_parser(
        "phylo",
        help="Run multiple sequence alignment and IQ-TREE phylogenetic inference",
        description=(
            "Aligns sequences with hmmalign (profile-guided, default) or MAFFT\n"
            "(de-novo), optionally trims poorly-aligned columns with ClipKit,\n"
            "and infers maximum-likelihood trees with IQ-TREE.  Produces one\n"
            "tree per HMM family (and optionally a combined tree).\n\n"
            "ALIGNMENT STRATEGIES\n"
            "--------------------\n"
            "  hmmalign (default)  — profile-guided alignment using the HMM\n"
            "                        profiles from inputs.hmm_input; fast and\n"
            "                        consistent across families.\n"
            "  MAFFT (--mafft)     — de-novo progressive alignment; no HMM\n"
            "                        profiles needed.  Use --mafft_mode to choose\n"
            "                        the MAFFT strategy (auto/linsi/ginsi/einsi).\n"
            "  ClipKit             — trims gappy/poorly-conserved columns from the\n"
            "                        alignment before tree inference.  Skip with\n"
            "                        --skip_clipkit.\n\n"
            "TYPICAL USE (after extract in a pipeline run)\n"
            "---------------------------------------------\n"
            "  phylofoundry phylo --outdir ./results --cpu 16\n\n"
            "  Expects per-HMM FASTA files in --outdir/fasta_per_hmm/.\n"
            "  Reads HMM profiles from inputs.hmm_input (needed for hmmalign).\n"
            "  Outputs written to --outdir:\n"
            "    alignments_hmm/          raw hmmalign/MAFFT alignments (.afa)\n"
            "    alignments_clipkit/      ClipKit-trimmed alignments (.clipkit.faa)\n"
            "    trees_iqtree/            IQ-TREE .treefile (Newick) outputs\n\n"
            "STANDALONE USE (skipping prep/hmmer/extract)\n"
            "--------------------------------------------\n"
            "  # With hmmalign (requires HMM profiles):\n"
            "  phylofoundry phylo --fasta_dir ./my_fastas \\\n"
            "                     --hmm_dir ./markers \\\n"
            "                     --outdir ./phylo_out\n\n"
            "  # With MAFFT (no HMM profiles needed):\n"
            "  phylofoundry phylo --fasta_dir ./my_fastas --mafft \\\n"
            "                     --outdir ./phylo_out\n\n"
            "  --fasta_dir should contain one .faa file per gene family.\n\n"
            "CHAINING FROM A PRIOR PIPELINE RUN\n"
            "-----------------------------------\n"
            "  phylofoundry phylo --input-run ./prior_pipeline_out \\\n"
            "                     --outdir ./phylo_out --cpu 16\n\n"
            "REQUIREMENTS\n"
            "------------\n"
            "  hmmalign  (HMMER ≥3.3, needed unless --mafft)\n"
            "  mafft     (≥7.0, needed only with --mafft)\n"
            "  clipkit   (Python package; skipped with --skip_clipkit)\n"
            "  iqtree3 / iqtree2 / iqtree  (auto-detected from PATH)\n\n"
            "USING A CONFIG FILE\n"
            "-------------------\n"
            "  phylofoundry phylo --config config/config.yaml --outdir ./results\n\n"
            "  Relevant config keys:\n"
            "    inputs.fasta_dir         per-HMM FASTA directory (standalone)\n"
            "    inputs.hmm_input         HMM profiles dir (for hmmalign)\n"
            "    phylo.mafft              use MAFFT instead of hmmalign\n"
            "    phylo.mafft_mode         MAFFT mode: auto|linsi|ginsi|einsi|\n"
            "                             fftns|fftnsi  (default: auto)\n"
            "    phylo.iqtree_bin         IQ-TREE binary (default: iqtree3)\n"
            "    phylo.iq_boot            ultrafast bootstrap replicates (default: 1000)\n"
            "    phylo.combined_tree      build one tree from all HMMs combined\n"
            "    phylo.skip_clipkit       skip ClipKit trimming step\n"
            "    phylo.no_trim_hmmalign   skip hmmalign column trimming\n"
            "    phylo.keep_all_hits      keep all hits (not just competitive best)\n"
            "    output.outdir, resources.cpu"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    _add_common_args(p_phylo)
    p_phylo.add_argument("--mafft", action="store_true",
                         help="Use MAFFT for multiple sequence alignment instead of "
                              "hmmalign.  Required when --hmm_dir is not available.  "
                              "Overrides phylo.mafft in the config.")
    p_phylo.add_argument("--mafft_mode", default=None,
                         choices=["auto", "linsi", "ginsi", "einsi", "fftns", "fftnsi"],
                         help="MAFFT alignment strategy when --mafft is set.  "
                              "'auto' (default) lets MAFFT choose automatically.  "
                              "'linsi'/'ginsi'/'einsi' are progressively more accurate "
                              "iterative strategies for locally/globally/locally-globally "
                              "similar sequences.  'fftns'/'fftnsi' are fast FFT-based "
                              "strategies suitable for large datasets.  "
                              "Overrides phylo.mafft_mode in the config.")
    p_phylo.add_argument("--skip_clipkit", action="store_true",
                         help="Skip the ClipKit alignment trimming step.  Use this if "
                              "you want to pass the raw alignment directly to IQ-TREE "
                              "or if ClipKit is not installed.  "
                              "Overrides phylo.skip_clipkit=true in the config.")
    p_phylo.add_argument("--combined", action="store_true",
                         help="Build a single combined tree from all HMM families "
                              "concatenated.  Overrides phylo.combined_tree in the "
                              "config.")
    p_phylo.add_argument("--iqtree_bin", default=None,
                         help="IQ-TREE executable name or full path (e.g. iqtree2, "
                              "iqtree3, /opt/iqtree/bin/iqtree).  Auto-detected from "
                              "PATH if not set.  Overrides phylo.iqtree_bin.")
    p_phylo.add_argument("--iq_boot", type=int, default=None,
                         help="Number of ultrafast bootstrap replicates for IQ-TREE "
                              "(e.g. 1000).  Set to 0 to disable bootstrapping.  "
                              "Overrides phylo.iq_boot in the config.")
    p_phylo.add_argument("--fasta_dir", default=None,
                         help="Directory containing per-HMM .faa files.  Use this to "
                              "run phylo independently without a prior pipeline run.  "
                              "Each file should be named <HMM_NAME>.faa.  Combine "
                              "with --mafft to skip the need for --hmm_dir.  "
                              "Overrides inputs.fasta_dir in the config.")

    # ── curate ─────────────────────────────────────────────────────────────
    p_curate = sub.add_parser(
        "curate",
        help="Curate sequences using TreeShrink outlier pruning and/or ESM filtering",
        description=(
            "Identifies and removes outlier sequences to improve downstream\n"
            "analyses.  Two complementary methods are available and can be\n"
            "combined:\n\n"
            "  TreeShrink  — detects sequences with abnormally long branch\n"
            "                lengths in the IQ-TREE phylogeny.  These often\n"
            "                represent misassembled, chimeric, or highly\n"
            "                divergent sequences that distort tree topology.\n"
            "                Requires run_treeshrink.py on PATH.\n\n"
            "  ESM filter  — clusters sequences in ESM-2 embedding space with\n"
            "                DBSCAN (cosine distance) and marks sequences assigned\n"
            "                to noise clusters (label=-1) as outliers.  Requires\n"
            "                embeddings from the embed step.\n\n"
            "Curated FASTAs and trees are written to curated/ as overlays;\n"
            "the original pipeline outputs are never modified.\n\n"
            "TYPICAL USE (after phylo)\n"
            "-------------------------\n"
            "  phylofoundry curate --outdir ./results\n\n"
            "  Expects in --outdir:\n"
            "    trees_iqtree/          IQ-TREE treefiles (for TreeShrink)\n"
            "    fasta_per_hmm/         per-HMM FASTA files\n"
            "    alignments_clipkit/    trimmed alignments (curated copies written)\n"
            "    embeddings/            (optional) ESM embedding TSV files\n\n"
            "  Outputs written to --outdir/curated/:\n"
            "    <HMM>.faa              curated protein FASTA\n"
            "    <HMM>.treefile         curated tree (outliers pruned)\n"
            "    <HMM>.clipkit.faa      curated trimmed alignment\n"
            "    <HMM>.curate.log       log of removed sequences and reasons\n\n"
            "REQUIREMENTS\n"
            "------------\n"
            "  treeshrink method: run_treeshrink.py (TreeShrink package)\n"
            "  esm filter method: scikit-learn Python package\n\n"
            "USING A CONFIG FILE\n"
            "-------------------\n"
            "  phylofoundry curate --config config/config.yaml --outdir ./results\n\n"
            "  Relevant config keys:\n"
            "    curate.enabled          must be true (auto-set standalone)\n"
            "    curate.use_treeshrink   enable TreeShrink outlier pruning (default: true)\n"
            "    curate.use_esm_filter   enable ESM-embedding DBSCAN filter (default: true)\n"
            "    output.outdir"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    _add_common_args(p_curate)
    p_curate.add_argument("--no_treeshrink", action="store_true",
                          help="Disable TreeShrink branch-length outlier pruning.  "
                               "Overrides curate.use_treeshrink=false in the config.")
    p_curate.add_argument("--no_esm_filter", action="store_true",
                          help="Disable the ESM-embedding DBSCAN outlier filter.  "
                               "Overrides curate.use_esm_filter=false in the config.")

    # ── taxonomy ───────────────────────────────────────────────────────────
    p_tax = sub.add_parser(
        "taxonomy",
        help="Integrate taxonomy from GTDB-Tk output or a custom genome→lineage TSV",
        description=(
            "Loads GTDB-style lineage information for each genome and annotates\n"
            "the best_hits table produced by the hmmer step with full taxonomic\n"
            "ranks (domain through species).  The annotated table is used by the\n"
            "detect-clades step when detect_method='taxonomy'.\n\n"
            "INPUT SOURCES (in order of precedence)\n"
            "---------------------------------------\n"
            "  --gtdb_dir           GTDB-Tk output directory containing a file\n"
            "                       named gtdbtk.bac120.summary.tsv (or ar53 /\n"
            "                       ar122 equivalent).  Taxonomy is parsed from\n"
            "                       the 'classification' column.\n"
            "  --taxonomy_file      Custom two-column TSV.  Must have a header row\n"
            "                       with columns 'genome' and 'lineage'.  The\n"
            "                       'lineage' field should be semicolon-separated\n"
            "                       GTDB-style ranks, e.g.:\n"
            "                         d__Bacteria;p__Proteobacteria;...\n"
            "  --globdb_taxonomy    Headerless two-column TSV (col1=genome_id,\n"
            "                       col2=GTDB taxonomy string).\n\n"
            "TYPICAL USE (after hmmer)\n"
            "-------------------------\n"
            "  phylofoundry taxonomy --outdir ./results \\\n"
            "                        --gtdb_dir ./gtdbtk_output\n\n"
            "  Outputs written to --outdir:\n"
            "    summary/best_hits.with_taxonomy.tsv\n"
            "      added columns: domain, phylum, class, order, family, genus, species\n\n"
            "CUSTOM TAXONOMY FILE\n"
            "--------------------\n"
            "  phylofoundry taxonomy --outdir ./results \\\n"
            "    --taxonomy_file ./my_taxonomy.tsv\n\n"
            "USING A CONFIG FILE\n"
            "-------------------\n"
            "  phylofoundry taxonomy --config config/config.yaml --outdir ./results\n\n"
            "  Relevant config keys:\n"
            "    inputs.gtdb_dir              GTDB-Tk output directory\n"
            "    inputs.taxonomy_file         custom genome→lineage TSV\n"
            "    inputs.globdb_taxonomy_file  GlobDB headerless TSV\n"
            "    taxonomy_integrate.enabled   must be true (auto-set standalone)\n"
            "    output.outdir"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    _add_common_args(p_tax)
    p_tax.add_argument("--gtdb_dir", default=None,
                       help="Path to a GTDB-Tk output directory (contains "
                            "gtdbtk.bac120.summary.tsv or similar).  Overrides "
                            "inputs.gtdb_dir in the config.")
    p_tax.add_argument("--taxonomy_file", default=None,
                       help="Path to a custom two-column TSV with columns 'genome' "
                            "and 'lineage' (semicolon-separated GTDB-style ranks).  "
                            "Overrides inputs.taxonomy_file in the config.")
    p_tax.add_argument("--globdb_taxonomy", default=None,
                       help="Path to a GlobDB-style headerless TSV (col1=genome_id, "
                            "col2=GTDB taxonomy string).  Overrides "
                            "inputs.globdb_taxonomy_file in the config.")

    # ── conservation ───────────────────────────────────────────────────────
    p_cons = sub.add_parser(
        "conservation",
        help="Compute per-site conservation scores and KL divergence",
        description=(
            "Computes per-column sequence conservation metrics from ClipKit-\n"
            "trimmed multiple sequence alignments and, optionally, KL/Jensen-\n"
            "Shannon divergence between user-specified clade amino-acid\n"
            "frequency distributions.\n\n"
            "CONSERVATION METRICS\n"
            "--------------------\n"
            "  inverse_shannon_uncertainty (default)\n"
            "      1 − H(column) / log2(20), where H is the Shannon entropy over\n"
            "      the 20 amino-acid frequencies at a column.  Ranges 0 (fully\n"
            "      variable) to 1 (fully conserved).  Ignores gap characters.\n"
            "  Other scikit-bio TabularMSA-compatible metric names can also be\n"
            "  specified via conservation_metrics.conservation_metric.\n\n"
            "KL / JENSEN-SHANNON DIVERGENCE\n"
            "-------------------------------\n"
            "  When a clade TSV is provided (--clades_tsv or via detect-clades\n"
            "  output), per-column KL divergence is computed between each pair\n"
            "  of clade groups specified in conservation_metrics.kl_pairs, e.g.:\n"
            "    kl_pairs: 'clade_A:clade_B,clade_A:background'\n\n"
            "TYPICAL USE (after phylo)\n"
            "-------------------------\n"
            "  phylofoundry conservation --outdir ./results\n\n"
            "  Expects ClipKit-trimmed alignments in --outdir/alignments_clipkit/.\n"
            "  Outputs written to --outdir/summary/post_scikitbio/:\n"
            "    <HMM>.conservation.tsv    per-site conservation scores\n"
            "    <HMM>.kl_divergence.tsv   per-site KL divergence (if clades set)\n\n"
            "USING A CONFIG FILE\n"
            "-------------------\n"
            "  phylofoundry conservation --config config/config.yaml \\\n"
            "                            --outdir ./results\n\n"
            "  Relevant config keys:\n"
            "    conservation_metrics.enabled             must be true (auto-set standalone)\n"
            "    conservation_metrics.compute_conservation run conservation scoring\n"
            "    conservation_metrics.conservation_metric  metric name (default:\n"
            "                                              inverse_shannon_uncertainty)\n"
            "    conservation_metrics.compute_kl          compute KL divergence\n"
            "    conservation_metrics.kl_pairs            clade pairs for KL, e.g.\n"
            "                                             'A:B,A:background'\n"
            "    conservation_metrics.clades_tsv          clade assignment TSV path\n"
            "    output.outdir"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    _add_common_args(p_cons)
    p_cons.add_argument("--clades_tsv", default=None,
                        help="Path to a clade assignment TSV (tab-separated, columns: "
                             "protein_id, clade).  Used to compute per-clade KL/JS "
                             "divergence.  Overrides conservation_metrics.clades_tsv.")
    p_cons.add_argument("--conservation_metric", default=None,
                        help="Conservation metric name to pass to scikit-bio "
                             "TabularMSA.conservation().  Default is "
                             "'inverse_shannon_uncertainty'.  Overrides "
                             "conservation_metrics.conservation_metric in the config.")
    p_cons.add_argument("--compute_kl", action="store_true",
                        help="Compute KL/Jensen-Shannon divergence between clade "
                             "amino-acid distributions (requires --clades_tsv or "
                             "detect-clades output).  Overrides "
                             "conservation_metrics.compute_kl=true.")
    p_cons.add_argument("--kl_pairs", default=None,
                        help="Comma-separated list of clade pairs for KL divergence "
                             "computation, e.g. 'clade_A:clade_B,clade_A:background'. "
                             "Overrides conservation_metrics.kl_pairs in the config.")

    # ── detect-clades ──────────────────────────────────────────────────────
    p_dc = sub.add_parser(
        "detect-clades",
        help="Detect clades via taxonomy, TreeCluster, or tree+embedding",
        description=(
            "Assigns protein sequences to named clades using one of three\n"
            "automated strategies, or from a user-supplied TSV.  Clade\n"
            "assignments are used downstream by conservation, hyphy, maape,\n"
            "and discover-motifs steps.\n\n"
            "DETECTION METHODS\n"
            "-----------------\n"
            "  taxonomy\n"
            "      Groups sequences by GTDB taxonomic rank (default: genus).\n"
            "      Requires prior taxonomy_integrate output.\n"
            "      Config keys: detect_clades.taxonomy_clade_level\n\n"
            "  treecluster\n"
            "      Clusters tree tips by branch-length threshold using\n"
            "      TreeCluster.  Sequences within the threshold are grouped.\n"
            "      Config keys: detect_clades.treecluster_threshold,\n"
            "                   detect_clades.treecluster_method\n\n"
            "  tree_embed\n"
            "      Joint clustering using IQ-TREE topology (monophyly constraint)\n"
            "      + ESM-embedding space (Euclidean or cosine distance).  Finds\n"
            "      clades that are both phylogenetically coherent and\n"
            "      sequence-space coherent.\n"
            "      Config keys: detect_clades.embedtree_*\n\n"
            "  (pre-built TSV)\n"
            "      Set detect_clades.clades_tsv to a two-column TSV with columns\n"
            "      'tip' and 'clade_name' to skip automatic detection.\n\n"
            "TYPICAL USE (after phylo and taxonomy)\n"
            "--------------------------------------\n"
            "  phylofoundry detect-clades --outdir ./results\n\n"
            "  Outputs written to --outdir:\n"
            "    clade_assignments/<HMM>.tsv   per-HMM clade assignment tables\n"
            "                                  columns: tip, clade_name\n"
            "    summary/detected_clades.tsv   combined clade table (all HMMs)\n"
            "    summary/node_scores.embedtree.tsv  (tree_embed method only)\n\n"
            "USING A CONFIG FILE\n"
            "-------------------\n"
            "  phylofoundry detect-clades --config config/config.yaml \\\n"
            "                             --outdir ./results\n\n"
            "  Relevant config keys:\n"
            "    detect_clades.enabled               must be true (auto-set standalone)\n"
            "    detect_clades.detect_method         taxonomy|treecluster|tree_embed\n"
            "    detect_clades.clades_tsv            path to pre-built clade TSV\n"
            "    detect_clades.taxonomy_clade_level  GTDB rank for 'taxonomy' method\n"
            "    detect_clades.treecluster_threshold branch-length threshold (default: 0.045)\n"
            "    detect_clades.treecluster_method    TreeCluster linkage method\n"
            "    detect_clades.embedtree_*           tree_embed method parameters\n"
            "    output.outdir"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    _add_common_args(p_dc)
    p_dc.add_argument("--detect_method", default=None,
                      choices=["taxonomy", "treecluster", "tree_embed"],
                      help="Clade detection method.  'taxonomy' requires prior "
                           "taxonomy_integrate output.  'treecluster' requires trees "
                           "in trees_iqtree/.  'tree_embed' requires both trees and "
                           "embeddings.  Overrides detect_clades.detect_method.")

    # ── aa-composition ─────────────────────────────────────────────────────
    p_aac = sub.add_parser(
        "aa-composition",
        help="Compute per-gene amino acid composition and comparative statistics",
        description=(
            "Computes per-gene amino acid fractional composition (all 20 standard\n"
            "amino acids) and a suite of biochemical metrics derived from the\n"
            "genomeSPOT framework (Barnum et al.) for all HMM-hit protein\n"
            "sequences.  Optionally runs clade-level statistical comparisons and\n"
            "generates visualisation plots.\n\n"
            "COMPUTED METRICS (per gene)\n"
            "---------------------------\n"
            "  aa_A_comp … aa_Y_comp   fractional composition of each of the 20\n"
            "                          standard amino acids (0–1)\n"
            "  len                     sequence length (residues)\n"
            "  Zc                      mean carbon oxidation state of side chains\n"
            "  gravy                   GRAVY score (Grand Average of Hydropathy);\n"
            "                          positive = hydrophobic, negative = hydrophilic\n"
            "  nH2O                    estimated water activity per residue\n"
            "  pi                      isoelectric point (requires BioPython)\n"
            "  S_content               fractional sulphur-containing residue content\n"
            "                          (Cys + Met)\n"
            "  N_content               fractional nitrogen-rich residue content\n"
            "  thermo_freq             thermostability-associated residue frequency\n"
            "                          (IVYWREL)\n\n"
            "STATISTICAL COMPARISONS\n"
            "-----------------------\n"
            "  When clade assignments are available (from detect-clades or\n"
            "  clade_assignments/), runs Kruskal-Wallis tests across clade groups\n"
            "  for each metric, followed by Benjamini-Hochberg FDR correction.\n"
            "  Results written to aa_comp_stats.tsv.\n\n"
            "TYPICAL USE (after hmmer/extract)\n"
            "---------------------------------\n"
            "  phylofoundry aa-composition --outdir ./results\n\n"
            "  Outputs written to --outdir/summary/aa_composition/:\n"
            "    aa_comp_per_gene.tsv   gene × metric matrix\n"
            "    aa_comp_stats.tsv      Kruskal-Wallis + BH-FDR table (if ≥2 clades)\n"
            "    plots/                 boxplots and heatmap PNGs\n\n"
            "USING A CONFIG FILE\n"
            "-------------------\n"
            "  phylofoundry aa-composition --config config/config.yaml \\\n"
            "                              --outdir ./results\n\n"
            "  Relevant config keys:\n"
            "    aa_composition.enabled             must be true (auto-set standalone)\n"
            "    aa_composition.remove_initial_met  strip leading Met (default: true)\n"
            "    aa_composition.compute_pi          compute isoelectric point (default: true)\n"
            "    aa_composition.generate_plots      generate PNG plots (default: true)\n"
            "    aa_composition.top_n_aas_heatmap   AAs shown in heatmap (default: 10)\n"
            "    output.outdir"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    _add_common_args(p_aac)
    p_aac.add_argument("--no_plots", action="store_true",
                       help="Suppress plot generation "
                            "(overrides aa_composition.generate_plots=false)")
    p_aac.add_argument("--no_pi", action="store_true",
                       help="Skip isoelectric-point computation "
                            "(overrides aa_composition.compute_pi=false)")

    # ── post ───────────────────────────────────────────────────────────────
    p_post = sub.add_parser(
        "post",
        help="Run legacy post-processing (conservation + clade detection combined)",
        description=(
            "Backward-compatibility shim that runs conservation metrics and\n"
            "clade detection in a single combined step.  This step is equivalent\n"
            "to running 'conservation' followed by 'detect-clades' with the same\n"
            "configuration.\n\n"
            "NOTE: For new workflows, prefer the dedicated 'conservation' and\n"
            "'detect-clades' subcommands which offer finer individual control.\n\n"
            "TYPICAL USE (after phylo)\n"
            "-------------------------\n"
            "  phylofoundry post --outdir ./results\n\n"
            "  Expects in --outdir: alignments_clipkit/, trees_iqtree/,\n"
            "  and (optionally) embeddings/ and summary/best_hits.with_taxonomy.tsv.\n\n"
            "  Outputs written to --outdir:\n"
            "    summary/post_scikitbio/    conservation metric TSVs\n"
            "    clade_assignments/         per-HMM clade assignment TSVs\n"
            "    summary/detected_clades.tsv  combined clade table\n\n"
            "USING A CONFIG FILE\n"
            "-------------------\n"
            "  phylofoundry post --config config/config.yaml --outdir ./results\n\n"
            "  Relevant config keys:\n"
            "    post.enabled                  must be true (auto-set standalone)\n"
            "    post.compute_conservation     run conservation scoring\n"
            "    post.conservation_metric      metric name (inverse_shannon_uncertainty)\n"
            "    post.compute_kl               compute KL divergence\n"
            "    post.kl_pairs                 clade pairs, e.g. 'A:B,A:background'\n"
            "    post.clades_tsv               pre-built clade TSV path\n"
            "    post.detect_clades_method     taxonomy|treecluster|tree_embed\n"
            "    post.taxonomy_clade_level     GTDB rank for taxonomy method\n"
            "    post.treecluster_threshold    branch-length threshold (default: 0.045)\n"
            "    post.embedtree_*              tree_embed method parameters\n"
            "    output.outdir"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    _add_common_args(p_post)

    # ── tree-viz ────────────────────────────────────────────────────────────
    p_tree_viz = sub.add_parser(
        "tree-viz",
        help="Generate annotated ggtree plots from IQ-TREE treefiles",
        description=(
            "Produces annotated phylogenetic tree plots (PNG/PDF/SVG) from\n"
            "IQ-TREE treefiles using ggtree (R).  Optionally overlays clade\n"
            "designations from the detect-clades step, bootstrap support values,\n"
            "and taxonomy annotations from the taxonomy step.\n\n"
            "R SCRIPTS\n"
            "---------\n"
            "  scripts/plot_iqtree_ggtree.R    standard per-HMM tree plots\n"
            "  scripts/plot_synteny_tree.R      combined synteny+tree panel\n"
            "  Both scripts are called automatically; no manual R invocation needed.\n\n"
            "REQUIRES\n"
            "--------\n"
            "  R ≥4.3 with the following packages installed:\n"
            "    ggtree, treeio, ggtreeExtra, ggplot2, dplyr, RColorBrewer\n"
            "  All packages are included in the phylofoundry conda environment.\n\n"
            "TYPICAL USE (after phylo)\n"
            "-------------------------\n"
            "  phylofoundry tree-viz --outdir ./results\n\n"
            "  Expects in --outdir: trees_iqtree/ (IQ-TREE treefiles).\n"
            "  Optional inputs: clade_assignments/, summary/best_hits.with_taxonomy.tsv\n\n"
            "  Outputs written to --outdir/tree_viz/:\n"
            "    <HMM>.tree.png / .pdf     annotated tree in requested format(s)\n\n"
            "USING A CONFIG FILE\n"
            "-------------------\n"
            "  phylofoundry tree-viz --config config/config.yaml --outdir ./results\n\n"
            "  Relevant config keys:\n"
            "    phylo.tree_viz.enabled          must be true (auto-set standalone)\n"
            "    phylo.tree_viz.formats          output formats: ['png','pdf'] (default)\n"
            "    phylo.tree_viz.bootstrap_min    min bootstrap value to display (default: 80)\n"
            "    phylo.tree_viz.show_tip_labels  'auto'|'true'|'false'\n"
            "                                    (auto hides labels when >100 tips)\n"
            "    phylo.tree_viz.width            plot width in inches (default: 10)\n"
            "    phylo.tree_viz.height           plot height (null = auto-scale)\n"
            "    phylo.tree_viz.color_palette    RColorBrewer palette (default: Set3)\n"
            "    phylo.tree_viz.annotate_taxonomy overlay taxonomy labels (default: true)\n"
            "    phylo.tree_viz.tax_level        taxonomy rank to display (default: genus)\n"
            "    output.outdir"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    _add_common_args(p_tree_viz)

    # ── synteny ────────────────────────────────────────────────────────────
    p_syn = sub.add_parser(
        "synteny",
        help="Visualise synteny context around HMM hits",
        description=(
            "Extracts gene neighbourhoods around HMM hits from GenBank (.gbk/.gbff)\n"
            "or GFF3 annotation files and generates synteny track plots showing the\n"
            "genomic context ±window_genes around each hit.  Also generates a\n"
            "combined synteny+tree panel when IQ-TREE treefiles are available.\n\n"
            "ANNOTATION INPUT FORMATS\n"
            "------------------------\n"
            "  --gbk_dir   directory of GenBank flat files, one per genome\n"
            "              (.gbk or .gbff).  Protein IDs are extracted from\n"
            "              CDS qualifier fields specified in\n"
            "              synteny.protein_id_field (default: ID, protein_id,\n"
            "              locus_tag).\n"
            "  --gff_dir   directory of GFF3 files (alternative to GenBank).\n"
            "              Requires genome FASTA files in synteny.genome_fasta_dir.\n\n"
            "TYPICAL USE (after hmmer)\n"
            "-------------------------\n"
            "  phylofoundry synteny --outdir ./results \\\n"
            "                       --gbk_dir ./genbank_files\n\n"
            "  Expects in --outdir: summary/best_hits.competitive.tsv.\n"
            "  Optional: trees_iqtree/ for the combined synteny+tree panel.\n\n"
            "  Outputs written to --outdir:\n"
            "    synteny/                synteny track plots per HMM family\n"
            "    synteny/<HMM>.html      interactive clinker HTML report\n"
            "    synteny/<HMM>.png       static pygenomeviz track plot\n"
            "    synteny/tree_panel/     combined synteny+tree panel plots\n\n"
            "REQUIREMENTS\n"
            "------------\n"
            "  clinker (pip install clinker-clustermap) — for HTML reports\n"
            "  pygenomeviz (pip install pygenomeviz) — for static track plots\n"
            "  R + ggtree (for tree panel; same as tree-viz step)\n\n"
            "USING A CONFIG FILE\n"
            "-------------------\n"
            "  phylofoundry synteny --config config/config.yaml --outdir ./results\n\n"
            "  Relevant config keys:\n"
            "    synteny.enabled              must be true (auto-set standalone)\n"
            "    synteny.gbk_dir              GenBank directory\n"
            "    synteny.gff_dir              GFF3 directory (alternative to gbk)\n"
            "    synteny.genome_fasta_dir     genome FASTA dir (required with gff)\n"
            "    synteny.window_genes         genes either side of hit (default: 10)\n"
            "    synteny.max_hits_per_hmm     max hits to plot per HMM (default: 50)\n"
            "    synteny.dedup_by_genome      one hit per genome (default: true)\n"
            "    synteny.protein_id_field     CDS qualifier(s) for protein ID lookup\n"
            "    synteny.gene_label_field     CDS qualifier(s) for gene labels\n"
            "    synteny.include_tree         include IQ-TREE ordered tracks (true)\n"
            "    synteny.generate_tree_panel  generate combined tree+synteny plot\n"
            "    output.outdir"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    _add_common_args(p_syn)
    p_syn.add_argument("--gbk_dir", default=None,
                       help="Directory of GenBank (.gbk/.gbff) annotation files, one "
                            "per genome.  Used to extract gene neighbourhood windows "
                            "around HMM hits.  Overrides synteny.gbk_dir in the config.")
    p_syn.add_argument("--gff_dir", default=None,
                       help="Directory of GFF3 annotation files.  Alternative to "
                            "--gbk_dir; requires --faa_dir (or genome FASTAs) to "
                            "resolve sequences.  Overrides synteny.gff_dir in the "
                            "config.")
    p_syn.add_argument("--window_genes", type=int, default=None,
                       help="Number of genes to include on each side of the focal "
                            "HMM hit in the synteny window (default: 10).  Larger "
                            "values provide more genomic context but produce wider "
                            "plots.  Overrides synteny.window_genes in the config.")

    # ── codon ──────────────────────────────────────────────────────────────
    p_codon = sub.add_parser(
        "codon",
        help="Build codon-aware alignments with pal2nal",
        description=(
            "Threads nucleotide CDS sequences onto the protein multiple sequence\n"
            "alignments from the phylo step using pal2nal to produce codon-aware\n"
            "nucleotide alignments.  These alignments are required by the hyphy\n"
            "step for molecular-evolution tests (dN/dS, RELAX, aBSREL, MEME).\n\n"
            "CDS ID MATCHING\n"
            "---------------\n"
            "  The CDS nucleotide FASTA headers must be reconcilable with the\n"
            "  protein IDs used in the alignments.  Three modes are supported via\n"
            "  codon.cds_id_mode:\n"
            "    same              IDs match exactly\n"
            "    strip_pipe        strip everything after the first '|'\n"
            "    after_last_pipe   use only the part after the last '|' (default)\n\n"
            "TYPICAL USE (after phylo)\n"
            "-------------------------\n"
            "  phylofoundry codon --outdir ./results\n\n"
            "  Expects in --outdir: alignments_clipkit/ (protein MSAs).\n"
            "  Requires inputs.cds_dir: directory of per-genome CDS nucleotide\n"
            "  FASTA files (one file per genome, matching .faa filenames).\n\n"
            "  Outputs written to --outdir/codon_alignments/:\n"
            "    <HMM>.codon.fna    codon-aware nucleotide alignment (PHYLIP or FASTA)\n\n"
            "CHAINING FROM A PRIOR PIPELINE RUN\n"
            "-----------------------------------\n"
            "  phylofoundry codon --input-run ./phylo_out --outdir ./codon_out\n\n"
            "REQUIREMENTS\n"
            "------------\n"
            "  pal2nal.pl  (Perl script; add to PATH or set --pal2nal_cmd)\n\n"
            "USING A CONFIG FILE\n"
            "-------------------\n"
            "  phylofoundry codon --config config/config.yaml --outdir ./results\n\n"
            "  Relevant config keys:\n"
            "    codon.enabled                must be true (auto-set standalone)\n"
            "    codon.build_codon_alignments enable codon alignment step (default: false)\n"
            "    codon.cds_id_mode            ID matching mode: same|strip_pipe|\n"
            "                                 after_last_pipe  (default: after_last_pipe)\n"
            "    codon.pal2nal_cmd            pal2nal.pl path (default: pal2nal.pl)\n"
            "    inputs.cds_dir               directory of per-genome CDS FASTA files\n"
            "    output.outdir"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    _add_common_args(p_codon)
    p_codon.add_argument("--pal2nal_cmd", default=None,
                         help="Name or full path to the pal2nal.pl script "
                              "(e.g. pal2nal.pl, /opt/pal2nal/pal2nal.pl).  "
                              "Overrides codon.pal2nal_cmd in the config.")
    p_codon.add_argument("--cds_id_mode", default=None,
                         choices=["same", "strip_pipe", "after_last_pipe"],
                         help="Strategy for matching CDS nucleotide FASTA headers to "
                              "protein IDs in the alignments.  'same': IDs must match "
                              "exactly.  'strip_pipe': remove everything after the "
                              "first '|'.  'after_last_pipe': use only the text after "
                              "the last '|' (default).  Overrides codon.cds_id_mode.")

    # ── hyphy ──────────────────────────────────────────────────────────────
    p_hyphy = sub.add_parser(
        "hyphy",
        help="Run HyPhy selection tests (RELAX, aBSREL, MEME)",
        description=(
            "Runs HyPhy molecular-evolution tests on codon alignments + labelled\n"
            "phylogenetic trees to detect sites or lineages under positive or\n"
            "relaxed selection.  Clade assignments from detect-clades are used\n"
            "to label test vs. reference branches.\n\n"
            "SUPPORTED HYPHY TESTS\n"
            "---------------------\n"
            "  RELAX   — Detects relaxation or intensification of selection on a\n"
            "            labelled set of branches (test vs. reference).  Outputs\n"
            "            a K parameter: K<1 = relaxation, K>1 = intensification.\n\n"
            "  aBSREL  — Adaptive Branch-Site Random Effects Likelihood.  Tests\n"
            "            each branch individually for episodic positive selection\n"
            "            (dN/dS > 1).  More sensitive than PAML branch-site model.\n\n"
            "  MEME    — Mixed Effects Model of Evolution.  Site-level test for\n"
            "            episodic diversifying selection at individual codons\n"
            "            across a subset of lineages.  More powerful than FEL/FUBAR\n"
            "            for detecting selection on a subset of branches.\n\n"
            "TYPICAL USE (after codon and detect-clades)\n"
            "-------------------------------------------\n"
            "  phylofoundry hyphy --outdir ./results\n\n"
            "  Expects in --outdir:\n"
            "    codon_alignments/     codon-aware nucleotide alignments\n"
            "    trees_iqtree/         IQ-TREE treefiles\n"
            "    clade_assignments/    (optional) per-HMM clade TSVs for labelling\n\n"
            "  Outputs written to --outdir/summary/hyphy/:\n"
            "    <HMM>.<TEST>.json     raw HyPhy JSON output\n"
            "    <HMM>.<TEST>.tsv      parsed per-site/per-branch results\n"
            "    hyphy_summary.tsv     cross-HMM summary table\n\n"
            "CHAINING FROM A PRIOR PIPELINE RUN\n"
            "-----------------------------------\n"
            "  phylofoundry hyphy --input-run ./prior_run --outdir ./hyphy_out\n\n"
            "REQUIREMENTS\n"
            "------------\n"
            "  hyphy ≥2.5  (binary on PATH, or specify with --hyphy_bin)\n\n"
            "USING A CONFIG FILE\n"
            "-------------------\n"
            "  phylofoundry hyphy --config config/config.yaml --outdir ./results\n\n"
            "  Relevant config keys:\n"
            "    hyphy.enabled              must be true (auto-set standalone)\n"
            "    hyphy.run_hyphy            enable HyPhy execution (default: false)\n"
            "    hyphy.hyphy_bin            HyPhy executable (default: hyphy)\n"
            "    hyphy.hyphy_tests          comma-separated tests: RELAX,aBSREL,MEME\n"
            "    hyphy.use_detected_clades  use clade_assignments/ to label branches\n"
            "    hyphy.min_clade_size       min sequences per clade to run test (default: 4)\n"
            "    hyphy.label_mode           branch labelling: 'crown' or 'stem'\n"
            "    hyphy.relax_label_reference label a reference clade for RELAX\n"
            "    output.outdir, resources.cpu"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    _add_common_args(p_hyphy)
    p_hyphy.add_argument("--hyphy_bin", default=None,
                         help="HyPhy executable name or full path (e.g. hyphy, "
                              "/usr/local/bin/hyphy).  Overrides hyphy.hyphy_bin.")
    p_hyphy.add_argument("--hyphy_tests", default=None,
                         help="Comma-separated list of HyPhy tests to run.  Supported "
                              "values: RELAX, aBSREL, MEME.  Example: 'RELAX,MEME'.  "
                              "Overrides hyphy.hyphy_tests in the config.")
    p_hyphy.add_argument("--min_clade_size", type=int, default=None,
                         help="Minimum number of sequences a clade must contain to "
                              "be included in a HyPhy test run (default: 4).  Clades "
                              "smaller than this are skipped.  Overrides "
                              "hyphy.min_clade_size in the config.")
    p_hyphy.add_argument("--label_mode", default=None,
                         choices=["crown", "stem"],
                         help="Branch-labelling strategy for test/reference clade "
                              "designation.  'crown': label only the internal branches "
                              "within the clade (default).  'stem': also label the "
                              "stem branch leading to the clade.  Overrides "
                              "hyphy.label_mode in the config.")

    # ── score-motifs ───────────────────────────────────────────────────────
    p_sm = sub.add_parser(
        "score-motifs",
        help="Score known motifs using ESM-2 attention weights",
        description=(
            "Quantifies the structural/functional importance of user-supplied\n"
            "sequence motifs by aggregating ESM-2 self-attention weights directed\n"
            "at the motif residue positions across all sequences in each HMM\n"
            "family.  High attention scores indicate that the model 'focuses' on\n"
            "the motif region when processing surrounding context, suggesting\n"
            "structural or functional significance.\n\n"
            "HOW IT WORKS\n"
            "------------\n"
            "  1. For each sequence, all occurrences of each motif are located.\n"
            "  2. ESM-2 attention tensors are extracted (last N layers, default: 4).\n"
            "  3. Mean attention weight directed at motif positions from all other\n"
            "     residues is computed per sequence.\n"
            "  4. Results are aggregated across the HMM family and written to TSV.\n\n"
            "TYPICAL USE (after embed)\n"
            "-------------------------\n"
            "  phylofoundry score-motifs --outdir ./results \\\n"
            "                            --motifs 'HPEVY,HPEVF'\n\n"
            "  Expects in --outdir:\n"
            "    fasta_per_hmm/    per-HMM FASTA files (protein sequences)\n"
            "  Requires: ESM-2 model loaded via embeddings.model config key.\n\n"
            "  Outputs written to --outdir/motifs/:\n"
            "    <HMM>.<MOTIF>.attention.tsv   per-sequence attention scores\n"
            "    <HMM>.<MOTIF>.summary.tsv     per-motif aggregate statistics\n"
            "    <HMM>.<MOTIF>.png             attention score distribution plot\n\n"
            "REQUIREMENTS\n"
            "------------\n"
            "  ESM-2 model ('esm' Python package or HuggingFace transformers)\n\n"
            "USING A CONFIG FILE\n"
            "-------------------\n"
            "  phylofoundry score-motifs --config config/config.yaml \\\n"
            "                            --outdir ./results\n\n"
            "  Relevant config keys:\n"
            "    motifs.enabled           must be true (auto-set when --motifs given)\n"
            "    motifs.motif_list        list of amino acid motif strings to score\n"
            "    motifs.attention_layers  last N ESM-2 layers to average (default: 4)\n"
            "    embeddings.model         ESM-2 model name (default: esm2_t33_650M_UR50D)\n"
            "    embeddings.device        'cuda' or 'cpu'\n"
            "    output.outdir"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    _add_common_args(p_sm)
    p_sm.add_argument("--motifs", default=None,
                      help="Comma-separated list of amino acid motifs to score "
                           "(e.g. 'HPEVY,HPEVF').  Overrides motifs.motif_list "
                           "in the config.")
    p_sm.add_argument("--attention_layers", type=int, default=None,
                      help="Number of ESM-2 transformer layers (counting from the "
                           "last) to average when computing attention scores (default: "
                           "4).  Higher values capture more global attention patterns; "
                           "lower values focus on local interactions.  Overrides "
                           "motifs.attention_layers in the config.")

    # ── discover-motifs ────────────────────────────────────────────────────
    p_dm = sub.add_parser(
        "discover-motifs",
        help="Discover novel motifs using ESM-2 attention and k-mer analysis",
        description=(
            "Discovers clade-specific sequence motifs by comparing ESM-2\n"
            "attention profiles and enriched k-mers between a 'standard' reference\n"
            "clade and a 'novel' focal clade.  Identifies positions that are\n"
            "differentially attended to by the language model and enriched for\n"
            "specific amino acid patterns in the novel clade.\n\n"
            "HOW IT WORKS\n"
            "------------\n"
            "  1. For each HMM family, sequences are split into standard and novel\n"
            "     clade groups (from detect-clades, or manually specified).\n"
            "  2. Per-residue ESM-2 attention profiles are computed for all sequences.\n"
            "  3. Positions with significantly higher attention in the novel clade\n"
            "     are identified (high-attention peaks).\n"
            "  4. k-mer enrichment analysis finds amino acid patterns that are\n"
            "     over-represented around high-attention peaks in the novel clade.\n"
            "  5. Candidate motifs are ranked by a composite score (attention delta,\n"
            "     amino-acid shift, Jensen-Shannon divergence).\n\n"
            "TYPICAL USE (after embed and detect-clades)\n"
            "-------------------------------------------\n"
            "  phylofoundry discover-motifs --outdir ./results\n\n"
            "  Expects in --outdir:\n"
            "    fasta_per_hmm/       per-HMM FASTA files\n"
            "    clade_assignments/   per-HMM clade TSVs (from detect-clades)\n\n"
            "  Outputs written to --outdir/discover/:\n"
            "    <HMM>.attention_peaks.tsv   high-attention positions per sequence\n"
            "    <HMM>.kmer_enrichment.tsv   enriched k-mers in novel clade\n"
            "    <HMM>.candidate_motifs.tsv  ranked candidate motif list\n"
            "    <HMM>.motif_plot.png        attention profile comparison plot\n\n"
            "SPECIFYING CLADES\n"
            "-----------------\n"
            "  By default, the two largest HDBSCAN clusters from detect-clades are\n"
            "  used as the standard and novel clades.  Override with:\n"
            "    --standard_clade  <clade name>  reference/background clade\n"
            "    --novel_clade     <clade name>  focal/foreground clade\n\n"
            "REQUIREMENTS\n"
            "------------\n"
            "  ESM-2 model ('esm' Python package or HuggingFace transformers)\n\n"
            "USING A CONFIG FILE\n"
            "-------------------\n"
            "  phylofoundry discover-motifs --config config/config.yaml \\\n"
            "                               --outdir ./results\n\n"
            "  Relevant config keys:\n"
            "    discover.enabled          must be true (auto-set standalone)\n"
            "    discover.standard_clade   reference clade ID/name (null = auto)\n"
            "    discover.novel_clade      focal clade ID/name (null = auto)\n"
            "    discover.kmer_size        k-mer length for enrichment (default: 5)\n"
            "    discover.top_n_peaks      top attention peaks to analyse (default: 20)\n"
            "    discover.attention_layers last N ESM-2 layers to average (default: 4)\n"
            "    discover.candidates.top_n top N candidate motifs to report (default: 200)\n"
            "    discover.cross_hmm_comparison  compare motifs across HMM families\n"
            "    embeddings.model          ESM-2 model name\n"
            "    embeddings.device         'cuda' or 'cpu'\n"
            "    output.outdir"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    _add_common_args(p_dm)
    p_dm.add_argument("--standard_clade", default=None,
                      help="Name or ID of the reference (standard/background) clade "
                           "to compare against.  When not set, the largest HDBSCAN "
                           "cluster is used automatically.  Overrides "
                           "discover.standard_clade in the config.")
    p_dm.add_argument("--novel_clade", default=None,
                      help="Name or ID of the focal (novel/foreground) clade to "
                           "discover motifs for.  When not set, the second-largest "
                           "HDBSCAN cluster is used automatically.  Overrides "
                           "discover.novel_clade in the config.")
    p_dm.add_argument("--kmer_size", type=int, default=None,
                      help="Length of k-mers used for amino acid enrichment analysis "
                           "(default: 5).  Larger k-mers capture longer sequence "
                           "patterns but require more sequences for reliable statistics. "
                           "Overrides discover.kmer_size in the config.")
    p_dm.add_argument("--attention_layers", type=int, default=None,
                      help="Number of ESM-2 transformer layers (counting from the "
                           "last) to average when computing attention profiles "
                           "(default: 4).  Overrides discover.attention_layers in the "
                           "config.")

    # ── Utility subcommands ────────────────────────────────────────────────

    sub.add_parser(
        "list-steps",
        help="List all registered workflow steps with descriptions and exit",
    )

    p_plan = sub.add_parser(
        "plan",
        help="Show the execution plan for a given config and exit",
        description=(
            "Resolves the configuration and prints which steps would run "
            "without actually executing any pipeline code."
        ),
    )
    _add_common_args(p_plan)
    p_plan.add_argument("--start_at", choices=STEPS, default=None,
                        help="Override workflow.start_at")
    p_plan.add_argument("--stop_after", choices=STEPS, default=None,
                        help="Override workflow.stop_after")

    p_val = sub.add_parser(
        "validate",
        help="Validate configuration without running the pipeline and exit",
        description=(
            "Checks that all required config fields are present and valid "
            "without launching any pipeline steps."
        ),
    )
    _add_common_args(p_val)

    p_doc = sub.add_parser(
        "doctor",
        help="Check tool availability and Python-package health, then exit",
        description=(
            "Scans PATH for required and optional external tools, checks "
            "optional Python packages, and optionally validates a config."
        ),
    )
    _add_common_args(p_doc)

    sub.add_parser(
        "dump-config",
        help="Print the default annotated YAML config to stdout and exit",
    )

    p_lr = sub.add_parser(
        "list-runs",
        help="List available prior runs in a directory and show their metadata",
        description=(
            "Scans a directory (default: current working directory) for "
            "PhyloFoundry output directories and prints a table with run IDs, "
            "timestamps, completed steps, and artifact paths.  Use this to "
            "identify which directories can be passed to --input-run."
        ),
    )
    p_lr.add_argument(
        "search_dir",
        nargs="?",
        default=".",
        metavar="DIR",
        help="Directory to scan for prior runs (default: current working directory)",
    )
    p_lr.add_argument(
        "--json",
        action="store_true",
        dest="list_runs_json",
        help="Output run metadata as JSON instead of a human-readable table",
    )

    return ap


# ── Step-specific config override application ──────────────────────────────────

def _apply_step_args(args: argparse.Namespace, cfg: dict) -> None:
    """Apply step-specific CLI overrides from *args* to *cfg* in-place."""
    subcmd = getattr(args, "subcommand", None)

    if subcmd == "embed":
        if getattr(args, "model", None):
            cfg["embeddings"]["model"] = args.model
        if getattr(args, "device", None):
            cfg["embeddings"]["device"] = args.device
        if getattr(args, "batch_size", None) is not None:
            cfg["embeddings"]["batch_size"] = args.batch_size
        if getattr(args, "backend", None):
            cfg["embeddings"]["backend"] = args.backend
        if getattr(args, "repr_layer", None) is not None:
            cfg["embeddings"]["repr_layer"] = args.repr_layer
        if getattr(args, "fasta_dir", None):
            cfg["inputs"]["fasta_dir"] = args.fasta_dir

    elif subcmd == "maape":
        if getattr(args, "knn_k", None) is not None:
            cfg.setdefault("maape", {})["knn_k"] = args.knn_k
        if getattr(args, "knn_threshold", None) is not None:
            cfg.setdefault("maape", {})["knn_threshold"] = args.knn_threshold
        if getattr(args, "window_sizes", None):
            ws = [int(w.strip()) for w in args.window_sizes.split(",") if w.strip()]
            if ws:
                cfg.setdefault("maape", {})["window_sizes"] = ws
        if getattr(args, "no_aggregated", False):
            cfg.setdefault("maape", {})["generate_aggregated"] = False

    elif subcmd == "phylo":
        if getattr(args, "mafft", False):
            cfg["phylo"]["mafft"] = True
        if getattr(args, "mafft_mode", None):
            cfg["phylo"]["mafft_mode"] = args.mafft_mode
        if getattr(args, "skip_clipkit", False):
            cfg["phylo"]["skip_clipkit"] = True
        if getattr(args, "combined", False):
            cfg["phylo"]["combined_tree"] = True
        if getattr(args, "iqtree_bin", None):
            cfg["phylo"]["iqtree_bin"] = args.iqtree_bin
        if getattr(args, "iq_boot", None) is not None:
            cfg["phylo"]["iq_boot"] = args.iq_boot
        if getattr(args, "fasta_dir", None):
            cfg["inputs"]["fasta_dir"] = args.fasta_dir

    elif subcmd == "curate":
        if getattr(args, "no_treeshrink", False):
            cfg.setdefault("curate", {})["use_treeshrink"] = False
        if getattr(args, "no_esm_filter", False):
            cfg.setdefault("curate", {})["use_esm_filter"] = False

    elif subcmd == "taxonomy":
        if getattr(args, "gtdb_dir", None):
            cfg["inputs"]["gtdb_dir"] = args.gtdb_dir
        if getattr(args, "taxonomy_file", None):
            cfg["inputs"]["taxonomy_file"] = args.taxonomy_file
        if getattr(args, "globdb_taxonomy", None):
            cfg["inputs"]["globdb_taxonomy_file"] = args.globdb_taxonomy

    elif subcmd == "conservation":
        if getattr(args, "clades_tsv", None):
            cfg["conservation_metrics"]["clades_tsv"] = args.clades_tsv
        if getattr(args, "conservation_metric", None):
            cfg.setdefault("conservation_metrics", {})["conservation_metric"] = args.conservation_metric
        if getattr(args, "compute_kl", False):
            cfg.setdefault("conservation_metrics", {})["compute_kl"] = True
        if getattr(args, "kl_pairs", None):
            cfg.setdefault("conservation_metrics", {})["kl_pairs"] = args.kl_pairs

    elif subcmd == "detect-clades":
        if getattr(args, "detect_method", None):
            cfg["detect_clades"]["detect_method"] = args.detect_method

    elif subcmd == "aa-composition":
        if getattr(args, "no_plots", False):
            cfg.setdefault("aa_composition", {})["generate_plots"] = False
        if getattr(args, "no_pi", False):
            cfg.setdefault("aa_composition", {})["compute_pi"] = False

    elif subcmd == "hmmer":
        if getattr(args, "diamond_query", None):
            cfg["inputs"]["diamond_query"] = args.diamond_query
        if getattr(args, "diamond_mode", False):
            cfg.setdefault("diamond", {})["enabled"] = True
        if getattr(args, "diamond_db", None):
            cfg["inputs"]["diamond_db"] = args.diamond_db
        if getattr(args, "combined_faa", None):
            cfg["inputs"]["combined_faa"] = args.combined_faa
        if getattr(args, "min_score", None) is not None:
            cfg.setdefault("filtering", {})["global_min_score"] = args.min_score
        if getattr(args, "min_coverage", None) is not None:
            cfg.setdefault("filtering", {})["min_coverage"] = args.min_coverage

    elif subcmd == "extract":
        if getattr(args, "min_score", None) is not None:
            cfg.setdefault("filtering", {})["global_min_score"] = args.min_score
        if getattr(args, "min_coverage", None) is not None:
            cfg.setdefault("filtering", {})["min_coverage"] = args.min_coverage

    elif subcmd == "synteny":
        if getattr(args, "gbk_dir", None):
            cfg["synteny"]["gbk_dir"] = args.gbk_dir
        if getattr(args, "gff_dir", None):
            cfg["synteny"]["gff_dir"] = args.gff_dir
        if getattr(args, "window_genes", None) is not None:
            cfg.setdefault("synteny", {})["window_genes"] = args.window_genes

    elif subcmd == "codon":
        if getattr(args, "pal2nal_cmd", None):
            cfg["codon"]["pal2nal_cmd"] = args.pal2nal_cmd
        if getattr(args, "cds_id_mode", None):
            cfg["codon"]["cds_id_mode"] = args.cds_id_mode

    elif subcmd == "hyphy":
        if getattr(args, "hyphy_bin", None):
            cfg["hyphy"]["hyphy_bin"] = args.hyphy_bin
        if getattr(args, "hyphy_tests", None):
            cfg["hyphy"]["hyphy_tests"] = args.hyphy_tests
        if getattr(args, "min_clade_size", None) is not None:
            cfg.setdefault("hyphy", {})["min_clade_size"] = args.min_clade_size
        if getattr(args, "label_mode", None):
            cfg.setdefault("hyphy", {})["label_mode"] = args.label_mode

    elif subcmd == "score-motifs":
        if getattr(args, "motifs", None):
            motif_list = [m.strip() for m in args.motifs.split(",") if m.strip()]
            if motif_list:
                cfg.setdefault("motifs", {})["enabled"] = True
                cfg["motifs"]["motif_list"] = motif_list
        if getattr(args, "attention_layers", None) is not None:
            cfg.setdefault("motifs", {})["attention_layers"] = args.attention_layers

    elif subcmd == "discover-motifs":
        if getattr(args, "standard_clade", None):
            cfg.setdefault("discover", {})["standard_clade"] = args.standard_clade
        if getattr(args, "novel_clade", None):
            cfg.setdefault("discover", {})["novel_clade"] = args.novel_clade
        if getattr(args, "kmer_size", None) is not None:
            cfg.setdefault("discover", {})["kmer_size"] = args.kmer_size
        if getattr(args, "attention_layers", None) is not None:
            cfg.setdefault("discover", {})["attention_layers"] = args.attention_layers

    elif subcmd == "run":
        if getattr(args, "combined", False):
            cfg["phylo"]["combined_tree"] = True
        if getattr(args, "motifs", None):
            motif_list = [m.strip() for m in args.motifs.split(",") if m.strip()]
            if motif_list:
                cfg.setdefault("motifs", {})["enabled"] = True
                cfg["motifs"]["motif_list"] = motif_list
        if getattr(args, "diamond_query", None):
            cfg["inputs"]["diamond_query"] = args.diamond_query
        if getattr(args, "diamond_mode", False):
            cfg.setdefault("diamond", {})["enabled"] = True
        if getattr(args, "diamond_db", None):
            cfg["inputs"]["diamond_db"] = args.diamond_db
        if getattr(args, "combined_faa", None):
            cfg["inputs"]["combined_faa"] = args.combined_faa
        if getattr(args, "globdb_taxonomy", None):
            cfg["inputs"]["globdb_taxonomy_file"] = args.globdb_taxonomy

    elif subcmd == "prep":
        if getattr(args, "combined_faa", None):
            cfg["inputs"]["combined_faa"] = args.combined_faa


# ── Dependency-list builder ────────────────────────────────────────────────────

def _build_deps(cfg: dict, step_internal: "str | None" = None) -> list[str]:
    """Return the list of external tool names required for the run.

    When *step_internal* is given (single-step run), only tools relevant to
    that step are returned, avoiding spurious dependency errors for tools that
    belong to steps that will not run.
    """
    # Steps with no external tool requirements
    _no_deps = {
        "prep", "extract", "embed", "maape", "score_motifs", "discover_motifs",
        "taxonomy_integrate", "conservation_metrics", "detect_clades",
        "post", "curate",
    }
    if step_internal in _no_deps:
        return []
    if step_internal == "hmmer":
        if cfg.get("diamond", {}).get("enabled", False):
            return ["diamond"]
        return ["hmmscan", "hmmsearch"]
    if step_internal == "phylo":
        return [cfg["phylo"].get("iqtree_bin", "iqtree"), "mafft", "clipkit"]
    if step_internal == "synteny":
        sim = cfg["synteny"].get("similarity", {}).get("method", "diamond")
        if sim == "mmseqs":
            return [cfg["synteny"].get("similarity", {}).get("mmseqs_bin", "mmseqs")]
        return [cfg["synteny"].get("similarity", {}).get("diamond_bin", "diamond")]
    if step_internal == "codon":
        return [cfg["codon"].get("pal2nal_cmd", "pal2nal.pl")]
    if step_internal == "hyphy":
        return [cfg["hyphy"].get("hyphy_bin", "hyphy")]

    # Full-pipeline dependency collection
    deps: list[str] = ["hmmscan", "hmmsearch"]
    if cfg.get("diamond", {}).get("enabled", False):
        deps = ["diamond"]
    if cfg.get("phylo"):
        deps += [cfg["phylo"].get("iqtree_bin", "iqtree"), "mafft", "clipkit"]
    if cfg.get("codon", {}).get("enabled", False):
        deps.append(cfg["codon"].get("pal2nal_cmd", "pal2nal.pl"))
    if cfg.get("hyphy", {}).get("enabled", False):
        deps.append(cfg["hyphy"].get("hyphy_bin", "hyphy"))
    if cfg.get("synteny", {}).get("enabled", False):
        sim = cfg["synteny"].get("similarity", {}).get("method", "diamond")
        if sim == "diamond":
            deps.append(
                cfg["synteny"].get("similarity", {}).get("diamond_bin", "diamond")
            )
        elif sim == "mmseqs":
            deps.append(
                cfg["synteny"].get("similarity", {}).get("mmseqs_bin", "mmseqs")
            )
    return deps


# ── Entry point ────────────────────────────────────────────────────────────────

def main() -> None:
    """PhyloFoundry command-line entry point.

    Supports both the new subcommand-style interface::

        phylofoundry embed --outdir ./results --cpu 8
        phylofoundry phylo --outdir ./results --combined

    and the legacy flag-style interface (automatically routed to the ``run``
    subcommand for backward compatibility)::

        phylofoundry --config config/config.yaml
        phylofoundry --faa_dir ./proteomes --hmm_dir ./markers --outdir ./results
    """
    argv = sys.argv[1:]

    # Backward compatibility: when no recognised subcommand is present (old-style
    # invocation such as ``phylofoundry --config foo.json``), prepend ``run`` so
    # the existing flag interface continues to work unchanged.
    if not argv or argv[0] not in _ALL_SUBCMDS:
        argv = ["run"] + argv

    ap = _build_parser()
    args = ap.parse_args(argv)

    # ── Utility commands that need no pipeline config ──────────────────────
    if args.subcommand == "list-steps":
        from .preflight import print_list_steps
        print_list_steps()
        sys.exit(0)

    if args.subcommand == "list-runs":
        from .provenance import list_runs
        list_runs(
            args.search_dir,
            as_json=getattr(args, "list_runs_json", False),
        )
        sys.exit(0)

    if args.subcommand == "dump-config":
        _print_default_config_yaml()
        sys.exit(0)

    # Legacy --list-steps flag on 'run' subcommand
    if args.subcommand == "run" and getattr(args, "list_steps", False):
        from .preflight import print_list_steps
        print_list_steps()
        sys.exit(0)

    # Ensure attributes used by resolve_config are always present
    for _attr, _default in (
        ("start_at", None),
        ("stop_after", None),
        ("dump_default_config", False),
    ):
        if not hasattr(args, _attr):
            setattr(args, _attr, _default)

    # For step subcommands, restrict execution to that single step
    step_internal: "str | None" = _STEP_SUBCMD_MAP.get(args.subcommand)
    if step_internal:
        args.start_at = step_internal
        args.stop_after = step_internal

    # ── Utility commands that need a resolved config ───────────────────────
    _is_plan = args.subcommand == "plan" or getattr(args, "plan", False)
    _is_validate = args.subcommand == "validate" or getattr(
        args, "validate_config", False
    )
    _is_doctor = args.subcommand == "doctor" or getattr(args, "doctor", False)

    if _is_plan or _is_validate or _is_doctor:
        cfg = resolve_config(args)
        if cfg is None:
            sys.exit(0)
        if _is_plan:
            from .preflight import print_plan
            print_plan(cfg)
            sys.exit(0)
        if _is_validate:
            from .preflight import validate_config_only
            ok = validate_config_only(cfg)
            sys.exit(0 if ok else 1)
        if _is_doctor:
            from .preflight import run_doctor
            ok = run_doctor(cfg)
            sys.exit(0 if ok else 1)

    # dump_default_config on the 'run' subcommand
    if getattr(args, "dump_default_config", False):
        _print_default_config_yaml()
        sys.exit(0)

    # ── Normal pipeline / single-step execution ────────────────────────────
    cfg = resolve_config(args)
    if cfg is None:
        sys.exit(0)

    # Apply subcommand-specific argument overrides
    _apply_step_args(args, cfg)

    # Auto-enable optional steps when invoked as a standalone subcommand.
    # Keys may use dot-notation to address nested config values, e.g.
    # "tree_viz.enabled" means cfg[section]["tree_viz"]["enabled"].
    if step_internal and step_internal in _STEP_ENABLE_MAP:
        section, key = _STEP_ENABLE_MAP[step_internal]
        if "." in key:
            # Walk the config hierarchy, creating intermediate dicts as needed.
            # All paths in _STEP_ENABLE_MAP are expected to point to dicts that
            # already exist after resolve_config() (e.g. cfg["phylo"]["tree_viz"]),
            # so the isinstance guard here is a safety net for unexpected configs.
            parts = key.split(".")
            target = cfg.setdefault(section, {})
            for part in parts[:-1]:
                if not isinstance(target.get(part), dict):
                    target[part] = {}
                target = target[part]
            target[parts[-1]] = True
        else:
            cfg.setdefault(section, {})[key] = True

    validate_config(cfg)

    from .utils.helpers import check_dependencies
    check_dependencies(_build_deps(cfg, step_internal))

    # ── Resume / checkpoint resolution ────────────────────────────────────
    resume_from = getattr(args, "resume_from", None)
    resume = getattr(args, "resume", False)
    no_resume = getattr(args, "no_resume", False)
    force = getattr(args, "force", False)

    # Embed resume flags into cfg so run_pipeline can use them
    cfg.setdefault("_checkpoint", {})
    cfg["_checkpoint"]["resume"] = resume
    cfg["_checkpoint"]["resume_from"] = resume_from
    cfg["_checkpoint"]["no_resume"] = no_resume

    # ── Input provenance (--input-run) ─────────────────────────────────────
    input_run = getattr(args, "input_run", None)
    if input_run is not None:
        cfg["_checkpoint"]["input_run"] = input_run

    from .pipeline import run_pipeline
    run_pipeline(cfg)


if __name__ == "__main__":
    main()

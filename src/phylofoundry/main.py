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
            "Combines all input FAA proteome files into a single\n"
            "combined_proteomes.faa and all HMM profiles into a combined.hmm\n"
            "database (including running hmmpress).\n\n"
            "TYPICAL USE\n"
            "-----------\n"
            "  phylofoundry prep --faa_dir ./proteomes --hmm_dir ./markers \\\n"
            "                    --outdir ./results\n\n"
            "  Outputs written to --outdir:\n"
            "    combined_proteomes.faa   concatenated protein FASTA\n"
            "    combined.hmm             pressed HMM database\n\n"
            "USING A CONFIG FILE\n"
            "-------------------\n"
            "  phylofoundry prep --config config/config.yaml --outdir ./results\n\n"
            "  Relevant config keys: inputs.faa_dir, inputs.hmm_input,\n"
            "  output.outdir, output.workdir (for large intermediate files),\n"
            "  prep.cleanup_combined_faa\n\n"
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
            "Runs hmmscan (per-genome) and hmmsearch (per-HMM) to identify\n"
            "competitive HMM hits across all input proteomes, then builds a\n"
            "filtered best-hits table.\n\n"
            "TYPICAL USE (after prep)\n"
            "------------------------\n"
            "  phylofoundry hmmer --outdir ./results\n\n"
            "  Expects in --outdir from a previous prep run:\n"
            "    combined_proteomes.faa\n"
            "    combined.hmm\n\n"
            "  Outputs written to --outdir:\n"
            "    hmmscan_tbl/           per-genome domain hit tables\n"
            "    hmmsearch_tbl/         per-HMM domain hit tables\n"
            "    summary/best_hits.competitive.tsv\n\n"
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
            "USING A CONFIG FILE\n"
            "-------------------\n"
            "  phylofoundry hmmer --config config/config.yaml --outdir ./results\n\n"
            "  Relevant config keys: inputs.faa_dir, inputs.hmm_input,\n"
            "  inputs.combined_faa, output.outdir, hmmer.run_scan,\n"
            "  hmmer.run_search, filtering.*, diamond.*"
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

    # ── extract ────────────────────────────────────────────────────────────
    p_extract = sub.add_parser(
        "extract",
        help="Extract sequences for HMM hits into per-HMM FASTA files",
        description=(
            "Reads the filtered best-hits table from the hmmer step and\n"
            "extracts the matching protein sequences from the combined proteome\n"
            "into individual per-HMM FASTA files under fasta_per_hmm/.\n\n"
            "TYPICAL USE (after hmmer)\n"
            "-------------------------\n"
            "  phylofoundry extract --outdir ./results\n\n"
            "  Expects in --outdir from a previous hmmer run:\n"
            "    combined_proteomes.faa\n"
            "    summary/best_hits.competitive.tsv\n\n"
            "  Outputs written to --outdir:\n"
            "    fasta_per_hmm/<HMM_NAME>.faa   one FASTA per HMM family\n\n"
            "USING A CONFIG FILE\n"
            "-------------------\n"
            "  phylofoundry extract --config config/config.yaml --outdir ./results\n\n"
            "  Relevant config keys: inputs.faa_dir, output.outdir,\n"
            "  filtering.global_min_score, filtering.min_coverage"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    _add_common_args(p_extract)

    # ── embed ──────────────────────────────────────────────────────────────
    p_embed = sub.add_parser(
        "embed",
        help="Generate protein language model embeddings (ESM-2 or HuggingFace)",
        description=(
            "Generates per-sequence embeddings using ESM-2 (default) or any\n"
            "HuggingFace transformer model.  Optionally runs PCA dimensionality\n"
            "reduction, UMAP projection, and HDBSCAN clustering.\n\n"
            "TYPICAL USE (after extract in a pipeline run)\n"
            "---------------------------------------------\n"
            "  phylofoundry embed --outdir ./results --cpu 8\n\n"
            "  Expects per-HMM FASTA files in --outdir/fasta_per_hmm/.\n"
            "  Outputs written to --outdir/embeddings/:\n"
            "    <HMM>.embeddings.npy     raw embedding vectors\n"
            "    <HMM>.pca.tsv            PCA coordinates + metadata\n"
            "    <HMM>.umap.tsv           UMAP coordinates\n"
            "    <HMM>.umap.png           UMAP scatter plot\n"
            "    <HMM>.knn.tsv            k-NN distance table\n"
            "    <HMM>.dispersion.tsv     per-cluster dispersion metrics\n\n"
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
            "USING A CONFIG FILE\n"
            "-------------------\n"
            "  phylofoundry embed --config config/config.yaml --outdir ./results\n\n"
            "  Relevant config keys:\n"
            "    embeddings.enabled        must be true (auto-set when run standalone)\n"
            "    embeddings.backend        'esm' (default) or 'transformers'\n"
            "    embeddings.model          model name/path (e.g. esm2_t33_650M_UR50D)\n"
            "    embeddings.device         'cuda' or 'cpu'\n"
            "    embeddings.batch_size     sequences per inference batch\n"
            "    embeddings.repr_layer     ESM-2 layer to extract (default: last)\n"
            "    embeddings.pca_components number of PCA components\n"
            "    embeddings.cluster_embeddings  run HDBSCAN clustering\n"
            "    embeddings.hdbscan_min_cluster_size\n"
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
            "language model embeddings using the MAAPE algorithm.  Reveals\n"
            "evolutionary relationships and directional flow between protein\n"
            "subfamilies in embedding space.\n\n"
            "TYPICAL USE (after embed in a pipeline run)\n"
            "-------------------------------------------\n"
            "  phylofoundry maape --outdir ./results --cpu 8\n\n"
            "  Expects embedding files in --outdir/embeddings/.\n"
            "  Outputs written to --outdir/maape/:\n"
            "    <HMM>.maape_network.pkl      serialised directed graph\n"
            "    <HMM>.edge_list.txt          weighted directed edge list\n"
            "    <HMM>.maape_network.png      evolutionary network plot\n"
            "    <HMM>.maape_aggregated.png   condensed cluster-level plot\n"
            "    <HMM>.paths.pkl              assembly paths (intermediate)\n"
            "    maape_summary.tsv            per-HMM network statistics\n\n"
            "USING A CONFIG FILE\n"
            "-------------------\n"
            "  phylofoundry maape --config config/config.yaml --outdir ./results\n\n"
            "  Relevant config keys:\n"
            "    maape.enabled                    must be true (auto-set standalone)\n"
            "    maape.window_sizes               sliding window sizes (default: [5,10,20,40,80])\n"
            "    maape.knn_k                      KNN graph neighbours (default: 20)\n"
            "    maape.knn_threshold              cosine similarity threshold (default: 0.5)\n"
            "    maape.pca_components             PCA dims (null = auto)\n"
            "    maape.similarity_threshold_base  base threshold at window_size=5\n"
            "    maape.generate_aggregated        enable Step 6 aggregated plot\n"
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
            "Aligns sequences with hmmalign (or MAFFT), trims with ClipKit,\n"
            "and infers maximum-likelihood trees with IQ-TREE.\n\n"
            "TYPICAL USE (after extract in a pipeline run)\n"
            "---------------------------------------------\n"
            "  phylofoundry phylo --outdir ./results --cpu 16\n\n"
            "  Expects per-HMM FASTA files in --outdir/fasta_per_hmm/.\n"
            "  Reads HMM profiles from inputs.hmm_input (needed for hmmalign).\n"
            "  Outputs written to --outdir:\n"
            "    alignments_hmm/          raw hmmalign/MAFFT alignments\n"
            "    alignments_clipkit/      ClipKit-trimmed alignments\n"
            "    trees_iqtree/            IQ-TREE .treefile outputs\n\n"
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
            "USING A CONFIG FILE\n"
            "-------------------\n"
            "  phylofoundry phylo --config config/config.yaml --outdir ./results\n\n"
            "  Relevant config keys:\n"
            "    inputs.fasta_dir    per-HMM FASTA directory (for standalone use)\n"
            "    inputs.hmm_input    HMM profiles directory (for hmmalign)\n"
            "    phylo.mafft         use MAFFT instead of hmmalign\n"
            "    phylo.mafft_mode    MAFFT alignment mode (auto, linsi, etc.)\n"
            "    phylo.iqtree_bin    IQ-TREE binary (auto-detected if not set)\n"
            "    phylo.iq_boot       bootstrap replicates (default: 1000)\n"
            "    phylo.combined_tree build one tree from all HMMs combined\n"
            "    phylo.skip_clipkit  skip ClipKit trimming\n"
            "    output.outdir, resources.cpu"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    _add_common_args(p_phylo)
    p_phylo.add_argument("--mafft", action="store_true",
                         help="Use MAFFT for multiple sequence alignment instead of "
                              "hmmalign.  Required when --hmm_dir is not available.  "
                              "Overrides phylo.mafft in the config.")
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
            "Identifies and removes outlier sequences by branch-length "
            "(TreeShrink) and/or ESM-embedding distance. Writes curated overlays "
            "to curated/ without overwriting the raw pipeline outputs.\n\n"
            "TYPICAL USE (after phylo)\n"
            "-------------------------\n"
            "  phylofoundry curate --outdir ./results\n\n"
            "  Expects in --outdir: trees_iqtree/, fasta_per_hmm/,\n"
            "  alignments_clipkit/, and (optionally) embeddings/.\n\n"
            "USING A CONFIG FILE\n"
            "-------------------\n"
            "  phylofoundry curate --config config/config.yaml --outdir ./results\n\n"
            "  Relevant config keys: curate.enabled, curate.use_treeshrink,\n"
            "  curate.use_esm_filter, output.outdir"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    _add_common_args(p_curate)

    # ── taxonomy ───────────────────────────────────────────────────────────
    p_tax = sub.add_parser(
        "taxonomy",
        help="Integrate taxonomy from GTDB-Tk output or a custom genome→lineage TSV",
        description=(
            "Loads taxonomy from a GTDB-Tk summary directory\n"
            "(inputs.gtdb_dir) or a custom TSV (inputs.taxonomy_file) and\n"
            "annotates the best_hits table with lineage information.\n\n"
            "TYPICAL USE (after hmmer)\n"
            "-------------------------\n"
            "  phylofoundry taxonomy --outdir ./results \\\n"
            "                        --gtdb_dir ./gtdbtk_output\n\n"
            "  Outputs written to --outdir:\n"
            "    summary/best_hits.with_taxonomy.tsv\n\n"
            "CUSTOM TAXONOMY FILE\n"
            "--------------------\n"
            "  phylofoundry taxonomy --outdir ./results \\\n"
            "    --taxonomy_file ./my_taxonomy.tsv\n\n"
            "  The TSV must have columns: genome (genome filename/ID) and\n"
            "  lineage (semicolon-separated rank string).\n\n"
            "USING A CONFIG FILE\n"
            "-------------------\n"
            "  phylofoundry taxonomy --config config/config.yaml --outdir ./results\n\n"
            "  Relevant config keys: inputs.gtdb_dir, inputs.taxonomy_file,\n"
            "  inputs.globdb_taxonomy_file, taxonomy_integrate.enabled,\n"
            "  output.outdir"
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
            "Uses scikit-bio to compute per-column conservation scores and\n"
            "Jensen-Shannon / KL divergence between user-specified clade groups.\n\n"
            "TYPICAL USE (after phylo)\n"
            "-------------------------\n"
            "  phylofoundry conservation --outdir ./results\n\n"
            "  Expects ClipKit-trimmed alignments in\n"
            "  --outdir/alignments_clipkit/.\n"
            "  Outputs written to --outdir/summary/post_scikitbio/.\n\n"
            "USING A CONFIG FILE\n"
            "-------------------\n"
            "  phylofoundry conservation --config config/config.yaml \\\n"
            "                            --outdir ./results\n\n"
            "  Relevant config keys: conservation_metrics.enabled,\n"
            "  conservation_metrics.clades_tsv, output.outdir"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    _add_common_args(p_cons)
    p_cons.add_argument("--clades_tsv", default=None,
                        help="Path to a clade assignment TSV (tab-separated, columns: "
                             "protein_id, clade).  Used to compute per-clade KL/JS "
                             "divergence.  Overrides conservation_metrics.clades_tsv.")

    # ── detect-clades ──────────────────────────────────────────────────────
    p_dc = sub.add_parser(
        "detect-clades",
        help="Detect clades via taxonomy, TreeCluster, or tree+embedding",
        description=(
            "Assigns protein sequences to clades using one of three strategies:\n"
            "  taxonomy    — GTDB rank-based grouping from taxonomy_integrate output\n"
            "  treecluster — branch-length thresholding via TreeCluster\n"
            "  tree_embed  — joint tree + embedding clustering\n\n"
            "TYPICAL USE (after phylo and taxonomy)\n"
            "--------------------------------------\n"
            "  phylofoundry detect-clades --outdir ./results\n\n"
            "  Outputs written to --outdir:\n"
            "    clade_assignments/<HMM>.tsv   per-HMM clade assignment tables\n"
            "    summary/detected_clades.tsv   combined clade table\n\n"
            "USING A CONFIG FILE\n"
            "-------------------\n"
            "  phylofoundry detect-clades --config config/config.yaml \\\n"
            "                             --outdir ./results\n\n"
            "  Relevant config keys: detect_clades.enabled,\n"
            "  detect_clades.detect_method, output.outdir"
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
            "Computes per-gene amino acid fractional composition and biochemical\n"
            "metrics (Zc, GRAVY, nH2O, pI, S/N ratio) for all HMM-hit proteins,\n"
            "optionally runs clade-level Kruskal-Wallis + Benjamini-Hochberg FDR\n"
            "tests, and generates boxplots and heatmaps.\n\n"
            "TYPICAL USE (after hmmer/extract)\n"
            "---------------------------------\n"
            "  phylofoundry aa-composition --outdir ./results\n\n"
            "  Outputs written to --outdir/summary/aa_composition/.\n\n"
            "USING A CONFIG FILE\n"
            "-------------------\n"
            "  phylofoundry aa-composition --config config/config.yaml \\\n"
            "                              --outdir ./results\n\n"
            "  Relevant config keys: aa_composition.enabled,\n"
            "  aa_composition.generate_plots, aa_composition.compute_pi,\n"
            "  output.outdir"
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
            "clade detection in a single step.\n\n"
            "NOTE: For new workflows, prefer the dedicated 'conservation' and\n"
            "'detect-clades' subcommands which offer finer control.\n\n"
            "TYPICAL USE (after phylo)\n"
            "-------------------------\n"
            "  phylofoundry post --outdir ./results\n\n"
            "USING A CONFIG FILE\n"
            "-------------------\n"
            "  phylofoundry post --config config/config.yaml --outdir ./results\n\n"
            "  Relevant config keys: post.enabled, post.clades_tsv, output.outdir"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    _add_common_args(p_post)

    # ── synteny ────────────────────────────────────────────────────────────
    p_syn = sub.add_parser(
        "synteny",
        help="Visualise synteny context around HMM hits",
        description=(
            "Extracts gene neighbourhoods around HMM hits from GenBank (.gbk)\n"
            "or GFF3 annotation files and generates synteny track plots using\n"
            "clinker or pygenomeviz.\n\n"
            "TYPICAL USE (after hmmer)\n"
            "-------------------------\n"
            "  phylofoundry synteny --outdir ./results \\\n"
            "                       --gbk_dir ./genbank_files\n\n"
            "  Outputs written to --outdir.\n\n"
            "USING A CONFIG FILE\n"
            "-------------------\n"
            "  phylofoundry synteny --config config/config.yaml --outdir ./results\n\n"
            "  Relevant config keys: synteny.enabled, synteny.gbk_dir,\n"
            "  synteny.gff_dir, synteny.window_genes, synteny.genome_fasta_dir,\n"
            "  synteny.protein_id_field, synteny.gene_label_field, output.outdir"
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

    # ── codon ──────────────────────────────────────────────────────────────
    p_codon = sub.add_parser(
        "codon",
        help="Build codon-aware alignments with pal2nal",
        description=(
            "Converts protein alignments to codon-aware nucleotide alignments\n"
            "using pal2nal.  Requires protein alignments from the phylo step and\n"
            "nucleotide CDS sequences.\n\n"
            "TYPICAL USE (after phylo)\n"
            "-------------------------\n"
            "  phylofoundry codon --outdir ./results\n\n"
            "  Expects in --outdir: alignments_clipkit/ (protein alignments).\n"
            "  Requires inputs.cds_dir (CDS nucleotide sequences in FASTA format).\n"
            "  Outputs written to --outdir/codon_alignments/.\n\n"
            "CHAINING FROM A PRIOR PIPELINE RUN\n"
            "-----------------------------------\n"
            "  phylofoundry codon --input-run ./phylo_out --outdir ./codon_out\n\n"
            "USING A CONFIG FILE\n"
            "-------------------\n"
            "  phylofoundry codon --config config/config.yaml --outdir ./results\n\n"
            "  Relevant config keys: codon.enabled, codon.build_codon_alignments,\n"
            "  codon.cds_id_mode, codon.pal2nal_cmd, inputs.cds_dir, output.outdir"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    _add_common_args(p_codon)
    p_codon.add_argument("--pal2nal_cmd", default=None,
                         help="Name or full path to the pal2nal.pl script "
                              "(e.g. pal2nal.pl, /opt/pal2nal/pal2nal.pl).  "
                              "Overrides codon.pal2nal_cmd in the config.")

    # ── hyphy ──────────────────────────────────────────────────────────────
    p_hyphy = sub.add_parser(
        "hyphy",
        help="Run HyPhy selection tests (RELAX, aBSREL, MEME)",
        description=(
            "Runs HyPhy molecular-evolution tests on codon alignments to detect\n"
            "sites or branches under positive/relaxed selection.\n\n"
            "TYPICAL USE (after codon and detect-clades)\n"
            "-------------------------------------------\n"
            "  phylofoundry hyphy --outdir ./results\n\n"
            "  Expects in --outdir: codon_alignments/, trees_iqtree/,\n"
            "  and (optionally) clade_assignments/ for clade-aware tests.\n"
            "  Outputs written to --outdir/summary/hyphy/.\n\n"
            "CHAINING FROM A PRIOR PIPELINE RUN\n"
            "-----------------------------------\n"
            "  phylofoundry hyphy --input-run ./prior_run --outdir ./hyphy_out\n\n"
            "USING A CONFIG FILE\n"
            "-------------------\n"
            "  phylofoundry hyphy --config config/config.yaml --outdir ./results\n\n"
            "  Relevant config keys: hyphy.enabled, hyphy.run_hyphy,\n"
            "  hyphy.hyphy_bin, hyphy.hyphy_tests (RELAX,aBSREL,MEME),\n"
            "  hyphy.use_detected_clades, output.outdir"
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

    # ── score-motifs ───────────────────────────────────────────────────────
    p_sm = sub.add_parser(
        "score-motifs",
        help="Score known motifs using ESM-2 attention weights",
        description=(
            "Quantifies the structural importance of user-supplied sequence\n"
            "motifs by aggregating ESM-2 attention scores at the motif positions\n"
            "across all sequences.\n\n"
            "TYPICAL USE (after embed)\n"
            "-------------------------\n"
            "  phylofoundry score-motifs --outdir ./results \\\n"
            "                            --motifs 'HPEVY,HPEVF'\n\n"
            "  Expects embeddings in --outdir/embeddings/ and per-HMM FASTAs in\n"
            "  --outdir/fasta_per_hmm/.  Outputs written to --outdir/motifs/.\n\n"
            "USING A CONFIG FILE\n"
            "-------------------\n"
            "  phylofoundry score-motifs --config config/config.yaml \\\n"
            "                            --outdir ./results\n\n"
            "  Relevant config keys: motifs.enabled, motifs.motif_list,\n"
            "  motifs.attention_layers, embeddings.model, output.outdir"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    _add_common_args(p_sm)
    p_sm.add_argument("--motifs", default=None,
                      help="Comma-separated list of amino acid motifs to score "
                           "(e.g. 'HPEVY,HPEVF').  Overrides motifs.motif_list "
                           "in the config.")

    # ── discover-motifs ────────────────────────────────────────────────────
    p_dm = sub.add_parser(
        "discover-motifs",
        help="Discover novel motifs using ESM-2 attention and k-mer analysis",
        description=(
            "Discovers clade-specific sequence motifs by comparing ESM-2\n"
            "attention profiles and enriched k-mers across HDBSCAN-defined\n"
            "clades.\n\n"
            "TYPICAL USE (after embed)\n"
            "-------------------------\n"
            "  phylofoundry discover-motifs --outdir ./results\n\n"
            "  Expects embeddings in --outdir/embeddings/ and per-HMM FASTAs in\n"
            "  --outdir/fasta_per_hmm/.  Outputs written to --outdir/discover/.\n\n"
            "USING A CONFIG FILE\n"
            "-------------------\n"
            "  phylofoundry discover-motifs --config config/config.yaml \\\n"
            "                               --outdir ./results\n\n"
            "  Relevant config keys: discover.enabled, discover.standard_clade,\n"
            "  discover.novel_clade, discover.kmer_size, discover.top_n_peaks,\n"
            "  embeddings.model, output.outdir"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    _add_common_args(p_dm)

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
        if getattr(args, "combined", False):
            cfg["phylo"]["combined_tree"] = True
        if getattr(args, "iqtree_bin", None):
            cfg["phylo"]["iqtree_bin"] = args.iqtree_bin
        if getattr(args, "iq_boot", None) is not None:
            cfg["phylo"]["iq_boot"] = args.iq_boot
        if getattr(args, "fasta_dir", None):
            cfg["inputs"]["fasta_dir"] = args.fasta_dir

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

    elif subcmd == "synteny":
        if getattr(args, "gbk_dir", None):
            cfg["synteny"]["gbk_dir"] = args.gbk_dir
        if getattr(args, "gff_dir", None):
            cfg["synteny"]["gff_dir"] = args.gff_dir

    elif subcmd == "codon":
        if getattr(args, "pal2nal_cmd", None):
            cfg["codon"]["pal2nal_cmd"] = args.pal2nal_cmd

    elif subcmd == "hyphy":
        if getattr(args, "hyphy_bin", None):
            cfg["hyphy"]["hyphy_bin"] = args.hyphy_bin
        if getattr(args, "hyphy_tests", None):
            cfg["hyphy"]["hyphy_tests"] = args.hyphy_tests

    elif subcmd == "score-motifs":
        if getattr(args, "motifs", None):
            motif_list = [m.strip() for m in args.motifs.split(",") if m.strip()]
            if motif_list:
                cfg.setdefault("motifs", {})["enabled"] = True
                cfg["motifs"]["motif_list"] = motif_list

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

    # Auto-enable optional steps when invoked as a standalone subcommand
    if step_internal and step_internal in _STEP_ENABLE_MAP:
        section, key = _STEP_ENABLE_MAP[step_internal]
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

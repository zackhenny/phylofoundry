import argparse
import sys
from .config import resolve_config, validate_config, STEPS

# ── Constants ──────────────────────────────────────────────────────────────────

# Maps CLI subcommand names to internal pipeline step names.
# Hyphenated names are preferred for multi-word steps.
_STEP_SUBCMD_MAP: dict[str, str] = {
    "prep": "prep",
    "hmmer": "hmmer",
    "extract": "extract",
    "embed": "embed",
    "phylo": "phylo",
    "curate": "curate",
    "taxonomy": "taxonomy_integrate",
    "conservation": "conservation_metrics",
    "detect-clades": "detect_clades",
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
    "curate": ("curate", "enabled"),
    "taxonomy_integrate": ("taxonomy_integrate", "enabled"),
    "conservation_metrics": ("conservation_metrics", "enabled"),
    "detect_clades": ("detect_clades", "enabled"),
    "post": ("post", "enabled"),
    "synteny": ("synteny", "enabled"),
    "codon": ("codon", "enabled"),
    "hyphy": ("hyphy", "enabled"),
    "score_motifs": ("motifs", "enabled"),
    "discover_motifs": ("discover", "enabled"),
}

# All recognised top-level subcommand tokens (used for legacy routing).
_ALL_SUBCMDS: frozenset[str] = frozenset(
    {"run", "list-steps", "plan", "validate", "doctor", "dump-config"}
    | set(_STEP_SUBCMD_MAP.keys())
)


# ── Argument-parser helpers ────────────────────────────────────────────────────

def _add_common_args(parser: argparse.ArgumentParser) -> None:
    """Add shared arguments used by most subparsers."""
    parser.add_argument("--config", default=None,
                        help="JSON config file")
    parser.add_argument("--faa_dir", default=None,
                        help="Override inputs.faa_dir")
    parser.add_argument("--hmm_dir", default=None,
                        help="Override inputs.hmm_input (dir or single .hmm)")
    parser.add_argument("--outdir", default=None,
                        help="Override output.outdir")
    parser.add_argument("--cpu", type=int, default=None,
                        help="Override resources.cpu")
    parser.add_argument("--force", action="store_true",
                        help="Override workflow.force=True")


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
            "  phylofoundry plan  --config config.json\n"
            "  phylofoundry validate --config config.json\n"
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
            "Combines all input FAA proteome files into a single combined_proteomes.faa "
            "and all HMM profiles into a combined.hmm database."
        ),
    )
    _add_common_args(p_prep)
    p_prep.add_argument("--combined_faa", default=None,
                        help="Path to a prebuilt combined proteomes FASTA "
                             "(overrides inputs.combined_faa; skips prep build)")
    p_hmmer = sub.add_parser(
        "hmmer",
        help="Run HMM scan/search (or DIAMOND blastp) against proteomes",
        description=(
            "Runs hmmscan and/or hmmsearch to find competitive HMM hits. "
            "Use --diamond_mode to run DIAMOND blastp instead."
        ),
    )
    _add_common_args(p_hmmer)
    p_hmmer.add_argument("--diamond_query", default=None,
                         help="Override inputs.diamond_query "
                              "(FASTA file/dir for DIAMOND blastp)")
    p_hmmer.add_argument("--diamond_mode", action="store_true",
                         help="Enable DIAMOND search mode")
    p_hmmer.add_argument("--diamond_db", default=None,
                         help="Path to a prebuilt DIAMOND .dmnd database "
                              "(overrides inputs.diamond_db; skips makedb step)")
    p_hmmer.add_argument("--combined_faa", default=None,
                         help="Path to a prebuilt combined proteomes FASTA "
                              "(overrides inputs.combined_faa)")

    # ── extract ────────────────────────────────────────────────────────────
    p_extract = sub.add_parser(
        "extract",
        help="Extract sequences for HMM hits into per-HMM FASTA files",
        description=(
            "Reads HMM hit tables and extracts the matching protein sequences "
            "from the combined proteome into individual per-HMM FASTA files."
        ),
    )
    _add_common_args(p_extract)

    # ── embed ──────────────────────────────────────────────────────────────
    p_embed = sub.add_parser(
        "embed",
        help="Generate protein language model embeddings (ESM-2 or HuggingFace)",
        description=(
            "Generates per-sequence embeddings using ESM-2 (default) or any "
            "HuggingFace transformer model, then runs PCA/UMAP and optional "
            "HDBSCAN clustering."
        ),
    )
    _add_common_args(p_embed)
    p_embed.add_argument("--model", default=None,
                         help="Model name or path "
                              "(overrides embeddings.model, "
                              "e.g. esm2_t6_8M_UR50D)")
    p_embed.add_argument("--device", default=None,
                         help="Compute device: 'cuda' or 'cpu' "
                              "(overrides embeddings.device)")
    p_embed.add_argument("--batch_size", type=int, default=None,
                         help="Batch size for embedding inference "
                              "(overrides embeddings.batch_size)")
    p_embed.add_argument("--backend", default=None,
                         choices=["esm", "transformers"],
                         help="Embedding backend "
                              "(overrides embeddings.backend)")

    # ── phylo ──────────────────────────────────────────────────────────────
    p_phylo = sub.add_parser(
        "phylo",
        help="Run multiple sequence alignment and IQ-TREE phylogenetic inference",
        description=(
            "Aligns sequences with hmmalign (or MAFFT), trims with ClipKit, "
            "and infers maximum-likelihood trees with IQ-TREE."
        ),
    )
    _add_common_args(p_phylo)
    p_phylo.add_argument("--mafft", action="store_true",
                         help="Use MAFFT alignment instead of hmmalign "
                              "(overrides phylo.mafft)")
    p_phylo.add_argument("--combined", action="store_true",
                         help="Build a single combined tree from all HMMs "
                              "(overrides phylo.combined_tree)")
    p_phylo.add_argument("--iqtree_bin", default=None,
                         help="IQ-TREE binary name or path "
                              "(overrides phylo.iqtree_bin)")
    p_phylo.add_argument("--iq_boot", type=int, default=None,
                         help="IQ-TREE bootstrap replicates "
                              "(overrides phylo.iq_boot)")

    # ── curate ─────────────────────────────────────────────────────────────
    p_curate = sub.add_parser(
        "curate",
        help="Curate sequences using TreeShrink outlier pruning and/or ESM filtering",
        description=(
            "Identifies and removes outlier sequences by branch-length "
            "(TreeShrink) and/or ESM-embedding distance. Writes to curated/ "
            "overlay without overwriting raw pipeline outputs."
        ),
    )
    _add_common_args(p_curate)

    # ── taxonomy ───────────────────────────────────────────────────────────
    p_tax = sub.add_parser(
        "taxonomy",
        help="Integrate taxonomy from GTDB-Tk output or a custom genome→lineage TSV",
        description=(
            "Loads taxonomy from a GTDB-Tk summary directory "
            "(inputs.gtdb_dir) or a custom TSV (inputs.taxonomy_file) "
            "and annotates best_hits tables."
        ),
    )
    _add_common_args(p_tax)
    p_tax.add_argument("--gtdb_dir", default=None,
                       help="Override inputs.gtdb_dir")
    p_tax.add_argument("--taxonomy_file", default=None,
                       help="Override inputs.taxonomy_file "
                            "(TSV with columns: genome, lineage)")
    p_tax.add_argument("--globdb_taxonomy", default=None,
                       help="Path to a GlobDB-style headerless taxonomy TSV "
                            "(col1=genome_id, col2=GTDB taxonomy; overrides "
                            "inputs.globdb_taxonomy_file)")

    # ── conservation ───────────────────────────────────────────────────────
    p_cons = sub.add_parser(
        "conservation",
        help="Compute per-site conservation scores and KL divergence",
        description=(
            "Uses scikit-bio to compute per-column conservation metrics and "
            "Jensen-Shannon / KL divergence between user-specified clade groups."
        ),
    )
    _add_common_args(p_cons)
    p_cons.add_argument("--clades_tsv", default=None,
                        help="Path to clade assignment TSV "
                             "(overrides conservation_metrics.clades_tsv)")

    # ── detect-clades ──────────────────────────────────────────────────────
    p_dc = sub.add_parser(
        "detect-clades",
        help="Detect clades via taxonomy, TreeCluster, or tree+embedding",
        description=(
            "Assigns tip labels to clades using one of three strategies: "
            "GTDB taxonomy rank, TreeCluster branch-length thresholding, "
            "or tree+embedding joint analysis."
        ),
    )
    _add_common_args(p_dc)
    p_dc.add_argument("--detect_method", default=None,
                      choices=["taxonomy", "treecluster", "tree_embed"],
                      help="Clade detection method "
                           "(overrides detect_clades.detect_method)")

    # ── post ───────────────────────────────────────────────────────────────
    p_post = sub.add_parser(
        "post",
        help="Run legacy post-processing (conservation + clade detection combined)",
        description=(
            "Backward-compatibility shim that runs conservation metrics and "
            "clade detection in one step. Prefer 'conservation' and "
            "'detect-clades' for new workflows."
        ),
    )
    _add_common_args(p_post)

    # ── synteny ────────────────────────────────────────────────────────────
    p_syn = sub.add_parser(
        "synteny",
        help="Visualise synteny context around HMM hits",
        description=(
            "Extracts gene neighbourhoods around HMM hits from GenBank/GFF "
            "files and plots synteny tracks using clinker or pygenomeviz."
        ),
    )
    _add_common_args(p_syn)
    p_syn.add_argument("--gbk_dir", default=None,
                       help="Override synteny.gbk_dir "
                            "(directory of GenBank files)")
    p_syn.add_argument("--gff_dir", default=None,
                       help="Override synteny.gff_dir "
                            "(directory of GFF3 files)")

    # ── codon ──────────────────────────────────────────────────────────────
    p_codon = sub.add_parser(
        "codon",
        help="Build codon-aware alignments with pal2nal",
        description=(
            "Converts protein alignments to codon (nucleotide) alignments "
            "using pal2nal. Requires nucleotide CDS sequences in inputs.cds_dir."
        ),
    )
    _add_common_args(p_codon)
    p_codon.add_argument("--pal2nal_cmd", default=None,
                         help="pal2nal command name or path "
                              "(overrides codon.pal2nal_cmd)")

    # ── hyphy ──────────────────────────────────────────────────────────────
    p_hyphy = sub.add_parser(
        "hyphy",
        help="Run HyPhy selection tests (RELAX, aBSREL, MEME)",
        description=(
            "Runs HyPhy molecular-evolution tests on codon alignments. "
            "Supports clade-aware labelling driven by detect-clades output."
        ),
    )
    _add_common_args(p_hyphy)
    p_hyphy.add_argument("--hyphy_bin", default=None,
                         help="HyPhy binary name or path "
                              "(overrides hyphy.hyphy_bin)")
    p_hyphy.add_argument("--hyphy_tests", default=None,
                         help="Comma-separated HyPhy tests to run "
                              "(overrides hyphy.hyphy_tests, "
                              "e.g. 'RELAX,aBSREL,MEME')")

    # ── score-motifs ───────────────────────────────────────────────────────
    p_sm = sub.add_parser(
        "score-motifs",
        help="Score known motifs using ESM-2 attention weights",
        description=(
            "Quantifies structural importance of user-supplied sequence motifs "
            "by aggregating ESM-2 attention scores at the motif positions."
        ),
    )
    _add_common_args(p_sm)
    p_sm.add_argument("--motifs", default=None,
                      help="Comma-separated motif list to score "
                           "(e.g. 'HPEVY,HPEVF')")

    # ── discover-motifs ────────────────────────────────────────────────────
    p_dm = sub.add_parser(
        "discover-motifs",
        help="Discover novel motifs using ESM-2 attention and k-mer analysis",
        description=(
            "Compares attention profiles across HDBSCAN clades to identify "
            "clade-specific high-attention sequence regions (novel motifs)."
        ),
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
        help="Print the default config JSON to stdout and exit",
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

    elif subcmd == "phylo":
        if getattr(args, "mafft", False):
            cfg["phylo"]["mafft"] = True
        if getattr(args, "combined", False):
            cfg["phylo"]["combined_tree"] = True
        if getattr(args, "iqtree_bin", None):
            cfg["phylo"]["iqtree_bin"] = args.iqtree_bin
        if getattr(args, "iq_boot", None) is not None:
            cfg["phylo"]["iq_boot"] = args.iq_boot

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
        "prep", "extract", "embed", "score_motifs", "discover_motifs",
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

        phylofoundry --config config.json
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

    if args.subcommand == "dump-config":
        import json
        from .constants import DEFAULT_CONFIG
        print(json.dumps(DEFAULT_CONFIG, indent=2, sort_keys=True))
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
        import json
        from .constants import DEFAULT_CONFIG
        print(json.dumps(DEFAULT_CONFIG, indent=2, sort_keys=True))
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

    from .pipeline import run_pipeline
    run_pipeline(cfg)


if __name__ == "__main__":
    main()

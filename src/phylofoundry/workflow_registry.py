"""PhyloFoundry workflow registry.

Registers all current workflow steps with their metadata contracts. The
registry is validated at module import time so that dependency errors are
caught early.

Architecture notes
------------------
The former monolithic ``post`` step has been split into three focused steps:

- ``taxonomy_integrate`` — load taxonomy and annotate best_hits
- ``conservation_metrics`` — scikit-bio conservation / KL divergence metrics
- ``detect_clades`` — clade detection (taxonomy / TreeCluster / tree+embedding)

The legacy ``post`` step is preserved as a backward-compatibility shim.

Other planned renames (not yet implemented):

- ``tasks/motifs.py`` → ``tasks/score_motifs.py``
- ``tasks/discover.py`` → ``tasks/discover_motifs.py``
- ``curate`` will write only to a ``curated/`` overlay directory, never
  overwriting the raw ``trees_iqtree/`` outputs.
"""

from .workflow_schema import (
    ArtifactKind,
    ArtifactSpec,
    CheckpointStrategy,
    StepDefinition,
    ToolRequirement,
    WorkflowRegistry,
)

# ---------------------------------------------------------------------------
# Registry instance
# ---------------------------------------------------------------------------

REGISTRY = WorkflowRegistry()

# ---------------------------------------------------------------------------
# Step definitions (listed in default pipeline execution order)
# ---------------------------------------------------------------------------

REGISTRY.register(StepDefinition(
    name="prep",
    description="Combine input FAA files and HMM files into pipeline inputs.",
    outputs=[
        ArtifactSpec("combined_proteomes", ArtifactKind.FASTA, "combined_proteomes.faa"),
        ArtifactSpec("combined_hmm", ArtifactKind.FASTA, "combined.hmm"),
    ],
    tool_requirements=[
        ToolRequirement(
            "hmmpress",
            optional=True,
            description="Presses the combined HMM database for hmmscan.",
        ),
    ],
    dependencies=[],
    optional=False,
))

REGISTRY.register(StepDefinition(
    name="hmmer",
    description="Run hmmscan and/or hmmsearch; filter and aggregate hits.",
    outputs=[
        ArtifactSpec(
            "hmmscan_hits",
            ArtifactKind.TABLE,
            "summary/hmmscan_hits.filtered.tsv",
        ),
        ArtifactSpec(
            "hmmsearch_hits",
            ArtifactKind.TABLE,
            "summary/hmmsearch_hits.filtered.tsv",
        ),
    ],
    tool_requirements=[
        ToolRequirement("hmmscan", optional=True),
        ToolRequirement("hmmsearch", optional=True),
    ],
    dependencies=["prep"],
    optional=False,
))

REGISTRY.register(StepDefinition(
    name="extract",
    description="Extract sequences from FASTA inputs that matched HMM hits.",
    outputs=[
        ArtifactSpec("fasta_per_hmm", ArtifactKind.DIRECTORY, "fasta_per_hmm/"),
    ],
    dependencies=["hmmer"],
    optional=False,
))

REGISTRY.register(StepDefinition(
    name="embed",
    description=(
        "Generate protein language model embeddings (ESM-2 or HuggingFace transformers)."
    ),
    outputs=[
        ArtifactSpec("embeddings", ArtifactKind.DIRECTORY, "embeddings/"),
    ],
    tool_requirements=[
        ToolRequirement(
            "torch",
            optional=True,
            description="Required for ESM-2 or HuggingFace embeddings.",
        ),
    ],
    dependencies=["extract"],
    optional=True,
    enabled_config_key="embeddings.enabled",
))

REGISTRY.register(StepDefinition(
    name="phylo",
    description=(
        "Run multiple sequence alignment (hmmalign/mafft) and phylogenetic tree "
        "inference (IQ-TREE), optionally followed by ASR parsing."
    ),
    outputs=[
        ArtifactSpec("alignments_hmm", ArtifactKind.DIRECTORY, "alignments_hmm/"),
        ArtifactSpec(
            "alignments_clipkit",
            ArtifactKind.DIRECTORY,
            "alignments_clipkit/",
        ),
        ArtifactSpec("trees_iqtree", ArtifactKind.DIRECTORY, "trees_iqtree/"),
    ],
    tool_requirements=[
        ToolRequirement(
            "hmmalign",
            description="Default multiple sequence alignment method.",
        ),
        ToolRequirement(
            "mafft",
            optional=True,
            description="Used when mafft alignment mode is selected.",
        ),
        ToolRequirement("clipkit", description="Alignment trimming."),
        ToolRequirement(
            "iqtree",
            config_key="phylo.iqtree_bin",
            description="Maximum-likelihood tree inference.",
        ),
    ],
    dependencies=["extract"],
    optional=False,
))

REGISTRY.register(StepDefinition(
    name="curate",
    description=(
        "Curate sequences using TreeShrink outlier pruning and/or ESM-based "
        "filtering."
    ),
    outputs=[
        ArtifactSpec(
            "curated",
            ArtifactKind.DIRECTORY,
            "curated/",
            optional=True,
        ),
    ],
    tool_requirements=[
        ToolRequirement(
            "treeshrink",
            optional=True,
            description="Optional branch-length outlier pruning.",
        ),
    ],
    dependencies=["phylo"],
    optional=True,
    enabled_config_key="curate.enabled",
    notes=(
        "Future target: curate will write only to a curated/ overlay directory "
        "and will never overwrite raw trees_iqtree/ outputs."
    ),
))

REGISTRY.register(StepDefinition(
    name="taxonomy_integrate",
    description=(
        "Load taxonomy from a GTDB-Tk summary directory (inputs.gtdb_dir) or "
        "a custom genome→lineage TSV (inputs.taxonomy_file) and write "
        "summary/genome_taxonomy.tsv and summary/best_hits.with_taxonomy.tsv."
    ),
    outputs=[
        ArtifactSpec(
            "genome_taxonomy",
            ArtifactKind.TABLE,
            "summary/genome_taxonomy.tsv",
            optional=True,
        ),
        ArtifactSpec(
            "best_hits_with_taxonomy",
            ArtifactKind.TABLE,
            "summary/best_hits.with_taxonomy.tsv",
            optional=True,
        ),
    ],
    dependencies=["phylo"],
    optional=True,
    enabled_config_key="taxonomy_integrate.enabled",
))

REGISTRY.register(StepDefinition(
    name="conservation_metrics",
    description=(
        "Compute per-site conservation scores and KL divergence using "
        "scikit-bio; write outputs under summary/post_scikitbio/."
    ),
    outputs=[
        ArtifactSpec(
            "conservation_metrics_dir",
            ArtifactKind.DIRECTORY,
            "summary/post_scikitbio/",
        ),
    ],
    dependencies=["phylo"],
    optional=True,
    enabled_config_key="conservation_metrics.enabled",
))

REGISTRY.register(StepDefinition(
    name="detect_clades",
    description=(
        "Detect clades via taxonomy, TreeCluster, or tree+embedding strategies; "
        "write summary/detected_clades.tsv and per-HMM clade_assignments/."
    ),
    outputs=[
        ArtifactSpec(
            "detected_clades",
            ArtifactKind.TABLE,
            "summary/detected_clades.tsv",
            optional=True,
        ),
        ArtifactSpec(
            "clade_assignments",
            ArtifactKind.DIRECTORY,
            "clade_assignments/",
        ),
    ],
    dependencies=["phylo"],
    optional=True,
    enabled_config_key="detect_clades.enabled",
))

REGISTRY.register(StepDefinition(
    name="post",
    description=(
        "Backward-compatibility shim: conservation metrics, KL divergence, and "
        "clade detection in a single step. Prefer the dedicated steps "
        "taxonomy_integrate, conservation_metrics, and detect_clades for new "
        "workflows."
    ),
    outputs=[
        ArtifactSpec(
            "post_scikitbio",
            ArtifactKind.DIRECTORY,
            "summary/post_scikitbio/",
        ),
        ArtifactSpec(
            "clade_assignments",
            ArtifactKind.DIRECTORY,
            "clade_assignments/",
        ),
    ],
    dependencies=["phylo"],
    optional=True,
    enabled_config_key="post.enabled",
    notes=(
        "Backward-compatibility shim. New workflows should use "
        "taxonomy_integrate, conservation_metrics, and detect_clades instead. "
        "Each of those steps has an independent enabled flag so that failures "
        "in one do not block the others."
    ),
))

REGISTRY.register(StepDefinition(
    name="synteny",
    description=(
        "Visualise synteny context around HMM hits using clinker / pygenomeviz."
    ),
    outputs=[
        ArtifactSpec("synteny", ArtifactKind.DIRECTORY, "synteny/"),
    ],
    tool_requirements=[
        ToolRequirement("clinker", optional=True),
    ],
    dependencies=["hmmer"],
    optional=True,
    enabled_config_key="synteny.enabled",
))

REGISTRY.register(StepDefinition(
    name="codon",
    description="Build codon-aware alignments with pal2nal.",
    outputs=[
        ArtifactSpec(
            "codon_alignments",
            ArtifactKind.DIRECTORY,
            "codon_alignments/",
        ),
    ],
    tool_requirements=[
        ToolRequirement(
            "pal2nal.pl",
            config_key="codon.pal2nal_cmd",
            description="Protein-to-nucleotide codon alignment.",
        ),
    ],
    dependencies=["phylo"],
    optional=True,
    enabled_config_key="codon.enabled",
))

REGISTRY.register(StepDefinition(
    name="hyphy",
    description="Run HyPhy selection tests (RELAX, aBSREL, MEME).",
    outputs=[
        ArtifactSpec("hyphy", ArtifactKind.DIRECTORY, "summary/hyphy/"),
    ],
    tool_requirements=[
        ToolRequirement("hyphy", config_key="hyphy.hyphy_bin"),
    ],
    dependencies=["codon", "detect_clades"],
    optional=True,
    enabled_config_key="hyphy.enabled",
    notes=(
        "Depends on detect_clades for clade assignments used in RELAX test "
        "labelling. Also accepts clade assignments written by the legacy post "
        "step when detect_clades is not enabled."
    ),
))

REGISTRY.register(StepDefinition(
    name="score_motifs",
    description=(
        "Score known motifs using ESM-2 attention weights. "
        "Current module: tasks/motifs.py "
        "(future rename target: tasks/score_motifs.py)."
    ),
    outputs=[
        ArtifactSpec(
            "motif_scores",
            ArtifactKind.TABLE,
            "summary/motif_scores.tsv",
            optional=True,
        ),
    ],
    tool_requirements=[
        ToolRequirement("torch", optional=True),
    ],
    dependencies=["extract"],
    optional=True,
    enabled_config_key="motifs.enabled",
    notes="Future: tasks/motifs.py will be renamed to tasks/score_motifs.py.",
))

REGISTRY.register(StepDefinition(
    name="discover_motifs",
    description=(
        "Discover novel motifs using ESM-2 attention and k-mer analysis. "
        "Current module: tasks/discover.py "
        "(future rename target: tasks/discover_motifs.py)."
    ),
    outputs=[
        ArtifactSpec(
            "discovered_motifs",
            ArtifactKind.TABLE,
            "summary/discovered_motifs.tsv",
            optional=True,
        ),
    ],
    tool_requirements=[
        ToolRequirement("torch", optional=True),
    ],
    dependencies=["extract", "detect_clades"],
    optional=True,
    enabled_config_key="discover.enabled",
    notes="Future: tasks/discover.py will be renamed to tasks/discover_motifs.py.",
))

# ---------------------------------------------------------------------------
# Validate the full registry at import time
# ---------------------------------------------------------------------------

REGISTRY.validate()

"""Artifact path resolution helpers for PhyloFoundry.

Centralizes canonical artifact path definitions and provides a single place
where path construction logic lives.  Task modules and the pipeline orchestrator
should import from here rather than constructing paths ad hoc.

Usage
-----
::

    from phylofoundry.artifact_paths import ArtifactPaths

    paths = ArtifactPaths(outdir)
    print(paths.best_hits_tsv)
    print(paths.treefile_for_hmm("MyHMM"))
    safe_mkdir(paths.trees_dir)

The module-level constants expose the *relative* path fragments so that
callers who already hold an absolute path can still introspect the canonical
names without constructing an :class:`ArtifactPaths` instance.
"""

from __future__ import annotations

import os
from dataclasses import dataclass

# ---------------------------------------------------------------------------
# Canonical artifact names / relative path fragments
# ---------------------------------------------------------------------------

# Top-level pipeline inputs reconstructed during prep
COMBINED_FAA = "combined_proteomes.faa"
COMBINED_HMM = "combined.hmm"

# Intermediate HMMER tables
HMMSCAN_TBL_DIR = "hmmscan_tbl"
HMMSEARCH_TBL_DIR = "hmmsearch_tbl"

# Per-HMM FASTA sequences extracted from hits
FASTA_PER_HMM_DIR = "fasta_per_hmm"

# Multiple-sequence alignments
ALIGNMENTS_HMM_DIR = "alignments_hmm"
ALIGNMENTS_CLIPKIT_DIR = "alignments_clipkit"

# IQ-TREE phylogenetic trees
TREES_IQTREE_DIR = "trees_iqtree"

# Protein language model embeddings
EMBEDDINGS_DIR = "embeddings"

# Codon-aware nucleotide alignments
CODON_ALIGNMENTS_DIR = "codon_alignments"

# Clade assignment tables (one per HMM)
CLADE_ASSIGNMENTS_DIR = "clade_assignments"

# Curated overlay root and sub-directories
CURATED_DIR = "curated"
CURATED_TREES_DIR = "curated/trees"
CURATED_FASTA_DIR = "curated/fasta_per_hmm"
CURATED_ALIGNMENTS_DIR = "curated/alignments_clipkit"
CURATED_MANIFESTS_DIR = "curated/manifests"
CURATED_AUDIT_DIR = "curated/audit"

# Summary directory and its key files
SUMMARY_DIR = "summary"
BEST_HITS_TSV = "summary/best_hits.competitive.tsv"
BEST_HITS_WITH_TAXONOMY_TSV = "summary/best_hits.with_taxonomy.tsv"
DETECTED_CLADES_TSV = "summary/detected_clades.tsv"
CONSERVATION_METRICS_DIR = "summary/post_scikitbio"
HYPHY_DIR = "summary/hyphy"
RESOLVED_CONFIG_JSON = "summary/resolved_config.json"

# Logs directory
LOGS_DIR = "logs"


# ---------------------------------------------------------------------------
# ArtifactPaths helper class
# ---------------------------------------------------------------------------


@dataclass
class ArtifactPaths:
    """Resolve canonical artifact paths relative to a pipeline output directory.

    All properties return absolute paths constructed from *outdir* and the
    canonical relative fragments defined at module level.

    Parameters
    ----------
    outdir:
        Absolute path to the pipeline output directory.
    """

    outdir: str

    # ── Combined pipeline inputs ────────────────────────────────────────────

    @property
    def combined_faa(self) -> str:
        """Combined proteome FASTA used for HMMER search/scan."""
        return os.path.join(self.outdir, COMBINED_FAA)

    @property
    def combined_hmm(self) -> str:
        """Combined HMM database built from all input HMM files."""
        return os.path.join(self.outdir, COMBINED_HMM)

    # ── HMMER hit tables ────────────────────────────────────────────────────

    @property
    def hmmscan_dir(self) -> str:
        """Directory containing per-genome hmmscan domain-table files."""
        return os.path.join(self.outdir, HMMSCAN_TBL_DIR)

    @property
    def hmmsearch_dir(self) -> str:
        """Directory containing per-HMM hmmsearch domain-table files."""
        return os.path.join(self.outdir, HMMSEARCH_TBL_DIR)

    # ── Per-HMM extracted sequences ─────────────────────────────────────────

    @property
    def fasta_dir(self) -> str:
        """Directory containing per-HMM FASTA files of extracted hits."""
        return os.path.join(self.outdir, FASTA_PER_HMM_DIR)

    # ── Alignments ──────────────────────────────────────────────────────────

    @property
    def alignments_hmm_dir(self) -> str:
        """Directory containing raw hmmalign (or mafft) alignment files."""
        return os.path.join(self.outdir, ALIGNMENTS_HMM_DIR)

    @property
    def alignments_clipkit_dir(self) -> str:
        """Directory containing ClipKIT-trimmed alignment files."""
        return os.path.join(self.outdir, ALIGNMENTS_CLIPKIT_DIR)

    # ── Phylogenetic trees ──────────────────────────────────────────────────

    @property
    def trees_dir(self) -> str:
        """Directory containing IQ-TREE .treefile outputs."""
        return os.path.join(self.outdir, TREES_IQTREE_DIR)

    # ── Embeddings ──────────────────────────────────────────────────────────

    @property
    def embeddings_dir(self) -> str:
        """Directory containing protein language model embedding files."""
        return os.path.join(self.outdir, EMBEDDINGS_DIR)

    # ── Codon alignments ────────────────────────────────────────────────────

    @property
    def codon_dir(self) -> str:
        """Directory containing codon-aware nucleotide alignment files."""
        return os.path.join(self.outdir, CODON_ALIGNMENTS_DIR)

    # ── Clade assignments ───────────────────────────────────────────────────

    @property
    def clade_assign_dir(self) -> str:
        """Directory containing per-HMM clade assignment TSV files."""
        return os.path.join(self.outdir, CLADE_ASSIGNMENTS_DIR)

    # ── Summary artifacts ───────────────────────────────────────────────────

    @property
    def summary_dir(self) -> str:
        """Top-level summary directory."""
        return os.path.join(self.outdir, SUMMARY_DIR)

    @property
    def best_hits_tsv(self) -> str:
        """Competitive best-hits TSV produced by the HMMER step."""
        return os.path.join(self.outdir, BEST_HITS_TSV)

    @property
    def best_hits_with_taxonomy_tsv(self) -> str:
        """Best-hits TSV annotated with taxonomy by the taxonomy_integrate step."""
        return os.path.join(self.outdir, BEST_HITS_WITH_TAXONOMY_TSV)

    @property
    def detected_clades_tsv(self) -> str:
        """Detected-clades TSV produced by the detect_clades step."""
        return os.path.join(self.outdir, DETECTED_CLADES_TSV)

    @property
    def conservation_metrics_dir(self) -> str:
        """Directory containing scikit-bio conservation metric outputs."""
        return os.path.join(self.outdir, CONSERVATION_METRICS_DIR)

    @property
    def hyphy_dir(self) -> str:
        """Directory containing HyPhy JSON result files."""
        return os.path.join(self.outdir, HYPHY_DIR)

    @property
    def resolved_config_json(self) -> str:
        """Resolved pipeline configuration snapshot written at startup."""
        return os.path.join(self.outdir, RESOLVED_CONFIG_JSON)

    # ── Logs ────────────────────────────────────────────────────────────────

    @property
    def logs_dir(self) -> str:
        """Directory containing pipeline and per-step log files."""
        return os.path.join(self.outdir, LOGS_DIR)

    # ── Curated overlay ─────────────────────────────────────────────────────

    @property
    def curated_dir(self) -> str:
        """Root directory for curated overlay artifacts."""
        return os.path.join(self.outdir, CURATED_DIR)

    @property
    def curated_trees_dir(self) -> str:
        """Directory containing curated per-HMM treefile overlays."""
        return os.path.join(self.outdir, CURATED_TREES_DIR)

    @property
    def curated_fasta_dir(self) -> str:
        """Directory containing curated per-HMM FASTA overlays."""
        return os.path.join(self.outdir, CURATED_FASTA_DIR)

    @property
    def curated_alignments_dir(self) -> str:
        """Directory containing curated per-HMM alignment overlays."""
        return os.path.join(self.outdir, CURATED_ALIGNMENTS_DIR)

    @property
    def curated_manifests_dir(self) -> str:
        """Directory containing per-HMM curation manifest JSON files."""
        return os.path.join(self.outdir, CURATED_MANIFESTS_DIR)

    @property
    def curated_audit_dir(self) -> str:
        """Directory containing append-only JSONL curation audit logs."""
        return os.path.join(self.outdir, CURATED_AUDIT_DIR)

    # ── Per-HMM artifact helpers ─────────────────────────────────────────────

    def fasta_for_hmm(self, hmm: str) -> str:
        """Return the path to the extracted per-HMM FASTA file."""
        return os.path.join(self.fasta_dir, f"{hmm}.faa")

    def treefile_for_hmm(self, hmm: str) -> str:
        """Return the path to the IQ-TREE treefile for *hmm*."""
        return os.path.join(self.trees_dir, f"{hmm}.treefile")

    def alignment_hmm_for_hmm(self, hmm: str) -> str:
        """Return the path to the raw hmmalign alignment for *hmm*."""
        return os.path.join(self.alignments_hmm_dir, f"{hmm}.afa")

    def alignment_clipkit_for_hmm(self, hmm: str) -> str:
        """Return the path to the ClipKIT-trimmed alignment for *hmm*."""
        return os.path.join(self.alignments_clipkit_dir, f"{hmm}.clipkit.faa")

    def codon_alignment_for_hmm(self, hmm: str) -> str:
        """Return the path to the codon-aware nucleotide alignment for *hmm*."""
        return os.path.join(self.codon_dir, f"{hmm}.codon.fasta")

    def hyphy_output_for_hmm(self, hmm: str, test: str, clade_name: str) -> str:
        """Return the HyPhy JSON output path for a clade-level test run.

        Parameters
        ----------
        hmm:
            HMM (gene family) name.
        test:
            HyPhy test name; case-insensitive (e.g. ``"relax"``, ``"RELAX"``,
            ``"aBSREL"``, ``"MEME"``).  The name is normalized to upper-case
            in the output directory structure.
        clade_name:
            Clade label used during labelled-tree HyPhy runs.
        """
        # HyPhy test names are stored in upper-case directories to match the
        # convention established in tasks/hyphy.py.
        return os.path.join(self.hyphy_dir, hmm, test.upper(), f"{clade_name}.json")

    # ── Bulk directory list ──────────────────────────────────────────────────

    def all_output_dirs(self) -> list:
        """Return the ordered list of directories that should exist at pipeline startup.

        The pipeline orchestrator creates each of these at startup with
        :func:`~phylofoundry.utils.helpers.safe_mkdir`.
        """
        return [
            self.hmmscan_dir,
            self.hmmsearch_dir,
            self.fasta_dir,
            self.alignments_hmm_dir,
            self.alignments_clipkit_dir,
            self.trees_dir,
            self.summary_dir,
            self.conservation_metrics_dir,
            self.codon_dir,
            self.hyphy_dir,
            self.embeddings_dir,
            self.clade_assign_dir,
            self.curated_dir,
        ]

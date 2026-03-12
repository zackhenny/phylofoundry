"""Tests for the curate overlay-only refactor (Phase 7)."""

import glob
import json
import os

import pytest
from Bio import Phylo, SeqIO


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

SIMPLE_TREE = "(seqA:0.1,seqB:0.2,(seqC:0.3,seqD:0.4):0.5);\n"
SIMPLE_FASTA = ">seqA\nMAA\n>seqB\nMBB\n>seqC\nMCC\n>seqD\nMDD\n"
SIMPLE_ALN   = ">seqA\nMAA-\n>seqB\nMBB-\n>seqC\nMCC-\n>seqD\nMDD-\n"


def _make_raw_dirs(tmp_path, hmm="TestHMM"):
    tree_dir    = tmp_path / "trees_iqtree"
    fasta_dir   = tmp_path / "fasta_per_hmm"
    clipkit_dir = tmp_path / "alignments_clipkit"
    summary_dir = tmp_path / "summary"
    emb_dir     = tmp_path / "embeddings"

    for d in [tree_dir, fasta_dir, clipkit_dir, summary_dir, emb_dir]:
        d.mkdir(parents=True, exist_ok=True)

    (tree_dir    / f"{hmm}.treefile").write_text(SIMPLE_TREE)
    (fasta_dir   / f"{hmm}.fasta").write_text(SIMPLE_FASTA)
    (clipkit_dir / f"{hmm}.clipkit.faa").write_text(SIMPLE_ALN)

    return tree_dir, fasta_dir, clipkit_dir, emb_dir, summary_dir


# ---------------------------------------------------------------------------
# _prune_fasta
# ---------------------------------------------------------------------------

class TestPruneFasta:
    def test_prune_keeps_subset(self, tmp_path):
        from phylofoundry.tasks.curate import _prune_fasta

        in_fp  = str(tmp_path / "in.fasta")
        out_fp = str(tmp_path / "out.fasta")
        with open(in_fp, "w") as f:
            f.write(SIMPLE_FASTA)

        _prune_fasta(in_fp, out_fp, {"seqA", "seqC"})

        recs = list(SeqIO.parse(out_fp, "fasta"))
        ids = {r.id for r in recs}
        assert ids == {"seqA", "seqC"}

    def test_prune_does_not_modify_source(self, tmp_path):
        from phylofoundry.tasks.curate import _prune_fasta

        in_fp  = str(tmp_path / "in.fasta")
        out_fp = str(tmp_path / "out.fasta")
        with open(in_fp, "w") as f:
            f.write(SIMPLE_FASTA)

        original_content = open(in_fp).read()
        _prune_fasta(in_fp, out_fp, {"seqA"})
        assert open(in_fp).read() == original_content, "Source FASTA must not be modified"

    def test_prune_missing_source_is_noop(self, tmp_path):
        from phylofoundry.tasks.curate import _prune_fasta

        # Should not raise even if source doesn't exist
        _prune_fasta(str(tmp_path / "missing.fasta"), str(tmp_path / "out.fasta"), {"seqA"})


# ---------------------------------------------------------------------------
# run_curate — overlay-only behaviour
# ---------------------------------------------------------------------------

class TestRunCurate:
    def _cfg(self, use_treeshrink=False, use_esm_filter=False):
        return {
            "curate": {
                "enabled": True,
                "use_treeshrink": use_treeshrink,
                "use_esm_filter": use_esm_filter,
            }
        }

    def test_raw_files_are_not_modified(self, tmp_path):
        """Raw pipeline outputs must be immutable after curate runs."""
        from phylofoundry.tasks.curate import run_curate

        hmm = "TestHMM"
        tree_dir, fasta_dir, clipkit_dir, emb_dir, summary_dir = _make_raw_dirs(tmp_path, hmm)

        raw_tree_before  = (tree_dir    / f"{hmm}.treefile").read_text()
        raw_fasta_before = (fasta_dir   / f"{hmm}.fasta").read_text()
        raw_aln_before   = (clipkit_dir / f"{hmm}.clipkit.faa").read_text()

        curated_dir = str(tmp_path / "curated")
        run_curate(
            self._cfg(),
            str(tree_dir), str(fasta_dir), str(clipkit_dir),
            str(emb_dir), str(summary_dir), None, force=False,
            curated_dir=curated_dir,
        )

        assert (tree_dir    / f"{hmm}.treefile").read_text()    == raw_tree_before,  "Raw tree was overwritten!"
        assert (fasta_dir   / f"{hmm}.fasta").read_text()       == raw_fasta_before, "Raw FASTA was overwritten!"
        assert (clipkit_dir / f"{hmm}.clipkit.faa").read_text() == raw_aln_before,   "Raw alignment was overwritten!"

    def test_curated_overlay_is_created(self, tmp_path):
        """Curated overlay files must be written under curated/."""
        from phylofoundry.tasks.curate import run_curate

        hmm = "TestHMM"
        tree_dir, fasta_dir, clipkit_dir, emb_dir, summary_dir = _make_raw_dirs(tmp_path, hmm)

        curated_dir = str(tmp_path / "curated")
        run_curate(
            self._cfg(),
            str(tree_dir), str(fasta_dir), str(clipkit_dir),
            str(emb_dir), str(summary_dir), None, force=False,
            curated_dir=curated_dir,
        )

        assert os.path.exists(os.path.join(curated_dir, "trees", f"{hmm}.treefile"))
        assert os.path.exists(os.path.join(curated_dir, "fasta_per_hmm", f"{hmm}.fasta"))
        assert os.path.exists(os.path.join(curated_dir, "alignments_clipkit", f"{hmm}.clipkit.faa"))

    def test_manifest_is_written(self, tmp_path):
        """A per-HMM JSON manifest must be created in curated/manifests/."""
        from phylofoundry.tasks.curate import run_curate

        hmm = "TestHMM"
        tree_dir, fasta_dir, clipkit_dir, emb_dir, summary_dir = _make_raw_dirs(tmp_path, hmm)

        curated_dir = str(tmp_path / "curated")
        run_curate(
            self._cfg(),
            str(tree_dir), str(fasta_dir), str(clipkit_dir),
            str(emb_dir), str(summary_dir), None, force=False,
            curated_dir=curated_dir,
        )

        manifest_fp = os.path.join(curated_dir, "manifests", f"{hmm}.json")
        assert os.path.exists(manifest_fp), "Manifest JSON not written"

        with open(manifest_fp) as f:
            m = json.load(f)

        assert m["hmm"] == hmm
        assert "timestamp" in m
        assert m["source"] == "automated"
        assert "dropped_taxa" in m
        assert "retained_taxa_count" in m
        assert isinstance(m["was_modified"], bool)

    def test_audit_log_is_appended(self, tmp_path):
        """An audit JSONL record must be appended for each HMM."""
        from phylofoundry.tasks.curate import run_curate

        hmm = "TestHMM"
        tree_dir, fasta_dir, clipkit_dir, emb_dir, summary_dir = _make_raw_dirs(tmp_path, hmm)

        curated_dir = str(tmp_path / "curated")
        run_curate(
            self._cfg(),
            str(tree_dir), str(fasta_dir), str(clipkit_dir),
            str(emb_dir), str(summary_dir), None, force=False,
            curated_dir=curated_dir,
        )

        audit_dir = os.path.join(curated_dir, "audit")
        audit_files = [f for f in os.listdir(audit_dir) if f.endswith(".jsonl")]
        assert len(audit_files) >= 1, "No audit JSONL file found"

        audit_fp = os.path.join(audit_dir, audit_files[0])
        lines = [l for l in open(audit_fp).readlines() if l.strip()]
        assert len(lines) >= 1
        record = json.loads(lines[0])
        assert record["hmm"] == hmm

    def test_disabled_curate_returns_none(self, tmp_path):
        from phylofoundry.tasks.curate import run_curate

        hmm = "TestHMM"
        tree_dir, fasta_dir, clipkit_dir, emb_dir, summary_dir = _make_raw_dirs(tmp_path, hmm)

        cfg = {"curate": {"enabled": False}}
        result = run_curate(
            cfg,
            str(tree_dir), str(fasta_dir), str(clipkit_dir),
            str(emb_dir), str(summary_dir), None,
        )
        assert result is None

    def test_no_raw_backup_files_created(self, tmp_path):
        """No .raw.treefile / .raw.fasta legacy backup files must be created."""
        from phylofoundry.tasks.curate import run_curate

        hmm = "TestHMM"
        tree_dir, fasta_dir, clipkit_dir, emb_dir, summary_dir = _make_raw_dirs(tmp_path, hmm)

        curated_dir = str(tmp_path / "curated")
        run_curate(
            self._cfg(),
            str(tree_dir), str(fasta_dir), str(clipkit_dir),
            str(emb_dir), str(summary_dir), None, force=False,
            curated_dir=curated_dir,
        )

        # No .raw.* backup files should be created anywhere in the raw dirs
        raw_backups = (
            glob.glob(str(tree_dir    / "*.raw.treefile")) +
            glob.glob(str(fasta_dir   / "*.raw.fasta")) +
            glob.glob(str(clipkit_dir / "*.raw.clipkit.faa"))
        )
        assert raw_backups == [], f"Unexpected raw backup files: {raw_backups}"

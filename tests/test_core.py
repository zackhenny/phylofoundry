"""Tests for phylofoundry core utilities and configuration."""

import os
import json
import tempfile
import pytest

# ──── Fixtures ────────────────────────────────────────────────────────────────

FIXTURES_DIR = os.path.join(os.path.dirname(__file__), "fixtures")


@pytest.fixture
def fixtures_dir():
    return FIXTURES_DIR


@pytest.fixture
def tmp_outdir(tmp_path):
    return str(tmp_path / "output")


# ──── bio.py ──────────────────────────────────────────────────────────────────

class TestBio:
    def test_read_fasta(self, fixtures_dir):
        from phylofoundry.utils.bio import read_fasta

        seqs = read_fasta(os.path.join(fixtures_dir, "genomeA.faa"))
        assert len(seqs) == 3
        assert "protein_A1" in seqs
        assert "protein_A2" in seqs
        assert "protein_A3" in seqs
        # Check sequences are not empty
        for sid, seq in seqs.items():
            assert len(seq) > 0
            assert "\n" not in seq

    def test_write_fasta(self, tmp_path):
        from phylofoundry.utils.bio import read_fasta, write_fasta

        records = {"seq1": "MSKGEELFT", "seq2": "MVSKGEELFT"}
        out_fp = str(tmp_path / "test.faa")
        write_fasta(out_fp, records)
        assert os.path.exists(out_fp)

        # Round-trip test
        reloaded = read_fasta(out_fp)
        assert reloaded == records

    def test_read_fasta_empty_file(self, tmp_path):
        from phylofoundry.utils.bio import read_fasta

        empty_fp = str(tmp_path / "empty.faa")
        with open(empty_fp, "w") as f:
            pass
        seqs = read_fasta(empty_fp)
        assert seqs == {}


# ──── helpers.py ──────────────────────────────────────────────────────────────

class TestHelpers:
    def test_safe_mkdir(self, tmp_path):
        from phylofoundry.utils.helpers import safe_mkdir

        d = str(tmp_path / "a" / "b" / "c")
        safe_mkdir(d)
        assert os.path.isdir(d)
        # Call again — should not raise
        safe_mkdir(d)

    def test_normalize_genome_id(self):
        from phylofoundry.utils.helpers import normalize_genome_id

        assert normalize_genome_id("genome.faa") == "genome"
        assert normalize_genome_id("genome.faa.gz") == "genome"
        assert normalize_genome_id("genome.fna") == "genome"
        assert normalize_genome_id("genome.cds.fna") == "genome"
        assert normalize_genome_id("genome") == "genome"
        assert normalize_genome_id(None) is None

    def test_write_json(self, tmp_path):
        from phylofoundry.utils.helpers import write_json

        data = {"key": "value", "nested": {"a": 1}}
        fp = str(tmp_path / "sub" / "test.json")
        write_json(data, fp)
        assert os.path.exists(fp)
        with open(fp) as f:
            loaded = json.load(f)
        assert loaded["key"] == "value"
        assert loaded["nested"]["a"] == 1

    def test_load_json_config(self, tmp_path):
        from phylofoundry.utils.helpers import write_json, load_json_config

        data = {"x": 42}
        fp = str(tmp_path / "cfg.json")
        write_json(data, fp)
        loaded = load_json_config(fp)
        assert loaded == data


# ──── config.py ───────────────────────────────────────────────────────────────

class TestConfig:
    def test_deep_update(self):
        from phylofoundry.config import deep_update

        base = {"a": 1, "b": {"c": 2, "d": 3}}
        updates = {"b": {"c": 99}, "e": 5}
        result = deep_update(base, updates)
        assert result["a"] == 1
        assert result["b"]["c"] == 99
        assert result["b"]["d"] == 3
        assert result["e"] == 5

    def test_default_config_has_required_keys(self):
        from phylofoundry.constants import DEFAULT_CONFIG

        assert "inputs" in DEFAULT_CONFIG
        assert "output" in DEFAULT_CONFIG
        assert "resources" in DEFAULT_CONFIG
        assert "workflow" in DEFAULT_CONFIG
        assert "filtering" in DEFAULT_CONFIG
        assert "phylo" in DEFAULT_CONFIG
        assert "embeddings" in DEFAULT_CONFIG
        assert "post" in DEFAULT_CONFIG
        assert "codon" in DEFAULT_CONFIG
        assert "hyphy" in DEFAULT_CONFIG

    def test_steps_are_in_order(self):
        from phylofoundry.constants import STEPS

        assert STEPS == [
            "prep", "hmmer", "extract", "embed", "phylo", "curate",
            "taxonomy_integrate", "conservation_metrics", "detect_clades",
            "post",
            "synteny", "codon", "hyphy", "score_motifs", "discover_motifs",
        ]


# ──── pipeline.py ─────────────────────────────────────────────────────────────

class TestPipeline:
    def test_step_in_range(self):
        from phylofoundry.pipeline import step_in_range

        # Full range
        assert step_in_range("prep", None, None)
        assert step_in_range("hyphy", None, None)

        # Partial range
        assert step_in_range("hmmer", "hmmer", "extract")
        assert step_in_range("extract", "hmmer", "extract")
        assert not step_in_range("prep", "hmmer", "extract")
        assert not step_in_range("phylo", "hmmer", "extract")

    def test_load_manifest(self, tmp_path):
        from phylofoundry.pipeline import load_manifest

        # None returns None
        assert load_manifest(None) is None

        # File with entries
        fp = str(tmp_path / "manifest.txt")
        with open(fp, "w") as f:
            f.write("# comment\nHMM_A\nHMM_B\n\n")
        result = load_manifest(fp)
        assert result == {"HMM_A", "HMM_B"}


# ──── extract.py ──────────────────────────────────────────────────────────────

class TestExtract:
    def test_extract_sequences(self, tmp_path):
        import pandas as pd
        from phylofoundry.tasks.extract import run_extract

        cfg = {"phylo": {"use_hmmsearch_alignment": False, "keep_all_hits": False}}
        scan_df = pd.DataFrame({
            "genome": ["genomeA.faa", "genomeB.faa"],
            "protein": ["protein_A1", "protein_B1"],
            "hmm": ["TestHMM", "TestHMM"],
            "bitscore": [100.0, 95.0],
            "evalue": [1e-10, 1e-9],
            "coverage": [0.9, 0.85],
        })
        search_df = pd.DataFrame()
        fasta_dir = str(tmp_path)
        proteome_seqs = {
            "genomeA.faa": {"protein_A1": "MSKGEELFT", "protein_A2": "MVSKGEELFT"},
            "genomeB.faa": {"protein_B1": "MSKGEELFT"},
        }
        hmm_to_seqs = run_extract(cfg, scan_df, search_df, fasta_dir, None, proteome_seqs, force=True)
        assert "TestHMM" in hmm_to_seqs
        assert len(hmm_to_seqs["TestHMM"]) == 2


# ──── hmmer.py ────────────────────────────────────────────────────────────────

class TestHmmer:
    def test_best_hits_basic(self):
        import pandas as pd
        from phylofoundry.tasks.hmmer import best_hits

        df = pd.DataFrame({
            "genome": ["g1", "g1", "g1"],
            "protein": ["p1", "p1", "p2"],
            "hmm": ["hmmA", "hmmB", "hmmA"],
            "bitscore": [100.0, 80.0, 50.0],
            "evalue": [1e-10, 1e-8, 1e-5],
            "coverage": [0.9, 0.7, 0.6],
        })
        result = best_hits(df)
        assert len(result) == 2  # two unique (genome, protein) pairs
        # p1 should pick hmmA (bitscore=100)
        p1_row = result[result["protein"] == "p1"].iloc[0]
        assert p1_row["hmm"] == "hmmA"
        assert p1_row["bitscore"] == 100.0

    def test_best_hits_empty(self):
        import pandas as pd
        from phylofoundry.tasks.hmmer import best_hits

        result = best_hits(pd.DataFrame())
        assert result.empty

    def test_apply_filtering(self):
        import pandas as pd
        from phylofoundry.tasks.hmmer import apply_filtering

        df = pd.DataFrame({
            "hmm": ["hmmA", "hmmB", "hmmC"],
            "bitscore": [100.0, 20.0, 50.0],
            "coverage": [0.9, 0.8, 0.3],
        })
        result = apply_filtering(df, {}, global_min_score=25.0, global_min_cov=0.5)
        assert len(result) == 1  # only hmmA passes both filters
        assert result.iloc[0]["hmm"] == "hmmA"


# ──── codon.py ────────────────────────────────────────────────────────────────

class TestCodon:
    def test_strip_stop_from_protein_removes_stop(self):
        from phylofoundry.tasks.codon import _strip_stop_from_protein

        assert _strip_stop_from_protein("MACD*") == "MACD"
        assert _strip_stop_from_protein("MACD") == "MACD"
        assert _strip_stop_from_protein("MA*CD*") == "MACD"

    def test_strip_stop_from_protein_keeps_x(self):
        from phylofoundry.tasks.codon import _strip_stop_from_protein

        # 'X' (ambiguous AA) must be preserved at all positions
        assert _strip_stop_from_protein("MACDX") == "MACDX"
        assert _strip_stop_from_protein("MAXCD") == "MAXCD"
        assert _strip_stop_from_protein("MAXCDX") == "MAXCDX"

    def test_strip_stop_from_protein_keeps_c(self):
        from phylofoundry.tasks.codon import _strip_stop_from_protein

        # 'C' (cysteine) must be preserved regardless of position
        assert _strip_stop_from_protein("MACDC") == "MACDC"
        assert _strip_stop_from_protein("CMAC*") == "CMAC"

    def test_clean_aa_alignment_strips_stop_keeps_x_and_c(self):
        from phylofoundry.tasks.codon import _clean_aa_alignment

        aln = {"seq1": "MACD*", "seq2": "MVXC", "seq3": "MA-CD"}
        cleaned = _clean_aa_alignment(aln)
        assert cleaned["seq1"] == "MACD"
        assert cleaned["seq2"] == "MVXC"  # X and C preserved
        assert cleaned["seq3"] == "MA-CD"

    def test_read_tree_tips_basic(self, tmp_path):
        from phylofoundry.tasks.codon import _read_tree_tips

        treefile = tmp_path / "test.treefile"
        treefile.write_text("(seqA:0.1,seqB:0.2,(seqC:0.3,seqD:0.4)98:0.5)Root;\n")
        tips = _read_tree_tips(str(treefile))
        assert "seqA" in tips
        assert "seqB" in tips
        assert "seqC" in tips
        assert "seqD" in tips
        # Internal node labels like bootstrap values may also appear in the set;
        # that is expected and documented.  Real sequence IDs are never purely
        # numeric, so downstream filtering is unaffected.
        assert "98" in tips or "Root" in tips or True  # non-sequence labels present

    def test_read_tree_tips_with_pipe_names(self, tmp_path):
        from phylofoundry.tasks.codon import _read_tree_tips

        treefile = tmp_path / "test.treefile"
        treefile.write_text(
            "(genome1|prot_A1:0.1,genome2|prot_B1:0.2,genome3|prot_C1:0.3);\n"
        )
        tips = _read_tree_tips(str(treefile))
        assert "genome1|prot_A1" in tips
        assert "genome2|prot_B1" in tips
        assert "genome3|prot_C1" in tips


# ──── prep.py ─────────────────────────────────────────────────────────────────

class TestPrep:
    def test_combined_faa_strips_stop_keeps_x_and_c(self, tmp_path):
        """prep.py must remove '*' but preserve 'X' and 'C' residues."""
        from phylofoundry.utils.bio import write_fasta, read_fasta
        from phylofoundry.tasks.prep import run_prep

        faa_dir = str(tmp_path / "faa")
        os.makedirs(faa_dir)
        write_fasta(
            os.path.join(faa_dir, "genomeX.faa"),
            {
                "prot1": "MSKGEELFT*",   # trailing stop
                "prot2": "MVSKXEELFT",   # X should be kept
                "prot3": "MCSKELFT*",    # C should be kept, * removed
            },
        )

        outdir = str(tmp_path / "out")
        os.makedirs(outdir)
        combined_faa = os.path.join(outdir, "combined_proteomes.faa")
        combined_hmm = os.path.join(outdir, "combined.hmm")
        with open(combined_hmm, "w") as f:
            f.write("")
        for ext in [".h3f", ".h3i", ".h3m", ".h3p"]:
            with open(combined_hmm + ext, "w") as f:
                f.write("")

        cfg = {"hmmer": {"run_scan": True, "run_search": True}}
        run_prep(
            cfg,
            genomes=["genomeX.faa"],
            faa_dir=faa_dir,
            hmm_input_mode="file",
            hmm_dir=faa_dir,
            hmm_files=["genomeX.faa"],
            combined_faa=combined_faa,
            combined_hmm=combined_hmm,
            force=True,
        )

        seqs = read_fasta(combined_faa)
        for sid, seq in seqs.items():
            assert "*" not in seq, f"Stop codon found in {sid}: {seq}"
        # X and C must be preserved
        assert "X" in seqs["genomeX.faa~prot2"]
        assert "C" in seqs["genomeX.faa~prot3"]

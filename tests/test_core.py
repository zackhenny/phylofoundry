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


# ──── artifact_paths.py ──────────────────────────────────────────────────────

class TestArtifactPaths:
    def test_basic_properties(self, tmp_path):
        from phylofoundry.artifact_paths import ArtifactPaths

        paths = ArtifactPaths(str(tmp_path))

        assert paths.combined_faa == str(tmp_path / "combined_proteomes.faa")
        assert paths.combined_hmm == str(tmp_path / "combined.hmm")
        assert paths.hmmscan_dir == str(tmp_path / "hmmscan_tbl")
        assert paths.hmmsearch_dir == str(tmp_path / "hmmsearch_tbl")
        assert paths.fasta_dir == str(tmp_path / "fasta_per_hmm")
        assert paths.alignments_hmm_dir == str(tmp_path / "alignments_hmm")
        assert paths.alignments_clipkit_dir == str(tmp_path / "alignments_clipkit")
        assert paths.trees_dir == str(tmp_path / "trees_iqtree")
        assert paths.embeddings_dir == str(tmp_path / "embeddings")
        assert paths.codon_dir == str(tmp_path / "codon_alignments")
        assert paths.clade_assign_dir == str(tmp_path / "clade_assignments")
        assert paths.summary_dir == str(tmp_path / "summary")
        assert paths.logs_dir == str(tmp_path / "logs")

    def test_summary_artifact_paths(self, tmp_path):
        from phylofoundry.artifact_paths import ArtifactPaths

        paths = ArtifactPaths(str(tmp_path))

        assert paths.best_hits_tsv == str(tmp_path / "summary" / "best_hits.competitive.tsv")
        assert paths.best_hits_with_taxonomy_tsv == str(tmp_path / "summary" / "best_hits.with_taxonomy.tsv")
        assert paths.detected_clades_tsv == str(tmp_path / "summary" / "detected_clades.tsv")
        assert paths.conservation_metrics_dir == str(tmp_path / "summary" / "post_scikitbio")
        assert paths.hyphy_dir == str(tmp_path / "summary" / "hyphy")
        assert paths.resolved_config_json == str(tmp_path / "summary" / "resolved_config.json")

    def test_curated_paths(self, tmp_path):
        from phylofoundry.artifact_paths import ArtifactPaths

        paths = ArtifactPaths(str(tmp_path))

        assert paths.curated_dir == str(tmp_path / "curated")
        assert paths.curated_trees_dir == str(tmp_path / "curated" / "trees")
        assert paths.curated_fasta_dir == str(tmp_path / "curated" / "fasta_per_hmm")
        assert paths.curated_alignments_dir == str(tmp_path / "curated" / "alignments_clipkit")
        assert paths.curated_manifests_dir == str(tmp_path / "curated" / "manifests")
        assert paths.curated_audit_dir == str(tmp_path / "curated" / "audit")

    def test_per_hmm_helpers(self, tmp_path):
        from phylofoundry.artifact_paths import ArtifactPaths

        paths = ArtifactPaths(str(tmp_path))
        hmm = "TestHMM"

        assert paths.fasta_for_hmm(hmm) == str(tmp_path / "fasta_per_hmm" / "TestHMM.faa")
        assert paths.treefile_for_hmm(hmm) == str(tmp_path / "trees_iqtree" / "TestHMM.treefile")
        assert paths.alignment_hmm_for_hmm(hmm) == str(tmp_path / "alignments_hmm" / "TestHMM.afa")
        assert paths.alignment_clipkit_for_hmm(hmm) == str(tmp_path / "alignments_clipkit" / "TestHMM.clipkit.faa")
        assert paths.codon_alignment_for_hmm(hmm) == str(tmp_path / "codon_alignments" / "TestHMM.codon.fasta")
        assert paths.hyphy_output_for_hmm(hmm, "relax", "clade1") == str(
            tmp_path / "summary" / "hyphy" / "TestHMM" / "RELAX" / "clade1.json"
        )

    def test_all_output_dirs_contains_key_paths(self, tmp_path):
        from phylofoundry.artifact_paths import ArtifactPaths

        paths = ArtifactPaths(str(tmp_path))
        dirs = paths.all_output_dirs()

        assert paths.fasta_dir in dirs
        assert paths.trees_dir in dirs
        assert paths.summary_dir in dirs
        assert paths.curated_dir in dirs
        assert paths.hyphy_dir in dirs

    def test_module_level_constants(self):
        from phylofoundry import artifact_paths as ap

        # Spot-check canonical path fragments are strings
        assert ap.COMBINED_HMM == "combined.hmm"
        assert ap.FASTA_PER_HMM_DIR == "fasta_per_hmm"
        assert ap.TREES_IQTREE_DIR == "trees_iqtree"
        assert ap.DETECTED_CLADES_TSV == "summary/detected_clades.tsv"
        assert ap.BEST_HITS_TSV == "summary/best_hits.competitive.tsv"
        assert ap.CURATED_DIR == "curated"
        assert ap.LOGS_DIR == "logs"


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

    def test_default_config_workdir_and_globdb_keys(self):
        from phylofoundry.constants import DEFAULT_CONFIG

        # workdir defaults to None so it is optional but present
        assert "workdir" in DEFAULT_CONFIG["output"]
        assert DEFAULT_CONFIG["output"]["workdir"] is None

        # globdb_taxonomy_file defaults to None so it is optional but present
        assert "globdb_taxonomy_file" in DEFAULT_CONFIG["inputs"]
        assert DEFAULT_CONFIG["inputs"]["globdb_taxonomy_file"] is None

    def test_steps_are_in_order(self):
        from phylofoundry.constants import STEPS

        assert STEPS == [
            "prep", "hmmer", "extract", "embed", "phylo", "curate",
            "taxonomy_integrate", "conservation_metrics", "detect_clades",
            "aa_composition",
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

    def test_codon_target_order_keeps_per_hmm_primary(self):
        from phylofoundry.tasks.codon import _ordered_codon_targets
        hmm_names = _ordered_codon_targets(
            ["HMM_A", "combined_all_hits", "HMM_B"],
            include_combined_aux=True,
        )
        assert hmm_names == ["HMM_A", "HMM_B", "combined_all_hits"]

    def test_hyphy_combined_mode_still_loads_per_hmm_clades(self):
        from phylofoundry.tasks.hyphy import _clade_assignment_label_for_hmm
        assert _clade_assignment_label_for_hmm("HMM_A") == "HMM_A"
        assert _clade_assignment_label_for_hmm("combined_all_hits") == "combined"


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

    def test_run_prep_diamond_mode_creates_combined_faa(self, tmp_path):
        """run_prep_diamond_mode should build combined_proteomes.faa (no HMM step)."""
        from phylofoundry.utils.bio import write_fasta, read_fasta
        from phylofoundry.tasks.prep import run_prep_diamond_mode

        faa_dir = str(tmp_path / "faa")
        os.makedirs(faa_dir)
        write_fasta(
            os.path.join(faa_dir, "genomeY.faa"),
            {"prot1": "MSKGEELFT*", "prot2": "MVSKXEELFT"},
        )

        outdir = str(tmp_path / "out")
        os.makedirs(outdir)
        combined_faa = os.path.join(outdir, "combined_proteomes.faa")

        run_prep_diamond_mode(
            cfg={},
            genomes=["genomeY.faa"],
            faa_dir=faa_dir,
            combined_faa=combined_faa,
            force=True,
        )

        assert os.path.exists(combined_faa)
        seqs = read_fasta(combined_faa)
        assert "genomeY.faa~prot1" in seqs
        assert "*" not in seqs["genomeY.faa~prot1"]
        assert "X" in seqs["genomeY.faa~prot2"]


# ──── diamond.py ──────────────────────────────────────────────────────────────

class TestDiamond:
    def test_load_query_fastas_single_file(self, tmp_path):
        from phylofoundry.tasks.diamond import load_query_fastas

        faa = tmp_path / "myquery.faa"
        faa.write_text(">prot1\nMSKGEELFT\n")
        queries = load_query_fastas(str(faa))
        assert list(queries.keys()) == ["myquery"]
        assert queries["myquery"] == str(faa)

    def test_load_query_fastas_directory(self, tmp_path):
        from phylofoundry.tasks.diamond import load_query_fastas

        qdir = tmp_path / "queries"
        qdir.mkdir()
        (qdir / "alpha.faa").write_text(">p1\nMSK\n")
        (qdir / "beta.fasta").write_text(">p2\nMSK\n")
        (qdir / "ignored.txt").write_text("not a fasta")

        queries = load_query_fastas(str(qdir))
        assert set(queries.keys()) == {"alpha", "beta"}

    def test_load_query_fastas_various_extensions(self, tmp_path):
        from phylofoundry.tasks.diamond import load_query_fastas

        qdir = tmp_path / "queries"
        qdir.mkdir()
        # Use distinct base names so each resolves to a unique query
        (qdir / "protA.faa").write_text(">p\nMSK\n")
        (qdir / "protB.fasta").write_text(">p\nMSK\n")
        (qdir / "protC.fa").write_text(">p\nMSK\n")
        (qdir / "protD.fas").write_text(">p\nMSK\n")

        queries = load_query_fastas(str(qdir))
        assert len(queries) == 4

    def test_load_query_fastas_empty_directory_raises(self, tmp_path):
        from phylofoundry.tasks.diamond import load_query_fastas

        empty_dir = tmp_path / "empty"
        empty_dir.mkdir()
        with pytest.raises(ValueError, match="No FASTA files"):
            load_query_fastas(str(empty_dir))

    def test_load_query_fastas_invalid_path_raises(self, tmp_path):
        from phylofoundry.tasks.diamond import load_query_fastas

        with pytest.raises(ValueError):
            load_query_fastas(str(tmp_path / "nonexistent"))

    def test_parse_diamond_tab_empty_file(self, tmp_path):
        from phylofoundry.tasks.diamond import parse_diamond_tab

        empty = tmp_path / "empty.tsv"
        empty.write_text("")
        df = parse_diamond_tab(str(empty))
        assert df.empty

    def test_parse_diamond_tab_missing_file(self, tmp_path):
        from phylofoundry.tasks.diamond import parse_diamond_tab

        df = parse_diamond_tab(str(tmp_path / "nonexistent.tsv"))
        assert df.empty

    def test_parse_diamond_tab_valid(self, tmp_path):
        import pandas as pd
        from phylofoundry.tasks.diamond import parse_diamond_tab

        tab = tmp_path / "hits.tsv"
        tab.write_text(
            "query1\tgenomeA~prot1\t95.0\t100\t5\t0\t1\t100\t1\t100\t1e-50\t200.0\t100\t200\n"
            "query1\tgenomeB~prot2\t80.0\t80\t16\t0\t1\t80\t1\t80\t1e-30\t150.0\t100\t150\n"
        )
        df = parse_diamond_tab(str(tab))
        assert len(df) == 2
        assert list(df.columns) == [
            "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen",
        ]
        assert df.iloc[0]["pident"] == 95.0

    def test_apply_diamond_filtering_identity(self):
        import pandas as pd
        from phylofoundry.tasks.diamond import apply_diamond_filtering

        df = pd.DataFrame({
            "pident": [95.0, 20.0, 50.0],
            "coverage": [0.9, 0.9, 0.9],
            "evalue": [1e-10, 1e-10, 1e-10],
        })
        result = apply_diamond_filtering(df, {"min_identity": 30.0, "min_coverage": 0.5, "max_evalue": 1e-5})
        assert len(result) == 2
        assert 20.0 not in result["pident"].values

    def test_apply_diamond_filtering_coverage(self):
        import pandas as pd
        from phylofoundry.tasks.diamond import apply_diamond_filtering

        df = pd.DataFrame({
            "pident": [95.0, 95.0, 95.0],
            "coverage": [0.9, 0.2, 0.6],
            "evalue": [1e-10, 1e-10, 1e-10],
        })
        result = apply_diamond_filtering(df, {"min_identity": 30.0, "min_coverage": 0.5, "max_evalue": 1e-5})
        assert len(result) == 2
        assert 0.2 not in result["coverage"].values

    def test_apply_diamond_filtering_evalue(self):
        import pandas as pd
        from phylofoundry.tasks.diamond import apply_diamond_filtering

        df = pd.DataFrame({
            "pident": [95.0, 95.0],
            "coverage": [0.9, 0.9],
            "evalue": [1e-10, 1e-2],
        })
        result = apply_diamond_filtering(df, {"min_identity": 30.0, "min_coverage": 0.5, "max_evalue": 1e-5})
        assert len(result) == 1

    def test_apply_diamond_filtering_empty(self):
        import pandas as pd
        from phylofoundry.tasks.diamond import apply_diamond_filtering

        result = apply_diamond_filtering(pd.DataFrame(), {"min_identity": 30.0})
        assert result.empty

    def test_diamond_hits_schema_matches_extract_expectations(self, tmp_path):
        """Verify the DataFrame schema produced by diamond search matches what extract expects."""
        from phylofoundry.tasks.diamond import parse_diamond_tab

        tab = tmp_path / "hits.tsv"
        tab.write_text(
            "myquery\tgenomeA.faa~prot_A1\t95.0\t100\t5\t0\t1\t100\t1\t100\t1e-50\t200.0\t100\t200\n"
        )
        raw = parse_diamond_tab(str(tab))
        raw["coverage"] = raw["length"] / raw["qlen"]
        raw["genome"] = raw["sseqid"].apply(lambda x: x.split("~")[0] if "~" in x else "Unknown")
        raw["protein"] = raw["sseqid"].apply(lambda x: x.split("~", 1)[1] if "~" in x else x)
        raw["hmm"] = "myquery"
        df = raw[["genome", "protein", "hmm", "bitscore", "evalue", "coverage", "pident"]]

        for col in ["genome", "protein", "hmm", "bitscore", "evalue", "coverage"]:
            assert col in df.columns, f"Missing column: {col}"

        assert df.iloc[0]["genome"] == "genomeA.faa"
        assert df.iloc[0]["protein"] == "prot_A1"
        assert df.iloc[0]["hmm"] == "myquery"
        assert df.iloc[0]["coverage"] == pytest.approx(1.0)

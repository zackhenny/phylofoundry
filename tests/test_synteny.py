import sys
import os
import pytest
from unittest.mock import MagicMock, patch

# Mock pygenomeviz before importing synteny
sys.modules["pygenomeviz"] = MagicMock()
sys.modules["pygenomeviz.GenomeViz"] = MagicMock()

from phylofoundry.tasks import synteny

# ──── Fixtures ────────────────────────────────────────────────────────────────

@pytest.fixture
def dummy_gff_content():
    return """##gff-version 3
contig1	Prodigal	CDS	100	200	.	+	0	ID=prot1;locus_tag=LOC1;product=hypothetical
contig1	Prodigal	CDS	300	400	.	-	0	ID=prot2;locus_tag=LOC2;product=target_protein
contig1	Prodigal	CDS	500	600	.	+	0	ID=prot3;locus_tag=LOC3;product=downstream
"""

@pytest.fixture
def dummy_gbk_content():
    # Minimal GBK for Bio.SeqIO
    return """LOCUS       contig1                  600 bp    DNA     linear   UNK 01-JAN-1980
FEATURES             Location/Qualifiers
     CDS             100..200
                     /locus_tag="LOC1"
                     /protein_id="prot1"
                     /translation="MTEST"
     CDS             complement(300..400)
                     /locus_tag="LOC2"
                     /protein_id="prot2"
                     /translation="MTARGET"
     CDS             500..600
                     /locus_tag="LOC3"
                     /protein_id="prot3"
                     /translation="MDOWN"
ORIGIN
        1 aaaaaa
//
"""

# ──── Tests ───────────────────────────────────────────────────────────────────

class TestSyntenyParsing:
    def test_parse_gff_line(self):
        line = "contig1\tProdigal\tCDS\t100\t200\t.\t+\t0\tID=prot1;locus_tag=LOC1"
        res = synteny.parse_gff_line(line)
        assert res["contig"] == "contig1"
        assert res["start"] == 100
        assert res["end"] == 200
        assert res["strand"] == "+"
        assert res["attrs"]["ID"] == "prot1"

    def test_load_gff_neighborhood(self, tmp_path, dummy_gff_content):
        gff_fp = tmp_path / "test.gff"
        gff_fp.write_text(dummy_gff_content)
        
        # Target LOC2 (prot2)
        # Window=1 should get LOC1, LOC2, LOC3
        feats = synteny.load_gff_neighborhood(
            str(gff_fp), "prot2", window_genes=1,
            protein_id_fields=["ID"], gene_label_fields=["locus_tag"]
        )
        assert len(feats) == 3
        # Middle one should be focal
        assert feats[0]["label"] == "LOC1"
        assert not feats[0]["is_focal"]
        assert feats[1]["label"] == "LOC2"
        assert feats[1]["is_focal"]
        assert feats[2]["label"] == "LOC3"

    @patch("phylofoundry.tasks.synteny.SeqIO.parse")
    def test_load_genbank_neighborhood(self, mock_parse, dummy_gbk_content):
        # We need to mock SeqIO.parse to return a Record from string
        from Bio import SeqIO
        from io import StringIO
        
        # Actually standard Bio.SeqIO.parse works on file handles, so let's write to file
        # and not mock inner parsing if we can avoid it. But Gbk parsing is complex.
        # Let's rely on standard Bio installed behavior if possible.
        # If Bio is not installed, this test will fail, which is fine (module needs Bio).
        
        # Re-write passing a real file
        pass 
        
    def test_load_genbank_real_bio(self, tmp_path, dummy_gbk_content):
        gbk_fp = tmp_path / "test.gbk"
        gbk_fp.write_text(dummy_gbk_content)
        
        try:
            feats = synteny.load_genbank_neighborhood(
                str(gbk_fp), "prot2", 1, ["protein_id"], ["locus_tag"]
            )
        except ImportError:
            pytest.skip("Biopython not installed")
            return

        assert len(feats) == 3
        assert feats[1]["protein_id"] == "prot2"
        assert feats[1]["is_focal"]
        assert feats[1]["strand"] == -1 # complement

class TestSyntenyRun:
    @patch("phylofoundry.tasks.synteny.run_cmd")
    def test_run_synteny_flow(self, mock_run, tmp_path):
        import pandas as pd
        
        # We need to ensure the lazy import picks up our mock
        # sys.modules is already mocked at top of file
        # But run_synteny does: from pygenomeviz import GenomeViz
        # So we need to ensure pygenomeviz.GenomeViz is our mock class
        mock_gv_cls = sys.modules["pygenomeviz"].GenomeViz

        cfg = {
            "output": {"outdir": str(tmp_path)},
            "resources": {"cpu": 1},
            "synteny": {
                "enabled": True,
                "gbk_dir": None, # Skip parsing for this flow test
                "gff_dir": None,
                "window_genes": 1
            }
        }
        
        # Make a dummy hits file
        summary_dir = tmp_path / "summary"
        summary_dir.mkdir(parents=True)
        hits_tsv = summary_dir / "best_hits.competitive.tsv"
        pd.DataFrame({
            "genome": ["g1"], "protein": ["p1"], "hmm": ["HMM1"], "bitscore": [100.0]
        }).to_csv(hits_tsv, sep="\t", index=False)
        
        # Test 1: No gbk/gff provided -> should log warning and return
        synteny.run_synteny(cfg, str(tmp_path/"synteny"), str(tmp_path/"trees"), pd.DataFrame(), pd.DataFrame(), None)
        # Should not have tried to plot because no input dirs
        mock_gv_cls.assert_not_called()


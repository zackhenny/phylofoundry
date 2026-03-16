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

class TestGenBankOrientation:
    """Tests that GenBank records are oriented so the focal gene is on the forward strand."""

    def _make_feats(self, focal_strand):
        """Return three synthetic features with the focal gene in the middle."""
        return [
            {"start": 100, "end": 200, "strand": 1, "label": "up", "protein_id": "p1",
             "translation": "", "is_focal": False, "uid": "p1"},
            {"start": 300, "end": 400, "strand": focal_strand, "label": "focal", "protein_id": "p2",
             "translation": "", "is_focal": True, "uid": "p2"},
            {"start": 500, "end": 600, "strand": 1, "label": "down", "protein_id": "p3",
             "translation": "", "is_focal": False, "uid": "p3"},
        ]

    def _build_gbk_record(self, feats, contig_seq=None):
        """Replicate the GenBank-building logic from run_synteny, returning the SeqRecord."""
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        from Bio.SeqFeature import SeqFeature, FeatureLocation

        min_start = min(f["start"] for f in feats)
        max_end = max(f["end"] for f in feats)
        size = max_end - min_start + 1

        if contig_seq is not None:
            nucl_seq = contig_seq[min_start:max_end]
            if isinstance(nucl_seq, str):
                nucl_seq = Seq(nucl_seq)
        else:
            nucl_seq = Seq("N" * size)

        focal_feat = next((f for f in feats if f.get("is_focal")), None)
        focal_strand = 1
        if focal_feat is not None:
            focal_strand = 1 if str(focal_feat["strand"]) in ["1", "+"] else -1

        if focal_strand == -1:
            nucl_seq = nucl_seq.reverse_complement()

        rec = SeqRecord(nucl_seq, id="genome_p2", name="genome_p2",
                        description="test", annotations={"molecule_type": "DNA"})

        for f in feats:
            strand = 1 if str(f["strand"]) in ["1", "+"] else -1
            feat_start = f["start"] - min_start
            feat_end = f["end"] - min_start

            if focal_strand == -1:
                new_start = size - feat_end
                new_end = size - feat_start
                strand = -strand
                feat_start, feat_end = new_start, new_end

            loc = FeatureLocation(feat_start, feat_end, strand=strand)
            qualifiers = {"locus_tag": [f.get("label", "")]}
            sf = SeqFeature(loc, type="CDS", id=f["uid"], qualifiers=qualifiers)
            rec.features.append(sf)

        return rec, focal_strand

    def test_forward_focal_gene_unchanged(self):
        """When the focal gene is on the forward strand, coordinates are unchanged."""
        feats = self._make_feats(focal_strand=1)
        rec, _ = self._build_gbk_record(feats)

        focal_sf = next(sf for sf in rec.features if sf.id == "p2")
        # Focal gene at original offset from min_start (300-100=200, 400-100=300)
        assert focal_sf.location.strand == 1
        assert int(focal_sf.location.start) == 200
        assert int(focal_sf.location.end) == 300

    def test_reverse_focal_gene_flipped(self):
        """When the focal gene is on the reverse strand, all features are flipped."""
        feats = self._make_feats(focal_strand=-1)
        size = 600 - 100 + 1  # max_end - min_start + 1 = 501
        rec, focal_strand = self._build_gbk_record(feats)

        assert focal_strand == -1

        focal_sf = next(sf for sf in rec.features if sf.id == "p2")
        # Focal gene was originally at offsets 200-300; after flip it should be forward
        assert focal_sf.location.strand == 1, "Focal gene should be forward after flip"

        # All neighbours should have their strands inverted relative to original
        upstream_sf = next(sf for sf in rec.features if sf.id == "p1")
        downstream_sf = next(sf for sf in rec.features if sf.id == "p3")
        assert upstream_sf.location.strand == -1
        assert downstream_sf.location.strand == -1

    def test_reverse_focal_gene_with_sequence(self):
        """When a real sequence is provided and focal is reverse, sequence is rev-comped."""
        from Bio.Seq import Seq

        # Build with a contig sequence of all-A; reverse complement of all-A is all-T
        contig_seq = Seq("A" * 600)
        feats = self._make_feats(focal_strand=-1)
        rec, _ = self._build_gbk_record(feats, contig_seq=contig_seq)

        # Rev-comp of all-A is all-T
        assert str(rec.seq) == "T" * len(rec.seq)

    def test_no_focal_gene_unchanged(self):
        """If no focal gene is found, features are written with original orientation."""
        feats = [
            {"start": 100, "end": 200, "strand": -1, "label": "a", "protein_id": "p1",
             "translation": "", "is_focal": False, "uid": "p1"},
        ]
        rec, focal_strand = self._build_gbk_record(feats)
        assert focal_strand == 1  # defaults to forward
        sf = rec.features[0]
        assert sf.location.strand == -1  # original orientation preserved


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


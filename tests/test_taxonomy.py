import os
import pandas as pd
import pytest
from unittest.mock import MagicMock, patch
from phylofoundry.tasks import post

def test_load_taxonomy_gtdb(tmp_path):
    # Setup GTDB dir
    gtdb_dir = tmp_path / "gtdb_out"
    gtdb_dir.mkdir()
    
    # Create dummy summary
    summary_fp = gtdb_dir / "gtdbtk.bac120.summary.tsv"
    with open(summary_fp, "w") as f:
        f.write("user_genome\tclassification\n")
        f.write("genomeA\td_Bacteria;p_Proteobacteria;...\n")
        f.write("genomeB.faa\td_Archaea;p_Euryarchaeota;...\n") # Check if normalization handles .faa if user put it there

    # Call helper
    # We need to access the private helper or test via run_post
    # Accessing via post._load_taxonomy works if it's there
    tax_map = post._load_taxonomy(str(gtdb_dir), None)
    
    assert "genomeA" in tax_map
    assert tax_map["genomeA"].startswith("d_Bacteria")
    
    # Check normalization logic: "genomeB.faa" -> "genomeB"
    assert "genomeB" in tax_map
    assert tax_map["genomeB"].startswith("d_Archaea")

def test_load_taxonomy_custom_tax(tmp_path):
    tax_fp = tmp_path / "custom.tsv"
    with open(tax_fp, "w") as f:
        f.write("genome\tlineage\n")
        f.write("genomeC\tRoot;Life\n")
        
    tax_map = post._load_taxonomy(None, str(tax_fp))
    assert tax_map["genomeC"] == "Root;Life"

def test_run_post_merges_taxonomy(tmp_path):
    # Setup directories
    outdir = tmp_path / "results"
    summary_dir = outdir / "summary"
    post_dir = summary_dir / "post_scikitbio"
    tree_dir = outdir / "trees"
    os.makedirs(summary_dir)
    os.makedirs(post_dir)
    os.makedirs(tree_dir)
    
    # Create best_hits
    best_hits_fp = summary_dir / "best_hits.competitive.tsv"
    pd.DataFrame({
        "genome": ["genomeA.faa", "genomeB", "genomeC.fna"],
        "protein": ["p1", "p2", "p3"]
    }).to_csv(best_hits_fp, sep="\t", index=False)
    
    # Create GTDB output
    gtdb_dir = tmp_path / "gtdb_out"
    gtdb_dir.mkdir()
    with open(gtdb_dir / "gtdbtk.bac120.summary.tsv", "w") as f:
        f.write("user_genome\tclassification\n")
        f.write("genomeA\tTaxA\n")
        f.write("genomeB\tTaxB\n")
        f.write("genomeC\tTaxC\n")

    cfg = {
        "inputs": {"gtdb_dir": str(gtdb_dir)},
        "post": {"enabled": True}
    }
    
    # Run
    # We need to prevent it from running "real" post steps (scikit-bio).
    # Either mock them or ensure no trees exist so loop skips.
    # The taxonomy part happens BEFORE the per-HMM loop in my implementation?
    # No, looking at my `post.py` edit...
    # I put it at the very top of `run_post`.
    # So it should run even if hmm_keep is empty or no trees.
    
    post.run_post(cfg, str(tree_dir), "", "", str(post_dir), str(summary_dir), None, force=True)
    
    # Check output
    out_fp = summary_dir / "best_hits.with_taxonomy.tsv"
    assert os.path.exists(out_fp)
    
    df = pd.read_csv(out_fp, sep="\t")
    assert "taxonomy" in df.columns
    
    # Verify mapping
    # genomeA.faa -> normalize -> genomeA -> TaxA
    rowA = df[df["genome"] == "genomeA.faa"].iloc[0]
    assert rowA["taxonomy"] == "TaxA"
    
    # genomeC.fna -> normalize -> genomeC -> TaxC
    rowC = df[df["genome"] == "genomeC.fna"].iloc[0]
    assert rowC["taxonomy"] == "TaxC"


def test_detect_taxonomy_clades(tmp_path):
    summary_dir = tmp_path / "summary"
    summary_dir.mkdir()
    pd.DataFrame({
        "genome": ["genomeA.faa", "genomeB.faa"],
        "protein": ["p1", "p2"],
    }).to_csv(summary_dir / "best_hits.competitive.tsv", sep="\t", index=False)
    tax_map = {
        "genomeA": "d__Bacteria;p__Firmicutes;g__Bacillus",
        "genomeB": "d__Bacteria;p__Firmicutes;g__Bacillus",
    }

    clades = post._detect_taxonomy_clades(str(summary_dir), tax_map, "genus")
    assert "g__Bacillus" in clades
    assert sorted(clades["g__Bacillus"]) == ["genomeA.faa|p1", "genomeB.faa|p2"]


def test_compute_kl_defaults_to_clade_vs_others(tmp_path, monkeypatch):
    outdir = tmp_path / "results"
    tree_dir = outdir / "trees"
    clipkit_dir = outdir / "clipkit"
    summary_dir = outdir / "summary"
    post_dir = summary_dir / "post_scikitbio"
    os.makedirs(tree_dir)
    os.makedirs(clipkit_dir)
    os.makedirs(summary_dir)
    os.makedirs(post_dir)

    with open(tree_dir / "HMM1.treefile", "w") as f:
        f.write("(g1|p1:0.1,g2|p2:0.1,g3|p3:0.1);")
    with open(clipkit_dir / "HMM1.clipkit.faa", "w") as f:
        f.write(">g1|p1\nAAAA\n>g2|p2\nAAAA\n>g3|p3\nAAAA\n")
    clades_fp = tmp_path / "clades.tsv"
    with open(clades_fp, "w") as f:
        f.write("clade_name\ttip\n")
        f.write("CladeA\tg1|p1\n")
        f.write("CladeB\tg2|p2\n")

    def fake_counts(_aln_seqs, subset_tips=None):
        n = len(_aln_seqs) if subset_tips is None else len(subset_tips)
        return [{"A": n}]

    monkeypatch.setattr(post, "site_counts_from_subset", fake_counts)
    monkeypatch.setattr(post, "kl_divergence", lambda *_args, **_kwargs: 0.0)

    cfg = {
        "inputs": {"gtdb_dir": None, "taxonomy_file": None},
        "post": {"enabled": True, "compute_kl": True, "clades_tsv": str(clades_fp), "kl_pairs": None}
    }
    post.run_post(cfg, str(tree_dir), str(clipkit_dir), "", str(post_dir), str(summary_dir), None, force=True)

    kl_fp = post_dir / "kl_divergence.tsv"
    assert os.path.exists(kl_fp)
    kl_df = pd.read_csv(kl_fp, sep="\t")
    assert set(kl_df["pair"]) == {"CladeA:others::CladeA", "CladeB:others::CladeB"}
    assert sorted(kl_df["nB"].tolist()) == [2, 2]


def test_compute_kl_uses_detected_clades_when_clades_tsv_missing(tmp_path, monkeypatch):
    outdir = tmp_path / "results"
    tree_dir = outdir / "trees"
    clipkit_dir = outdir / "clipkit"
    summary_dir = outdir / "summary"
    post_dir = summary_dir / "post_scikitbio"
    os.makedirs(tree_dir)
    os.makedirs(clipkit_dir)
    os.makedirs(summary_dir)
    os.makedirs(post_dir)

    with open(tree_dir / "HMM1.treefile", "w") as f:
        f.write("(g1|p1:0.1,g2|p2:0.1,g3|p3:0.1);")
    with open(clipkit_dir / "HMM1.clipkit.faa", "w") as f:
        f.write(">g1|p1\nAAAA\n>g2|p2\nAAAA\n>g3|p3\nAAAA\n")

    with open(summary_dir / "detected_clades.tsv", "w") as f:
        f.write("clade_name\ttip\n")
        f.write("CladeA\tg1|p1\n")
        f.write("CladeB\tg2|p2\n")

    def fake_counts(_aln_seqs, subset_tips=None):
        n = len(_aln_seqs) if subset_tips is None else len(subset_tips)
        return [{"A": n}]

    monkeypatch.setattr(post, "site_counts_from_subset", fake_counts)
    monkeypatch.setattr(post, "kl_divergence", lambda *_args, **_kwargs: 0.0)

    cfg = {
        "inputs": {"gtdb_dir": None, "taxonomy_file": None},
        "post": {"enabled": True, "compute_kl": True, "clades_tsv": None, "kl_pairs": None}
    }
    post.run_post(cfg, str(tree_dir), str(clipkit_dir), "", str(post_dir), str(summary_dir), None, force=True)

    kl_fp = post_dir / "kl_divergence.tsv"
    assert os.path.exists(kl_fp)

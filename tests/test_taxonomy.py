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

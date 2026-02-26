from pathlib import Path

# Paths
ROOT_DIR = Path(__file__).parent
DEFAULT_CONFIG_FILE = "config.json"

# Bio
AA_ALPHABET = set(list("ACDEFGHIKLMNPQRSTVWY"))
GAP_CHARS = set(["-", ".", "X", "x", "?", "*"])

# Workflow
STEPS = ["prep", "hmmer", "extract", "embed", "phylo", "curate", "post",
         "synteny", "codon", "hyphy", "score_motifs", "discover_motifs"]

# Defaults
DEFAULT_CONFIG = {
    "inputs": {
        "faa_dir": None,     # directory of *.faa OR a single .faa
        "hmm_input": None,   # directory of *.hmm OR a single .hmm
        "cds_dir": None,      # optional directory of CDS nucleotide FASTAs (per genome)
        "gtdb_dir": None,     # optional directory of GTA-Tk output
        "taxonomy_file": None, # optional TSV with genome->lineage mapping
        "pfam_dir": None      # optional directory of Pfam .dat and .hmm files
    },
    "output": {
        "outdir": None
    },
    "resources": {
        "cpu": 8
    },
    "workflow": {
        "start_at": None,
        "stop_after": None,
        "force": False,
        "hmm_manifest": None
    },
    "hmmer": {
        "run_scan": True,
        "run_search": True
    },
    "filtering": {
        "scores_tsv": None,
        "global_min_score": 25.0,
        "min_coverage": 0.5,
        "keep_tbl": False
    },
    "phylo": {
        "mafft": False,
        "also_mafft": False,
        "mafft_for_tree": False,
        "mafft_mode": "auto", # auto, ginsi, linsi, einsi, fftns, fftnsi
        "no_trim_hmmalign": False,
        "skip_clipkit": False,
        "no_asr": False,
        "iq_boot": 1000,
        "use_hmmsearch_alignment": False,
        "keep_all_hits": False,
        "combined_tree": False, # If true, combine all hits into a single tree
        "iqtree_bin": "iqtree" # Default to iqtree, falling back to iqtree2 or iqtree3 if needed
    },
    "embeddings": {
        "enabled": False,
        "backend": "esm",            # "esm" or "transformers"
        "model": "esm2_t33_650M_UR50D",  # for backend=esm; or HF model id/path for transformers
        "device": "cuda",            # "cuda" or "cpu"
        "batch_size": 8,
        "repr_layer": None,          # if None, choose last layer (esm) or last_hidden_state (transformers)
        "pooling": "mean",           # "mean" (implemented)
        "pca_components": 3,
        "write_full_vectors": False,  # if True, write TSV with all dims (can be huge); always writes .npy
        "cluster_embeddings": True,   # run HDBSCAN on raw embeddings
        "hdbscan_min_cluster_size": 5  # HDBSCAN min_cluster_size param
    },
    "curate": {
        "enabled": False,
        "use_treeshrink": True,
        "use_esm_filter": True
    },
    "post": {
        "enabled": False,
        "compute_conservation": False,
        "conservation_metric": "inverse_shannon_uncertainty",
        "compute_kl": False,
        "clades_tsv": None,  # TSV columns: clade_name, tip (tip label must match alignment tip labels)
        "kl_pairs": None     # "A:B,A:background"
    },
    "synteny": {
        "enabled": False,
        "gbk_dir": None,
        "gff_dir": None,
        "genome_fasta_dir": None,
        "window_genes": 10,
        "max_hits_per_hmm": 50,
        "dedup_by_genome": True,
        "prefer_best_hit": True,
        "protein_id_field": ["ID", "protein_id", "locus_tag"],
        "gene_label_field": ["gene", "product", "Name", "locus_tag"],
        "plot_width": 14,
        "include_tree": True,
    },
    "codon": {
        "enabled": False,
        "build_codon_alignments": False,
        "cds_id_mode": "after_last_pipe",  # "same"|"strip_pipe"|"after_last_pipe"
        "pal2nal_cmd": "pal2nal.pl"
    },
    "hyphy": {
        "enabled": False,
        "run_hyphy": False,
        "hyphy_bin": "hyphy",
        "hyphy_tests": "RELAX,aBSREL,MEME",
        "hyphy_args": {
            "MEME": ["--branches", "All"],
            "aBSREL": ["--branches", "All"],
            "RELAX": ["--test", "Test"]
        }
    },
    "motifs": {
        "enabled": False,
        "motif_list": [],           # e.g. ["HPEVY", "HPEVF"]
        "attention_layers": 4,      # last N ESM-2 layers to average
    },
    "discover": {
        "enabled": False,
        "kmer_size": 5,
        "top_n_peaks": 20,
        "attention_layers": 4,
    },
}

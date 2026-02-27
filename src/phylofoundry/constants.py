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
    "prep": {
        "cleanup_combined_faa": False
    },
    "filtering": {
        "scores_tsv": None,
        "global_min_score": 25.0,
        "min_coverage": 0.5,
        "keep_tbl": False,
        "max_evalue": None,
        "use_evalue": False,
        "disable_bitscore_filter": False,
        "disable_coverage_filter": False
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
        "kl_pairs": None,    # "A:B,A:background"
        "detect_clades_method": None,  # None|"taxonomy"|"treecluster"|"tree_embed"
        "taxonomy_clade_level": "genus",
        "treecluster_threshold": 0.045,
        "treecluster_method": "max_clade",
        "embedtree_support_min": 80,
        "embedtree_min_size": 5,
        "embedtree_max_size": 5000,
        "embedtree_top_k": 10,
        "embedtree_pcs": 10,
        "embedtree_distance": "euclidean",  # "euclidean"|"cosine"
        "embedtree_allow_nested": False,
        "embedtree_require_monophyly": True,
        "embedtree_emit_all": False
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
        "use_detected_clades": True,
        "min_clade_size": 4,
        "label_mode": "crown",  # "crown"|"stem"
        "relax_label_reference": True,
        "hyphy_args": {
            "MEME": ["--branches", "All"],
            "aBSREL": ["--branches", "All"],
            "RELAX": ["--test", "test", "--reference", "reference"]
        }
    },
    "ha": {
        "enabled": False,
        "layer_mode": "middle",  # "middle"|"range"
        "layer_start": None,
        "layer_end": None,
        "agg": "median",  # "median"|"mean"
        "call_mode": "percentile",  # "percentile"|"topk"
        "percentile": 0.05,
        "topk": 20,
        "min_sites": 8,
        "max_sites": 60
    },
    "motifs": {
        "enabled": False,
        "motif_list": [],           # e.g. ["HPEVY", "HPEVF"]
        "attention_layers": 4,      # last N ESM-2 layers to average
        "use_ha": False
    },
    "discover": {
        "enabled": False,
        "standard_clade": None,   # optional: manually specify the reference clade ID/name
        "novel_clade": None,       # optional: manually specify the focal clade ID/name
        "kmer_size": 5,
        "top_n_peaks": 20,
        "attention_layers": 4,
        "use_ha": False,
        "ha_window": 9,
        "ha_delta_min": 0.2,
        "ha_gap_frac_max": 0.6
    },
}

import os
import glob
from collections import defaultdict
from ..utils.bio import write_fasta

def run_extract(cfg, scan_df, search_df, fasta_dir, hmm_keep, proteome_seqs, force=False):
    print("\n[extract] Building per-HMM FASTAs...")
    
    phy_cfg = cfg.get("phylo", {})
    hmm_to_seqs = defaultdict(dict)

    # choose base hits for membership
    if phy_cfg.get("use_hmmsearch_alignment", False) and not search_df.empty:
        base_df = search_df
    else:
        base_df = scan_df if not scan_df.empty else search_df

    if base_df.empty:
        print("WARNING: No hits available to extract sequences.")
        return hmm_to_seqs

    if not phy_cfg.get("keep_all_hits", False):
        base_df = base_df.sort_values("bitscore", ascending=False).groupby(["genome", "hmm"], as_index=False).head(1)

    for _, r in base_df.iterrows():
        hmm = r["hmm"]
        # manifest filtering is on hmm label as used in hits tables
        if hmm_keep is not None and hmm not in hmm_keep:
            continue
        
        genome = r["genome"]
        prot = r["protein"]
        seq = proteome_seqs.get(genome, {}).get(prot)
        if not seq:
            continue
        tip = f"{genome}|{prot}"
        hmm_to_seqs[hmm][tip] = seq

    for hmm, seqs in hmm_to_seqs.items():
        out_fp = os.path.join(fasta_dir, f"{hmm}.faa")
        if (not os.path.exists(out_fp)) or force:
            write_fasta(out_fp, seqs)
            
    return hmm_to_seqs

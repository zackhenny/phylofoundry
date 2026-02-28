import os
import shutil
from ..utils.bio import read_fasta
from ..utils.helpers import run_cmd

def run_prep(cfg, genomes, faa_dir, hmm_input_mode, hmm_dir, hmm_files, combined_faa, combined_hmm, force=False):
    print("\n[prep] Building combined proteomes FASTA and combined HMM DB...")

    hmmer_cfg = cfg.get("hmmer", {})
    run_search = hmmer_cfg.get("run_search", True)

    if run_search:
        if (not os.path.exists(combined_faa)) or force:
            with open(combined_faa, "w") as out_faa:
                for g in genomes:
                    seqs = read_fasta(os.path.join(faa_dir, g))
                    for pid, s in seqs.items():
                        # Strip stop codons ('*') from protein sequences before
                        # they enter the pipeline.  'X' and other ambiguous
                        # residues are intentionally kept to preserve original
                        # AA positions.
                        s_clean = s.replace("*", "")
                        out_faa.write(f">{g}~{pid}\n{s_clean}\n")
    else:
        print("[prep] Skipping combined_proteomes.faa (hmmsearch disabled).")

    index_exts = [".h3f", ".h3i", ".h3m", ".h3p"]

    if (not os.path.exists(combined_hmm)) or force:
        if hmm_input_mode == "file":
            src = os.path.join(hmm_dir, hmm_files[0])
            run_cmd(f"cp {src} {combined_hmm}", quiet=True, shell=True)
            # Copy sibling index files if they exist
            for ext in index_exts:
                src_idx = src + ext
                if os.path.exists(src_idx):
                    shutil.copy2(src_idx, combined_hmm + ext)
        else:
            run_cmd(f"cat {hmm_dir}/*.hmm > {combined_hmm}", quiet=True, shell=True)

    all_indices_exist = all(os.path.exists(combined_hmm + ext) for ext in index_exts)
    if not all_indices_exist:
        run_cmd(["hmmpress", "-f", combined_hmm], quiet=True)
    else:
        print("[prep] HMM indices already exist, skipping hmmpress.")

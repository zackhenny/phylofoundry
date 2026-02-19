import os
from ..utils.bio import read_fasta
from ..utils.helpers import run_cmd

def run_prep(cfg, genomes, faa_dir, hmm_input_mode, hmm_dir, hmm_files, combined_faa, combined_hmm, force=False):
    print("\n[prep] Building combined proteomes FASTA and combined HMM DB...")

    if (not os.path.exists(combined_faa)) or force:
        with open(combined_faa, "w") as out_faa:
            for g in genomes:
                seqs = read_fasta(os.path.join(faa_dir, g))
                for pid, s in seqs.items():
                    out_faa.write(f">{g}~{pid}\n{s}\n")

    if (not os.path.exists(combined_hmm)) or force:
        if hmm_input_mode == "file":
            src = os.path.join(hmm_dir, hmm_files[0])
            run_cmd(f"cp {src} {combined_hmm}", quiet=True, shell=True)
        else:
            run_cmd(f"cat {hmm_dir}/*.hmm > {combined_hmm}", quiet=True, shell=True)
        run_cmd(["hmmpress", "-f", combined_hmm], quiet=True)

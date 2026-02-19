import os
import glob
from collections import defaultdict
from ..utils.bio import read_fasta, write_fasta
from ..utils.helpers import run_cmd, find_cds_fasta_for_genome

def map_aa_id_to_cds_id(aa_id, mode="after_last_pipe"):
    if mode == "same":
        return aa_id
    if mode == "strip_pipe":
        return aa_id.split("|", 1)[-1]
    return aa_id.split("|")[-1]

def build_codon_alignment_pal2nal(aa_aln_fp, cds_subset_fp, out_codon_fp, pal2nal_cmd="pal2nal.pl", codon_format="fasta"):
    cmd = f"{pal2nal_cmd} {aa_aln_fp} {cds_subset_fp} -output {codon_format} > {out_codon_fp}"
    run_cmd(cmd, quiet=True, shell=True)

def run_codon(cfg, tree_dir, clipkit_dir, aln_dir, codon_dir, hmm_keep, force=False):
    print("\n[codon] Building per-HMM codon alignments (PAL2NAL)...")
    
    codon_cfg = cfg.get("codon", {})
    cds_dir = cfg["inputs"].get("cds_dir", None)
    if not cds_dir:
        raise SystemExit("codon.build_codon_alignments requires inputs.cds_dir")

    hmm_names = sorted([os.path.basename(x).split(".")[0] for x in glob.glob(os.path.join(tree_dir, "*.treefile"))])
    if hmm_keep is not None:
        hmm_names = [h for h in hmm_names if h in hmm_keep]

    for hmm in hmm_names:
        tree_fp = os.path.join(tree_dir, f"{hmm}.treefile")
        clip_aln = os.path.join(clipkit_dir, f"{hmm}.clipkit.faa")
        if os.path.exists(clip_aln):
            aa_aln_fp = clip_aln
        else:
            cand1 = os.path.join(aln_dir, f"{hmm}.afa")
            cand2 = os.path.join(aln_dir, f"{hmm}.mafft.fasta")
            aa_aln_fp = cand2 if os.path.exists(cand2) else cand1 if os.path.exists(cand1) else None
        if not aa_aln_fp or not os.path.exists(tree_fp):
            continue

        codon_aln_fp = os.path.join(codon_dir, f"{hmm}.codon.fasta")
        if os.path.exists(codon_aln_fp) and not force:
            continue

        aln_seqs = read_fasta(aa_aln_fp)
        cds_subset = {}
        tips_by_genome = defaultdict(list)
        for tip in aln_seqs.keys():
            genome = tip.split("|", 1)[0] if "|" in tip else "Unknown"
            tips_by_genome[genome].append(tip)

        for genome, tips in tips_by_genome.items():
            cds_fp = find_cds_fasta_for_genome(cds_dir, genome)
            if not cds_fp:
                continue
            cds_all = read_fasta(cds_fp)

            for tip in tips:
                cds_id = map_aa_id_to_cds_id(tip, mode=codon_cfg.get("cds_id_mode", "after_last_pipe"))
                candidates = [cds_id, cds_id.split("|")[-1], cds_id.replace("|", "_")]
                found = None
                for c in candidates:
                    if c in cds_all:
                        found = c
                        break
                if found:
                    cds_subset[tip] = cds_all[found]

        if len(cds_subset) < 3:
            continue

        cds_subset_fp = os.path.join(codon_dir, f"{hmm}.cds.fna")
        write_fasta(cds_subset_fp, cds_subset)

        try:
            build_codon_alignment_pal2nal(
                aa_aln_fp, cds_subset_fp, codon_aln_fp,
                pal2nal_cmd=codon_cfg.get("pal2nal_cmd", "pal2nal.pl")
            )
        except Exception:
            continue

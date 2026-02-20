import os
import sys
import glob
from collections import defaultdict
from ..utils.bio import read_fasta, write_fasta
from ..utils.helpers import run_cmd, find_cds_fasta_for_genome, normalize_genome_id

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

        missing_cds_files = []
        missing_cds_ids = []

        for genome, tips in tips_by_genome.items():
            cds_fp = find_cds_fasta_for_genome(cds_dir, genome)
            if not cds_fp:
                missing_cds_files.append(genome)
                continue
            cds_all = read_fasta(cds_fp)

            for tip in tips:
                cds_id = map_aa_id_to_cds_id(tip, mode=codon_cfg.get("cds_id_mode", "after_last_pipe"))
                # Build candidate list for fuzzy matching
                candidates = [
                    cds_id,                                # exact protein ID
                    cds_id.split("|")[-1],                 # last pipe segment
                    cds_id.replace("|", "_"),               # pipes to underscores
                    tip,                                    # full tip label
                    tip.split("|", 1)[-1] if "|" in tip else tip, # everything after first pipe
                ]
                # Also try with/without common suffixes
                for sfx in ["_cds", "_CDS", ".cds", ".p01", "_1"]:
                    candidates.append(cds_id + sfx)
                    if cds_id.endswith(sfx):
                        candidates.append(cds_id[:-len(sfx)])

                found = None
                for c in candidates:
                    if c in cds_all:
                        found = c
                        break

                # Last-resort: substring match (protein ID anywhere in CDS header)
                if found is None:
                    for cds_header in cds_all:
                        if cds_id in cds_header or cds_header in cds_id:
                            found = cds_header
                            break

                if found:
                    # Key in CDS subset must match the AA alignment header exactly
                    # so pal2nal can pair them by name
                    cds_subset[tip] = cds_all[found]
                else:
                    missing_cds_ids.append((genome, tip, cds_id))

        if missing_cds_files:
            print(f"[codon] Warning ({hmm}): No CDS file found for genomes: "
                  f"{', '.join(missing_cds_files[:5])}{'...' if len(missing_cds_files) > 5 else ''}",
                  file=sys.stderr)
        if missing_cds_ids:
            print(f"[codon] Warning ({hmm}): {len(missing_cds_ids)} protein(s) could not be matched to CDS. "
                  f"First few: {[(g, cid) for g, _, cid in missing_cds_ids[:3]]}",
                  file=sys.stderr)

        if len(cds_subset) < 3:
            print(f"[codon] Skipping {hmm}: only {len(cds_subset)} CDS sequences matched (need ≥3).", file=sys.stderr)
            continue

        cds_subset_fp = os.path.join(codon_dir, f"{hmm}.cds.fna")
        write_fasta(cds_subset_fp, cds_subset)

        try:
            build_codon_alignment_pal2nal(
                aa_aln_fp, cds_subset_fp, codon_aln_fp,
                pal2nal_cmd=codon_cfg.get("pal2nal_cmd", "pal2nal.pl")
            )
            # Validate output is non-empty
            if os.path.exists(codon_aln_fp):
                if os.path.getsize(codon_aln_fp) == 0:
                    print(f"[codon] Warning: pal2nal produced empty output for {hmm}. "
                          f"Check AA alignment and CDS name correspondence.", file=sys.stderr)
                    os.remove(codon_aln_fp)
                else:
                    print(f"[codon] Built codon alignment for {hmm} ({len(cds_subset)} seqs)")
        except Exception as e:
            print(f"[codon] pal2nal FAILED for {hmm}: {e}", file=sys.stderr)
            # Clean up partial output
            if os.path.exists(codon_aln_fp):
                os.remove(codon_aln_fp)
            continue

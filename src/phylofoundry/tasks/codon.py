import os
import sys
import glob
import re
from collections import defaultdict
from ..utils.bio import read_fasta, write_fasta
from ..utils.helpers import run_cmd, find_cds_fasta_for_genome, normalize_genome_id

# ---------------------------------------------------------------------------
# Stop-codon handling utilities
# ---------------------------------------------------------------------------

STOP_CODONS = {"TAA", "TAG", "TGA", "taa", "tag", "tga"}


def _strip_stop_from_protein(seq: str) -> str:
    """Remove trailing '*' (stop) from a protein sequence, including gap-padded."""
    seq = seq.rstrip("-").rstrip("*").rstrip("-")
    # Also remove any internal '*' that might exist (e.g. mid-sequence stops)
    return seq.replace("*", "")


def _trim_stop_codon_from_cds(seq: str) -> str:
    """Remove terminal stop codon (3 nt) from a CDS nucleotide sequence."""
    seq = seq.rstrip()
    if len(seq) >= 3 and seq[-3:].upper() in {"TAA", "TAG", "TGA"}:
        return seq[:-3]
    return seq


def _clean_cds_for_pal2nal(cds_seq: str) -> str:
    """Clean CDS sequence: remove terminal stop codon, strip whitespace."""
    return _trim_stop_codon_from_cds(cds_seq.replace("\n", "").replace(" ", ""))


def _clean_aa_alignment(aln_seqs: dict) -> dict:
    """Strip stop codons from all sequences in an AA alignment dict."""
    cleaned = {}
    for k, v in aln_seqs.items():
        cleaned[k] = _strip_stop_from_protein(v)
    return cleaned


# ---------------------------------------------------------------------------
# ID mapping
# ---------------------------------------------------------------------------

def map_aa_id_to_cds_id(aa_id, mode="after_last_pipe"):
    if mode == "same":
        return aa_id
    if mode == "strip_pipe":
        return aa_id.split("|", 1)[-1]
    return aa_id.split("|")[-1]


# ---------------------------------------------------------------------------
# pal2nal wrapper
# ---------------------------------------------------------------------------

def build_codon_alignment_pal2nal(aa_aln_fp, cds_subset_fp, out_codon_fp,
                                   pal2nal_cmd="pal2nal.pl",
                                   codon_format="fasta"):
    """Run pal2nal with -nogap -nomismatch to tolerate minor differences."""
    cmd = (f"{pal2nal_cmd} {aa_aln_fp} {cds_subset_fp} "
           f"-output {codon_format} -nogap -nomismatch > {out_codon_fp}")
    run_cmd(cmd, quiet=True, shell=True)


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def run_codon(cfg, tree_dir, clipkit_dir, aln_dir, codon_dir, hmm_keep, force=False):
    print("\n[codon] Building per-HMM codon alignments (PAL2NAL)...")

    codon_cfg = cfg.get("codon", {})
    cds_dir = cfg["inputs"].get("cds_dir", None)
    if not cds_dir:
        raise SystemExit("codon.build_codon_alignments requires inputs.cds_dir")

    hmm_names = sorted([
        os.path.basename(x).split(".")[0]
        for x in glob.glob(os.path.join(tree_dir, "*.treefile"))
    ])
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

        # ── Read and clean AA alignment (strip stop codons) ──────────────
        aln_seqs_raw = read_fasta(aa_aln_fp)
        aln_seqs = _clean_aa_alignment(aln_seqs_raw)

        # Write cleaned AA alignment to temp file for pal2nal
        cleaned_aa_fp = os.path.join(codon_dir, f"{hmm}.cleaned_aa.faa")
        write_fasta(cleaned_aa_fp, aln_seqs)

        # ── Find matching CDS for each tip ───────────────────────────────
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
                cds_id = map_aa_id_to_cds_id(
                    tip, mode=codon_cfg.get("cds_id_mode", "after_last_pipe")
                )
                # Build candidate list for fuzzy matching
                candidates = [
                    cds_id,
                    cds_id.split("|")[-1],
                    cds_id.replace("|", "_"),
                    tip,
                    tip.split("|", 1)[-1] if "|" in tip else tip,
                ]
                for sfx in ["_cds", "_CDS", ".cds", ".p01", "_1"]:
                    candidates.append(cds_id + sfx)
                    if cds_id.endswith(sfx):
                        candidates.append(cds_id[:-len(sfx)])

                found = None
                for c in candidates:
                    if c in cds_all:
                        found = c
                        break

                # Last-resort: substring match
                if found is None:
                    for cds_header in cds_all:
                        if cds_id in cds_header or cds_header in cds_id:
                            found = cds_header
                            break

                if found:
                    # Clean CDS: strip terminal stop codon so pal2nal
                    # can verify CDS translates to the (stop-stripped) protein
                    cleaned_cds = _clean_cds_for_pal2nal(cds_all[found])
                    # Key must match the AA alignment header exactly
                    cds_subset[tip] = cleaned_cds
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
            print(f"[codon] Skipping {hmm}: only {len(cds_subset)} CDS sequences matched (need ≥3).",
                  file=sys.stderr)
            continue

        cds_subset_fp = os.path.join(codon_dir, f"{hmm}.cds.fna")
        write_fasta(cds_subset_fp, cds_subset)

        # ── Run pal2nal ──────────────────────────────────────────────────
        try:
            build_codon_alignment_pal2nal(
                cleaned_aa_fp, cds_subset_fp, codon_aln_fp,
                pal2nal_cmd=codon_cfg.get("pal2nal_cmd", "pal2nal.pl")
            )
            # Validate output is non-empty
            if os.path.exists(codon_aln_fp):
                if os.path.getsize(codon_aln_fp) == 0:
                    print(f"[codon] Warning: pal2nal produced empty output for {hmm}. "
                          f"Check AA alignment and CDS name correspondence.",
                          file=sys.stderr)
                    os.remove(codon_aln_fp)
                else:
                    print(f"[codon] Built codon alignment for {hmm} ({len(cds_subset)} seqs)")
        except Exception as e:
            print(f"[codon] pal2nal FAILED for {hmm}: {e}", file=sys.stderr)
            if os.path.exists(codon_aln_fp):
                os.remove(codon_aln_fp)
            continue
        finally:
            # Clean up temp cleaned AA file
            if os.path.exists(cleaned_aa_fp):
                try:
                    os.remove(cleaned_aa_fp)
                except OSError:
                    pass

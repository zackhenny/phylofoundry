import os
import sys
import glob
import re
import subprocess
from collections import defaultdict
from ..utils.bio import read_fasta, write_fasta
from ..utils.helpers import run_cmd, find_cds_fasta_for_genome, normalize_genome_id

# ---------------------------------------------------------------------------
# Stop-codon handling utilities
# ---------------------------------------------------------------------------

STOP_CODONS = {"TAA", "TAG", "TGA", "taa", "tag", "tga"}


def _strip_stop_from_protein(seq: str) -> str:
    """Remove '*' (stop codon) from a protein sequence.

    Trailing stops are stripped first (including gap-padded ends), then any
    remaining internal '*' are removed.  All other residues—including 'X'
    (ambiguous) and 'C' (cysteine)—are preserved to avoid changing original
    AA positions.
    """
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
        cleaned[k] = _strip_stop_from_protein(v).upper()
    return cleaned


# ---------------------------------------------------------------------------
# ID mapping
# ---------------------------------------------------------------------------

def _read_tree_tips(treefile: str) -> set:
    """Return the set of leaf/tip labels from a Newick treefile.

    IQ-TREE writes tip names as bare tokens before branch-length fields.
    We extract every token that precedes a `:`/`,`/`)`, which includes both
    tip labels and internal node labels (e.g. bootstrap support values).
    Internal labels such as numeric bootstrap values (e.g. '98') are also
    captured, but because real sequence IDs are never purely numeric, the
    downstream ``k in tree_tips`` filter remains correct in practice.

    Note: if sequence IDs happen to be purely numeric they may be confused
    with bootstrap support values; avoid such IDs or use a full Newick
    parser (e.g. ``Bio.Phylo``) if this is a concern.
    """
    with open(treefile) as f:
        tree_str = f.read()
    tips = set(re.findall(r'([^(),;:\s]+)(?=:|[,)])', tree_str))
    return tips


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
    """Run pal2nal to thread CDS onto a protein alignment.

    The protein alignment gaps drive codon placement; do NOT use -nogap
    (which strips all gap-containing columns, destroying the alignment)
    or -nomismatch (which silently drops mismatched codons).

    Stderr from pal2nal is captured and printed so the user can see diagnostic
    messages (e.g. "ERROR: CDS does not match protein sequence").
    """
    cmd = (f"{pal2nal_cmd} {aa_aln_fp} {cds_subset_fp} "
           f"-output {codon_format} > {out_codon_fp}")
    result = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE, text=True)
    if result.stderr:
        print(f"[codon] pal2nal stderr:\n{result.stderr}", file=sys.stderr)


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
        os.path.basename(x).replace(".treefile", "")
        for x in glob.glob(os.path.join(tree_dir, "*.treefile"))
    ])
    if hmm_keep is not None:
        hmm_names = [h for h in hmm_names if h in hmm_keep]

    for hmm in hmm_names:
        tree_fp = os.path.join(tree_dir, f"{hmm}.treefile")
        # pal2nal needs full-length (untrimmed) protein sequences so that
        # len(CDS) == len(ungapped_AA) * 3.  ClipKit-trimmed alignments have
        # columns removed, which shrinks the ungapped AA count and causes
        # length mismatches.  Prefer the untrimmed alignment; fall back to
        # ClipKit only if no untrimmed version is available.
        cand1 = os.path.join(aln_dir, f"{hmm}.afa")
        cand2 = os.path.join(aln_dir, f"{hmm}.mafft.fasta")
        clip_aln = os.path.join(clipkit_dir, f"{hmm}.clipkit.faa")
        if os.path.exists(cand2):
            aa_aln_fp = cand2
        elif os.path.exists(cand1):
            aa_aln_fp = cand1
        elif os.path.exists(clip_aln):
            aa_aln_fp = clip_aln
            print(f"[codon] Warning ({hmm}): Using ClipKit-trimmed alignment "
                  f"(no untrimmed alignment found). Length mismatches may occur.",
                  file=sys.stderr)
        else:
            aa_aln_fp = None
        if not aa_aln_fp or not os.path.exists(tree_fp):
            continue

        codon_aln_fp = os.path.join(codon_dir, f"{hmm}.codon.fasta")
        if os.path.exists(codon_aln_fp) and not force:
            continue

        # ── Read and clean AA alignment (strip stop codons) ──────────────
        aln_seqs_raw = read_fasta(aa_aln_fp)
        aln_seqs = _clean_aa_alignment(aln_seqs_raw)

        # ── Filter to sequences present in the tree ───────────────────────
        # IQ-TREE may prune identical or problematic sequences; only retain
        # sequences that actually appear as tips in the treefile so the
        # resulting codon alignment is consistent with the tree.
        tree_tips = _read_tree_tips(tree_fp)
        n_before = len(aln_seqs)
        aln_seqs = {k: v for k, v in aln_seqs.items() if k in tree_tips}
        if len(aln_seqs) < n_before:
            print(f"[codon] {hmm}: Filtered alignment from {n_before} to "
                  f"{len(aln_seqs)} sequences to match tree tips.",
                  file=sys.stderr)

        cleaned_aa_fp = os.path.join(codon_dir, f"{hmm}.cleaned_aa.faa")

        # ── Find matching CDS for each tip ───────────────────────────────
        cds_subset = {}
        tips_by_genome = defaultdict(list)
        for tip in aln_seqs.keys():
            raw_genome = tip.split("|", 1)[0] if "|" in tip else "Unknown"
            genome = normalize_genome_id(raw_genome)
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

        # Reorder cds_subset to EXACTLY match aln_seqs order, required by pal2nal
        ordered_cds_subset = {}
        for tip in aln_seqs.keys():
            if tip in cds_subset:
                ordered_cds_subset[tip] = cds_subset[tip]

        # ── Debug: show first few AA and CDS headers ─────────────────────
        aa_headers = list(aln_seqs.keys())[:3]
        cds_headers = list(ordered_cds_subset.keys())[:3]
        print(f"[codon] {hmm}: AA headers (first 3): {aa_headers}", file=sys.stderr)
        print(f"[codon] {hmm}: CDS headers (first 3): {cds_headers}", file=sys.stderr)

        # ── Validate CDS/AA length compatibility ─────────────────────────
        length_ok = {}
        length_bad = []
        for tip, cds_seq in ordered_cds_subset.items():
            aa_seq = aln_seqs[tip]
            non_gap_aa = aa_seq.replace("-", "").replace(".", "")
            expected_cds_len = len(non_gap_aa) * 3
            if len(cds_seq) != expected_cds_len:
                length_bad.append((tip, len(non_gap_aa), len(cds_seq)))
            else:
                length_ok[tip] = cds_seq
        if length_bad:
            print(f"[codon] Warning ({hmm}): {len(length_bad)} sequence(s) excluded due to "
                  f"CDS/AA length mismatch (expected CDS=AA*3). "
                  f"First few: {[(t, aa, cds) for t, aa, cds in length_bad[:3]]}",
                  file=sys.stderr)
            ordered_cds_subset = length_ok

        if len(ordered_cds_subset) < 3:
            print(f"[codon] Skipping {hmm}: only {len(ordered_cds_subset)} sequences pass "
                  f"length validation (need ≥3).", file=sys.stderr)
            continue

        cds_subset_fp = os.path.join(codon_dir, f"{hmm}.cds.fna")
        write_fasta(cds_subset_fp, ordered_cds_subset)

        # Write cleaned AA alignment with only the sequences that have a matching CDS entry,
        # so that pal2nal receives files with the exact same set of sequences in the same order.
        matched_aa_seqs = {tip: aln_seqs[tip] for tip in ordered_cds_subset}
        write_fasta(cleaned_aa_fp, matched_aa_seqs)

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

"""
asr.py — Ancestral Sequence Reconstruction parser for IQ-TREE .state files.

Parses the marginal ancestral reconstruction output from IQ-TREE
(produced with the -asr flag) and reconstructs ancestral protein sequences
for each internal tree node.
"""

import os
import sys
from collections import defaultdict
from ..utils.bio import write_fasta
from ..utils.helpers import safe_mkdir


def parse_iqtree_state(state_fp: str) -> dict:
    """Parse IQ-TREE .state file → dict of {node_id: reconstructed_sequence}.

    The .state file format:
        # Ancestral state reconstruction by IQ-TREE
        # Columns: Node  Site  State  p_A  p_R  p_N  p_D ...
        Node1   1   A   0.999   0.000   ...
        Node1   2   R   0.002   0.995   ...
        Node2   1   M   0.800   0.100   ...
        ...

    For each node, the most-probable state per site is taken to reconstruct
    the ancestral sequence. Sites are 1-indexed.

    Returns
    -------
    dict
        {node_id: amino_acid_sequence_string}
    """
    if not os.path.exists(state_fp):
        print(f"[asr] .state file not found: {state_fp}", file=sys.stderr)
        return {}

    node_sites = defaultdict(dict)  # node -> {site_int: state_char}

    with open(state_fp) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                parts = line.split()
            if len(parts) < 3:
                continue

            node_id = parts[0]
            try:
                site = int(parts[1])
            except ValueError:
                continue
            state = parts[2]

            node_sites[node_id][site] = state

    # Reconstruct sequences: join states in site order
    ancestral_seqs = {}
    for node_id, sites in node_sites.items():
        if not sites:
            continue
        max_site = max(sites.keys())
        seq_chars = []
        for i in range(1, max_site + 1):
            seq_chars.append(sites.get(i, "X"))  # X for any missing site
        seq = "".join(seq_chars)
        # Skip gap-only or very short sequences
        if len(seq.replace("-", "").replace("X", "")) < 10:
            continue
        ancestral_seqs[node_id] = seq

    return ancestral_seqs


def run_asr_parse(cfg, tree_dir, fasta_dir, hmm_keep, force=False):
    """Parse IQ-TREE .state files for all HMMs and write ancestral FASTAs.

    Parameters
    ----------
    cfg : dict
        Pipeline config
    tree_dir : str
        Directory containing .treefile and .state files
    fasta_dir : str
        Directory to write ancestral_nodes FASTAs
    hmm_keep : set or None
        HMM names to process (None = all)
    force : bool
        Overwrite existing output

    Returns
    -------
    dict
        {hmm_name: {node_id: sequence}} for all HMMs with ASR output
    """
    import glob

    phy_cfg = cfg.get("phylo", {})
    if phy_cfg.get("no_asr", False):
        print("[asr] ASR disabled (phylo.no_asr = true). Skipping.")
        return {}

    all_ancestral = {}

    state_files = sorted(glob.glob(os.path.join(tree_dir, "*.state")))
    if not state_files:
        print("[asr] No .state files found in tree directory. "
              "Run IQ-TREE with -asr flag to generate them.")
        return {}

    for state_fp in state_files:
        hmm = os.path.basename(state_fp).replace(".state", "")
        if hmm_keep is not None and hmm not in hmm_keep:
            continue

        out_fasta = os.path.join(fasta_dir, f"{hmm}.ancestral_nodes.fasta")
        if os.path.exists(out_fasta) and not force:
            # Load existing
            from ..utils.bio import read_fasta
            all_ancestral[hmm] = read_fasta(out_fasta)
            print(f"[asr] Loaded existing ancestral sequences for {hmm}: "
                  f"{len(all_ancestral[hmm])} nodes")
            continue

        print(f"[asr] Parsing ancestral states for {hmm}...")
        anc_seqs = parse_iqtree_state(state_fp)

        if not anc_seqs:
            print(f"[asr] No ancestral sequences reconstructed for {hmm}.",
                  file=sys.stderr)
            continue

        write_fasta(out_fasta, anc_seqs)
        all_ancestral[hmm] = anc_seqs
        print(f"[asr] Wrote {len(anc_seqs)} ancestral sequences for {hmm} → {out_fasta}")

    # Also handle combined tree if present
    combined_state = os.path.join(tree_dir, "combined_all_hits.state")
    if os.path.exists(combined_state):
        out_combined = os.path.join(fasta_dir, "combined_all_hits.ancestral_nodes.fasta")
        if not os.path.exists(out_combined) or force:
            print("[asr] Parsing ancestral states for combined tree...")
            anc_seqs = parse_iqtree_state(combined_state)
            if anc_seqs:
                write_fasta(out_combined, anc_seqs)
                all_ancestral["combined_all_hits"] = anc_seqs
                print(f"[asr] Wrote {len(anc_seqs)} ancestral sequences "
                      f"for combined tree → {out_combined}")

    return all_ancestral

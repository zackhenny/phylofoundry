"""
motifs.py — ESM-2 attention-based targeted motif scoring.

For each protein sequence, extracts attention weights directed at
user-specified motif positions (e.g., HPEVY, HPEVF) to quantify
how structurally important the model considers that motif.
"""

import os
import sys
import re
import numpy as np
import pandas as pd
from ..utils.bio import read_fasta


def _find_motif_positions(sequence: str, motif: str) -> list:
    """Find all start positions (0-indexed) of `motif` in `sequence`."""
    positions = []
    start = 0
    while True:
        idx = sequence.find(motif, start)
        if idx == -1:
            break
        positions.append(idx)
        start = idx + 1
    return positions


def _extract_attention_at_positions(model, alphabet, sequence: str,
                                    positions: list, motif_len: int,
                                    device: str, n_layers: int = 4):
    """Extract attention weights directed at motif positions.

    Uses the last `n_layers` layers, averages across all attention heads.
    For each motif position i..i+motif_len, sums the attention column
    (how much the rest of the sequence "looks at" the motif).

    Returns
    -------
    float
        Mean attention score for the motif positions, normalized.
    """
    import torch

    batch_converter = alphabet.get_batch_converter()
    _, _, toks = batch_converter([("seq", sequence)])
    toks = toks.to(device)

    with torch.no_grad():
        results = model(toks, repr_layers=[], return_contacts=False)

    # results["attentions"] shape: (1, n_layers, n_heads, seq_len, seq_len)
    # Note: toks includes BOS/EOS tokens, so position offset is +1
    attentions = results["attentions"]  # (1, L, H, T, T)
    total_layers = attentions.shape[1]
    start_layer = max(0, total_layers - n_layers)

    # Average over last n_layers and all heads
    attn_avg = attentions[0, start_layer:, :, :, :].mean(dim=(0, 1))  # (T, T)

    scores = []
    for pos in positions:
        # +1 offset for BOS token
        motif_indices = list(range(pos + 1, pos + 1 + motif_len))
        # Clamp to valid range
        motif_indices = [i for i in motif_indices if i < attn_avg.shape[0]]
        if not motif_indices:
            continue
        # Sum attention FROM all positions TO motif positions
        attn_to_motif = attn_avg[:, motif_indices].sum().item()
        # Normalize by sequence length
        seq_len = attn_avg.shape[0]
        normalized_score = attn_to_motif / seq_len if seq_len > 0 else 0.0
        scores.append(normalized_score)

    return np.mean(scores) if scores else 0.0


def score_motifs(cfg, fasta_dir, summary_dir, hmm_keep, force=False):
    """Score user-defined motifs using ESM-2 attention weights.

    Parameters
    ----------
    cfg : dict
        Pipeline config with cfg["motifs"]["motif_list"]
    fasta_dir : str
        Directory containing per-HMM FASTA files
    summary_dir : str
        Output directory for scores TSV
    hmm_keep : set or None
    force : bool

    Returns
    -------
    pd.DataFrame or None
    """
    import glob
    import torch

    motif_cfg = cfg.get("motifs", {})
    if not motif_cfg.get("enabled", False):
        return None

    motif_list = motif_cfg.get("motif_list", [])
    if isinstance(motif_list, str):
        motif_list = [m.strip() for m in motif_list.split(",") if m.strip()]
    if not motif_list:
        print("[motifs] No motifs specified. Set motifs.motif_list in config.",
              file=sys.stderr)
        return None

    out_fp = os.path.join(summary_dir, "motif_attention_scores.tsv")
    if os.path.exists(out_fp) and not force:
        print(f"[motifs] Output exists: {out_fp}. Use --force to override.")
        return pd.read_csv(out_fp, sep="\t")

    # Load clade assignments if available
    clade_fp = os.path.join(summary_dir, "clade_assignment.tsv")
    clade_map = {}
    if os.path.exists(clade_fp):
        clade_df = pd.read_csv(clade_fp, sep="\t")
        for _, row in clade_df.iterrows():
            key = row.get("protein", row.get("seq_id", ""))
            clade_map[key] = str(row.get("cluster_id", ""))

    # Load ESM model
    emb_cfg = cfg.get("embeddings", {})
    model_name = emb_cfg.get("model", "esm2_t33_650M_UR50D")
    device = emb_cfg.get("device", "cuda")
    n_layers = motif_cfg.get("attention_layers", 4)

    print(f"[motifs] Loading ESM model {model_name} for attention extraction...")
    try:
        import esm
        model, alphabet = esm.pretrained.load_model_and_alphabet(model_name)
        model.eval()
        model = model.to(device)
    except Exception as e:
        print(f"[motifs] Failed to load ESM model: {e}", file=sys.stderr)
        return None

    # Process all HMM FASTAs
    fasta_files = sorted(glob.glob(os.path.join(fasta_dir, "*.faa")))
    all_rows = []

    for fasta_fp in fasta_files:
        hmm = os.path.basename(fasta_fp).replace(".faa", "")
        if hmm_keep is not None and hmm not in hmm_keep:
            continue

        seqs = read_fasta(fasta_fp)
        print(f"[motifs] Scoring {len(seqs)} sequences for {hmm} "
              f"({len(motif_list)} motifs)...")

        # Also check for ancestral sequences
        anc_fp = os.path.join(fasta_dir, f"{hmm}.ancestral_nodes.fasta")
        if os.path.exists(anc_fp):
            anc_seqs = read_fasta(anc_fp)
            for k, v in anc_seqs.items():
                seqs[f"ANC|{k}"] = v

        for seq_id, seq in seqs.items():
            # Clean sequence
            clean_seq = seq.replace("*", "").replace("-", "").replace("X", "")
            if len(clean_seq) < 10:
                continue

            seq_type = "ancestral" if seq_id.startswith("ANC|") else "modern"

            for motif in motif_list:
                positions = _find_motif_positions(clean_seq, motif)
                if not positions:
                    # Record zero score for sequences lacking the motif
                    all_rows.append({
                        "hmm": hmm,
                        "seq_id": seq_id,
                        "motif": motif,
                        "start_pos": -1,
                        "end_pos": -1,
                        "attention_score": 0.0,
                        "motif_present": False,
                        "clade_id": clade_map.get(seq_id, ""),
                        "type": seq_type,
                    })
                    continue

                try:
                    score = _extract_attention_at_positions(
                        model, alphabet, clean_seq, positions,
                        len(motif), device, n_layers=n_layers
                    )
                except Exception as e:
                    print(f"[motifs] Attention extraction failed for {seq_id}: {e}",
                          file=sys.stderr)
                    score = float("nan")

                for pos in positions:
                    all_rows.append({
                        "hmm": hmm,
                        "seq_id": seq_id,
                        "motif": motif,
                        "start_pos": pos,
                        "end_pos": pos + len(motif),
                        "attention_score": score,
                        "motif_present": True,
                        "clade_id": clade_map.get(seq_id, ""),
                        "type": seq_type,
                    })

    if not all_rows:
        print("[motifs] No motif scores generated.", file=sys.stderr)
        return None

    df = pd.DataFrame(all_rows)
    os.makedirs(summary_dir, exist_ok=True)
    df.to_csv(out_fp, sep="\t", index=False)
    print(f"[motifs] Wrote motif attention scores: {out_fp} "
          f"({len(df)} rows)")
    return df

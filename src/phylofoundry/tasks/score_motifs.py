"""
score_motifs.py — ESM-2 attention-based targeted motif scoring.

For each protein sequence, extracts attention weights directed at
user-specified motif positions (e.g., HPEVY, HPEVF) to quantify
how structurally important the model considers that motif.
"""

import os
import sys
import numpy as np
import pandas as pd
from ..utils.bio import read_fasta
from ..utils.ha import compute_ha_from_layer_profiles, write_ha_outputs


def _find_motif_positions(sequence: str, motif: str) -> list:
    positions = []
    start = 0
    while True:
        idx = sequence.find(motif, start)
        if idx == -1:
            break
        positions.append(idx)
        start = idx + 1
    return positions


def _get_attention_tensor(model, alphabet, sequence: str, device: str):
    import torch

    batch_converter = alphabet.get_batch_converter()
    _, _, toks = batch_converter([("seq", sequence)])
    toks = toks.to(device)
    with torch.no_grad():
        results = model(toks, repr_layers=[], return_contacts=False)
    return results["attentions"][0]  # (L, H, T, T)


def _attention_received_by_layer(attentions, seq_len: int) -> np.ndarray:
    # keep amino-acid positions only (skip BOS/EOS)
    # per layer+position: mean_h(sum_i A[l, h, i, j]) -> (L, S)
    by_layer = attentions[:, :, 1 : seq_len + 1, 1 : seq_len + 1].sum(dim=2).mean(dim=1)
    return by_layer.cpu().numpy()


def _extract_attention_at_positions(attentions, positions: list, motif_len: int, n_layers: int = 4):
    total_layers = attentions.shape[0]
    start_layer = max(0, total_layers - n_layers)
    attn_avg = attentions[start_layer:, :, :, :].mean(dim=(0, 1))

    scores = []
    for pos in positions:
        motif_indices = list(range(pos + 1, pos + 1 + motif_len))
        motif_indices = [i for i in motif_indices if i < attn_avg.shape[0]]
        if not motif_indices:
            continue
        attn_to_motif = attn_avg[:, motif_indices].sum().item()
        seq_tok_len = attn_avg.shape[0]
        scores.append(attn_to_motif / seq_tok_len if seq_tok_len > 0 else 0.0)
    return float(np.mean(scores)) if scores else 0.0


def _motif_ha_metrics(positions, motif_len, ha_scores, ha_mask):
    if not positions:
        return 0, np.nan, np.nan, np.nan, np.nan

    overlaps = []
    score_means = []
    score_sums = []
    for pos in positions:
        idx = np.arange(pos, min(pos + motif_len, len(ha_scores)))
        if len(idx) == 0:
            continue
        overlaps.append(float(np.mean(ha_mask[idx])))
        score_means.append(float(np.mean(ha_scores[idx])))
        score_sums.append(float(np.sum(ha_scores[idx])))

    if not overlaps:
        return len(positions), np.nan, np.nan, np.nan, np.nan
    return len(positions), float(np.mean(overlaps)), float(np.max(overlaps)), float(np.mean(score_means)), float(np.mean(score_sums))


def score_motifs(cfg, fasta_dir, summary_dir, motifs_dir, hmm_keep, force=False):
    import glob

    motif_cfg = cfg.get("motifs", {})
    if not motif_cfg.get("enabled", False):
        return None

    motif_list = motif_cfg.get("motif_list", [])
    if isinstance(motif_list, str):
        motif_list = [m.strip() for m in motif_list.split(",") if m.strip()]
    if not motif_list:
        print("[motifs] No motifs specified. Set motifs.motif_list in config.", file=sys.stderr)
        return None

    ha_cfg = cfg.get("ha", {})
    use_ha = bool(ha_cfg.get("enabled", False) and motif_cfg.get("use_ha", False))
    if use_ha and not cfg.get("embeddings", {}).get("write_full_vectors", False):
        raise SystemExit(
            "HA motif scoring requires embeddings.write_full_vectors: true to ensure attention-derived "
            "per-residue vectors are available/reproducible."
        )

    os.makedirs(motifs_dir, exist_ok=True)
    out_fp = os.path.join(motifs_dir, "motif_attention_scores.tsv")
    out_ha_fp = os.path.join(motifs_dir, "motif_ha_scores.tsv")
    if os.path.exists(out_fp) and not force and (not use_ha or os.path.exists(out_ha_fp)):
        print(f"[motifs] Output exists: {out_fp}. Use --force to override.")
        return pd.read_csv(out_fp, sep="\t")

    clade_fp = os.path.join(summary_dir, "clade_assignment.tsv")
    clade_map = {}
    if os.path.exists(clade_fp):
        clade_df = pd.read_csv(clade_fp, sep="\t")
        for _, row in clade_df.iterrows():
            key = row.get("protein", row.get("seq_id", ""))
            clade_map[key] = str(row.get("cluster_id", ""))

    # Load detected clades from post step if available
    detected_fp = os.path.join(summary_dir, "detected_clades.tsv")
    detected_clade_map = {}
    if os.path.exists(detected_fp):
        detected_df = pd.read_csv(detected_fp, sep="\t")
        for _, row in detected_df.iterrows():
            detected_clade_map[row["tip"]] = row["clade_name"]

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

    fasta_files = sorted(glob.glob(os.path.join(fasta_dir, "*.faa")))
    all_rows = []
    ha_rows = []
    attention_dir = os.path.join(os.path.dirname(summary_dir), "attention")

    for fasta_fp in fasta_files:
        hmm = os.path.basename(fasta_fp).replace(".faa", "")
        if hmm_keep is not None and hmm not in hmm_keep:
            continue

        seqs = read_fasta(fasta_fp)
        print(f"[motifs] Scoring {len(seqs)} sequences for {hmm} ({len(motif_list)} motifs)...")

        anc_fp = os.path.join(fasta_dir, f"{hmm}.ancestral_nodes.fasta")
        if os.path.exists(anc_fp):
            anc_seqs = read_fasta(anc_fp)
            for k, v in anc_seqs.items():
                seqs[f"ANC|{k}"] = v

        hmm_ha_scores = {}
        hmm_ha_masks = {}
        hmm_loc_layers = {}
        used_layer_start = None
        used_layer_end = None

        for seq_id, seq in seqs.items():
            clean_seq = seq.replace("*", "").replace("-", "").replace("X", "")
            if len(clean_seq) < 10:
                continue

            seq_type = "ancestral" if seq_id.startswith("ANC|") else "modern"
            try:
                attentions = _get_attention_tensor(model, alphabet, clean_seq, device)
            except Exception as e:
                print(f"[motifs] Attention extraction failed for {seq_id}: {e}", file=sys.stderr)
                continue

            if use_ha:
                layer_by_pos = _attention_received_by_layer(attentions, len(clean_seq))
                ha_scores, ha_mask, ls, le, loc_layer = compute_ha_from_layer_profiles(layer_by_pos, ha_cfg)
                hmm_ha_scores[seq_id] = ha_scores
                hmm_ha_masks[seq_id] = ha_mask
                if loc_layer is not None:
                    hmm_loc_layers[seq_id] = loc_layer
                used_layer_start, used_layer_end = ls, le

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
                        "detected_clade": detected_clade_map.get(seq_id, ""),
                        "type": seq_type,
                    })
                    continue

                score = _extract_attention_at_positions(attentions, positions, len(motif), n_layers=n_layers)
                for pos in positions:
                    all_rows.append(
                        {
                            "hmm": hmm,
                            "seq_id": seq_id,
                            "motif": motif,
                            "start_pos": pos,
                            "end_pos": pos + len(motif),
                            "attention_score": score,
                            "motif_present": True,
                            "clade_id": clade_map.get(seq_id, ""),
                            "type": seq_type,
                        }
                    )

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
                        "detected_clade": detected_clade_map.get(seq_id, ""),
                        "type": seq_type,
                    })

        if use_ha and hmm_ha_scores:
            write_ha_outputs(
                attention_dir,
                summary_dir,
                hmm,
                {k: v for k, v in seqs.items() if k in hmm_ha_scores},
                hmm_ha_scores,
                hmm_ha_masks,
                used_layer_start,
                used_layer_end,
                ha_cfg,
                per_seq_loc_layer=hmm_loc_layers,
            )

    if not all_rows:
        print("[motifs] No motif scores generated.", file=sys.stderr)
        return None

    df = pd.DataFrame(all_rows)
    os.makedirs(summary_dir, exist_ok=True)
    df.to_csv(out_fp, sep="\t", index=False)
    print(f"[motifs] Wrote motif attention scores: {out_fp} ({len(df)} rows)")

    if use_ha and ha_rows:
        pd.DataFrame(ha_rows).to_csv(out_ha_fp, sep="\t", index=False)
        print(f"[motifs] Wrote motif HA scores: {out_ha_fp} ({len(ha_rows)} rows)")

    return df

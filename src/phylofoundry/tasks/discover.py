"""
discover.py — Unsupervised motif discovery via comparative ESM-2 attention profiles.

Compares per-residue 1D attention profiles between two HDBSCAN clades
(e.g., standard vs novel) to find attention peaks unique to the novel clade.
Extracts k-mers at those peak positions as candidate structural hubs.
"""

import os
import sys
import numpy as np
import pandas as pd
from ..utils.bio import read_fasta


def _compute_attention_profile(model, alphabet, sequence: str, device: str,
                                n_layers: int = 4) -> np.ndarray:
    """Compute per-residue 1D attention profile for a sequence.

    For each residue position, sums how much attention it receives
    from ALL other positions (column-sum of the attention matrix).

    Returns
    -------
    np.ndarray
        Shape (seq_len,) attention profile.
    """
    import torch

    batch_converter = alphabet.get_batch_converter()
    _, _, toks = batch_converter([("seq", sequence)])
    toks = toks.to(device)

    with torch.no_grad():
        results = model(toks, repr_layers=[], return_contacts=False)

    attentions = results["attentions"]  # (1, L, H, T, T)
    total_layers = attentions.shape[1]
    start_layer = max(0, total_layers - n_layers)

    # Average over last n_layers and all heads → (T, T)
    attn_avg = attentions[0, start_layer:, :, :, :].mean(dim=(0, 1))

    # Column sum: how much attention each position RECEIVES
    # Exclude BOS (0) and EOS (last) tokens
    seq_len = len(sequence)
    profile = attn_avg[1:seq_len + 1, 1:seq_len + 1].sum(dim=0).cpu().numpy()

    return profile


def _align_profiles_to_fixed_length(profiles: list, target_len: int = 500) -> np.ndarray:
    """Interpolate variable-length profiles to a fixed-length grid.

    This allows averaging profiles across sequences of different lengths.
    """
    from scipy.interpolate import interp1d

    aligned = []
    for p in profiles:
        if len(p) == 0:
            continue
        x_old = np.linspace(0, 1, len(p))
        x_new = np.linspace(0, 1, target_len)
        f = interp1d(x_old, p, kind="linear", fill_value="extrapolate")
        aligned.append(f(x_new))

    return np.array(aligned) if aligned else np.zeros((0, target_len))


def _find_peaks_in_delta(delta: np.ndarray, top_n: int = 20,
                          min_distance: int = 5) -> list:
    """Find top-N peaks in the attention delta signal.

    Returns list of (position, delta_value) tuples, sorted by delta descending.
    """
    try:
        from scipy.signal import find_peaks
        peaks, properties = find_peaks(delta, distance=min_distance)
        if len(peaks) == 0:
            # Fallback: just take top positions
            top_idx = np.argsort(delta)[::-1][:top_n]
            return [(int(i), float(delta[i])) for i in top_idx]

        # Sort peaks by delta value
        peak_vals = delta[peaks]
        sorted_idx = np.argsort(peak_vals)[::-1][:top_n]
        return [(int(peaks[i]), float(peak_vals[i])) for i in sorted_idx]
    except ImportError:
        # Without scipy, just take top-N positions
        top_idx = np.argsort(delta)[::-1][:top_n]
        return [(int(i), float(delta[i])) for i in top_idx]


def _extract_kmer_at_position(sequence: str, position: int, k: int = 5) -> str:
    """Extract a k-mer centered at the given position."""
    half = k // 2
    start = max(0, position - half)
    end = min(len(sequence), start + k)
    start = max(0, end - k)  # adjust if near the end
    return sequence[start:end]


def discover_motifs(cfg, fasta_dir, summary_dir, hmm_keep, force=False):
    """Discover novel motifs by comparing attention profiles between HDBSCAN clades.

    Parameters
    ----------
    cfg : dict
        Pipeline config with cfg["discover"] section
    fasta_dir : str
        Directory containing per-HMM FASTA files
    summary_dir : str
        Directory containing clade_assignment.tsv and for output
    hmm_keep : set or None
    force : bool

    Returns
    -------
    pd.DataFrame or None
    """
    import glob
    import torch

    disc_cfg = cfg.get("discover", {})
    if not disc_cfg.get("enabled", False):
        return None

    standard_clade = disc_cfg.get("standard_clade", None)
    novel_clade = disc_cfg.get("novel_clade", None)
    if standard_clade is None or novel_clade is None:
        print("[discover] standard_clade and novel_clade must be set. Skipping.",
              file=sys.stderr)
        return None

    standard_clade = int(standard_clade)
    novel_clade = int(novel_clade)
    kmer_size = disc_cfg.get("kmer_size", 5)
    top_n = disc_cfg.get("top_n_peaks", 20)
    n_layers = disc_cfg.get("attention_layers", 4)

    out_fp = os.path.join(summary_dir, "discovered_motifs.tsv")
    if os.path.exists(out_fp) and not force:
        print(f"[discover] Output exists: {out_fp}. Use --force to override.")
        return pd.read_csv(out_fp, sep="\t")

    # Load clade assignments
    clade_fp = os.path.join(summary_dir, "clade_assignment.tsv")
    if not os.path.exists(clade_fp):
        print(f"[discover] clade_assignment.tsv not found at {clade_fp}. "
              f"Run embed step with clustering first.", file=sys.stderr)
        return None

    clade_df = pd.read_csv(clade_fp, sep="\t")
    std_ids = set(clade_df[clade_df["cluster_id"] == standard_clade]["protein"].values)
    nov_ids = set(clade_df[clade_df["cluster_id"] == novel_clade]["protein"].values)

    if not std_ids or not nov_ids:
        print(f"[discover] Standard clade {standard_clade} has {len(std_ids)} seqs, "
              f"novel clade {novel_clade} has {len(nov_ids)} seqs. "
              f"Need both to be non-empty.", file=sys.stderr)
        return None

    print(f"[discover] Standard clade {standard_clade}: {len(std_ids)} seqs, "
          f"Novel clade {novel_clade}: {len(nov_ids)} seqs")

    # Load all sequences
    all_seqs = {}
    fasta_files = sorted(glob.glob(os.path.join(fasta_dir, "*.faa")))
    for fp in fasta_files:
        hmm = os.path.basename(fp).replace(".faa", "")
        if hmm_keep is not None and hmm not in hmm_keep:
            continue
        seqs = read_fasta(fp)
        all_seqs.update(seqs)

        # Also load ancestral sequences if available
        anc_fp = os.path.join(fasta_dir, f"{hmm}.ancestral_nodes.fasta")
        if os.path.exists(anc_fp):
            anc_seqs = read_fasta(anc_fp)
            for k, v in anc_seqs.items():
                all_seqs[f"ANC|{k}"] = v

    # Filter sequences by clade
    std_seqs = {k: v.replace("*", "").replace("-", "")
                for k, v in all_seqs.items() if k in std_ids}
    nov_seqs = {k: v.replace("*", "").replace("-", "")
                for k, v in all_seqs.items() if k in nov_ids}

    if not std_seqs or not nov_seqs:
        print("[discover] Could not find sequences for the specified clades.",
              file=sys.stderr)
        return None

    # Load ESM model
    emb_cfg = cfg.get("embeddings", {})
    model_name = emb_cfg.get("model", "esm2_t33_650M_UR50D")
    device = emb_cfg.get("device", "cuda")

    print(f"[discover] Loading ESM model {model_name} for attention profile extraction...")
    try:
        import esm
        model, alphabet = esm.pretrained.load_model_and_alphabet(model_name)
        model.eval()
        model = model.to(device)
    except Exception as e:
        print(f"[discover] Failed to load ESM model: {e}", file=sys.stderr)
        return None

    # Compute attention profiles for each clade
    print(f"[discover] Computing attention profiles for standard clade "
          f"({len(std_seqs)} seqs)...")
    std_profiles = []
    for seq_id, seq in std_seqs.items():
        if len(seq) < 10:
            continue
        try:
            profile = _compute_attention_profile(model, alphabet, seq, device,
                                                  n_layers=n_layers)
            std_profiles.append(profile)
        except Exception as e:
            print(f"[discover] Profile extraction failed for {seq_id}: {e}",
                  file=sys.stderr)

    print(f"[discover] Computing attention profiles for novel clade "
          f"({len(nov_seqs)} seqs)...")
    nov_profiles = []
    nov_profile_seqs = []  # Track which sequences produced profiles
    for seq_id, seq in nov_seqs.items():
        if len(seq) < 10:
            continue
        try:
            profile = _compute_attention_profile(model, alphabet, seq, device,
                                                  n_layers=n_layers)
            nov_profiles.append(profile)
            nov_profile_seqs.append((seq_id, seq))
        except Exception as e:
            print(f"[discover] Profile extraction failed for {seq_id}: {e}",
                  file=sys.stderr)

    if not std_profiles or not nov_profiles:
        print("[discover] Insufficient profiles computed.", file=sys.stderr)
        return None

    # Align profiles to common length and compute delta
    target_len = 500
    std_aligned = _align_profiles_to_fixed_length(std_profiles, target_len)
    nov_aligned = _align_profiles_to_fixed_length(nov_profiles, target_len)

    std_mean = std_aligned.mean(axis=0) if len(std_aligned) > 0 else np.zeros(target_len)
    nov_mean = nov_aligned.mean(axis=0) if len(nov_aligned) > 0 else np.zeros(target_len)

    attention_delta = nov_mean - std_mean

    # Find peaks in the delta (where novel clade has higher attention)
    peaks = _find_peaks_in_delta(attention_delta, top_n=top_n)

    # Map peaks back to raw sequences and extract k-mers
    all_rows = []
    for peak_pos_norm, delta_val in peaks:
        if delta_val <= 0:
            continue  # Only interested in positive delta (novel > standard)

        # Map normalized position back to raw sequence positions
        for seq_id, seq in nov_profile_seqs:
            raw_pos = int(peak_pos_norm / target_len * len(seq))
            raw_pos = min(raw_pos, len(seq) - 1)
            kmer = _extract_kmer_at_position(seq, raw_pos, k=kmer_size)

            all_rows.append({
                "kmer": kmer,
                "position": raw_pos,
                "normalized_position": peak_pos_norm,
                "attention_delta": delta_val,
                "source_clade": novel_clade,
                "representative_seq_id": seq_id,
            })

    if not all_rows:
        print("[discover] No significant attention peaks found.", file=sys.stderr)
        return None

    # Aggregate: count unique k-mers and their frequencies
    df = pd.DataFrame(all_rows)

    # Group by kmer and compute summary stats
    kmer_summary = (
        df.groupby("kmer")
        .agg(
            n_sequences=("representative_seq_id", "nunique"),
            mean_attention_delta=("attention_delta", "mean"),
            median_position=("position", "median"),
            representative_seq_id=("representative_seq_id", "first"),
        )
        .reset_index()
        .sort_values("mean_attention_delta", ascending=False)
        .head(top_n)
    )
    kmer_summary["source_clade"] = novel_clade
    kmer_summary["reference_clade"] = standard_clade

    os.makedirs(summary_dir, exist_ok=True)
    kmer_summary.to_csv(out_fp, sep="\t", index=False)
    print(f"[discover] Wrote discovered motifs: {out_fp} ({len(kmer_summary)} k-mers)")

    return kmer_summary

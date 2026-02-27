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

    # Build clade_to_proteins mapping from HDBSCAN clusters (protein column)
    # Keys are cluster IDs (numeric); values are sets of protein identifiers
    clade_to_proteins = {}
    for cluster_id, grp in clade_df.groupby("cluster_id"):
        if pd.isna(cluster_id) or cluster_id == -1:
            continue
        clade_to_proteins[cluster_id] = set(grp["protein"].values)

    # Load detected clades from post step if available
    detected_fp = os.path.join(summary_dir, "detected_clades.tsv")
    if os.path.exists(detected_fp):
        detected_df = pd.read_csv(detected_fp, sep="\t")
        n_detected_proteins = 0
        for clade_name, grp in detected_df.groupby("clade_name"):
            tips = set(grp["tip"].values)
            clade_to_proteins[clade_name] = tips
            n_detected_proteins += len(tips)
        n_detected_clades = detected_df["clade_name"].nunique()
        print(f"[discover] Loaded {n_detected_clades} detected clades from "
              f"detected_clades.tsv ({n_detected_proteins} total proteins).")

    all_clades = list(clade_to_proteins.keys())

    if len(all_clades) < 2:
        print("[discover] Not enough clusters (< 2) for comparative motif discovery.", file=sys.stderr)
        return None

    print(f"[discover] Found {len(all_clades)} clades for motif discovery.")

    # Load all sequences
    all_seqs = {}
    fasta_files = sorted(glob.glob(os.path.join(fasta_dir, "*.faa")))
    for fp in fasta_files:
        hmm = os.path.basename(fp).replace(".faa", "")
        if hmm_keep is not None and hmm not in hmm_keep:
            continue
        seqs = read_fasta(fp)
        all_seqs.update({k: v.replace("*", "").replace("-", "") for k, v in seqs.items()})

    # Load ESM model globally
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

    # Pre-calculate attention profiles for all sequences to avoid re-computing
    print("[discover] Pre-computing attention profiles for all sequences...")
    all_profiles = {}
    for seq_id, seq in all_seqs.items():
        if len(seq) < 10: continue
        try:
            profile = _compute_attention_profile(model, alphabet, seq, device, n_layers=n_layers)
            all_profiles[seq_id] = profile
        except Exception as e:
            pass

    # Helper: resolve protein identifiers for a clade against available profiles.
    # HDBSCAN clusters store short protein IDs; sequence keys are genome|protein format.
    # Detected clades store tip as genome|protein directly.
    # Build a suffix lookup map once: protein_id -> seq_id (genome|protein)
    suffix_lookup = {}
    for seq_id in all_profiles:
        if "|" in seq_id:
            protein_part = seq_id.split("|", 1)[1]
            suffix_lookup.setdefault(protein_part, seq_id)

    def _resolve_seq_ids(protein_ids):
        """Return seq_ids present in all_profiles, trying direct match first,
        then looking up via the suffix map (genome|protein format)."""
        resolved = set()
        for pid in protein_ids:
            if pid in all_profiles:
                resolved.add(pid)
            elif pid in suffix_lookup:
                resolved.add(suffix_lookup[pid])
        return resolved

    # 1-vs-All clade comparison
    target_len = 500
    all_discovered_motifs = []

    # Determine which clades to iterate over
    if standard_clade is not None and novel_clade is not None:
        print(f"[discover] Manual override: comparing novel_clade={novel_clade} "
              f"vs standard_clade={standard_clade}.")
        comparison_clades = [novel_clade]
    else:
        comparison_clades = all_clades

    for focal_clade in comparison_clades:
        if focal_clade not in clade_to_proteins:
            print(f"[discover] Clade '{focal_clade}' not found in clade assignments. "
                  f"Skipping.", file=sys.stderr)
            continue

        nov_protein_ids = clade_to_proteins[focal_clade]

        if standard_clade is not None and novel_clade is not None:
            # Manual override: compare against specified standard_clade only
            if standard_clade not in clade_to_proteins:
                print(f"[discover] standard_clade '{standard_clade}' not found in "
                      f"clade assignments. Skipping.", file=sys.stderr)
                continue
            std_protein_ids = clade_to_proteins[standard_clade]
        else:
            # 1-vs-All: compare focal clade against all other clades combined
            std_protein_ids = set()
            for other_clade, other_ids in clade_to_proteins.items():
                if other_clade != focal_clade:
                    std_protein_ids.update(other_ids)

        nov_ids = _resolve_seq_ids(nov_protein_ids)
        std_ids = _resolve_seq_ids(std_protein_ids)

        nov_profiles = [all_profiles[s] for s in nov_ids if s in all_profiles]
        std_profiles = [all_profiles[s] for s in std_ids if s in all_profiles]

        if not nov_profiles or not std_profiles:
            continue

        nov_aligned = _align_profiles_to_fixed_length(nov_profiles, target_len)
        std_aligned = _align_profiles_to_fixed_length(std_profiles, target_len)

        nov_mean = nov_aligned.mean(axis=0)
        std_mean = std_aligned.mean(axis=0)

        attention_delta = nov_mean - std_mean
        peaks = _find_peaks_in_delta(attention_delta, top_n=top_n)

        clade_rows = []
        for peak_pos_norm, delta_val in peaks:
            if delta_val <= 0: continue

            # Map peaks back to raw sequences of the active target clade
            for seq_id in nov_ids:
                if seq_id not in all_seqs: continue
                seq = all_seqs[seq_id]
                if len(seq) < 10: continue

                raw_pos = int(peak_pos_norm / target_len * len(seq))
                raw_pos = min(raw_pos, len(seq) - 1)
                kmer = _extract_kmer_at_position(seq, raw_pos, k=kmer_size)

                clade_rows.append({
                    "kmer": kmer,
                    "position": raw_pos,
                    "normalized_position": peak_pos_norm,
                    "attention_delta": delta_val,
                    "source_clade": focal_clade,
                    "representative_seq_id": seq_id,
                })

        if clade_rows:
            df_clade = pd.DataFrame(clade_rows)
            kmer_summary = (
                df_clade.groupby("kmer")
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
            kmer_summary["source_clade"] = focal_clade
            if standard_clade is not None and novel_clade is not None:
                kmer_summary["reference_clade"] = standard_clade
            else:
                kmer_summary["reference_clade"] = "ALL_OTHERS"
            all_discovered_motifs.append(kmer_summary)

    if not all_discovered_motifs:
        print("[discover] No significant attention peaks found across all clades.", file=sys.stderr)
        return None

    # Combine all clade findings and export
    kmer_summary = pd.concat(all_discovered_motifs, ignore_index=True)
    
    os.makedirs(summary_dir, exist_ok=True)
    kmer_summary.to_csv(out_fp, sep="\t", index=False)
    print(f"[discover] Wrote discovered motifs: {out_fp} ({len(kmer_summary)} k-mers)")

    return kmer_summary

"""
discover.py — Unsupervised motif discovery via comparative ESM-2 attention profiles.
"""

import os
import sys
import numpy as np
import pandas as pd
from ..utils.bio import read_fasta
from ..utils.ha import (
    bh_fdr,
    compute_ha_from_layer_profiles,
    ungapped_to_msa_column,
    write_ha_outputs,
)


def _compute_attention_profile(model, alphabet, sequence: str, device: str, n_layers: int = 4) -> np.ndarray:
    import torch

    batch_converter = alphabet.get_batch_converter()
    _, _, toks = batch_converter([("seq", sequence)])
    toks = toks.to(device)

    with torch.no_grad():
        results = model(toks, repr_layers=[], return_contacts=False)

    attentions = results["attentions"]
    total_layers = attentions.shape[1]
    start_layer = max(0, total_layers - n_layers)
    attn_avg = attentions[0, start_layer:, :, :, :].mean(dim=(0, 1))
    seq_len = len(sequence)
    profile = attn_avg[1 : seq_len + 1, 1 : seq_len + 1].sum(dim=0).cpu().numpy()

    return profile


def _attention_received_by_layer(attentions, seq_len: int) -> np.ndarray:
    return attentions[:, :, 1 : seq_len + 1, 1 : seq_len + 1].sum(dim=2).cpu().numpy()


def _align_profiles_to_fixed_length(profiles: list, target_len: int = 500) -> np.ndarray:
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


def _find_peaks_in_delta(delta: np.ndarray, top_n: int = 20, min_distance: int = 5) -> list:
    try:
        from scipy.signal import find_peaks

        peaks, _ = find_peaks(delta, distance=min_distance)
        if len(peaks) == 0:
            top_idx = np.argsort(delta)[::-1][:top_n]
            return [(int(i), float(delta[i])) for i in top_idx]
        peak_vals = delta[peaks]
        sorted_idx = np.argsort(peak_vals)[::-1][:top_n]
        return [(int(peaks[i]), float(peak_vals[i])) for i in sorted_idx]
    except ImportError:
        top_idx = np.argsort(delta)[::-1][:top_n]
        return [(int(i), float(delta[i])) for i in top_idx]


def _extract_kmer_at_position(sequence: str, position: int, k: int = 5) -> str:
    half = k // 2
    start = max(0, position - half)
    end = min(len(sequence), start + k)
    start = max(0, end - k)
    return sequence[start:end]


def _compute_ha_enrichment_for_hmm(hmm, aln_seqs, clade_df_hmm, ha_cfg, disc_cfg, attention_dir, discover_dir):
    use_clusters = sorted([c for c in clade_df_hmm["cluster_id"].dropna().unique() if c != -1])
    if len(use_clusters) < 2:
        return

    ha_fp = os.path.join(attention_dir, f"{hmm}.ha_sites.tsv")
    if not os.path.exists(ha_fp):
        return
    ha_df = pd.read_csv(ha_fp, sep="\t")
    if ha_df.empty:
        return

    by_seq = {}
    for seq_id, sub in ha_df.groupby("seq_id"):
        pos = set(sub.loc[sub["is_ha"] == 1, "pos_ungapped"].astype(int).tolist())
        by_seq[seq_id] = pos

    maps = {sid: ungapped_to_msa_column(aln) for sid, aln in aln_seqs.items()}
    n_cols = len(next(iter(aln_seqs.values()))) if aln_seqs else 0

    enrich_rows = []
    hub_rows = []
    for clade_id in use_clusters:
        clade_ids = set(clade_df_hmm.loc[clade_df_hmm["cluster_id"] == clade_id, "protein"].astype(str))
        rest_ids = set(clade_df_hmm.loc[clade_df_hmm["cluster_id"] != clade_id, "protein"].astype(str))
        clade_ids = {s for s in clade_ids if s in maps}
        rest_ids = {s for s in rest_ids if s in maps}
        if not clade_ids or not rest_ids:
            continue

        pvals = []
        local_rows = []
        for col in range(n_cols):
            clade_has = 0
            rest_has = 0
            for sid in clade_ids:
                inv = maps[sid]
                if any((u in inv and inv[u] == col) for u in by_seq.get(sid, set())):
                    clade_has += 1
            for sid in rest_ids:
                inv = maps[sid]
                if any((u in inv and inv[u] == col) for u in by_seq.get(sid, set())):
                    rest_has += 1

            f_clade = clade_has / max(1, len(clade_ids))
            f_rest = rest_has / max(1, len(rest_ids))
            delta = f_clade - f_rest

            pval = np.nan
            try:
                from scipy.stats import fisher_exact

                table = np.array(
                    [
                        [clade_has, len(clade_ids) - clade_has],
                        [rest_has, len(rest_ids) - rest_has],
                    ]
                )
                _, pval = fisher_exact(table, alternative="greater")
            except Exception:
                pass

            local_rows.append(
                {
                    "clade_id": clade_id,
                    "msa_col": col,
                    "f_clade": f_clade,
                    "f_rest": f_rest,
                    "delta": delta,
                    "pval": pval,
                }
            )
            pvals.append(pval)

        qvals = bh_fdr(pvals)
        for r, q in zip(local_rows, qvals):
            r["qval"] = q
        enrich_rows.extend(local_rows)

        delta = np.array([r["delta"] for r in local_rows], dtype=float)
        win = int(disc_cfg.get("ha_window", 9))
        if win < 1:
            win = 1
        kernel = np.ones(win) / win
        smooth = np.convolve(delta, kernel, mode="same")
        thresh = float(disc_cfg.get("ha_delta_min", 0.2))
        above = smooth >= thresh

        i = 0
        while i < len(above):
            if not above[i]:
                i += 1
                continue
            j = i
            while j + 1 < len(above) and above[j + 1]:
                j += 1
            region = list(range(i, j + 1))

            cols = local_rows[i : j + 1]
            qvals_region = [r["qval"] for r in cols if not np.isnan(r.get("qval", np.nan))]
            qmin = float(min(qvals_region)) if qvals_region else np.nan
            del_mean = float(np.mean([r["delta"] for r in cols])) if cols else np.nan

            consensus_chars = []
            gap_fracs = []
            for c in region:
                aas = [aln_seqs[sid][c] for sid in clade_ids if c < len(aln_seqs[sid])]
                nongap = [aa for aa in aas if aa != "-"]
                gap_frac = 1.0 - (len(nongap) / max(1, len(aas)))
                gap_fracs.append(gap_frac)
                if nongap:
                    vals, counts = np.unique(nongap, return_counts=True)
                    consensus_chars.append(vals[np.argmax(counts)])
                else:
                    consensus_chars.append("-")
            gap_mean = float(np.mean(gap_fracs)) if gap_fracs else 1.0
            if gap_mean <= float(disc_cfg.get("ha_gap_frac_max", 0.6)):
                hub_rows.append(
                    {
                        "clade_id": clade_id,
                        "hub_start": i,
                        "hub_end": j,
                        "width": int(j - i + 1),
                        "delta_mean": del_mean,
                        "qval_min": qmin,
                        "consensus": "".join(consensus_chars),
                        "gap_frac_mean": gap_mean,
                    }
                )
            i = j + 1

    os.makedirs(discover_dir, exist_ok=True)
    pd.DataFrame(enrich_rows).to_csv(os.path.join(discover_dir, f"{hmm}.ha_enrichment.tsv"), sep="\t", index=False)
    pd.DataFrame(hub_rows).to_csv(os.path.join(discover_dir, f"{hmm}.ha_hubs.tsv"), sep="\t", index=False)


def discover_motifs(cfg, fasta_dir, summary_dir, hmm_keep, force=False):
    import glob

    disc_cfg = cfg.get("discover", {})
    if not disc_cfg.get("enabled", False):
        return None

    kmer_size = disc_cfg.get("kmer_size", 5)
    top_n = disc_cfg.get("top_n_peaks", 20)
    n_layers = disc_cfg.get("attention_layers", 4)

    ha_cfg = cfg.get("ha", {})
    use_ha = bool(ha_cfg.get("enabled", False) and disc_cfg.get("use_ha", False))
    if use_ha and not cfg.get("embeddings", {}).get("write_full_vectors", False):
        raise SystemExit(
            "HA motif discovery requires embeddings.write_full_vectors: true to ensure attention-derived "
            "per-residue vectors are available/reproducible."
        )

    out_fp = os.path.join(summary_dir, "discovered_motifs.tsv")
    if os.path.exists(out_fp) and not force:
        print(f"[discover] Output exists: {out_fp}. Use --force to override.")
        return pd.read_csv(out_fp, sep="\t")

    clade_fp = os.path.join(summary_dir, "clade_assignment.tsv")
    if not os.path.exists(clade_fp):
        print(f"[discover] clade_assignment.tsv not found at {clade_fp}. Run embed step with clustering first.", file=sys.stderr)
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
    hmm_to_seqs = {}
    fasta_files = sorted(glob.glob(os.path.join(fasta_dir, "*.faa")))
    for fp in fasta_files:
        hmm = os.path.basename(fp).replace(".faa", "")
        if hmm_keep is not None and hmm not in hmm_keep:
            continue
        seqs = read_fasta(fp)
        seqs = {k: v.replace("*", "").replace("-", "") for k, v in seqs.items()}
        hmm_to_seqs[hmm] = seqs
        all_seqs.update(seqs)

    emb_cfg = cfg.get("embeddings", {})
    model_name = emb_cfg.get("model", "esm2_t33_650M_UR50D")
    device = emb_cfg.get("device", "cuda")

    print(f"[discover] Loading ESM model {model_name} for attention profile extraction...")
    try:
        import torch
        import esm

        model, alphabet = esm.pretrained.load_model_and_alphabet(model_name)
        model.eval()
        model = model.to(device)
    except Exception as e:
        print(f"[discover] Failed to load ESM model: {e}", file=sys.stderr)
        return None

    print("[discover] Pre-computing attention profiles for all sequences...")
    all_profiles = {}
    all_layer_profiles = {}
    seq_ha_scores = {}
    seq_ha_mask = {}
    layer_start = None
    layer_end = None

    for seq_id, seq in all_seqs.items():
        if len(seq) < 10:
            continue
        try:
            batch_converter = alphabet.get_batch_converter()
            _, _, toks = batch_converter([("seq", seq)])
            toks = toks.to(device)
            with torch.no_grad():
                res = model(toks, repr_layers=[], return_contacts=False)
            attentions = res["attentions"][0]
            total_layers = attentions.shape[0]
            start_layer = max(0, total_layers - n_layers)
            attn_avg = attentions[start_layer:, :, :, :].mean(dim=(0, 1))
            profile = attn_avg[1 : len(seq) + 1, 1 : len(seq) + 1].sum(dim=0).cpu().numpy()
            all_profiles[seq_id] = profile

            if use_ha:
                layer_by_pos = _attention_received_by_layer(attentions, len(seq))
                all_layer_profiles[seq_id] = layer_by_pos
                h_scores, h_mask, ls, le = compute_ha_from_layer_profiles(layer_by_pos, ha_cfg)
                seq_ha_scores[seq_id] = h_scores
                seq_ha_mask[seq_id] = h_mask
                layer_start, layer_end = ls, le
        except Exception:
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
                if seq_id not in all_seqs:
                    continue
                seq = all_seqs[seq_id]
                if len(seq) < 10: continue

                raw_pos = int(peak_pos_norm / target_len * len(seq))
                raw_pos = min(raw_pos, len(seq) - 1)
                kmer = _extract_kmer_at_position(seq, raw_pos, k=kmer_size)
                clade_rows.append(
                    {
                        "kmer": kmer,
                        "position": raw_pos,
                        "normalized_position": peak_pos_norm,
                        "attention_delta": delta_val,
                        "source_clade": novel_clade,
                        "representative_seq_id": seq_id,
                    }
                )

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

    kmer_summary = pd.concat(all_discovered_motifs, ignore_index=True)
    os.makedirs(summary_dir, exist_ok=True)
    kmer_summary.to_csv(out_fp, sep="\t", index=False)
    print(f"[discover] Wrote discovered motifs: {out_fp} ({len(kmer_summary)} k-mers)")

    if use_ha:
        attention_dir = os.path.join(os.path.dirname(summary_dir), "attention")
        discover_dir = os.path.join(os.path.dirname(summary_dir), "discover")
        align_clipkit_dir = os.path.join(os.path.dirname(summary_dir), "alignments_clipkit")
        align_dir = os.path.join(os.path.dirname(summary_dir), "alignments")

        for hmm, seqs in hmm_to_seqs.items():
            hmm_scores = {sid: seq_ha_scores[sid] for sid in seqs if sid in seq_ha_scores}
            hmm_masks = {sid: seq_ha_mask[sid] for sid in seqs if sid in seq_ha_mask}
            if hmm_scores:
                write_ha_outputs(
                    attention_dir,
                    summary_dir,
                    hmm,
                    {k: v for k, v in seqs.items() if k in hmm_scores},
                    hmm_scores,
                    hmm_masks,
                    layer_start,
                    layer_end,
                    ha_cfg,
                )

            clip_fp = os.path.join(align_clipkit_dir, f"{hmm}.clipkit.faa")
            afa_fp = os.path.join(align_dir, f"{hmm}.afa")
            mafft_fp = os.path.join(align_dir, f"{hmm}.mafft.fasta")
            aln_fp = clip_fp if os.path.exists(clip_fp) else mafft_fp if os.path.exists(mafft_fp) else afa_fp
            if not os.path.exists(aln_fp):
                continue
            aln_seqs = read_fasta(aln_fp)

            clade_sub = clade_df[clade_df["hmm"] == hmm] if "hmm" in clade_df.columns else clade_df
            _compute_ha_enrichment_for_hmm(
                hmm,
                aln_seqs,
                clade_sub,
                ha_cfg,
                disc_cfg,
                attention_dir,
                discover_dir,
            )

    return kmer_summary

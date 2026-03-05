"""
discover.py — Unsupervised motif discovery via comparative ESM-2 attention profiles.
"""

import os
import sys
from collections import Counter
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
    by_layer = attentions[:, :, 1 : seq_len + 1, 1 : seq_len + 1].sum(dim=2).mean(dim=1)
    return by_layer.cpu().numpy()


def _resolve_seq_ids(protein_ids, available_ids):
    suffix_lookup = {}
    for seq_id in available_ids:
        if "|" in seq_id:
            suffix_lookup.setdefault(seq_id.split("|", 1)[1], seq_id)
    resolved = set()
    for pid in protein_ids:
        if pid in available_ids:
            resolved.add(pid)
        elif pid in suffix_lookup:
            resolved.add(suffix_lookup[pid])
    return resolved


def _aa_distribution(chars):
    aa = [c for c in chars if c != "-"]
    if not aa:
        return {}, "-", 0.0
    counts = Counter(aa)
    total = float(sum(counts.values()))
    consensus, cons_n = counts.most_common(1)[0]
    return {k: v / total for k, v in counts.items()}, consensus, cons_n / total


def _js_divergence_from_distributions(p, q):
    keys = sorted(set(p.keys()) | set(q.keys()))
    if not keys:
        return 0.0
    p_vec = np.array([p.get(k, 0.0) for k in keys], dtype=float)
    q_vec = np.array([q.get(k, 0.0) for k in keys], dtype=float)
    m_vec = 0.5 * (p_vec + q_vec)

    def _kl(a, b):
        mask = (a > 0) & (b > 0)
        if not np.any(mask):
            return 0.0
        return float(np.sum(a[mask] * np.log2(a[mask] / b[mask])))

    return 0.5 * _kl(p_vec, m_vec) + 0.5 * _kl(q_vec, m_vec)


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


def _compute_candidate_residues_for_hmm(hmm, aln_seqs, clade_df_hmm, disc_cfg, attention_dir, discover_dir):
    cand_cfg = disc_cfg.get("candidates", {})
    if not cand_cfg.get("enabled", True):
        return

    use_clusters = sorted([c for c in clade_df_hmm["cluster_id"].dropna().unique() if c != -1])
    if not use_clusters:
        return

    ha_fp = os.path.join(attention_dir, f"{hmm}.ha_sites.tsv")
    if not os.path.exists(ha_fp):
        return
    ha_df = pd.read_csv(ha_fp, sep="\t")
    if ha_df.empty:
        return

    by_seq = {}
    for seq_id, sub in ha_df.groupby("seq_id"):
        by_seq[seq_id] = set(sub.loc[sub["is_ha"] == 1, "pos_ungapped"].astype(int).tolist())

    maps = {sid: ungapped_to_msa_column(aln) for sid, aln in aln_seqs.items()}
    n_cols = len(next(iter(aln_seqs.values()))) if aln_seqs else 0

    min_delta_ha = float(cand_cfg.get("min_delta_ha", 0.1))
    min_cons_frac = float(cand_cfg.get("min_cons_frac", 0.6))
    max_gap_frac = float(cand_cfg.get("max_gap_frac", 0.6))
    w_delta = float(cand_cfg.get("w_delta_ha", 1.0))
    w_shift = float(cand_cfg.get("w_aa_shift", 1.0))
    w_js = float(cand_cfg.get("w_js", 0.5))
    top_n = int(cand_cfg.get("top_n", 200))

    all_rows = []
    region_rows = []
    for clade_id in use_clusters:
        raw_clade = set(clade_df_hmm.loc[clade_df_hmm["cluster_id"] == clade_id, "protein"].astype(str))
        raw_rest = set(clade_df_hmm.loc[clade_df_hmm["cluster_id"] != clade_id, "protein"].astype(str))
        clade_ids = _resolve_seq_ids(raw_clade, set(aln_seqs.keys()))
        rest_ids = _resolve_seq_ids(raw_rest, set(aln_seqs.keys()))
        if not clade_ids or not rest_ids:
            continue

        pvals = []
        clade_rows = []
        for col in range(n_cols):
            clade_has = 0
            rest_has = 0
            for sid in clade_ids:
                inv = maps.get(sid, {})
                if any((u in inv and inv[u] == col) for u in by_seq.get(sid, set())):
                    clade_has += 1
            for sid in rest_ids:
                inv = maps.get(sid, {})
                if any((u in inv and inv[u] == col) for u in by_seq.get(sid, set())):
                    rest_has += 1

            f_clade = clade_has / max(1, len(clade_ids))
            f_rest = rest_has / max(1, len(rest_ids))
            delta_ha = f_clade - f_rest

            pval = np.nan
            try:
                from scipy.stats import fisher_exact

                table = np.array(
                    [[clade_has, len(clade_ids) - clade_has], [rest_has, len(rest_ids) - rest_has]]
                )
                _, pval = fisher_exact(table, alternative="greater")
            except Exception:
                pass

            clade_chars = [aln_seqs[sid][col] for sid in clade_ids if col < len(aln_seqs[sid])]
            rest_chars = [aln_seqs[sid][col] for sid in rest_ids if col < len(aln_seqs[sid])]
            all_chars = clade_chars + rest_chars
            gap_frac_all = float(sum(c == "-" for c in all_chars) / max(1, len(all_chars)))

            clade_dist, clade_consensus, clade_cons_frac = _aa_distribution(clade_chars)
            rest_dist, rest_consensus, rest_cons_frac = _aa_distribution(rest_chars)
            aa_shift = int(
                clade_consensus != "-"
                and rest_consensus != "-"
                and clade_consensus != rest_consensus
                and clade_cons_frac >= min_cons_frac
                and rest_cons_frac >= min_cons_frac
            )
            js_div = _js_divergence_from_distributions(clade_dist, rest_dist)
            aa_shift_strength = max(clade_cons_frac, rest_cons_frac) if aa_shift else 0.0
            score = (w_delta * max(delta_ha, 0.0)) + (w_shift * aa_shift_strength) + (w_js * js_div)

            row = {
                "hmm": hmm,
                "clade_id": clade_id,
                "msa_col": col,
                "f_clade_ha": f_clade,
                "f_rest_ha": f_rest,
                "delta_ha": delta_ha,
                "pval": pval,
                "clade_consensus_aa": clade_consensus,
                "clade_cons_frac": clade_cons_frac,
                "rest_consensus_aa": rest_consensus,
                "rest_cons_frac": rest_cons_frac,
                "aa_shift": aa_shift,
                "js_div": js_div,
                "gap_frac_all": gap_frac_all,
                "candidate_score": score,
            }
            clade_rows.append(row)
            pvals.append(pval)

        qvals = bh_fdr(pvals)
        for row, qval in zip(clade_rows, qvals):
            row["qval"] = qval

        gated = [
            r
            for r in clade_rows
            if r["gap_frac_all"] <= max_gap_frac
            and r["delta_ha"] >= min_delta_ha
            and (r["clade_cons_frac"] >= min_cons_frac or r["rest_cons_frac"] >= min_cons_frac)
        ]
        gated.sort(key=lambda r: r["candidate_score"], reverse=True)
        all_rows.extend(gated[:top_n])

        region_threshold = min_delta_ha
        if gated:
            by_col = {r["msa_col"]: r for r in gated}
            cols = sorted(by_col)
            start = cols[0]
            prev = cols[0]
            for col in cols[1:] + [None]:
                if col is not None and col == prev + 1 and by_col[col]["candidate_score"] >= region_threshold:
                    prev = col
                    continue
                region_cols = list(range(start, prev + 1))
                region_scores = [by_col[c]["candidate_score"] for c in region_cols if c in by_col]
                if region_scores:
                    top_sites = sorted(region_cols, key=lambda c: by_col[c]["candidate_score"], reverse=True)[:5]
                    region_qvals = [by_col[c]["qval"] for c in region_cols if c in by_col and not np.isnan(by_col[c]["qval"])]
                    region_rows.append(
                        {
                            "hmm": hmm,
                            "clade_id": clade_id,
                            "start_col": start,
                            "end_col": prev,
                            "width": prev - start + 1,
                            "delta_ha_mean": float(np.mean([by_col[c]["delta_ha"] for c in region_cols if c in by_col])),
                            "qval_min": float(min(region_qvals)) if region_qvals else np.nan,
                            "js_div_mean": float(np.mean([by_col[c]["js_div"] for c in region_cols if c in by_col])),
                            "consensus_region": "".join(by_col[c]["clade_consensus_aa"] for c in region_cols if c in by_col),
                            "gap_frac_mean": float(np.mean([by_col[c]["gap_frac_all"] for c in region_cols if c in by_col])),
                            "top_sites": ",".join(str(c) for c in top_sites),
                        }
                    )
                if col is None:
                    break
                start = col
                prev = col

    os.makedirs(discover_dir, exist_ok=True)
    pd.DataFrame(all_rows).to_csv(os.path.join(discover_dir, f"{hmm}.candidate_residues.tsv"), sep="\t", index=False)
    pd.DataFrame(region_rows).to_csv(os.path.join(discover_dir, f"{hmm}.candidate_regions.tsv"), sep="\t", index=False)


def discover_motifs(cfg, fasta_dir, summary_dir, hmm_keep, force=False, clade_assign_dir=None):
    import glob

    disc_cfg = cfg.get("discover", {})
    if not disc_cfg.get("enabled", False):
        return None

    kmer_size = disc_cfg.get("kmer_size", 5)
    top_n = disc_cfg.get("top_n_peaks", 20)
    n_layers = disc_cfg.get("attention_layers", 4)
    standard_clade = disc_cfg.get("standard_clade")
    novel_clade = disc_cfg.get("novel_clade")
    cross_hmm = bool(disc_cfg.get("cross_hmm_comparison", True))

    ha_cfg = cfg.get("ha", {})
    candidates_enabled = bool(disc_cfg.get("candidates", {}).get("enabled", True))
    needs_ha = bool(ha_cfg.get("enabled", False) and (disc_cfg.get("use_ha", False) or candidates_enabled))
    use_ha = bool(ha_cfg.get("enabled", False) and disc_cfg.get("use_ha", False))
    if needs_ha and not cfg.get("embeddings", {}).get("write_full_vectors", False):
        raise SystemExit(
            "Attention-driven discovery outputs (HA/LoC enrichment, motif hubs, candidate residues) require "
            "embeddings.write_full_vectors: true so attention extraction is reproducible."
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

    # Load detected clades from per-HMM clade_assignments or global detected_clades.tsv
    detected_fp = os.path.join(summary_dir, "detected_clades.tsv")
    if clade_assign_dir and os.path.isdir(clade_assign_dir):
        from .post import load_per_hmm_clades
        clade_files = glob.glob(os.path.join(clade_assign_dir, "*.clades.tsv"))
        n_detected_clades = 0
        n_detected_proteins = 0
        for cf in clade_files:
            hmm_name = os.path.basename(cf).replace(".clades.tsv", "")
            per_hmm = load_per_hmm_clades(clade_assign_dir, hmm_name)
            if per_hmm:
                for clade_name, tips in per_hmm.items():
                    clade_to_proteins[clade_name] = set(tips)
                    n_detected_proteins += len(tips)
                    n_detected_clades += 1
        if n_detected_clades:
            print(f"[discover] Loaded {n_detected_clades} detected clades from "
                  f"clade_assignments/ ({n_detected_proteins} total proteins).")
    elif os.path.exists(detected_fp):
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
    seq_ha_scores = {}
    seq_ha_mask = {}
    seq_layer_bounds = {}
    seq_loc_layer = {}

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

            if needs_ha:
                layer_by_pos = _attention_received_by_layer(attentions, len(seq))
                h_scores, h_mask, ls, le, loc_layer = compute_ha_from_layer_profiles(layer_by_pos, ha_cfg)
                seq_ha_scores[seq_id] = h_scores
                seq_ha_mask[seq_id] = h_mask
                seq_layer_bounds[seq_id] = (ls, le)
                if loc_layer is not None:
                    seq_loc_layer[seq_id] = loc_layer
        except Exception:
            pass

    # 1-vs-All clade comparison
    target_len = 500
    all_discovered_motifs = []

    def _run_comparison(focal_clade, available_clades, clade_map, label_ref):
        """Compare one focal clade against others using attention delta peaks."""
        if focal_clade not in clade_map:
            print(f"[discover] Clade '{focal_clade}' not found in clade assignments. "
                  f"Skipping.", file=sys.stderr)
            return None

        nov_protein_ids = clade_map[focal_clade]
        std_protein_ids = set()
        for other_clade in available_clades:
            if other_clade != focal_clade:
                std_protein_ids.update(clade_map[other_clade])

        nov_ids = _resolve_seq_ids(nov_protein_ids, set(all_profiles.keys()))
        std_ids = _resolve_seq_ids(std_protein_ids, set(all_profiles.keys()))

        nov_profiles_list = [all_profiles[s] for s in nov_ids if s in all_profiles]
        std_profiles_list = [all_profiles[s] for s in std_ids if s in all_profiles]

        if not nov_profiles_list or not std_profiles_list:
            return None

        nov_aligned = _align_profiles_to_fixed_length(nov_profiles_list, target_len)
        std_aligned = _align_profiles_to_fixed_length(std_profiles_list, target_len)

        attention_delta = nov_aligned.mean(axis=0) - std_aligned.mean(axis=0)
        peaks = _find_peaks_in_delta(attention_delta, top_n=top_n)

        clade_rows = []
        for peak_pos_norm, delta_val in peaks:
            if delta_val <= 0:
                continue
            for seq_id in nov_ids:
                if seq_id not in all_seqs:
                    continue
                seq = all_seqs[seq_id]
                if len(seq) < 10:
                    continue
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

        if not clade_rows:
            return None

        df_clade = pd.DataFrame(clade_rows)
        kmer_summ = (
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
        kmer_summ["source_clade"] = focal_clade
        kmer_summ["reference_clade"] = label_ref
        return kmer_summ

    # Manual override: specific clade pair
    if standard_clade is not None and novel_clade is not None:
        print(f"[discover] Manual override: comparing novel_clade={novel_clade} "
              f"vs standard_clade={standard_clade}.")
        if standard_clade not in clade_to_proteins:
            print(f"[discover] standard_clade '{standard_clade}' not found.", file=sys.stderr)
        else:
            result = _run_comparison(novel_clade, [standard_clade, novel_clade],
                                     clade_to_proteins, standard_clade)
            if result is not None:
                all_discovered_motifs.append(result)

    elif cross_hmm:
        # Cross-HMM mode: pool all clades, compare each against all others
        print("[discover] cross_hmm_comparison=true: pooling clades across all HMMs.")
        for focal_clade in all_clades:
            result = _run_comparison(focal_clade, all_clades,
                                     clade_to_proteins, "ALL_OTHERS")
            if result is not None:
                all_discovered_motifs.append(result)

    else:
        # Per-HMM mode: compare clades only within each HMM
        print("[discover] cross_hmm_comparison=false: comparing clades within each HMM separately.")
        hmm_clade_groups = {}
        for cname, proteins in clade_to_proteins.items():
            # Determine which HMM a clade belongs to by checking prefix or clade_df
            if "|" in str(cname):
                hmm_prefix = str(cname).split("|", 1)[0]
            elif isinstance(cname, (int, float)):
                # Numeric HDBSCAN cluster — assign to HMM via clade_df
                hmm_prefix = "__hdbscan__"
            else:
                hmm_prefix = str(cname)
            hmm_clade_groups.setdefault(hmm_prefix, []).append(cname)

        # Also group HDBSCAN clusters per HMM from clade_df
        if "hmm" in clade_df.columns:
            for hmm_name, sub in clade_df.groupby("hmm"):
                for cid in sub["cluster_id"].dropna().unique():
                    if cid == -1:
                        continue
                    if cid in clade_to_proteins:
                        hmm_clade_groups.setdefault(hmm_name, []).append(cid)
            # Remove the catch-all hdbscan group since we've assigned properly
            hmm_clade_groups.pop("__hdbscan__", None)

        for hmm_prefix, hmm_clades in hmm_clade_groups.items():
            hmm_clades_unique = list(dict.fromkeys(hmm_clades))  # deduplicate, preserve order
            if len(hmm_clades_unique) < 2:
                continue
            print(f"[discover]   HMM '{hmm_prefix}': {len(hmm_clades_unique)} clades")
            for focal_clade in hmm_clades_unique:
                result = _run_comparison(
                    focal_clade, hmm_clades_unique,
                    clade_to_proteins, f"OTHERS_IN_{hmm_prefix}"
                )
                if result is not None:
                    result["hmm"] = hmm_prefix
                    all_discovered_motifs.append(result)

    if not all_discovered_motifs:
        print("[discover] No significant attention peaks found across all clades.", file=sys.stderr)
        return None

    kmer_summary = pd.concat(all_discovered_motifs, ignore_index=True)
    os.makedirs(summary_dir, exist_ok=True)
    kmer_summary.to_csv(out_fp, sep="\t", index=False)
    print(f"[discover] Wrote discovered motifs: {out_fp} ({len(kmer_summary)} k-mers)")

    if needs_ha:
        attention_dir = os.path.join(os.path.dirname(summary_dir), "attention")
        discover_dir = os.path.join(os.path.dirname(summary_dir), "discover")
        align_clipkit_dir = os.path.join(os.path.dirname(summary_dir), "alignments_clipkit")
        align_dir = os.path.join(os.path.dirname(summary_dir), "alignments")

        for hmm, seqs in hmm_to_seqs.items():
            hmm_scores = {sid: seq_ha_scores[sid] for sid in seqs if sid in seq_ha_scores}
            hmm_masks = {sid: seq_ha_mask[sid] for sid in seqs if sid in seq_ha_mask}
            if hmm_scores:
                bounds = [seq_layer_bounds.get(sid) for sid in hmm_scores if sid in seq_layer_bounds]
                layer_start = bounds[0][0] if bounds else 0
                layer_end = bounds[0][1] if bounds else 1
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
                    per_seq_loc_layer={sid: seq_loc_layer[sid] for sid in hmm_scores if sid in seq_loc_layer},
                )

            clip_fp = os.path.join(align_clipkit_dir, f"{hmm}.clipkit.faa")
            afa_fp = os.path.join(align_dir, f"{hmm}.afa")
            mafft_fp = os.path.join(align_dir, f"{hmm}.mafft.fasta")
            aln_fp = clip_fp if os.path.exists(clip_fp) else mafft_fp if os.path.exists(mafft_fp) else afa_fp
            if not os.path.exists(aln_fp):
                continue
            aln_seqs = read_fasta(aln_fp)

            clade_sub = clade_df[clade_df["hmm"] == hmm] if "hmm" in clade_df.columns else clade_df
            if use_ha:
                _compute_ha_enrichment_for_hmm(
                    hmm,
                    aln_seqs,
                    clade_sub,
                    ha_cfg,
                    disc_cfg,
                    attention_dir,
                    discover_dir,
                )
            if candidates_enabled:
                _compute_candidate_residues_for_hmm(
                    hmm,
                    aln_seqs,
                    clade_sub,
                    disc_cfg,
                    attention_dir,
                    discover_dir,
                )

    return kmer_summary

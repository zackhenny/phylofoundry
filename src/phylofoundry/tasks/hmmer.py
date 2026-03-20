import os
import pandas as pd
import pyhmmer

def _run_pyhmmscan(genomes, faa_dir, combined_hmm, cpu):
    """Run hmmscan using PyHMMER for all genomes against the combined HMM database.

    For each genome's FASTA, every sequence is queried against the combined HMM
    database.  Results mirror the columns produced by the legacy ``--domtblout``
    approach: full-sequence bitscore/e-value paired with per-domain alignment
    coordinates used to compute per-protein coverage.
    """
    rows = []

    if not os.path.exists(combined_hmm):
        return pd.DataFrame()

    with pyhmmer.plan7.HMMFile(combined_hmm) as hmm_file:
        profiles = list(hmm_file)

    if not profiles:
        return pd.DataFrame()

    for g in genomes:
        faa_path = os.path.join(faa_dir, g)
        if not os.path.exists(faa_path):
            continue

        with pyhmmer.easel.SequenceFile(faa_path, digital=True) as seq_file:
            seqs = list(seq_file)

        if not seqs:
            continue

        for tophits in pyhmmer.hmmer.hmmscan(seqs, profiles, cpus=cpu):
            protein_name = tophits.query.name
            seq_len = len(tophits.query.sequence)
            for hit in tophits:
                if not hit.included:
                    continue
                hmm_name = hit.name
                bitscore = float(hit.score)
                evalue = float(hit.evalue)
                for domain in hit.domains:
                    if domain.included:
                        ali_from = domain.alignment.target_from
                        ali_to = domain.alignment.target_to
                        coverage = (ali_to - ali_from + 1) / seq_len
                        rows.append({
                            "genome": g,
                            "protein": protein_name,
                            "hmm": hmm_name,
                            "bitscore": bitscore,
                            "evalue": evalue,
                            "coverage": coverage,
                        })

    return pd.DataFrame(rows) if rows else pd.DataFrame()


def _run_pyhmmsearch(hmm_files, hmm_dir, combined_faa, hmm_keep, cpu):
    """Run hmmsearch using PyHMMER for all HMMs against the combined FASTA.

    All HMM profiles are searched against *combined_faa* in a single PyHMMER
    call, using the requested number of CPU threads.  Protein names in
    *combined_faa* are expected to follow the ``genome~protein`` convention
    produced by the prep step.
    """
    if not os.path.exists(combined_faa):
        return pd.DataFrame()

    with pyhmmer.easel.SequenceFile(combined_faa, digital=True) as seq_file:
        seqs = list(seq_file)

    if not seqs:
        return pd.DataFrame()

    hmms = []
    for hf in hmm_files:
        hmm_name = os.path.splitext(hf)[0]
        if hmm_keep is not None and hmm_name not in hmm_keep:
            continue
        hmm_path = os.path.join(hmm_dir, hf)
        if not os.path.exists(hmm_path):
            continue
        with pyhmmer.plan7.HMMFile(hmm_path) as hmm_file:
            hmms.extend(list(hmm_file))

    if not hmms:
        return pd.DataFrame()

    rows = []
    for tophits in pyhmmer.hmmer.hmmsearch(hmms, seqs, cpus=cpu):
        query_hmm_name = tophits.query.name
        for hit in tophits:
            if not hit.included:
                continue
            raw_protein = hit.name
            bitscore = float(hit.score)
            evalue = float(hit.evalue)
            seq_len = hit.length
            for domain in hit.domains:
                if domain.included:
                    ali_from = domain.alignment.target_from
                    ali_to = domain.alignment.target_to
                    coverage = (ali_to - ali_from + 1) / seq_len
                    if "~" in raw_protein:
                        genome, protein = raw_protein.split("~", 1)
                    else:
                        genome, protein = "Unknown", raw_protein
                    rows.append({
                        "genome": genome,
                        "protein": protein,
                        "hmm": query_hmm_name,
                        "bitscore": bitscore,
                        "evalue": evalue,
                        "coverage": coverage,
                    })

    return pd.DataFrame(rows) if rows else pd.DataFrame()


def apply_filtering(df, thresholds_map, global_min_score, global_min_cov,
                    max_evalue=None, disable_bitscore_filter=False,
                    disable_coverage_filter=False):
    if df.empty:
        return df
    df = df.copy()
    if not disable_coverage_filter and global_min_cov > 0:
        df = df[df["coverage"] >= global_min_cov].copy()
    if not disable_bitscore_filter:
        df["min_score_required"] = global_min_score
        if thresholds_map:
            custom_scores = df["hmm"].map(thresholds_map)
            df["min_score_required"] = custom_scores.fillna(df["min_score_required"])
        df = df[df["bitscore"] >= df["min_score_required"]].copy()
    if max_evalue is not None:
        df = df[df["evalue"] <= max_evalue].copy()
    return df

def best_hits(df, use_evalue=False):
    """Find the single best HMM for each (genome, protein) pair.

    When use_evalue=True, rank by e-value ascending (lower is better).
    Otherwise rank by bitscore descending.
    """
    if df.empty:
        return pd.DataFrame()
    if use_evalue:
        df_sorted = df.sort_values(["genome", "protein", "evalue"], ascending=[True, True, True])
    else:
        df_sorted = df.sort_values(["genome", "protein", "bitscore"], ascending=[True, True, False])
    rows = []
    for (genome, prot), chunk in df_sorted.groupby(["genome", "protein"], sort=False):
        best = chunk.iloc[0]
        delta = float(best.bitscore)
        if len(chunk) > 1:
            delta -= float(chunk.iloc[1].bitscore)
        rows.append({
            "genome": genome,
            "protein": prot,
            "hmm": best.hmm,           # keep as "hmm" to match extract.py expectations
            "bitscore": float(best.bitscore),
            "evalue": float(best.evalue),
            "coverage": float(best.coverage),
            "delta_bitscore": float(delta)
        })
    return pd.DataFrame(rows)

def run_hmmer(cfg, genomes, faa_dir, hmm_files, hmm_dir, combined_hmm, combined_faa, outdir, summary_dir, hmmscan_dir, hmmsearch_dir, hmm_keep, force=False):
    print("\n[hmmer] Running hmmscan + hmmsearch with PyHMMER (with caching)...")

    cpu = int(cfg["resources"]["cpu"])
    filt_cfg = cfg["filtering"]
    hmmer_cfg = cfg.get("hmmer", {})
    run_scan = hmmer_cfg.get("run_scan", True)
    run_search = hmmer_cfg.get("run_search", True)

    max_evalue = filt_cfg.get("max_evalue", None)
    use_evalue = bool(filt_cfg.get("use_evalue", False))
    disable_bitscore_filter = bool(filt_cfg.get("disable_bitscore_filter", False))
    disable_coverage_filter = bool(filt_cfg.get("disable_coverage_filter", False))

    hits_scan_tsv = os.path.join(summary_dir, "hmmscan_hits.filtered.tsv")
    hits_search_tsv = os.path.join(summary_dir, "hmmsearch_hits.filtered.tsv")
    best_hits_tsv = os.path.join(summary_dir, "best_hits.competitive.tsv")

    # Load custom score thresholds if provided
    thresholds_map = {}
    scores_tsv = filt_cfg.get("scores_tsv", None)
    if scores_tsv and os.path.exists(scores_tsv):
        print(f"  Loading custom thresholds from {scores_tsv}...")
        with open(scores_tsv) as f:
            for line in f:
                if line.strip().startswith("#") or not line.strip():
                    continue
                parts = line.strip().split()
                if len(parts) >= 2:
                    thresholds_map[parts[0]] = float(parts[1])

    scan_df = pd.DataFrame()
    search_df = pd.DataFrame()
    best_df = pd.DataFrame()

    # hmmscan
    if os.path.exists(hits_scan_tsv) and not force:
        scan_df = pd.read_csv(hits_scan_tsv, sep="\t")
    elif run_scan:
        print("  [hmmscan] Running with PyHMMER...")
        scan_df = _run_pyhmmscan(genomes, faa_dir, combined_hmm, cpu)
        scan_df = apply_filtering(
            scan_df, thresholds_map,
            float(filt_cfg["global_min_score"]),
            float(filt_cfg["min_coverage"]),
            max_evalue=max_evalue,
            disable_bitscore_filter=disable_bitscore_filter,
            disable_coverage_filter=disable_coverage_filter,
        )
        if not scan_df.empty:
            scan_df.to_csv(hits_scan_tsv, sep="\t", index=False)

    # hmmsearch
    if os.path.exists(hits_search_tsv) and not force:
        search_df = pd.read_csv(hits_search_tsv, sep="\t")
    elif run_search:
        print("  [hmmsearch] Running with PyHMMER...")
        search_df = _run_pyhmmsearch(hmm_files, hmm_dir, combined_faa, hmm_keep, cpu)
        search_df = apply_filtering(
            search_df, thresholds_map,
            float(filt_cfg["global_min_score"]),
            float(filt_cfg["min_coverage"]),
            max_evalue=max_evalue,
            disable_bitscore_filter=disable_bitscore_filter,
            disable_coverage_filter=disable_coverage_filter,
        )
        if not search_df.empty:
            search_df.to_csv(hits_search_tsv, sep="\t", index=False)

    # Competitive best per (genome, protein): prefer hmmscan if available else hmmsearch
    if os.path.exists(best_hits_tsv) and not force:
        best_df = pd.read_csv(best_hits_tsv, sep="\t")
    else:
        if not scan_df.empty:
            best_df = best_hits(scan_df, use_evalue=use_evalue)
            best_df["source"] = "hmmscan"
        elif not search_df.empty:
            best_df = best_hits(search_df, use_evalue=use_evalue)
            best_df["source"] = "hmmsearch"
        else:
            best_df = pd.DataFrame()

    # Optional Taxonomy Integration
    gtdb_dir = cfg["inputs"].get("gtdb_dir")
    tax_file = cfg["inputs"].get("taxonomy_file")
    if (gtdb_dir or tax_file) and not best_df.empty:
        from .post import _load_taxonomy
        tax_map = _load_taxonomy(gtdb_dir, tax_file)
        if tax_map:
            from ..utils.helpers import normalize_genome_id
            
            def lookup_tax(g_filename):
                return tax_map.get(normalize_genome_id(str(g_filename)), "Unknown")
            
            # Map to all generated output tables
            if not best_df.empty and "genome" in best_df.columns:
                best_df["taxonomy"] = best_df["genome"].apply(lookup_tax)
            if not scan_df.empty and "genome" in scan_df.columns:
                scan_df["taxonomy"] = scan_df["genome"].apply(lookup_tax)
                scan_df.to_csv(hits_scan_tsv, sep="\t", index=False)
            if not search_df.empty and "genome" in search_df.columns:
                search_df["taxonomy"] = search_df["genome"].apply(lookup_tax)
                search_df.to_csv(hits_search_tsv, sep="\t", index=False)
    
    if not best_df.empty:
        best_df.to_csv(best_hits_tsv, sep="\t", index=False)
    
    return scan_df, search_df, best_df

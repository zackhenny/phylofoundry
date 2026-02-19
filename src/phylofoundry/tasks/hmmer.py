import os
import subprocess
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
from ..utils.helpers import run_cmd

def parse_tblout_fast(tbl_file):
    """Parse HMMER --domtblout format. Uses fixed-width aware parsing."""
    if not os.path.exists(tbl_file) or os.path.getsize(tbl_file) == 0:
        return pd.DataFrame()
    try:
        df = pd.read_csv(
            tbl_file,
            sep=r"\s+",  # replaces deprecated delim_whitespace
            comment="#",
            header=None,
            usecols=[0, 2, 3, 5, 6, 7, 17, 18],
            names=["target", "tlen", "query", "qlen", "evalue", "bitscore", "ali_from", "ali_to"],
            engine="python",  # python engine required for regex sep
        )
        return df
    except Exception:
        return pd.DataFrame()

def apply_filtering(df, thresholds_map, global_min_score, global_min_cov):
    if df.empty:
        return df
    df = df.copy()
    if global_min_cov > 0:
        df = df[df["coverage"] >= global_min_cov].copy()
    df["min_score_required"] = global_min_score
    if thresholds_map:
        custom_scores = df["hmm"].map(thresholds_map)
        df["min_score_required"] = custom_scores.fillna(df["min_score_required"])
    df = df[df["bitscore"] >= df["min_score_required"]].copy()
    return df

def best_hits(df):
    """Find the single best HMM for each (genome, protein) pair by bitscore."""
    if df.empty:
        return pd.DataFrame()
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

def worker_hmmscan(args_pack):
    cmd, tbl_path, g_name, keep_tbl = args_pack
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        df = parse_tblout_fast(tbl_path)

        if (not keep_tbl) and os.path.exists(tbl_path):
            try:
                os.remove(tbl_path)
            except OSError:
                pass

        if df.empty:
            return None

        df = df.rename(columns={"target": "hmm", "query": "protein"})
        df["coverage"] = (df["ali_to"] - df["ali_from"] + 1) / df["qlen"]
        df["genome"] = g_name
        return df[["genome", "protein", "hmm", "bitscore", "evalue", "coverage"]]
    except Exception:
        return None

def worker_hmmsearch(args_pack):
    cmd, tbl_path, hmm_name, keep_tbl = args_pack
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        df = parse_tblout_fast(tbl_path)

        if (not keep_tbl) and os.path.exists(tbl_path):
            try:
                os.remove(tbl_path)
            except OSError:
                pass

        if df.empty:
            return None

        df = df.rename(columns={"target": "protein", "query": "hmm"})
        df["coverage"] = (df["ali_to"] - df["ali_from"] + 1) / df["tlen"]

        df["genome"] = df["protein"].apply(lambda x: x.split("~")[0] if "~" in x else "Unknown")
        df["protein"] = df["protein"].apply(lambda x: x.split("~", 1)[1] if "~" in x else x)

        return df[["genome", "protein", "hmm", "bitscore", "evalue", "coverage"]]
    except Exception:
        return None

def run_hmmer(cfg, genomes, faa_dir, hmm_files, hmm_dir, combined_hmm, combined_faa, outdir, summary_dir, hmmscan_dir, hmmsearch_dir, hmm_keep, force=False):
    print("\n[hmmer] Running hmmscan + hmmsearch (with caching)...")
    
    cpu = int(cfg["resources"]["cpu"])
    filt_cfg = cfg["filtering"]
    hmmer_cfg = cfg.get("hmmer", {})
    run_scan = hmmer_cfg.get("run_scan", True)
    run_search = hmmer_cfg.get("run_search", True)
    
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

    # hmmscan per genome
    scan_tasks = []
    if run_scan:
        for g in genomes:
            tbl = os.path.join(hmmscan_dir, g + ".tbl")
            if os.path.exists(tbl) and not force:
                continue
            cmd = ["hmmscan", "--cpu", "1", "--domtblout", tbl, combined_hmm, os.path.join(faa_dir, g)]
            scan_tasks.append((cmd, tbl, g, bool(filt_cfg.get("keep_tbl", False))))

    if scan_tasks:
        all_hits_scan = []
        with ProcessPoolExecutor(max_workers=cpu) as exe:
            for res in exe.map(worker_hmmscan, scan_tasks):
                if res is not None:
                    all_hits_scan.append(res)
        if all_hits_scan:
            scan_df = pd.concat(all_hits_scan, ignore_index=True)

    if os.path.exists(hits_scan_tsv) and not force:
        scan_df = pd.read_csv(hits_scan_tsv, sep="\t")
    else:
        # parse existing tables if present
        if scan_df.empty:
            dfs = []
            for g in genomes:
                tbl = os.path.join(hmmscan_dir, g + ".tbl")
                if os.path.exists(tbl):
                    df = parse_tblout_fast(tbl)
                    if not df.empty:
                        df = df.rename(columns={"target": "hmm", "query": "protein"})
                        df["coverage"] = (df["ali_to"] - df["ali_from"] + 1) / df["qlen"]
                        df["genome"] = g
                        dfs.append(df[["genome", "protein", "hmm", "bitscore", "evalue", "coverage"]])
            if dfs:
                scan_df = pd.concat(dfs, ignore_index=True)

        scan_df = apply_filtering(
            scan_df, thresholds_map,
            float(filt_cfg["global_min_score"]),
            float(filt_cfg["min_coverage"])
        )
        if not scan_df.empty:
            scan_df.to_csv(hits_scan_tsv, sep="\t", index=False)

    # hmmsearch per HMM file
    search_tasks = []
    if run_search:
        for hf in hmm_files:
            hmm_name = os.path.splitext(hf)[0]
            if hmm_keep is not None and hmm_name not in hmm_keep:
                continue
            tbl = os.path.join(hmmsearch_dir, f"{hmm_name}_combined.tbl")
            if os.path.exists(tbl) and not force:
                continue
            cmd = ["hmmsearch", "--cpu", "1", "--domtblout", tbl, os.path.join(hmm_dir, hf), combined_faa]
            search_tasks.append((cmd, tbl, hmm_name, bool(filt_cfg.get("keep_tbl", False))))

    if search_tasks:
        all_hits_search = []
        with ProcessPoolExecutor(max_workers=cpu) as exe:
            for res in exe.map(worker_hmmsearch, search_tasks):
                if res is not None:
                    all_hits_search.append(res)
        if all_hits_search:
            search_df = pd.concat(all_hits_search, ignore_index=True)

    if os.path.exists(hits_search_tsv) and not force:
        search_df = pd.read_csv(hits_search_tsv, sep="\t")
    else:
        if search_df.empty:
            dfs = []
            for hf in hmm_files:
                hmm_name = os.path.splitext(hf)[0]
                if hmm_keep is not None and hmm_name not in hmm_keep:
                    continue
                tbl = os.path.join(hmmsearch_dir, f"{hmm_name}_combined.tbl")
                if os.path.exists(tbl):
                    df = parse_tblout_fast(tbl)
                    if not df.empty:
                        df = df.rename(columns={"target": "protein", "query": "hmm"})
                        df["coverage"] = (df["ali_to"] - df["ali_from"] + 1) / df["tlen"]
                        df["genome"] = df["protein"].apply(lambda x: x.split("~")[0] if "~" in x else "Unknown")
                        df["protein"] = df["protein"].apply(lambda x: x.split("~", 1)[1] if "~" in x else x)
                        dfs.append(df[["genome", "protein", "hmm", "bitscore", "evalue", "coverage"]])
            if dfs:
                search_df = pd.concat(dfs, ignore_index=True)

        search_df = apply_filtering(
            search_df, thresholds_map,
            float(filt_cfg["global_min_score"]),
            float(filt_cfg["min_coverage"])
        )
        if not search_df.empty:
            search_df.to_csv(hits_search_tsv, sep="\t", index=False)

    # Competitive best per (genome,protein): prefer hmmscan if available else hmmsearch
    if os.path.exists(best_hits_tsv) and not force:
        best_df = pd.read_csv(best_hits_tsv, sep="\t")
    else:
        if not scan_df.empty:
            best_df = best_hits(scan_df)
            best_df["source"] = "hmmscan"
        elif not search_df.empty:
            best_df = best_hits(search_df)
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

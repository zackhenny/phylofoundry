import os
import shutil
import subprocess
import glob
from collections import defaultdict
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from concurrent.futures import ProcessPoolExecutor

from ..utils.bio import read_fasta, write_fasta
from ..utils.helpers import run_cmd, normalize_genome_id, safe_mkdir

def load_genbank_neighborhood(gbk_path, protein_id, window_genes, protein_id_fields, gene_label_fields):
    """
    Scan a GenBank file for a specific protein_id (or locus_tag) and extract
    the neighborhood of +/- window_genes.
    Returns: list of dicts (feature info), focal_feature
    """
    neighborhood = []
    focal_feature = None
    
    try:
        # We assume one record per file usually, or iterate all
        for record in SeqIO.parse(gbk_path, "genbank"):
            # First pass: find the focal CDS index
            target_idx = -1
            cds_features = [f for f in record.features if f.type == "CDS"]
            
            for i, f in enumerate(cds_features):
                qualifiers = f.qualifiers
                # Check all possible ID fields
                found = False
                for field in protein_id_fields:
                    if field in qualifiers:
                        # Some fields are lists
                        vals = qualifiers[field]
                        if isinstance(vals, str): vals = [vals]
                        if protein_id in vals:
                            found = True
                            break
                if found:
                    target_idx = i
                    focal_feature = f
                    break
            
            if target_idx != -1:
                # Extract window
                start_i = max(0, target_idx - window_genes)
                end_i = min(len(cds_features), target_idx + window_genes + 1)
                
                for f in cds_features[start_i:end_i]:
                    # Extract robust label
                    label = "NA"
                    for field in gene_label_fields:
                        if field in f.qualifiers:
                            label = f.qualifiers[field][0]
                            break
                    
                    # Extract protein ID
                    pid = "NA"
                    for field in protein_id_fields:
                        if field in f.qualifiers:
                            pid = f.qualifiers[field][0]
                            break

                    # Extract translation
                    translation = ""
                    if "translation" in f.qualifiers:
                        translation = f.qualifiers["translation"][0]

                    neighborhood.append({
                        "contig": record.id,
                        "start": int(f.location.start),
                        "end": int(f.location.end),
                        "strand": f.location.strand,
                        "label": label,
                        "protein_id": pid,
                        "translation": translation,
                        "is_focal": (f is focal_feature)
                    })
                
                return neighborhood
                
    except Exception:
        pass
    return []

def parse_gff_line(line):
    parts = line.strip().split("\t")
    if len(parts) < 9: return None
    attrs = {}
    for x in parts[8].split(";"):
        if "=" in x:
            k, v = x.split("=", 1)
            attrs[k] = v
    return {
        "contig": parts[0],
        "type": parts[2],
        "start": int(parts[3]),
        "end": int(parts[4]),
        "strand": parts[6],
        "attrs": attrs
    }

def load_gff_neighborhood(gff_path, protein_id, window_genes, protein_id_fields, gene_label_fields, fasta_path=None):
    """
    Scan a GFF file for protein_id. If fasta_path provided, extract protein sequences.
    """
    # 1. Load all CDS from GFF
    cds_list = []
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#"): continue
            feat = parse_gff_line(line)
            if feat and feat["type"] == "CDS":
                cds_list.append(feat)

    # 2. Find focal
    target_idx = -1
    for i, f in enumerate(cds_list):
        for field in protein_id_fields:
            if field in f["attrs"] and f["attrs"][field] == protein_id:
                target_idx = i
                break
        if target_idx != -1: break
    
    if target_idx == -1:
        return []

    # 3. Extract window
    start_i = max(0, target_idx - window_genes)
    end_i = min(len(cds_list), target_idx + window_genes + 1)
    
    # Optional: Load genome sequence if available to translate?
    # For now, we assume user might provide protein FASTA separately or we rely on 'translation' in GFF? 
    # Standard GFF doesn't have translation. 
    # If fasta_path corresponds to protein fasta, we can try to look up.
    # If fasta_path is genome, we'd need to translate.
    # To keep it simple and robust matching User Request "GFF + corresponding genome FASTA(s) OR protein FASTA(s)":
    # We will assume if protein sequences are needed, we look them up in a provided dict or just skip translation for plotting.
    
    neighborhood = []
    for current_idx, f in enumerate(cds_list[start_i:end_i], start=start_i):
        label = "NA"
        for field in gene_label_fields:
            if field in f["attrs"]:
                label = f["attrs"][field]
                break
        
        pid = "NA"
        for field in protein_id_fields:
            if field in f["attrs"]:
                pid = f["attrs"][field]
                break

        neighborhood.append({
            "contig": f["contig"],
            "start": f["start"],
            "end": f["end"],
            "strand": 1 if f["strand"] == "+" else -1,
            "label": label,
            "protein_id": pid,
            "translation": "", # Filled later if possible
            "is_focal": (current_idx == target_idx)
        })
    return neighborhood

def make_diamond_db(input_faa, db_prefix, diamond_bin="diamond"):
    cmd = [diamond_bin, "makedb", "--in", input_faa, "--db", db_prefix, "--quiet"]
    run_cmd(cmd, quiet=True)

def run_diamond_blastp(query_faa, db_prefix, out_tsv, diamond_bin="diamond", threads=1):
    # outfmt 6: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
    cmd = [
        diamond_bin, "blastp",
        "--query", query_faa,
        "--db", db_prefix,
        "--out", out_tsv,
        "--outfmt", "6", "qseqid", "sseqid", "pident", "bitscore", "evalue",
        "--threads", str(threads),
        "--quiet"
    ]
    run_cmd(cmd, quiet=True)

def run_synteny(cfg, synteny_dir, tree_dir, scan_df, search_df, hmm_keep, force=False):
    print("\n[synteny] Extracting neighborhoods and plotting synteny...")
    
    syn_cfg = cfg.get("synteny", {})
    if not syn_cfg.get("enabled", False):
        return

    # Check dependencies
    try:
        import pygenomeviz
        from pygenomeviz import GenomeViz
    except ImportError:
        print("WARNING: pygenomeviz not installed. Skipping synteny plots.")
        print("Install with: pip install pygenomeviz")
        return

    gbk_dir = syn_cfg.get("gbk_dir")
    gff_dir = syn_cfg.get("gff_dir")
    
    # We need at least one source
    if not gbk_dir and not gff_dir:
        print("WARNING: No gbk_dir or gff_dir provided for synteny. Skipping.")
        return

    # Aggregate hits
    # We use competitive best hits usually, but user might want all
    # "use pipeline hit tables"
    # To be robust, let's load best_hits.competitive.tsv if strictly preferring best
    # Or merge the passed scan_df/search_df
    
    # Re-load best hits if prefer_best_hit is True
    summary_dir = os.path.join(cfg["output"]["outdir"], "summary")
    best_hits_tsv = os.path.join(summary_dir, "best_hits.competitive.tsv")
    
    hits = pd.DataFrame()
    if syn_cfg.get("prefer_best_hit", True) and os.path.exists(best_hits_tsv):
        hits = pd.read_csv(best_hits_tsv, sep="\t")
    else:
        # Fallback to whatever was passed (e.g. all hits)
        dfs = []
        if not scan_df.empty: dfs.append(scan_df)
        if not search_df.empty: dfs.append(search_df)
        if dfs:
            hits = pd.concat(dfs, ignore_index=True)
            # Dedup if needed
            hits = hits.drop_duplicates(subset=["genome", "protein", "hmm"])

    if hits.empty:
        print("No hits found for synteny analysis.")
        return

    # Filter hmm_keep
    if hmm_keep:
        hits = hits[hits["hmm"].isin(hmm_keep)]

    # Group by HMM
    for hmm, group in hits.groupby("hmm"):
        hmm_out_dir = os.path.join(synteny_dir, hmm)
        safe_mkdir(hmm_out_dir)
        
        pdf_out = os.path.join(hmm_out_dir, f"synteny.{hmm}.pdf")
        if os.path.exists(pdf_out) and not force:
            continue

        print(f"  Processing HMM: {hmm} ({len(group)} hits)...")
        
        # Deduplicate per genome if requested
        if syn_cfg.get("dedup_by_genome", True):
            # Sort by bitscore desc to keep best
            group = group.sort_values("bitscore", ascending=False).groupby("genome").head(1)
            # Limit total hits
            max_hits = int(syn_cfg.get("max_hits_per_hmm", 50))
            if len(group) > max_hits:
                group = group.head(max_hits)

        neighborhoods = [] # list of (genome, track_name, [features])
        all_prots = {}     # uid -> sequence

        # Extract neighborhoods
        window = int(syn_cfg.get("window_genes", 10))
        fields_id = syn_cfg.get("protein_id_field", ["ID", "protein_id", "locus_tag"])
        fields_label = syn_cfg.get("gene_label_field", ["gene", "product", "Name", "locus_tag"])

        for _, r in group.iterrows():
            genome = r["genome"]
            protein = r["protein"]
            
            # Locate input file
            # Try GBK
            found_feats = []
            if gbk_dir:
                # Try genome.gbk, genome.gbff, etc.
                base = normalize_genome_id(genome)
                cand = glob.glob(os.path.join(gbk_dir, f"{base}*.gb*"))
                if cand:
                    found_feats = load_genbank_neighborhood(cand[0], protein, window, fields_id, fields_label)
            
            # Try GFF if no GBK or not found
            if not found_feats and gff_dir:
                base = normalize_genome_id(genome)
                cand = glob.glob(os.path.join(gff_dir, f"{base}*.gff*"))
                if cand:
                    found_feats = load_gff_neighborhood(cand[0], protein, window, fields_id, fields_label)
                    # If GFF, we might need external protein sequences for similarity
                    # We can try looking up in the global proteome inputs if available
                    # For now, if translation is missing, similarity check will just skip or fail gracefully
            
            if not found_feats:
                continue

            # Assign unique IDs for all proteins in neighborhood
            # {genome}|{contig}|{locus}|{idx}
            track_feats = []
            for i, f in enumerate(found_feats):
                # Construct a stable UID
                # Note: protein_id might not be unique if paralogs, so use index
                uid = f"{genome}|{f['contig']}|{f['protein_id']}|{i}"
                f["uid"] = uid
                if f["translation"]:
                    all_prots[uid] = f["translation"]
                track_feats.append(f)
            
            neighborhoods.append((genome, track_feats))

        if not neighborhoods:
            print(f"    No neighborhoods extracted for {hmm}.")
            continue

        # Write neighborhood proteins
        prot_faa = os.path.join(hmm_out_dir, "neighborhood_proteins.faa")
        write_fasta(prot_faa, all_prots)

        # compute similarity
        sim_cfg = syn_cfg.get("similarity", {})
        sim_tsv = os.path.join(hmm_out_dir, "gene_similarity.tsv")
        
        if all_prots and (not os.path.exists(sim_tsv) or force):
            diamond_bin = sim_cfg.get("diamond_bin", "diamond")
            db_prefix = os.path.join(hmm_out_dir, "diamond_db")
            make_diamond_db(prot_faa, db_prefix, diamond_bin)
            run_diamond_blastp(prot_faa, db_prefix, sim_tsv, diamond_bin, threads=int(cfg["resources"]["cpu"]))
            
            # Cleanup DB
            for f in glob.glob(f"{db_prefix}.*"):
                try: os.remove(f)
                except: pass

        # Load links
        links = []
        if os.path.exists(sim_tsv):
            try:
                sdf = pd.read_csv(sim_tsv, sep="\t", names=["q", "s", "pident", "bitscore", "evalue"])
                # Filter
                sdf = sdf[sdf["q"] != sdf["s"]]
                sdf = sdf[sdf["pident"] >= float(sim_cfg.get("min_identity", 30))]
                sdf = sdf[sdf["bitscore"] >= float(sim_cfg.get("min_bitscore", 50))]
                
                # Make simple list of (q, s, pident)
                for _, row in sdf.iterrows():
                    links.append((row["q"], row["s"], row["pident"]))
            except Exception:
                pass

        # Tree ordering
        ordered_genomes = [n[0] for n in neighborhoods]
        if syn_cfg.get("include_tree", True):
             tree_fp = os.path.join(tree_dir, f"{hmm}.treefile")
             if os.path.exists(tree_fp):
                 try:
                     from Bio import Phylo
                     t = Phylo.read(tree_fp, "newick")
                     if syn_cfg.get("tree_order") == "ladderize":
                         t.ladderize()
                     # Get leaf order
                     leaves = [term.name for term in t.get_terminals()]
                     # Map leaves to genomes
                     # Leaf names might be "Genome|Protein" or just "Genome"
                     # We need to match to our genome keys
                     leaf_to_genome = {}
                     for leak in leaves:
                         # simple heuristic: splits
                         g_cand = leak.split("|")[0]
                         leaf_to_genome[leak] = g_cand
                     
                     # Reorder
                     ordered_temp = []
                     seen = set()
                     for leaf in leaves:
                         g = leaf_to_genome.get(leaf)
                         # Search if this g is in our neighborhoods
                         # This is O(N*M), but N is small (50)
                         for i, (ng, nfeats) in enumerate(neighborhoods):
                             # Fuzzy match or exact?
                             # normalize_genome_id is useful here usually
                             if ng == g and i not in seen:
                                 ordered_temp.append(neighborhoods[i])
                                 seen.add(i)
                                 break
                     
                     # Append leftovers
                     for i, n in enumerate(neighborhoods):
                         if i not in seen:
                             ordered_temp.append(n)
                     
                     if ordered_temp:
                         ordered_genomes = [n[0] for n in ordered_temp]
                         neighborhoods = ordered_temp
                 except Exception as e:
                     print(f"    Tree ordering failed: {e}")


        
        # Refactor slightly to map UIDs to feature objects
        # Re-build the GenomeViz loop to store feature objects
        gv = GenomeViz(
            fig_width=syn_cfg.get("plot_width", 14),
            fig_track_height=syn_cfg.get("plot_height_per_track", 0.35),
        )
        
        uid_to_feat = {}
        
        for genome, feats in neighborhoods:
            if not feats: continue
            min_start = min(f["start"] for f in feats)
            max_end = max(f["end"] for f in feats)
            size = max_end - min_start + 1
            
            track = gv.add_feature_track(genome, size)
            for f in feats:
                color = "tomato" if f["is_focal"] else "skyblue"
                s_val = 1 if str(f["strand"]) in ["1", "+"] else -1
                
                # Use label as display, but store UID for linking
                feature_obj = track.add_feature(
                    start=f["start"] - min_start,
                    end=f["end"] - min_start,
                    strand=s_val,
                    label=f["label"],
                    facecolor=color,
                    plotstyle="arrow"
                )
                if "uid" in f:
                    uid_to_feat[f["uid"]] = feature_obj
        
        # Now add links
        count_links = 0
        for q, s, ident in links:
            if q in uid_to_feat and s in uid_to_feat:
                # color by identity?
                # simple grey for now
                if q != s: # self links already filtered but check again
                    gv.add_link(uid_to_feat[q], uid_to_feat[s], color="grey", alpha=0.3)
                    count_links += 1

        print(f"    Plotting {len(neighborhoods)} tracks with {count_links} links...")
        try:
            gv.savefig(pdf_out)
        except Exception as e:
            print(f"    Plotting failed: {e}")

    return

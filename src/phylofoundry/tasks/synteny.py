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
    Returns: list of dicts (feature info), contig_seq (Bio.Seq object)
    """
    neighborhood = []
    focal_feature = None
    contig_seq = None
    
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
                
                return neighborhood, record.seq
                
    except Exception as e:
        import sys
        print(f"[synteny] Warning: Failed to parse GenBank {gbk_path}: {e}", file=sys.stderr)
    return [], None

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
        return [], None

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
        
    return neighborhood, None

# Removed custom diamond functions - using pygenomeviz builtin


def run_synteny(cfg, synteny_dir, tree_dir, scan_df, search_df, hmm_keep, force=False):
    print("\n[synteny] Extracting neighborhoods and plotting synteny...")
    
    syn_cfg = cfg.get("synteny", {})
    if not syn_cfg.get("enabled", False):
        return

    # Check for clinker availability
    if not shutil.which("clinker"):
        import sys
        print("[synteny] FAILED: clinker command line tool is not installed or not in PATH.", file=sys.stderr)
        print("[synteny] Install with: pip install clinker-py", file=sys.stderr)
        return

    # Check for pgv-mmseqs availability
    if not shutil.which("pgv-mmseqs"):
        import sys
        print("[synteny] WARNING: pgv-mmseqs command line tool is not installed or not in PATH. pyGenomeViz plot will be skipped.", file=sys.stderr)
        print("[synteny] Install with: pip install pygenomeviz", file=sys.stderr)

    gbk_dir = syn_cfg.get("gbk_dir")
    gff_dir = syn_cfg.get("gff_dir")
    pfam_dir = cfg["inputs"].get("pfam_dir")
    
    if not gbk_dir and not gff_dir:
        print("WARNING: No gbk_dir or gff_dir provided for synteny. Skipping.")
        return

    # Aggregate hits
    summary_dir = os.path.join(cfg["output"]["outdir"], "summary")
    best_hits_tsv = os.path.join(summary_dir, "best_hits.competitive.tsv")
    
    hits = pd.DataFrame()
    if syn_cfg.get("prefer_best_hit", True) and os.path.exists(best_hits_tsv):
        hits = pd.read_csv(best_hits_tsv, sep="\t")
    else:
        dfs = []
        if not scan_df.empty: dfs.append(scan_df)
        if not search_df.empty: dfs.append(search_df)
        if dfs:
            hits = pd.concat(dfs, ignore_index=True)
            hits = hits.drop_duplicates(subset=["genome", "protein", "hmm"])

    if hits.empty:
        print("No hits found for synteny analysis.")
        return

    if hmm_keep:
        hits = hits[hits["hmm"].isin(hmm_keep)]

    for hmm, group in hits.groupby("hmm"):
        hmm_out_dir = os.path.join(synteny_dir, hmm)
        anno_dir = os.path.join(hmm_out_dir, "annotations")
        seqs_dir = os.path.join(hmm_out_dir, "seqs")
        
        safe_mkdir(hmm_out_dir)
        safe_mkdir(anno_dir)
        safe_mkdir(seqs_dir)
        
        html_out = os.path.join(hmm_out_dir, f"synteny.{hmm}.html")
        if os.path.exists(html_out) and not force:
            continue

        print(f"  Processing HMM: {hmm} ({len(group)} hits)...")
        
        if syn_cfg.get("dedup_by_genome", True):
            group = group.sort_values("bitscore", ascending=False).groupby("genome").head(1)
            max_hits = int(syn_cfg.get("max_hits_per_hmm", 50))
            if len(group) > max_hits:
                group = group.head(max_hits)

        neighborhoods = [] # list of (genome, track_feats)
        all_prots = {}     # uid -> sequence

        window = int(syn_cfg.get("window_genes", 10))
        fields_id = syn_cfg.get("protein_id_field", ["ID", "protein_id", "locus_tag"])
        fields_label = syn_cfg.get("gene_label_field", ["gene", "product", "Name", "locus_tag"])

        for _, r in group.iterrows():
            genome = r["genome"]
            protein = r["protein"]
            
            found_feats = []
            contig_seq = None
            if gbk_dir:
                base = normalize_genome_id(genome)
                cand = glob.glob(os.path.join(gbk_dir, f"{base}*.gb*"))
                if cand:
                    found_feats, contig_seq = load_genbank_neighborhood(cand[0], protein, window, fields_id, fields_label)
            
            if not found_feats and gff_dir:
                base = normalize_genome_id(genome)
                cand = glob.glob(os.path.join(gff_dir, f"{base}*.gff*"))
                if cand:
                    found_feats, contig_seq = load_gff_neighborhood(cand[0], protein, window, fields_id, fields_label)
            
            if not found_feats:
                continue

            track_feats = []
            for i, f in enumerate(found_feats):
                uid = f"{genome}|{f['contig']}|{f['protein_id']}|{i}"
                f["uid"] = uid
                if f.get("translation"):
                    all_prots[uid] = f["translation"]
                track_feats.append(f)
            
            neighborhoods.append((genome, track_feats, contig_seq))

        if not neighborhoods:
            print(f"    No neighborhoods extracted for {hmm}.")
            continue

        # Optional Pfam annotation
        pfam_mapping = {}
        if pfam_dir and all_prots:
            print("    Running Pfam annotation on neighborhood proteins...")
            import tempfile, json
            with tempfile.TemporaryDirectory() as td:
                faa_path = os.path.join(td, "proteins.faa")
                write_fasta(faa_path, all_prots)
                
                pfam_out = os.path.join(td, "pfam.json")
                # pfam_scan.py was copied into src/phylofoundry/utils/pfam_scan.py
                script_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "utils", "pfam_scan.py")
                cmd = ["python", script_path, faa_path, pfam_dir, "-outfmt", "json", "-out", pfam_out]
                try:
                    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
                    with open(pfam_out) as f:
                        pfam_res = json.load(f)
                    for uid, domains in pfam_res.items():
                        dom_names = [d["hmm_name"] for d in domains]
                        if dom_names:
                            pfam_mapping[uid] = "Pfam: " + ", ".join(dom_names)
                except Exception as e:
                    import sys
                    print(f"    Warning: Pfam scan failed: {e}", file=sys.stderr)

        # Save local structural proteins
        if all_prots:
            seq_out = os.path.join(seqs_dir, "neighborhood_proteins.faa")
            write_fasta(seq_out, all_prots)

        # Build GenBank files for clinker
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        
        gbk_files = []
        for genome, feats, contig_seq in neighborhoods:
            if not feats: continue
            min_start = min(f["start"] for f in feats)
            max_end = max(f["end"] for f in feats)
            size = max_end - min_start + 1
            
            if contig_seq is not None:
                nucl_seq = contig_seq[min_start:max_end]
                if isinstance(nucl_seq, str):
                    nucl_seq = Seq(nucl_seq)
            else:
                nucl_seq = Seq("N" * size)
            
            # Using molecule_type="DNA" allows writing out to genbank
            rec = SeqRecord(nucl_seq, id=genome, name=genome[:16], description=f"Neighborhood for {hmm}", annotations={"molecule_type": "DNA"})
            
            for f in feats:
                strand = 1 if str(f["strand"]) in ["1", "+"] else -1
                loc = FeatureLocation(f["start"] - min_start, f["end"] - min_start, strand=strand)
                
                qualifiers = {
                    "locus_tag": [f.get("label", "")],
                }
                if f.get("translation"):
                    qualifiers["translation"] = [f["translation"]]

                # Add Pfam annotation if available
                uid = f.get("uid", "")
                product_str = f.get("label", "")
                if uid in pfam_mapping:
                    if product_str and product_str != "NA":
                        product_str += f" ({pfam_mapping[uid]})"
                    else:
                        product_str = pfam_mapping[uid]
                
                if product_str and product_str != "NA":
                    qualifiers["product"] = [product_str]
                    qualifiers["gene"] = [product_str]  # Add as gene so clinker puts it on the figure
                    qualifiers["name"] = [product_str]

                sf = SeqFeature(loc, type="CDS", id=uid, qualifiers=qualifiers)
                rec.features.append(sf)
            
            gbk_path = os.path.join(anno_dir, f"{genome}.gbk")
            with open(gbk_path, "w") as outf:
                SeqIO.write(rec, outf, "genbank")
            gbk_files.append(gbk_path)
            
        # Run clinker
        print("    Generating Clinker synteny plot...")
        try:
            cmd = ["clinker", *gbk_files, "-p", html_out]
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
            print(f"    Saved Clinker plot: {html_out}")
            
        except Exception as e:
            import sys
            import traceback
            print(f"    Clinker plotting failed for {hmm}: {e}", file=sys.stderr)
            traceback.print_exc()

        # Run pgv-mmseqs
        if shutil.which("pgv-mmseqs"):
            print("    Generating pyGenomeViz synteny plot (pgv-mmseqs)...")
            try:
                pgv_out_dir = os.path.join(hmm_out_dir, "pgv_out")
                pgv_cmd = ["pgv-mmseqs", *gbk_files, "-o", pgv_out_dir]
                # Default behavior is to launch the image to outdir/result.png (or pdf)
                subprocess.run(pgv_cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
                print(f"    Saved pyGenomeViz plot: {pgv_out_dir}")
            except Exception as e:
                import sys
                import traceback
                print(f"    pyGenomeViz plotting failed for {hmm}: {e}", file=sys.stderr)
                traceback.print_exc()
            
            # Optionally clean up intermediate genbanks here, uncomment if preferred
            # for f in gbk_files:
            #     os.remove(f)

    return

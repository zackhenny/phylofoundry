import os
import sys
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


def _annotate_neighborhood_proteins(proteins_faa, combined_hmm, annotation_evalue, cpu=1):
    """Run hmmsearch on neighborhood proteins against the pipeline HMM set.

    Returns dict {protein_uid: {hmm_hit, evalue, bitscore}} for best hits.
    """
    import tempfile
    annotations = {}
    if not os.path.exists(combined_hmm) or not os.path.exists(proteins_faa):
        return annotations

    with tempfile.NamedTemporaryFile(suffix=".tbl", delete=False) as tmp:
        tbl_path = tmp.name

    try:
        cmd = [
            "hmmsearch", "--cpu", str(cpu),
            "--domtblout", tbl_path,
            "-E", str(annotation_evalue),
            "--domE", str(annotation_evalue),
            combined_hmm, proteins_faa
        ]
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        # Parse domtblout
        if os.path.exists(tbl_path) and os.path.getsize(tbl_path) > 0:
            with open(tbl_path) as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    parts = line.split()
                    if len(parts) < 13:
                        continue
                    protein_id = parts[0]
                    hmm_name = parts[3]
                    evalue = float(parts[6])
                    bitscore = float(parts[7])
                    # Keep best hit per protein
                    if protein_id not in annotations or bitscore > annotations[protein_id]["bitscore"]:
                        annotations[protein_id] = {
                            "hmm_annotation": hmm_name,
                            "annotation_evalue": evalue,
                            "annotation_bitscore": bitscore,
                        }
    except Exception as e:
        print(f"[synteny] Warning: HMM annotation of neighborhood proteins failed: {e}", file=sys.stderr)
    finally:
        try:
            os.remove(tbl_path)
        except OSError:
            pass

    return annotations

def _get_tree_tip_order(tree_fp):
    """Return an ordered list of tip labels as rendered by a rectangular tree.

    Uses scikit-bio's TreeNode to traverse the tree and extract tips in
    the order they would appear in a left-to-right rectangular layout
    (post-order traversal of the Newick tree as written by IQ-TREE).

    Returns an empty list if the treefile cannot be read or scikit-bio
    is unavailable.
    """
    try:
        from skbio import TreeNode
    except ImportError:
        return []

    try:
        tree = TreeNode.read(tree_fp, format="newick")
        return [n.name for n in tree.tips() if n.name]
    except Exception:
        return []


def _call_synteny_tree_r(
    rscript_bin,
    r_script,
    hmm_name,
    tree_fp,
    synteny_tsv,
    tree_panel_dir,
    clades_tsv=None,
    tax_tsv=None,
    formats="png",
    width=20,
    height=0,
    bootstrap_min=80,
    show_tip_labels="auto",
    color_palette="Set3",
    tax_level="genus",
):
    """Invoke ``scripts/plot_synteny_tree.R`` for one HMM."""
    cmd = [
        rscript_bin, "--vanilla", r_script,
        "--treefile", tree_fp,
        "--synteny_tsv", synteny_tsv,
        "--outdir", tree_panel_dir,
        "--hmm_name", hmm_name,
        "--format", formats,
        "--width", str(width),
        "--height", str(height),
        "--bootstrap_min", str(bootstrap_min),
        "--show_tip_labels", show_tip_labels,
        "--color_palette", color_palette,
        "--tax_level", tax_level,
    ]
    if clades_tsv:
        cmd += ["--clades_tsv", clades_tsv]
    if tax_tsv:
        cmd += ["--taxonomy_tsv", tax_tsv]

    try:
        result = subprocess.run(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        if result.returncode != 0:
            print(
                f"    [synteny] WARNING: plot_synteny_tree.R failed for {hmm_name}:\n"
                f"{result.stderr[-2000:]}",
                file=sys.stderr,
            )
        else:
            print(f"    [synteny] Saved synteny+tree panel: {tree_panel_dir}")
    except Exception as exc:
        print(
            f"    [synteny] WARNING: R subprocess error for {hmm_name}: {exc}",
            file=sys.stderr,
        )


def run_synteny(cfg, synteny_dir, tree_dir, scan_df, search_df, hmm_keep, force=False, clade_assign_dir=None, summary_dir=None):
    print("\n[synteny] Extracting neighborhoods and plotting synteny...")
    
    syn_cfg = cfg.get("synteny", {})
    if not syn_cfg.get("enabled", False):
        return

    # Check for clinker availability
    if not shutil.which("clinker"):
        print("[synteny] FAILED: clinker command line tool is not installed or not in PATH.", file=sys.stderr)
        print("[synteny] Install with: pip install clinker-py", file=sys.stderr)
        return

    # Check for pgv-mmseqs availability
    if not shutil.which("pgv-mmseqs"):
        print("[synteny] WARNING: pgv-mmseqs command line tool is not installed or not in PATH. pyGenomeViz plot will be skipped.", file=sys.stderr)
        print("[synteny] Install with: pip install pygenomeviz", file=sys.stderr)

    gbk_dir = syn_cfg.get("gbk_dir")
    gff_dir = syn_cfg.get("gff_dir")
    pfam_dir = cfg["inputs"].get("pfam_dir")
    annotation_evalue = float(syn_cfg.get("annotation_evalue", 1e-5))
    cpu = int(cfg.get("resources", {}).get("cpu", 1))
    combined_hmm = os.path.join(cfg["output"]["outdir"], "combined.hmm")

    if not gbk_dir and not gff_dir:
        print("WARNING: No gbk_dir or gff_dir provided for synteny. Skipping.")
        return

    # Resolve summary_dir: prefer passed-in argument, fall back to cfg-derived path
    if summary_dir is None:
        summary_dir = os.path.join(cfg["output"]["outdir"], "summary")
    best_hits_tsv = os.path.join(summary_dir, "best_hits.competitive.tsv")

    # Locate taxonomy annotation file for the R tree-panel script
    tax_tsv_for_r = os.path.join(summary_dir, "genome_taxonomy.tsv")
    if not os.path.exists(tax_tsv_for_r):
        tax_tsv_for_r = os.path.join(summary_dir, "best_hits.with_taxonomy.tsv")
    if not os.path.exists(tax_tsv_for_r):
        tax_tsv_for_r = None

    # Resolve R script paths for the tree+synteny panel
    rscript_bin = shutil.which("Rscript")
    _tasks_dir = os.path.dirname(__file__)
    _repo_root = os.path.abspath(os.path.join(_tasks_dir, "..", "..", "..", ".."))
    synteny_r_script = os.path.join(_repo_root, "scripts", "plot_synteny_tree.R")
    if not os.path.exists(synteny_r_script):
        synteny_r_script = None

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

    all_cluster_rows = []  # aggregate across all HMMs

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

        print(f"  [synteny] Calculating synteny for {hmm} ({len(group)} hits)...")
        
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
            
            neighborhoods.append((genome, protein, track_feats, contig_seq))

        if not neighborhoods:
            print(f"    No neighborhoods extracted for {hmm}.")
            continue

        # ── Tree-ordered neighborhood sorting ─────────────────────────────
        # When tree_ordered=true (and include_tree=true), sort neighborhoods
        # by the IQ-TREE leaf order so that closely related organisms are
        # adjacent in the synteny plot.
        if syn_cfg.get("tree_ordered", True) and syn_cfg.get("include_tree", True):
            tree_fp_for_order = os.path.join(tree_dir, f"{hmm}.treefile")
            tip_order = _get_tree_tip_order(tree_fp_for_order)
            if tip_order:
                # Build a rank dict: tip_label → index in tree order
                tip_rank = {tip: i for i, tip in enumerate(tip_order)}

                def _neighborhood_rank(nbr):
                    genome_n, protein_n, _, _ = nbr
                    # Tips in tree are formatted as "genome|protein"
                    full_key = f"{genome_n}|{protein_n}"
                    if full_key in tip_rank:
                        return tip_rank[full_key]
                    # Fallback: match by genome prefix
                    for tip, rank in tip_rank.items():
                        if tip.startswith(f"{genome_n}|"):
                            return rank
                    return len(tip_rank)  # unmatched → end

                neighborhoods = sorted(neighborhoods, key=_neighborhood_rank)


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
                    print(f"    Warning: Pfam scan failed: {e}", file=sys.stderr)

        # Save local structural proteins
        if all_prots:
            seq_out = os.path.join(seqs_dir, "neighborhood_proteins.faa")
            write_fasta(seq_out, all_prots)

        # ── HMM Annotation of neighborhood proteins ──────────────────────
        hmm_annotations = {}
        if all_prots and os.path.exists(combined_hmm):
            seq_out = os.path.join(seqs_dir, "neighborhood_proteins.faa")
            print(f"    Annotating neighborhood proteins (e-value={annotation_evalue})...")
            hmm_annotations = _annotate_neighborhood_proteins(
                seq_out, combined_hmm, annotation_evalue, cpu=cpu
            )
            if hmm_annotations:
                print(f"    Annotated {len(hmm_annotations)} proteins with HMM hits.")

        # ── Gene Cluster Annotations Table ────────────────────────────────
        cluster_rows = []
        for genome, protein, feats, contig_seq in neighborhoods:
            for f in feats:
                uid = f.get("uid", "")
                anno = hmm_annotations.get(uid, {})
                # Compute relative position to focal gene
                focal_idx = next((i for i, feat in enumerate(feats) if feat.get("is_focal")), 0)
                current_idx = feats.index(f)
                position_relative = current_idx - focal_idx
                cluster_rows.append({
                    "hmm": hmm,
                    "genome": genome,
                    "focal_protein": protein,
                    "cluster_size": len(feats),
                    "protein_id": f.get("protein_id", "NA"),
                    "position_relative": position_relative,
                    "strand": f.get("strand", "NA"),
                    "gene_label": f.get("label", "NA"),
                    "hmm_annotation": anno.get("hmm_annotation", ""),
                    "annotation_evalue": anno.get("annotation_evalue", ""),
                    "annotation_bitscore": anno.get("annotation_bitscore", ""),
                    "is_focal": f.get("is_focal", False),
                })

        if cluster_rows:
            cluster_df = pd.DataFrame(cluster_rows)
            cluster_fp = os.path.join(anno_dir, f"{hmm}_gene_clusters.tsv")
            cluster_df.to_csv(cluster_fp, sep="\t", index=False)
            print(f"    Wrote gene cluster table: {cluster_fp} ({len(cluster_rows)} entries)")
            all_cluster_rows.extend(cluster_rows)

        # Build GenBank files for clinker
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        
        gbk_files = []
        for genome, protein, feats, contig_seq in neighborhoods:
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

            # Detect focal gene strand and orient everything consistently
            focal_feat = next((f for f in feats if f.get("is_focal")), None)
            focal_strand = 1
            if focal_feat is not None:
                focal_strand = 1 if str(focal_feat["strand"]) in ["1", "+"] else -1

            if focal_strand == -1:
                # Reverse-complement the sequence so the focal gene points forward
                nucl_seq = nucl_seq.reverse_complement()

            # Make track name unique to support multiple contigs/proteins per genome
            safe_prot = protein.replace("|", "_").split("~")[-1]
            track_name = f"{genome}_{safe_prot}"
            if len(track_name) > 16:
                track_name = track_name[:16]

            # Using molecule_type="DNA" allows writing out to genbank
            rec = SeqRecord(nucl_seq, id=f"{genome}_{safe_prot}", name=track_name, description=f"Neighborhood for {hmm}", annotations={"molecule_type": "DNA"})
            
            for f in feats:
                strand = 1 if str(f["strand"]) in ["1", "+"] else -1
                feat_start = f["start"] - min_start
                feat_end = f["end"] - min_start

                if focal_strand == -1:
                    # Flip coordinates relative to the reversed sequence
                    new_start = size - feat_end
                    new_end = size - feat_start
                    strand = -strand
                    feat_start, feat_end = new_start, new_end

                loc = FeatureLocation(feat_start, feat_end, strand=strand)
                
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
            
            gbk_path = os.path.join(anno_dir, f"{genome}_{safe_prot}.gbk")
            with open(gbk_path, "w") as outf:
                SeqIO.write(rec, outf, "genbank")
            gbk_files.append(gbk_path)
            
        # Run clinker
        print("    Generating Clinker synteny plot...")
        try:
            cmd = ["clinker", *gbk_files, "-p", html_out, "--no_plot"]
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
            print(f"    Saved Clinker plot: {html_out}")

        except Exception as e:
            import traceback
            print(f"    Clinker plotting failed for {hmm}: {e}", file=sys.stderr)
            traceback.print_exc()

        # Run pgv-mmseqs
        if shutil.which("pgv-mmseqs"):
            print("    Generating pyGenomeViz synteny plot (pgv-mmseqs)...")
            try:
                pgv_out_dir = os.path.join(hmm_out_dir, "pgv_out")
                pgv_cmd = ["pgv-mmseqs", *gbk_files, "-o", pgv_out_dir, "--no-show"]
                # Default behavior is to launch the image to outdir/result.png (or pdf)
                subprocess.run(pgv_cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
                print(f"    Saved pyGenomeViz plot: {pgv_out_dir}")
            except Exception as e:
                import traceback
                print(f"    pyGenomeViz plotting failed for {hmm}: {e}", file=sys.stderr)
                traceback.print_exc()

            # Optionally clean up intermediate genbanks here, uncomment if preferred
            # for f in gbk_files:
            #     os.remove(f)

        # ── ggtree-based combined tree + synteny panel (R) ────────────────
        # Requires: a .treefile for the HMM, the gene-cluster TSV, and Rscript.
        if syn_cfg.get("generate_tree_panel", True) and rscript_bin and synteny_r_script:
            tree_fp_for_panel = os.path.join(tree_dir, f"{hmm}.treefile")
            cluster_tsv_for_panel = os.path.join(anno_dir, f"{hmm}_gene_clusters.tsv")
            if os.path.exists(tree_fp_for_panel) and os.path.exists(cluster_tsv_for_panel):
                print("    Generating ggtree synteny+tree panel...")
                tree_panel_dir = os.path.join(hmm_out_dir, "tree_panel")
                safe_mkdir(tree_panel_dir)

                panel_formats = syn_cfg.get("tree_panel_format", ["png"])
                panel_formats_str = ",".join(panel_formats) if isinstance(panel_formats, list) else str(panel_formats)

                clades_tsv_for_panel = None
                if clade_assign_dir:
                    _ctsv = os.path.join(clade_assign_dir, f"{hmm}.clades.tsv")
                    if os.path.exists(_ctsv):
                        clades_tsv_for_panel = _ctsv

                viz_cfg = cfg.get("phylo", {}).get("tree_viz", {})
                _call_synteny_tree_r(
                    rscript_bin=rscript_bin,
                    r_script=synteny_r_script,
                    hmm_name=hmm,
                    tree_fp=tree_fp_for_panel,
                    synteny_tsv=cluster_tsv_for_panel,
                    tree_panel_dir=tree_panel_dir,
                    clades_tsv=clades_tsv_for_panel,
                    tax_tsv=tax_tsv_for_r,
                    formats=panel_formats_str,
                    width=float(syn_cfg.get("tree_panel_width", 20)),
                    height=0,
                    bootstrap_min=float(viz_cfg.get("bootstrap_min", 80)),
                    show_tip_labels=str(viz_cfg.get("show_tip_labels", "auto")),
                    color_palette=str(viz_cfg.get("color_palette", "Set3")),
                    tax_level=str(viz_cfg.get("tax_level", "genus")),
                )

    # Write aggregate gene cluster annotations table
    if all_cluster_rows:
        annotations_dir = os.path.join(synteny_dir, "annotations")
        safe_mkdir(annotations_dir)
        agg_fp = os.path.join(annotations_dir, "all_gene_clusters.tsv")
        pd.DataFrame(all_cluster_rows).to_csv(agg_fp, sep="\t", index=False)
        print(f"[synteny] Wrote aggregate gene cluster table: {agg_fp} ({len(all_cluster_rows)} entries)")

    return

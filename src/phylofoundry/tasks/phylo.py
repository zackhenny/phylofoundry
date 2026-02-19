import os
import subprocess
import glob
from concurrent.futures import ProcessPoolExecutor
from ..utils.bio import write_fasta, read_fasta
from ..utils.helpers import run_cmd

def worker_phylo(args_pack):
    (hmm_name, seqs, fasta_dir, aln_dir, clipkit_dir, tree_dir,
     hmm_path, threads, phy_cfg, force) = args_pack

    fa_out = os.path.join(fasta_dir, f"{hmm_name}.faa")
    if (not os.path.exists(fa_out)) or force:
        write_fasta(fa_out, seqs)

    mafft_aln = os.path.join(aln_dir, f"{hmm_name}.mafft.fasta")
    hmmalign_afa = os.path.join(aln_dir, f"{hmm_name}.afa")
    clip_out = os.path.join(clipkit_dir, f"{hmm_name}.clipkit.faa")
    tree_prefix = os.path.join(tree_dir, hmm_name)
    treefile = tree_prefix + ".treefile"

    if os.path.exists(treefile) and not force:
        return hmm_name

    # Align
    try:
        if phy_cfg["mafft"]:
            if (not os.path.exists(mafft_aln)) or force:
                mode = phy_cfg.get("mafft_mode", "auto")
                # Map common names to flags if necessary, or just pass through if user knows flags
                # Simple mapping for convenience:
                mode_flags = {
                    "auto": "--auto",
                    "ginsi": "--globalpair --maxiterate 1000",
                    "linsi": "--localpair --maxiterate 1000",
                    "einsi": "--genafpair --maxiterate 1000",
                    "fftns": "--fftns",
                    "fftnsi": "--fftnsi",
                }
                flags = mode_flags.get(mode, mode) # Default to passing raw string if not in map (e.g. "--auto")
                
                cmd = f"mafft {flags} --thread {threads} --quiet {fa_out} > {mafft_aln}"
                run_cmd(cmd, quiet=True, shell=True)
            alignment_to_use = mafft_aln
        else:
            sto_out = os.path.join(aln_dir, f"{hmm_name}.sto")
            if (not os.path.exists(hmmalign_afa)) or force:
                cmd = ["hmmalign", "--outformat", "stockholm", "-o", sto_out, hmm_path, fa_out]
                if not phy_cfg["no_trim_hmmalign"]:
                    cmd.insert(1, "--trim")
                subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                run_cmd(f"esl-reformat afa {sto_out} > {hmmalign_afa}", quiet=True, shell=True)

            if phy_cfg["also_mafft"]:
                if (not os.path.exists(mafft_aln)) or force:
                    cmd = f"mafft --thread {threads} --quiet {fa_out} > {mafft_aln}"
                    run_cmd(cmd, quiet=True, shell=True)

            if phy_cfg["mafft_for_tree"] and os.path.exists(mafft_aln):
                alignment_to_use = mafft_aln
            else:
                alignment_to_use = hmmalign_afa

    except Exception:
        return None

    # ClipKit
    if phy_cfg["skip_clipkit"]:
        tree_input = alignment_to_use
    else:
        if (not os.path.exists(clip_out)) or force:
            try:
                subprocess.run(["clipkit", alignment_to_use, "-o", clip_out],
                               check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            except Exception:
                return None
        tree_input = clip_out

    if len(seqs) < 3:
        return hmm_name

    # IQ-TREE
    iqtree_bin = phy_cfg.get("iqtree_bin", "iqtree")
    iq_cmd = [
        iqtree_bin, "-s", tree_input, "-m", "MFP",
        "-B", str(phy_cfg["iq_boot"]),
        "-T", str(threads),
        "-pre", tree_prefix, "-quiet"
    ]
    if not phy_cfg["no_asr"]:
        iq_cmd.insert(len(iq_cmd) - 1, "-asr")

    try:
        subprocess.run(iq_cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError as e:
        print(f"[phylo] warning: IQ-TREE failed for {hmm_name}: {e}", file=os.sys.stderr)

    return hmm_name

def run_phylo(cfg, hmm_to_seqs, fasta_dir, aln_dir, clipkit_dir, tree_dir, name_to_hmm_path, hmm_keep, force=False):
    print("\n[phylo] Per-HMM alignment + trimming + IQ-TREE...")

    phy_cfg = cfg.get("phylo", {})
    cpu = cfg["resources"]["cpu"]

    if not hmm_to_seqs:
        for fp in glob.glob(os.path.join(fasta_dir, "*.faa")):
            hmm = os.path.basename(fp).rsplit(".", 1)[0]
            if hmm_keep is not None and hmm not in hmm_keep:
                continue
            hmm_to_seqs[hmm] = read_fasta(fp)

    threads_per = max(1, min(4, cpu // 2))
    workers = max(1, cpu // threads_per)

    tasks = []
    for hmm, seqs in hmm_to_seqs.items():
        if hmm_keep is not None and hmm not in hmm_keep:
            continue
        
        hmm_path = None
        if not phy_cfg.get("mafft", False):
            if hmm not in name_to_hmm_path:
                print(f"[phylo] warning: HMM {hmm} not found in names map, skipping hmmalign.", file=os.sys.stderr)
                continue
            hmm_path = name_to_hmm_path[hmm]

        if phy_cfg.get("mafft", False):
            hmm_path = hmm_path or ""  # not used in MAFFT mode

        tasks.append((hmm, seqs, fasta_dir, aln_dir, clipkit_dir, tree_dir, hmm_path, threads_per, phy_cfg, force))

    with ProcessPoolExecutor(max_workers=workers) as exe:
        list(exe.map(worker_phylo, tasks))

    # Combined tree mapping
    if phy_cfg.get("combined_tree", False):
        print("[phylo] Building combined tree for all hits...")
        import sys
        
        combined_seqs = {}
        for hmm, seqs in hmm_to_seqs.items():
            if hmm_keep is not None and hmm not in hmm_keep:
                continue
            for tip, s in seqs.items():
                if tip not in combined_seqs:
                    combined_seqs[tip] = s
        
        if len(combined_seqs) >= 3:
            cmb_fasta = os.path.join(fasta_dir, "combined_all_hits.faa")
            cmb_aln = os.path.join(aln_dir, "combined_all_hits.mafft.fasta")
            cmb_clip = os.path.join(clipkit_dir, "combined_all_hits.clipkit.faa")
            tree_prefix = os.path.join(tree_dir, "combined_all_hits")
            treefile = tree_prefix + ".treefile"
            
            if not os.path.exists(treefile) or force:
                write_fasta(cmb_fasta, combined_seqs)
                
                # Mafft alignment
                if not os.path.exists(cmb_aln) or force:
                    print("  Running MAFFT on combined hits...")
                    flags = {"auto": "--auto", "ginsi": "--globalpair --maxiterate 1000"}.get(phy_cfg.get("mafft_mode", "auto"), "--auto")
                    cmd = f"mafft {flags} --thread {cpu} --quiet {cmb_fasta} > {cmb_aln}"
                    try:
                        run_cmd(cmd, quiet=True, shell=True)
                    except Exception as e:
                        print(f"[phylo] combined mafft failed: {e}", file=sys.stderr)
                
                # ClipKit
                if not os.path.exists(cmb_clip) or force:
                    print("  Running ClipKit on combined hits...")
                    try:
                        subprocess.run(["clipkit", cmb_aln, "-o", cmb_clip], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                    except Exception as e:
                        print(f"[phylo] combined clipkit failed: {e}", file=sys.stderr)
                
                # IQ-TREE
                if not os.path.exists(treefile) or force:
                    print("  Running IQ-TREE on combined hits...")
                    iqtree_bin = phy_cfg.get("iqtree_bin", "iqtree")
                    iq_cmd = [
                        iqtree_bin, "-s", cmb_clip, "-m", "MFP",
                        "-B", str(phy_cfg.get("iq_boot", 1000)),
                        "-T", str(cpu),
                        "-pre", tree_prefix, "-quiet"
                    ]
                    if not phy_cfg.get("no_asr", False):
                        iq_cmd.insert(len(iq_cmd) - 1, "-asr")
                    
                    try:
                        subprocess.run(iq_cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                    except subprocess.CalledProcessError as e:
                        print(f"[phylo] combined IQ-TREE failed: {e}", file=sys.stderr)
        else:
            print("[phylo] Less than 3 combined sequences, skipping combined tree.", file=os.sys.stderr)

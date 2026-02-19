#!/usr/bin/env python3
"""
Competitive HMM pipeline (hmmscan + hmmsearch + per-HMM phylogenies) with:
- JSON config support (defaults -> config.json -> CLI overrides) + resolved_config.json provenance
- HMM input can be a directory OR a single .hmm file
- Minimal "workflow construction" upgrades (shortest path):
  - caching by file existence (skip if outputs exist unless --force)
  - --start_at / --stop_after step slicing
  - --hmm_manifest (subset processing)
  - per-HMM isolation (one failure doesn't kill whole run)

Core steps:
  prep   : build combined_proteomes.faa + combined.hmm + hmmpress
  hmmer  : hmmscan per genome + hmmsearch per HMM + filter -> best competitive hits
  extract: per-HMM FASTA bins for phylogeny/embeddings
  embed  : OPTIONAL protein language model embeddings per-HMM (functional landscape proxy)
  phylo  : per-HMM align/trim/tree (IQ-TREE3) (+ optional ASR)
  post   : OPTIONAL scikit-bio conservation + KL divergence (clade comparisons)
  codon  : OPTIONAL PAL2NAL codon alignments (requires cds_dir)
  hyphy  : OPTIONAL HyPhy wrappers (generic; branch labeling not handled here)

Embedding notes (what you get here):
- Per-HMM embeddings computed on unaligned AA sequences
- Mean-pooled representations (excluding special/pad tokens)
- PCA coordinates (PC1..PC3) for quick plotting
- Optional clade dispersion metrics if clades_tsv provided (mean Euclidean distance to clade centroid)

Backends supported:
- "esm" (Meta FAIR ESM) via `esm` + `torch`
- "transformers" via `transformers` + `torch` (e.g., ProtT5/ESM2 HF checkpoints)
This script does NOT download models; you must have them available in your environment / cache.

Outputs:
- summary/resolved_config.json
- summary/hmmscan_hits.filtered.tsv
- summary/hmmsearch_hits.filtered.tsv
- summary/best_hits.competitive.tsv
- fasta_per_hmm/*.faa
- embeddings/*.pca.tsv + embeddings/*.npy (+ optional dispersion.tsv)
- trees_iqtree/* (per-HMM)
- summary/post_scikitbio/* (optional)
- codon_alignments/* (optional)
- summary/hyphy/* (optional)
"""

import os
import re
import glob
import math
import json
import argparse
import subprocess
from copy import deepcopy
from collections import defaultdict, Counter
from concurrent.futures import ProcessPoolExecutor

import pandas as pd

# scikit-bio (post-processing glue)
from skbio import TreeNode
from skbio.alignment import TabularMSA
from skbio.sequence import Protein
import warnings


###############################################################################
# Constants
###############################################################################

AA_ALPHABET = set(list("ACDEFGHIKLMNPQRSTVWY"))
GAP_CHARS = set(["-", ".", "X", "x", "?", "*"])

STEPS = ["prep", "hmmer", "extract", "embed", "phylo", "post", "codon", "hyphy"]


DEFAULT_CONFIG = {
    "inputs": {
        "faa_dir": None,     # directory of *.faa OR a single .faa
        "hmm_input": None,   # directory of *.hmm OR a single .hmm
        "cds_dir": None      # optional directory of CDS nucleotide FASTAs (per genome)
    },
    "output": {
        "outdir": None
    },
    "resources": {
        "cpu": 8
    },
    "workflow": {
        "start_at": None,
        "stop_after": None,
        "force": False,
        "hmm_manifest": None
    },
    "filtering": {
        "scores_tsv": None,
        "global_min_score": 25.0,
        "min_coverage": 0.5,
        "keep_tbl": False
    },
    "phylo": {
        "mafft": False,
        "also_mafft": False,
        "mafft_for_tree": False,
        "no_trim_hmmalign": False,
        "skip_clipkit": False,
        "no_asr": False,
        "iq_boot": 1000,
        "use_hmmsearch_alignment": False,
        "keep_all_hits": False
    },
    "embeddings": {
        "enabled": False,
        "backend": "esm",            # "esm" or "transformers"
        "model": "esm2_t33_650M_UR50D",  # for backend=esm; or HF model id/path for transformers
        "device": "cuda",            # "cuda" or "cpu"
        "batch_size": 8,
        "repr_layer": None,          # if None, choose last layer (esm) or last_hidden_state (transformers)
        "pooling": "mean",           # "mean" (implemented)
        "pca_components": 3,
        "write_full_vectors": False  # if True, write TSV with all dims (can be huge); always writes .npy
    },
    "post": {
        "enabled": False,
        "compute_conservation": False,
        "conservation_metric": "inverse_shannon_uncertainty",
        "compute_kl": False,
        "clades_tsv": None,  # TSV columns: clade_name, tip (tip label must match alignment tip labels)
        "kl_pairs": None     # "A:B,A:background"
    },
    "codon": {
        "enabled": False,
        "build_codon_alignments": False,
        "cds_id_mode": "after_last_pipe",  # "same"|"strip_pipe"|"after_last_pipe"
        "pal2nal_cmd": "pal2nal.pl"
    },
    "hyphy": {
        "enabled": False,
        "run_hyphy": False,
        "hyphy_bin": "hyphy",
        "hyphy_tests": "RELAX,aBSREL,MEME"
    }
}


###############################################################################
# Utility
###############################################################################

def safe_mkdir(p):
    os.makedirs(p, exist_ok=True)


def deep_update(base: dict, updates: dict) -> dict:
    for k, v in updates.items():
        if isinstance(v, dict) and isinstance(base.get(k), dict):
            deep_update(base[k], v)
        else:
            base[k] = v
    return base


def load_json_config(path: str) -> dict:
    with open(path) as f:
        return json.load(f)


def write_json(obj: dict, out_fp: str):
    safe_mkdir(os.path.dirname(out_fp))
    with open(out_fp, "w") as f:
        json.dump(obj, f, indent=2, sort_keys=True)


def run_cmd(cmd, quiet=False, shell=False):
    if not quiet:
        if shell and isinstance(cmd, str):
            print("Running:", cmd)
        else:
            print("Running:", " ".join(cmd))
    stdout_target = subprocess.DEVNULL if quiet else None
    stderr_target = subprocess.DEVNULL if quiet else None
    subprocess.run(cmd, check=True, stdout=stdout_target, stderr=stderr_target, shell=shell)


def normalize_genome_id(x: str) -> str:
    if x is None:
        return x
    x = str(x)
    x = re.sub(r"\.gz$", "", x, flags=re.IGNORECASE)
    while True:
        new = re.sub(r"\.(faa|fna|fa|fasta|ffn|cds)$", "", x, flags=re.IGNORECASE)
        if new == x:
            break
        x = new
    return x


def step_in_range(step, start_at, stop_after):
    i = STEPS.index(step)
    i0 = STEPS.index(start_at) if start_at else 0
    i1 = STEPS.index(stop_after) if stop_after else len(STEPS) - 1
    return i0 <= i <= i1


def read_fasta(fp):
    seqs = {}
    h = None
    parts = []
    with open(fp) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if h is not None:
                    seqs[h] = "".join(parts)
                h = line[1:].strip().split()[0]
                parts = []
            else:
                parts.append(line.strip())
        if h is not None:
            seqs[h] = "".join(parts)
    return seqs


def write_fasta(fp, records):
    with open(fp, "w") as out:
        for h, s in records.items():
            out.write(f">{h}\n{s}\n")


def load_manifest(hmm_manifest_fp):
    if not hmm_manifest_fp:
        return None
    keep = set()
    with open(hmm_manifest_fp) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            keep.add(s)
    return keep


def load_custom_thresholds(tsv_path):
    thresholds = {}
    if not tsv_path:
        return thresholds
    print(f"Loading custom thresholds from {tsv_path}...")
    with open(tsv_path) as f:
        for line in f:
            if line.strip().startswith("#") or not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) >= 2:
                thresholds[parts[0]] = float(parts[1])
    return thresholds


###############################################################################
# HMMER domtblout parsing
###############################################################################

def parse_tblout_fast(tbl_file):
    if not os.path.exists(tbl_file) or os.path.getsize(tbl_file) == 0:
        return pd.DataFrame()
    try:
        df = pd.read_csv(
            tbl_file,
            delim_whitespace=True,
            comment="#",
            header=None,
            usecols=[0, 2, 3, 5, 6, 7, 17, 18],
            names=["target", "tlen", "query", "qlen", "evalue", "bitscore", "ali_from", "ali_to"],
            engine="c",
        )
        return df
    except Exception:
        return pd.DataFrame()


def apply_filtering(df, thresholds_map, global_min_score, global_min_cov):
    if df.empty:
        return df
    if global_min_cov > 0:
        df = df[df["coverage"] >= global_min_cov].copy()
    df["min_score_required"] = global_min_score
    if thresholds_map:
        custom_scores = df["hmm"].map(thresholds_map)
        df["min_score_required"] = custom_scores.fillna(df["min_score_required"])
    df = df[df["bitscore"] >= df["min_score_required"]].copy()
    return df


def best_hits(df):
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
            "best_hmm": best.hmm,
            "best_bitscore": float(best.bitscore),
            "best_evalue": float(best.evalue),
            "best_coverage": float(best.coverage),
            "delta_bitscore": float(delta)
        })
    return pd.DataFrame(rows)


###############################################################################
# Parallel workers: hmmscan / hmmsearch
###############################################################################

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


###############################################################################
# Phylo worker
###############################################################################

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
                cmd = f"mafft --thread {threads} --quiet {fa_out} > {mafft_aln}"
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
    iq_cmd = [
        "iqtree3", "-s", tree_input, "-m", "MFP",
        "-B", str(phy_cfg["iq_boot"]),
        "-T", str(threads),
        "-pre", tree_prefix, "-quiet"
    ]
    if not phy_cfg["no_asr"]:
        iq_cmd.insert(len(iq_cmd) - 1, "-asr")

    try:
        subprocess.run(iq_cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError:
        pass

    return hmm_name


###############################################################################
# scikit-bio post-processing
###############################################################################

def msa_from_fasta_dict(aln_seqs: dict) -> TabularMSA:
    seqs = [Protein(seq, metadata={"id": tip}) for tip, seq in aln_seqs.items()]
    return TabularMSA(seqs)


def consurf_like_scores(aln_seqs: dict, metric="inverse_shannon_uncertainty"):
    msa = msa_from_fasta_dict(aln_seqs)
    scores = msa.conservation(metric=metric, gap_mode="exclude")
    return list(scores)


def site_counts_from_subset(aln_seqs: dict, subset_tips=None):
    if not aln_seqs:
        return []
    msa = msa_from_fasta_dict(aln_seqs)
    tip_set = set(subset_tips) if subset_tips is not None else None
    msa_ids = [seq.metadata["id"] for seq in msa]
    counts = []
    for pos in msa.iter_positions():
        c = Counter()
        s = str(pos)
        for idx, char in enumerate(s):
            tip = msa_ids[idx]
            if tip_set is not None and tip not in tip_set:
                continue
            aa = char
            if aa in GAP_CHARS:
                continue
            aa = aa.upper()
            if aa in AA_ALPHABET:
                c[aa] += 1
        counts.append(c)
    return counts


def kl_divergence(p_counts, q_counts, pseudocount=1e-6):
    alphabet = list(AA_ALPHABET)
    p_total = sum(p_counts.get(a, 0) for a in alphabet) + pseudocount * len(alphabet)
    q_total = sum(q_counts.get(a, 0) for a in alphabet) + pseudocount * len(alphabet)
    kl = 0.0
    for a in alphabet:
        p = (p_counts.get(a, 0) + pseudocount) / p_total
        q = (q_counts.get(a, 0) + pseudocount) / q_total
        kl += p * math.log(p / q, 2)
    return float(kl)


def load_clades_tsv(clades_tsv):
    df = pd.read_csv(clades_tsv, sep="\t", dtype=str)
    if "clade_name" not in df.columns or "tip" not in df.columns:
        raise SystemExit("clades_tsv must have columns: clade_name, tip")
    out = defaultdict(list)
    for _, r in df.iterrows():
        if pd.isna(r["clade_name"]) or pd.isna(r["tip"]):
            continue
        out[str(r["clade_name"])].append(str(r["tip"]))
    return dict(out)


def tree_load(tree_fp):
    t = TreeNode.read(tree_fp, format="newick")
    t.create_caches()
    t.assign_ids()
    return t


def validate_tips_in_tree(tree: TreeNode, tips):
    tip_names = set([n.name for n in tree.tips()])
    missing = [x for x in tips if x not in tip_names]
    return missing


###############################################################################
# Codon + HyPhy
###############################################################################

def find_cds_fasta_for_genome(cds_dir, genome):
    if not cds_dir:
        return None
    genome_norm = normalize_genome_id(genome)
    patterns = [
        f"{genome}*.fna", f"{genome}*.ffn", f"{genome}*.cds*.fna",
        f"{genome_norm}*.fna", f"{genome_norm}*.ffn", f"{genome_norm}*.cds*.fna",
        f"*{genome}*.fna", f"*{genome_norm}*.fna",
    ]
    for pat in patterns:
        hits = glob.glob(os.path.join(cds_dir, pat))
        if hits:
            return sorted(hits)[0]
    return None


def map_aa_id_to_cds_id(aa_id, mode="after_last_pipe"):
    if mode == "same":
        return aa_id
    if mode == "strip_pipe":
        return aa_id.split("|", 1)[-1]
    return aa_id.split("|")[-1]


def build_codon_alignment_pal2nal(aa_aln_fp, cds_subset_fp, out_codon_fp, pal2nal_cmd="pal2nal.pl", codon_format="fasta"):
    cmd = f"{pal2nal_cmd} {aa_aln_fp} {cds_subset_fp} -output {codon_format} > {out_codon_fp}"
    run_cmd(cmd, quiet=True, shell=True)


def run_hyphy(hyphy_bin, test_name, codon_aln_fp, tree_fp, out_json):
    cmd = [hyphy_bin, test_name, "--alignment", codon_aln_fp, "--tree", tree_fp, "--output", out_json]
    try:
        run_cmd(cmd, quiet=True, shell=False)
        return True
    except Exception:
        return False


###############################################################################
# Embeddings
###############################################################################

def _pca_fit_transform(X, n_components=3):
    from sklearn.decomposition import PCA
    pca = PCA(n_components=n_components, random_state=0)
    Z = pca.fit_transform(X)
    return Z, pca.explained_variance_ratio_.tolist()


def _embed_esm(seqs: dict, model_name: str, device: str, batch_size: int, repr_layer):
    import torch
    import esm

    model, alphabet = esm.pretrained.load_model_and_alphabet(model_name)
    model.eval()
    model = model.to(device)

    if repr_layer is None:
        repr_layer = model.num_layers

    batch_converter = alphabet.get_batch_converter()

    items = list(seqs.items())
    vectors = []
    ids = []

    with torch.no_grad():
        for i in range(0, len(items), batch_size):
            chunk = items[i:i + batch_size]
            # esm expects list of (label, seq)
            batch = [(k, v) for k, v in chunk]
            labels, strs, toks = batch_converter(batch)
            toks = toks.to(device)
            out = model(toks, repr_layers=[repr_layer], return_contacts=False)
            reps = out["representations"][repr_layer]  # (B, T, C)

            # mean pool excluding special tokens + padding
            # ESM2: tokens: [BOS] ... [EOS] plus padding
            # We'll mask padding and also drop first token (BOS) and any token == padding_idx
            pad_idx = alphabet.padding_idx
            # mask: toks != pad
            mask = (toks != pad_idx).float()
            # drop BOS (position 0)
            mask[:, 0] = 0.0
            # also drop EOS by masking positions where token == eos_idx
            eos_idx = alphabet.eos_idx
            mask = mask * (toks != eos_idx).float()

            denom = mask.sum(dim=1).clamp(min=1.0).unsqueeze(-1)  # (B,1)
            pooled = (reps * mask.unsqueeze(-1)).sum(dim=1) / denom  # (B,C)

            pooled = pooled.detach().cpu().numpy()
            for (k, _), v in zip(chunk, pooled):
                ids.append(k)
                vectors.append(v)

    import numpy as np
    X = np.vstack(vectors)
    return ids, X


def _embed_transformers(seqs: dict, model_id_or_path: str, device: str, batch_size: int):
    import torch
    from transformers import AutoTokenizer, AutoModel

    tok = AutoTokenizer.from_pretrained(model_id_or_path, do_lower_case=False)
    model = AutoModel.from_pretrained(model_id_or_path)
    model.eval()
    model = model.to(device)

    items = list(seqs.items())
    vectors = []
    ids = []

    with torch.no_grad():
        for i in range(0, len(items), batch_size):
            chunk = items[i:i + batch_size]
            labels = [k for k, _ in chunk]
            strings = [v for _, v in chunk]
            enc = tok(strings, return_tensors="pt", padding=True, truncation=True)
            enc = {k: v.to(device) for k, v in enc.items()}
            out = model(**enc)
            reps = out.last_hidden_state  # (B,T,C)

            # mean pool over non-pad tokens
            attn = enc.get("attention_mask", None)
            if attn is None:
                raise RuntimeError("transformers embedding requires attention_mask")
            mask = attn.float()  # (B,T)
            denom = mask.sum(dim=1).clamp(min=1.0).unsqueeze(-1)
            pooled = (reps * mask.unsqueeze(-1)).sum(dim=1) / denom  # (B,C)

            pooled = pooled.detach().cpu().numpy()
            for k, v in zip(labels, pooled):
                ids.append(k)
                vectors.append(v)

    import numpy as np
    X = np.vstack(vectors)
    return ids, X


def compute_embeddings_for_hmm(hmm_name: str, seqs: dict, emb_cfg: dict, outdir_embeddings: str, force: bool, clades: dict | None):
    """
    Writes:
      embeddings/<hmm>.embeddings.npy        (float32)
      embeddings/<hmm>.pca.tsv              (tip, genome, protein, PC1..)
      embeddings/<hmm>.umap.tsv             (tip, genome, protein, UMAP1, UMAP2)
      embeddings/<hmm>.pca.meta.json        (model, backend, var_explained)
      embeddings/<hmm>.dispersion.tsv       (optional; if clades provided)
      embeddings/<hmm>.vectors.tsv          (optional huge; if write_full_vectors)
    """
    safe_mkdir(outdir_embeddings)

    out_npy = os.path.join(outdir_embeddings, f"{hmm_name}.embeddings.npy")
    out_pca = os.path.join(outdir_embeddings, f"{hmm_name}.pca.tsv")
    out_umap = os.path.join(outdir_embeddings, f"{hmm_name}.umap.tsv")
    out_meta = os.path.join(outdir_embeddings, f"{hmm_name}.pca.meta.json")
    out_disp = os.path.join(outdir_embeddings, f"{hmm_name}.dispersion.tsv")
    out_vec_tsv = os.path.join(outdir_embeddings, f"{hmm_name}.vectors.tsv")

    if os.path.exists(out_pca) and os.path.exists(out_npy) and os.path.exists(out_umap) and not force:
        return

    # sanitize sequences: keep only plausible AA (PLMs can handle X; but keep as-is except whitespace)
    seqs = {k: v.replace(" ", "").replace("\n", "") for k, v in seqs.items()}
    if len(seqs) < 3:
        return

    backend = emb_cfg["backend"]
    device = emb_cfg["device"]
    batch_size = int(emb_cfg["batch_size"])
    model_name = emb_cfg["model"]
    repr_layer = emb_cfg.get("repr_layer", None)

    try:
        if backend == "esm":
            ids, X = _embed_esm(seqs, model_name=model_name, device=device, batch_size=batch_size, repr_layer=repr_layer)
        elif backend == "transformers":
            ids, X = _embed_transformers(seqs, model_id_or_path=model_name, device=device, batch_size=batch_size)
        else:
            raise SystemExit(f"Unknown embeddings.backend: {backend}")
    except Exception as e:
        print(f"[embed] FAILED {hmm_name}: {e}")
        return

    # write raw vectors
    import numpy as np
    X = X.astype(np.float32)
    np.save(out_npy, X)

    # PCA
    ncomp = int(emb_cfg["pca_components"])
    if ncomp < 2:
        ncomp = 2
    ncomp = min(ncomp, X.shape[1], max(2, X.shape[0] - 1))

    Z, var = _pca_fit_transform(X, n_components=ncomp)

    # pca tsv
    rows = []
    for tip, coords in zip(ids, Z):
        genome = tip.split("|", 1)[0] if "|" in tip else "Unknown"
        protein = tip.split("|", 1)[1] if "|" in tip else tip
        r = {"hmm": hmm_name, "tip": tip, "genome": genome, "protein": protein}
        for j in range(Z.shape[1]):
            r[f"PC{j+1}"] = float(coords[j])
        rows.append(r)
    pd.DataFrame(rows).to_csv(out_pca, sep="\t", index=False)

    # UMAP
    try:
        import umap
        # UMAP needs > min_samples, typically ok if len(ids) > 3 from earlier check
        reducer = umap.UMAP(n_components=2, random_state=42)
        # Suppress possible warnings if mostly small data
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            U = reducer.fit_transform(X)
        
        u_rows = []
        for tip, coords in zip(ids, U):
            genome = tip.split("|", 1)[0] if "|" in tip else "Unknown"
            protein = tip.split("|", 1)[1] if "|" in tip else tip
            u_rows.append({
                "hmm": hmm_name,
                "tip": tip,
                "genome": genome,
                "protein": protein,
                "UMAP1": float(coords[0]),
                "UMAP2": float(coords[1])
            })
        pd.DataFrame(u_rows).to_csv(out_umap, sep="\t", index=False)
    except ImportError:
        print(f"[embed] UMAP skip: umap-learn not installed.")
    except Exception as e:
        print(f"[embed] UMAP failed for {hmm_name}: {e}")

    # metadata
    meta = {
        "hmm": hmm_name,
        "backend": backend,
        "model": model_name,
        "device": device,
        "batch_size": batch_size,
        "repr_layer": repr_layer,
        "pooling": emb_cfg["pooling"],
        "pca_components": int(Z.shape[1]),
        "explained_variance_ratio": var,
        "n_sequences": int(len(ids)),
        "vector_dim": int(X.shape[1]),
    }
    write_json(meta, out_meta)

    # optional full vectors TSV (can be enormous)
    if emb_cfg.get("write_full_vectors", False):
        vec_rows = []
        for tip, v in zip(ids, X):
            r = {"tip": tip}
            for j, val in enumerate(v):
                r[f"d{j+1}"] = float(val)
            vec_rows.append(r)
        pd.DataFrame(vec_rows).to_csv(out_vec_tsv, sep="\t", index=False)

    # clade dispersion (functional “tightness” proxy)
    if clades:
        # clade dispersion = mean euclidean distance to clade centroid (in PCA space)
        Z_df = pd.DataFrame(Z, columns=[f"PC{j+1}" for j in range(Z.shape[1])])
        Z_df["tip"] = ids
        Z_df = Z_df.set_index("tip")

        disp_rows = []
        for cname, tips in clades.items():
            tips_in = [t for t in tips if t in Z_df.index]
            if len(tips_in) < 3:
                continue
            M = Z_df.loc[tips_in].to_numpy()
            centroid = M.mean(axis=0)
            d = np.sqrt(((M - centroid) ** 2).sum(axis=1))
            disp_rows.append({
                "hmm": hmm_name,
                "clade": cname,
                "n": int(len(tips_in)),
                "mean_dist_to_centroid": float(d.mean()),
                "median_dist_to_centroid": float(np.median(d)),
            })
        if disp_rows:
            pd.DataFrame(disp_rows).to_csv(out_disp, sep="\t", index=False)


###############################################################################
# Main
###############################################################################

def main():
    ap = argparse.ArgumentParser(description="Competitive HMM pipeline with JSON config + embeddings")
    ap.add_argument("--config", default=None, help="JSON config file")
    ap.add_argument("--dump_default_config", action="store_true", help="Print default config JSON and exit")

    # Minimal CLI overrides (everything else in JSON)
    ap.add_argument("--faa_dir", default=None, help="Override inputs.faa_dir")
    ap.add_argument("--hmm_dir", default=None, help="Override inputs.hmm_input (dir or single .hmm)")
    ap.add_argument("--outdir", default=None, help="Override output.outdir")
    ap.add_argument("--cpu", type=int, default=None, help="Override resources.cpu")
    ap.add_argument("--start_at", choices=STEPS, default=None, help="Override workflow.start_at")
    ap.add_argument("--stop_after", choices=STEPS, default=None, help="Override workflow.stop_after")
    ap.add_argument("--force", action="store_true", help="Override workflow.force=True")

    args = ap.parse_args()

    if args.dump_default_config:
        print(json.dumps(DEFAULT_CONFIG, indent=2, sort_keys=True))
        return

    cfg = deepcopy(DEFAULT_CONFIG)
    if args.config:
        deep_update(cfg, load_json_config(args.config))

    # CLI overrides
    if args.faa_dir is not None:
        cfg["inputs"]["faa_dir"] = args.faa_dir
    if args.hmm_dir is not None:
        cfg["inputs"]["hmm_input"] = args.hmm_dir
    if args.outdir is not None:
        cfg["output"]["outdir"] = args.outdir
    if args.cpu is not None:
        cfg["resources"]["cpu"] = int(args.cpu)
    if args.start_at is not None:
        cfg["workflow"]["start_at"] = args.start_at
    if args.stop_after is not None:
        cfg["workflow"]["stop_after"] = args.stop_after
    if args.force:
        cfg["workflow"]["force"] = True

    # Validate required
    faa_arg = cfg["inputs"]["faa_dir"]
    hmm_arg = cfg["inputs"]["hmm_input"]
    outdir = cfg["output"]["outdir"]
    if not faa_arg or not hmm_arg or not outdir:
        raise SystemExit("Config must specify inputs.faa_dir, inputs.hmm_input, output.outdir (or pass via CLI).")

    cpu = int(cfg["resources"]["cpu"])
    start_at = cfg["workflow"]["start_at"]
    stop_after = cfg["workflow"]["stop_after"]
    force = bool(cfg["workflow"]["force"])

    filt_cfg = cfg["filtering"]
    phy_cfg = cfg["phylo"]
    emb_cfg = cfg["embeddings"]
    post_cfg = cfg["post"]
    codon_cfg = cfg["codon"]
    hyphy_cfg = cfg["hyphy"]

    safe_mkdir(outdir)

    # output dirs
    hmmscan_dir = os.path.join(outdir, "hmmscan_tbl")
    hmmsearch_dir = os.path.join(outdir, "hmmsearch_tbl")
    fasta_dir = os.path.join(outdir, "fasta_per_hmm")
    aln_dir = os.path.join(outdir, "alignments_hmm")
    clipkit_dir = os.path.join(outdir, "alignments_clipkit")
    tree_dir = os.path.join(outdir, "trees_iqtree")
    summary_dir = os.path.join(outdir, "summary")
    post_dir = os.path.join(summary_dir, "post_scikitbio")
    codon_dir = os.path.join(outdir, "codon_alignments")
    hyphy_dir = os.path.join(summary_dir, "hyphy")
    emb_dir = os.path.join(outdir, "embeddings")

    for d in [hmmscan_dir, hmmsearch_dir, fasta_dir, aln_dir, clipkit_dir, tree_dir, summary_dir, post_dir, codon_dir, hyphy_dir, emb_dir]:
        safe_mkdir(d)

    # Provenance
    write_json(cfg, os.path.join(summary_dir, "resolved_config.json"))

    # thresholds / manifest
    thresholds_map = load_custom_thresholds(filt_cfg["scores_tsv"])
    hmm_keep = load_manifest(cfg["workflow"]["hmm_manifest"])

    ###########################################################################
    # FAA input discovery (dir or file)
    ###########################################################################
    faa_abs = os.path.abspath(faa_arg)
    if os.path.isfile(faa_abs):
        if not faa_abs.endswith(".faa"):
            raise SystemExit(f"inputs.faa_dir points to a file but not .faa: {faa_abs}")
        faa_dir = os.path.dirname(faa_abs) or "."
        genomes = [os.path.basename(faa_abs)]
    else:
        if not os.path.isdir(faa_abs):
            raise SystemExit(f"inputs.faa_dir must be a directory or a single .faa file; not found: {faa_abs}")
        faa_dir = faa_abs
        genomes = sorted(f for f in os.listdir(faa_dir) if f.endswith(".faa"))
    if not genomes:
        raise SystemExit("No .faa inputs found.")

    ###########################################################################
    # HMM input discovery (dir or file)
    ###########################################################################
    hmm_abs = os.path.abspath(hmm_arg)
    if os.path.isfile(hmm_abs):
        if not hmm_abs.endswith(".hmm"):
            raise SystemExit(f"inputs.hmm_input points to a file but not .hmm: {hmm_abs}")
        hmm_dir = os.path.dirname(hmm_abs) or "."
        hmm_files = [os.path.basename(hmm_abs)]
        hmm_input_mode = "file"
    else:
        if not os.path.isdir(hmm_abs):
            raise SystemExit(f"inputs.hmm_input must be a directory or a single .hmm file; not found: {hmm_abs}")
        hmm_dir = hmm_abs
        hmm_files = sorted(f for f in os.listdir(hmm_dir) if f.endswith(".hmm"))
        hmm_input_mode = "dir"
    if not hmm_files:
        raise SystemExit("No .hmm files found.")

    # checkpoint files
    combined_faa = os.path.join(outdir, "combined_proteomes.faa")
    combined_hmm = os.path.join(outdir, "combined.hmm")
    hits_scan_tsv = os.path.join(summary_dir, "hmmscan_hits.filtered.tsv")
    hits_search_tsv = os.path.join(summary_dir, "hmmsearch_hits.filtered.tsv")
    best_hits_tsv = os.path.join(summary_dir, "best_hits.competitive.tsv")

    ###########################################################################
    # STEP: prep
    ###########################################################################
    if step_in_range("prep", start_at, stop_after):
        print("\n[prep] Building combined proteomes FASTA and combined HMM DB...")

        if (not os.path.exists(combined_faa)) or force:
            with open(combined_faa, "w") as out_faa:
                for g in genomes:
                    seqs = read_fasta(os.path.join(faa_dir, g))
                    for pid, s in seqs.items():
                        out_faa.write(f">{g}~{pid}\n{s}\n")

        if (not os.path.exists(combined_hmm)) or force:
            if hmm_input_mode == "file":
                src = os.path.join(hmm_dir, hmm_files[0])
                run_cmd(["bash", "-c", f"cp {src} {combined_hmm}"], quiet=True)
            else:
                run_cmd(["bash", "-c", f"cat {hmm_dir}/*.hmm > {combined_hmm}"], quiet=True)
            run_cmd(["hmmpress", "-f", combined_hmm], quiet=True)

    if stop_after == "prep":
        print("\nStopped after prep.")
        return

    ###########################################################################
    # STEP: hmmer
    ###########################################################################
    scan_df = pd.DataFrame()
    search_df = pd.DataFrame()
    best_df = pd.DataFrame()

    if step_in_range("hmmer", start_at, stop_after):
        print("\n[hmmer] Running hmmscan + hmmsearch (with caching)...")

        # hmmscan per genome
        scan_tasks = []
        for g in genomes:
            tbl = os.path.join(hmmscan_dir, g + ".tbl")
            if os.path.exists(tbl) and not force:
                continue
            cmd = ["hmmscan", "--cpu", "1", "--domtblout", tbl, combined_hmm, os.path.join(faa_dir, g)]
            scan_tasks.append((cmd, tbl, g, bool(filt_cfg["keep_tbl"])))

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
            scan_df.to_csv(hits_scan_tsv, sep="\t", index=False)

        # hmmsearch per HMM file
        search_tasks = []
        for hf in hmm_files:
            hmm_name = os.path.splitext(hf)[0]
            if hmm_keep is not None and hmm_name not in hmm_keep:
                continue
            tbl = os.path.join(hmmsearch_dir, f"{hmm_name}_combined.tbl")
            if os.path.exists(tbl) and not force:
                continue
            cmd = ["hmmsearch", "--cpu", "1", "--domtblout", tbl, os.path.join(hmm_dir, hf), combined_faa]
            search_tasks.append((cmd, tbl, hmm_name, bool(filt_cfg["keep_tbl"])))

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
            best_df.to_csv(best_hits_tsv, sep="\t", index=False)

    if stop_after == "hmmer":
        print("\nStopped after hmmer.")
        return

    ###########################################################################
    # NAME -> hmm file mapping (needed for hmmalign)
    ###########################################################################
    name_to_hmm_path = {}
    for hf in hmm_files:
        fp = os.path.join(hmm_dir, hf)
        try:
            with open(fp) as f:
                for line in f:
                    if line.startswith("NAME"):
                        name_to_hmm_path[line.split()[1]] = fp
                        break
        except OSError:
            continue

    ###########################################################################
    # STEP: extract
    ###########################################################################
    hmm_to_seqs = defaultdict(dict)

    if step_in_range("extract", start_at, stop_after) or step_in_range("embed", start_at, stop_after) or step_in_range("phylo", start_at, stop_after):
        # load proteomes once (needed for extraction)
        proteome_seqs = {}
        for g in genomes:
            proteome_seqs[g] = read_fasta(os.path.join(faa_dir, g))

    if step_in_range("extract", start_at, stop_after):
        print("\n[extract] Building per-HMM FASTAs...")

        # choose base hits for membership
        if phy_cfg["use_hmmsearch_alignment"] and not search_df.empty:
            base_df = search_df
        else:
            base_df = scan_df if not scan_df.empty else search_df

        if base_df.empty:
            raise SystemExit("No hits available to extract sequences.")

        if not phy_cfg["keep_all_hits"]:
            base_df = base_df.sort_values("bitscore", ascending=False).groupby(["genome", "hmm"], as_index=False).head(1)

        for _, r in base_df.iterrows():
            hmm = r["hmm"]
            # manifest filtering is on hmm label as used in hits tables
            if hmm_keep is not None and hmm not in hmm_keep:
                continue
            genome = r["genome"]
            prot = r["protein"]
            seq = proteome_seqs.get(genome, {}).get(prot)
            if not seq:
                continue
            tip = f"{genome}|{prot}"
            hmm_to_seqs[hmm][tip] = seq

        for hmm, seqs in hmm_to_seqs.items():
            out_fp = os.path.join(fasta_dir, f"{hmm}.faa")
            if (not os.path.exists(out_fp)) or force:
                write_fasta(out_fp, seqs)

    if stop_after == "extract":
        print("\nStopped after extract.")
        return

    ###########################################################################
    # STEP: embed
    ###########################################################################
    if step_in_range("embed", start_at, stop_after) and emb_cfg.get("enabled", False):
        print("\n[embed] Computing per-HMM embeddings...")

        clades = None
        if post_cfg.get("clades_tsv", None):
            # same clades file used for post-processing; here we use it for dispersion
            try:
                clades = load_clades_tsv(post_cfg["clades_tsv"])
            except Exception:
                clades = None

        # ensure we have sequence bins (if extract step skipped in this invocation)
        if not hmm_to_seqs:
            for fp in glob.glob(os.path.join(fasta_dir, "*.faa")):
                hmm = os.path.basename(fp).rsplit(".", 1)[0]
                if hmm_keep is not None and hmm not in hmm_keep:
                    continue
                hmm_to_seqs[hmm] = read_fasta(fp)

        for hmm, seqs in hmm_to_seqs.items():
            if hmm_keep is not None and hmm not in hmm_keep:
                continue
            compute_embeddings_for_hmm(hmm, seqs, emb_cfg, emb_dir, force=force, clades=clades)

    if stop_after == "embed":
        print("\nStopped after embed.")
        return

    ###########################################################################
    # STEP: phylo
    ###########################################################################
    if step_in_range("phylo", start_at, stop_after):
        print("\n[phylo] Per-HMM alignment + trimming + IQ-TREE...")

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
            if hmm not in name_to_hmm_path:
                # If HMM NAME != hmm label used in hits table, hmmalign won't know which HMM to use.
                # In that case, either run with --phylo.mafft=true, or make your hit HMM labels match NAME fields.
                if not phy_cfg["mafft"]:
                    continue
                hmm_path = None
            else:
                hmm_path = name_to_hmm_path[hmm]

            if phy_cfg["mafft"]:
                hmm_path = hmm_path or ""  # not used in MAFFT mode

            tasks.append((hmm, seqs, fasta_dir, aln_dir, clipkit_dir, tree_dir, hmm_path, threads_per, phy_cfg, force))

        with ProcessPoolExecutor(max_workers=workers) as exe:
            list(exe.map(worker_phylo, tasks))

    if stop_after == "phylo":
        print("\nStopped after phylo.")
        return

    ###########################################################################
    # STEP: post
    ###########################################################################
    if step_in_range("post", start_at, stop_after) and post_cfg.get("enabled", False):
        if post_cfg.get("compute_conservation", False) or post_cfg.get("compute_kl", False):
            print("\n[post] scikit-bio post-processing...")

            clades = None
            if post_cfg.get("clades_tsv", None):
                clades = load_clades_tsv(post_cfg["clades_tsv"])

            if post_cfg.get("compute_kl", False) and not clades:
                raise SystemExit("post.compute_kl requires post.clades_tsv")

            kl_pairs = []
            if post_cfg.get("kl_pairs", None):
                for pair in str(post_cfg["kl_pairs"]).split(","):
                    pair = pair.strip()
                    if not pair:
                        continue
                    if ":" not in pair:
                        raise SystemExit("post.kl_pairs must look like 'A:B' or 'A:background'")
                    a, b = pair.split(":", 1)
                    kl_pairs.append((a.strip(), b.strip()))
            else:
                if post_cfg.get("compute_kl", False) and clades and len(clades) >= 2:
                    keys = list(clades.keys())
                    for i in range(len(keys)):
                        for j in range(i + 1, len(keys)):
                            kl_pairs.append((keys[i], keys[j]))
                if post_cfg.get("compute_kl", False) and clades and len(clades) == 1:
                    k = list(clades.keys())[0]
                    kl_pairs.append((k, "background"))

            cons_rows = []
            kl_rows = []
            mrca_rows = []

            hmm_names = sorted([os.path.basename(x).split(".")[0] for x in glob.glob(os.path.join(tree_dir, "*.treefile"))])
            if hmm_keep is not None:
                hmm_names = [h for h in hmm_names if h in hmm_keep]

            for hmm in hmm_names:
                tree_fp = os.path.join(tree_dir, f"{hmm}.treefile")

                clip_aln = os.path.join(clipkit_dir, f"{hmm}.clipkit.faa")
                if os.path.exists(clip_aln):
                    aln_fp = clip_aln
                else:
                    cand1 = os.path.join(aln_dir, f"{hmm}.afa")
                    cand2 = os.path.join(aln_dir, f"{hmm}.mafft.fasta")
                    aln_fp = cand2 if os.path.exists(cand2) else cand1 if os.path.exists(cand1) else None

                if not aln_fp or not os.path.exists(tree_fp):
                    continue

                aln_seqs = read_fasta(aln_fp)

                if post_cfg.get("compute_conservation", False):
                    try:
                        cons = consurf_like_scores(aln_seqs, metric=post_cfg.get("conservation_metric", "inverse_shannon_uncertainty"))
                        for i, v in enumerate(cons, start=1):
                            cons_rows.append({"hmm": hmm, "scope": "global", "site_1based": i, "conservation": float(v)})
                    except Exception:
                        pass

                if post_cfg.get("compute_kl", False) and clades:
                    counts = {}
                    counts["background"] = site_counts_from_subset(aln_seqs, subset_tips=None)
                    for cname, tips in clades.items():
                        counts[cname] = site_counts_from_subset(aln_seqs, subset_tips=tips)

                    L = len(counts["background"])
                    for a, b in kl_pairs:
                        if a not in counts or b not in counts:
                            continue
                        if len(counts[a]) != L or len(counts[b]) != L:
                            continue
                        for i in range(L):
                            kl = kl_divergence(counts[a][i], counts[b][i])
                            kl_rows.append({
                                "hmm": hmm,
                                "pair": f"{a}:{b}",
                                "site_1based": i + 1,
                                "kl_bits": kl,
                                "nA": sum(counts[a][i].values()),
                                "nB": sum(counts[b][i].values())
                            })

                # MRCA sanity only (no internal node mapping)
                if clades:
                    try:
                        t = tree_load(tree_fp)
                        for cname, tips in clades.items():
                            tips2 = [x for x in tips if x in aln_seqs]
                            missing = validate_tips_in_tree(t, tips2)
                            if missing:
                                mrca_rows.append({"hmm": hmm, "clade": cname, "status": "missing_tips", "detail": ",".join(missing[:10])})
                                continue
                            mrca = t.lca(tips2)
                            mrca_rows.append({"hmm": hmm, "clade": cname, "status": "ok", "mrca_id": mrca.id, "mrca_name": mrca.name or ""})
                    except Exception:
                        pass

            if cons_rows:
                pd.DataFrame(cons_rows).to_csv(os.path.join(post_dir, "conservation.tsv"), sep="\t", index=False)
            if kl_rows:
                pd.DataFrame(kl_rows).to_csv(os.path.join(post_dir, "kl_divergence.tsv"), sep="\t", index=False)
            if mrca_rows:
                pd.DataFrame(mrca_rows).to_csv(os.path.join(post_dir, "mrca_sanity.tsv"), sep="\t", index=False)

    if stop_after == "post":
        print("\nStopped after post.")
        return

    ###########################################################################
    # STEP: codon
    ###########################################################################
    if step_in_range("codon", start_at, stop_after) and codon_cfg.get("enabled", False):
        if codon_cfg.get("build_codon_alignments", False):
            cds_dir = cfg["inputs"].get("cds_dir", None)
            if not cds_dir:
                raise SystemExit("codon.build_codon_alignments requires inputs.cds_dir")
            print("\n[codon] Building per-HMM codon alignments (PAL2NAL)...")

            hmm_names = sorted([os.path.basename(x).split(".")[0] for x in glob.glob(os.path.join(tree_dir, "*.treefile"))])
            if hmm_keep is not None:
                hmm_names = [h for h in hmm_names if h in hmm_keep]

            for hmm in hmm_names:
                tree_fp = os.path.join(tree_dir, f"{hmm}.treefile")
                clip_aln = os.path.join(clipkit_dir, f"{hmm}.clipkit.faa")
                if os.path.exists(clip_aln):
                    aa_aln_fp = clip_aln
                else:
                    cand1 = os.path.join(aln_dir, f"{hmm}.afa")
                    cand2 = os.path.join(aln_dir, f"{hmm}.mafft.fasta")
                    aa_aln_fp = cand2 if os.path.exists(cand2) else cand1 if os.path.exists(cand1) else None
                if not aa_aln_fp or not os.path.exists(tree_fp):
                    continue

                codon_aln_fp = os.path.join(codon_dir, f"{hmm}.codon.fasta")
                if os.path.exists(codon_aln_fp) and not force:
                    continue

                aln_seqs = read_fasta(aa_aln_fp)
                cds_subset = {}
                tips_by_genome = defaultdict(list)
                for tip in aln_seqs.keys():
                    genome = tip.split("|", 1)[0] if "|" in tip else "Unknown"
                    tips_by_genome[genome].append(tip)

                for genome, tips in tips_by_genome.items():
                    cds_fp = find_cds_fasta_for_genome(cds_dir, genome)
                    if not cds_fp:
                        continue
                    cds_all = read_fasta(cds_fp)

                    for tip in tips:
                        cds_id = map_aa_id_to_cds_id(tip, mode=codon_cfg.get("cds_id_mode", "after_last_pipe"))
                        candidates = [cds_id, cds_id.split("|")[-1], cds_id.replace("|", "_")]
                        found = None
                        for c in candidates:
                            if c in cds_all:
                                found = c
                                break
                        if found:
                            cds_subset[tip] = cds_all[found]

                if len(cds_subset) < 3:
                    continue

                cds_subset_fp = os.path.join(codon_dir, f"{hmm}.cds.fna")
                write_fasta(cds_subset_fp, cds_subset)

                try:
                    build_codon_alignment_pal2nal(
                        aa_aln_fp, cds_subset_fp, codon_aln_fp,
                        pal2nal_cmd=codon_cfg.get("pal2nal_cmd", "pal2nal.pl")
                    )
                except Exception:
                    continue

    if stop_after == "codon":
        print("\nStopped after codon.")
        return

    ###########################################################################
    # STEP: hyphy
    ###########################################################################
    if step_in_range("hyphy", start_at, stop_after) and hyphy_cfg.get("enabled", False):
        if hyphy_cfg.get("run_hyphy", False):
            print("\n[hyphy] Running HyPhy tests (generic wrapper)...")
            tests = [x.strip() for x in str(hyphy_cfg.get("hyphy_tests", "")).split(",") if x.strip()]

            hmm_names = sorted([os.path.basename(x).split(".")[0] for x in glob.glob(os.path.join(codon_dir, "*.codon.fasta"))])
            if hmm_keep is not None:
                hmm_names = [h for h in hmm_names if h in hmm_keep]

            for hmm in hmm_names:
                codon_aln_fp = os.path.join(codon_dir, f"{hmm}.codon.fasta")
                tree_fp = os.path.join(tree_dir, f"{hmm}.treefile")
                if not os.path.exists(tree_fp):
                    continue

                for test in tests:
                    out_json = os.path.join(hyphy_dir, f"{hmm}.{test}.json")
                    if os.path.exists(out_json) and not force:
                        continue
                    _ = run_hyphy(hyphy_cfg.get("hyphy_bin", "hyphy"), test, codon_aln_fp, tree_fp, out_json)

    print("\nPipeline complete.")


if __name__ == "__main__":
    main()

import glob
import logging
import os
import time

import pandas as pd

from .constants import STEPS
from .qc import generate_qc_summary, write_run_manifest
from .tasks import asr, codon, curate, embed, extract, hmmer, hyphy, phylo, post, prep, synteny
from .methods import run_evidence_join, run_ha_sites, run_regime_shift
from .utils.bio import read_fasta
from .utils.helpers import safe_mkdir, write_json
from .utils.logging_utils import setup_logging, step_timer


def step_in_range(step, start_at, stop_after):
    i = STEPS.index(step)
    i0 = STEPS.index(start_at) if start_at else 0
    i1 = STEPS.index(stop_after) if stop_after else len(STEPS) - 1
    return i0 <= i <= i1


def load_manifest(hmm_manifest_fp):
    if not hmm_manifest_fp:
        return None
    keep = set()
    with open(hmm_manifest_fp) as f:
        for line in f:
            s = line.strip()
            if s and not s.startswith("#"):
                keep.add(s)
    return keep


def _load_proteomes_lazy(genomes, faa_dir):
    return {g: read_fasta(os.path.join(faa_dir, g)) for g in genomes}


def _mark_done(outdir: str, step: str, ok: bool = True):
    marker = os.path.join(outdir, "summary", f"{step}.{'done' if ok else 'failed'}")
    with open(marker, "w") as f:
        f.write(str(time.time()))


def run_pipeline(cfg):
    outdir = cfg["output"]["outdir"]
    safe_mkdir(outdir)
    logger = setup_logging(outdir, cfg.get("logging", {}).get("level", "INFO"))
    logger.info("Starting phylofoundry run")

    faa_arg = cfg["inputs"]["faa_dir"]
    hmm_arg = cfg["inputs"]["hmm_input"]
    start_at = cfg["workflow"]["start_at"]
    stop_after = cfg["workflow"]["stop_after"]
    force = bool(cfg["workflow"]["force"])

    phy_cfg = cfg["phylo"]
    emb_cfg = cfg["embeddings"]
    post_cfg = cfg["post"]
    synteny_cfg = cfg["synteny"]
    codon_cfg = cfg["codon"]
    hyphy_cfg = cfg["hyphy"]
    motif_cfg = cfg.get("motifs", {})
    discover_cfg = cfg.get("discover", {})
    regime_cfg = cfg.get("regime_shift", {})
    ha_cfg = cfg.get("ha", {})
    evidence_cfg = cfg.get("evidence_join", {})

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
    clade_assign_dir = os.path.join(outdir, "clade_assignments")

    for d in [hmmscan_dir, hmmsearch_dir, fasta_dir, aln_dir, clipkit_dir, tree_dir, summary_dir, post_dir, codon_dir, hyphy_dir, emb_dir, clade_assign_dir]:
        safe_mkdir(d)

    write_json(cfg, os.path.join(summary_dir, "resolved_config.json"))
    hmm_keep = load_manifest(cfg["workflow"].get("hmm_manifest"))

    faa_abs = os.path.abspath(faa_arg)
    if os.path.isfile(faa_abs):
        if not faa_abs.endswith(".faa"):
            raise SystemExit(f"inputs.faa_dir points to a file but not .faa: {faa_abs}")
        faa_dir = os.path.dirname(faa_abs) or "."
        genomes = [os.path.basename(faa_abs)]
    else:
        if not os.path.isdir(faa_abs):
            raise SystemExit(f"inputs.faa_dir must be a directory or a single .faa file: {faa_abs}")
        faa_dir = faa_abs
        genomes = sorted(f for f in os.listdir(faa_dir) if f.endswith(".faa"))
    if not genomes:
        raise SystemExit("No .faa inputs found.")

    hmm_abs = os.path.abspath(hmm_arg)
    if os.path.isfile(hmm_abs):
        if not hmm_abs.endswith(".hmm"):
            raise SystemExit(f"inputs.hmm_input points to a file but not .hmm: {hmm_abs}")
        hmm_dir = os.path.dirname(hmm_abs) or "."
        hmm_files = [os.path.basename(hmm_abs)]
        hmm_input_mode = "file"
    else:
        if not os.path.isdir(hmm_abs):
            raise SystemExit(f"inputs.hmm_input must be a directory or .hmm file: {hmm_abs}")
        hmm_dir = hmm_abs
        hmm_files = sorted(f for f in os.listdir(hmm_dir) if f.endswith(".hmm"))
        hmm_input_mode = "dir"
    if not hmm_files:
        raise SystemExit("No .hmm files found.")

    step_status = []

    def run_step(name, fn):
        if not step_in_range(name, start_at, stop_after):
            return None
        try:
            with step_timer(logger, name):
                res = fn()
            step_status.append({"step": name, "status": "success", "message": "ok"})
            _mark_done(outdir, name, ok=True)
            return res
        except Exception as e:
            logger.exception("Step failed: %s", name)
            step_status.append({"step": name, "status": "failed", "message": str(e)})
            _mark_done(outdir, name, ok=False)
            raise

    combined_faa = os.path.join(outdir, "combined_proteomes.faa")
    combined_hmm = os.path.join(outdir, "combined.hmm")

    run_step("prep", lambda: prep.run_prep(cfg, genomes, faa_dir, hmm_input_mode, hmm_dir, hmm_files, combined_faa, combined_hmm, force))
    if stop_after == "prep":
        return

    scan_df = search_df = best_df = None

    def _hmmer():
        return hmmer.run_hmmer(cfg, genomes, faa_dir, hmm_files, hmm_dir, combined_hmm, combined_faa, outdir, summary_dir, hmmscan_dir, hmmsearch_dir, hmm_keep, force)

    hmmer_res = run_step("hmmer", _hmmer)
    if hmmer_res:
        scan_df, search_df, best_df = hmmer_res
        if cfg.get("prep", {}).get("cleanup_combined_faa", False) and os.path.exists(combined_faa):
            os.remove(combined_faa)

    if stop_after == "hmmer":
        return

    if best_df is None:
        best_fp = os.path.join(summary_dir, "best_hits.tsv")
        if os.path.exists(best_fp):
            best_df = pd.read_csv(best_fp, sep="\t")
        else:
            raise SystemExit("Missing best_hits.tsv required for extract step")

    # ID audit for mapping consistency
    if best_df is not None and not best_df.empty and "protein" in best_df.columns:
        from .utils.id_normalization import parse_sequence_id

        rows = []
        mode = cfg.get("codon", {}).get("cds_id_mode", "after_last_pipe")
        for _, r in best_df.iterrows():
            nid = parse_sequence_id(str(r.get("protein")), mode=mode)
            rows.append({
                "original_header": nid.original,
                "normalized_id": nid.normalized,
                "genome_id": nid.genome_id,
                "protein_id": nid.protein_id,
                "where_used": "faa,embedding,cds",
            })
        pd.DataFrame(rows).drop_duplicates().to_csv(os.path.join(summary_dir, "id_audit.tsv"), sep="\t", index=False)

    hmm_to_seqs = run_step("extract", lambda: extract.run_extract(cfg, best_df, faa_dir, fasta_dir, hmm_keep, force, summary_dir=summary_dir))
    if stop_after == "extract":
        return

    tax_map = {}
    if cfg["inputs"].get("gtdb_dir") or cfg["inputs"].get("taxonomy_file"):
        tax_map = post._load_taxonomy(cfg["inputs"].get("gtdb_dir"), cfg["inputs"].get("taxonomy_file"))

    if step_in_range("embed", start_at, stop_after) and emb_cfg.get("enabled", False):
        run_step("embed", lambda: embed.run_embed(cfg, hmm_to_seqs, None, emb_dir, fasta_dir, hmm_keep, force, summary_dir=summary_dir, tax_map=tax_map))
    if stop_after == "embed":
        return

    run_step("phylo", lambda: phylo.run_phylo(cfg, hmm_to_seqs, fasta_dir, aln_dir, clipkit_dir, tree_dir, {}, hmm_keep, force))

    if not phy_cfg.get("no_asr", False):
        run_step("asr", lambda: asr.run_asr_parse(cfg, tree_dir, fasta_dir, hmm_keep, force))

    if stop_after == "phylo":
        return

    if step_in_range("curate", start_at, stop_after):
        run_step("curate", lambda: curate.run_curate(cfg, tree_dir, fasta_dir, clipkit_dir, emb_dir, summary_dir, hmm_keep, force))
    if step_in_range("post", start_at, stop_after) and post_cfg.get("enabled", False):
        run_step("post", lambda: post.run_post(cfg, tree_dir, clipkit_dir, aln_dir, post_dir, summary_dir, hmm_keep, force, clade_assign_dir=clade_assign_dir))
    if step_in_range("synteny", start_at, stop_after) and synteny_cfg.get("enabled", False):
        run_step("synteny", lambda: synteny.run_synteny(cfg, os.path.join(outdir, "synteny"), tree_dir, scan_df, search_df, hmm_keep, force, clade_assign_dir=clade_assign_dir))
    if step_in_range("codon", start_at, stop_after) and codon_cfg.get("enabled", False):
        run_step("codon", lambda: codon.run_codon(cfg, tree_dir, clipkit_dir, aln_dir, codon_dir, hmm_keep, force))
    if step_in_range("hyphy", start_at, stop_after) and hyphy_cfg.get("enabled", False):
        run_step("hyphy", lambda: hyphy.run_hyphy(cfg, codon_dir, tree_dir, hyphy_dir, hmm_keep, force, clade_assign_dir=clade_assign_dir))
    if step_in_range("score_motifs", start_at, stop_after) and motif_cfg.get("enabled", False):
        from .tasks import motifs
        run_step("score_motifs", lambda: motifs.score_motifs(cfg, fasta_dir, summary_dir, hmm_keep, force))
    if step_in_range("regime_shift", start_at, stop_after) and regime_cfg.get("enable", False):
        run_step("regime_shift", lambda: run_regime_shift(cfg, tree_dir, emb_dir, summary_dir, hmm_keep=hmm_keep))
    if step_in_range("ha_sites", start_at, stop_after) and ha_cfg.get("enabled", False):
        run_step("ha_sites", lambda: run_ha_sites(cfg, clipkit_dir if os.path.exists(clipkit_dir) else fasta_dir, emb_dir, summary_dir, os.path.join(outdir, "qc"), hmm_keep=hmm_keep))

    if step_in_range("discover_motifs", start_at, stop_after) and discover_cfg.get("enabled", False):
        from .tasks import discover
        ha_req = os.path.join(summary_dir, "ha_sites.tsv")
        if not os.path.exists(ha_req):
            raise SystemExit(f"DiscoverCandidates requires HA output: {ha_req}")
        run_step("discover_motifs", lambda: discover.discover_motifs(cfg, fasta_dir, summary_dir, hmm_keep, force, clade_assign_dir=clade_assign_dir))

    if step_in_range("evidence_join", start_at, stop_after) and evidence_cfg.get("enable", False):
        run_step("evidence_join", lambda: run_evidence_join(cfg, summary_dir))

    if step_in_range("qc_report", start_at, stop_after) and cfg.get("qc", {}).get("enable", True):
        generate_qc_summary(outdir, summary_dir, cfg=cfg)

    if not os.path.exists(os.path.join(summary_dir, "qc_manifest.tsv")):
        pd.DataFrame([]).to_csv(os.path.join(summary_dir, "qc_manifest.tsv"), sep="\t", index=False)

    write_run_manifest(outdir, summary_dir, step_status)
    logger.info("Pipeline complete")

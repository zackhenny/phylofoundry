import os
import glob
import json
import datetime
import shutil
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt
from Bio import Phylo, SeqIO
from .artifact_paths import ArtifactPaths

st.set_page_config(page_title="Phylofoundry Curation", layout="wide")
st.title("Phylofoundry Interactive Curation")
st.markdown("Interactively prune phylogenetic trees and label branches for HyPhy testing.")

config_path = st.sidebar.text_input("Path to project 'config.json'", value="")

if not config_path or not os.path.exists(config_path):
    st.info("Waiting for config.json path... Enter the absolute path to your config file in the sidebar.")
    st.stop()

with open(config_path) as f:
    cfg = json.load(f)

outdir = cfg.get("output", {}).get("outdir", "")
if not outdir or not os.path.exists(outdir):
    st.error(f"Invalid outdir from config: {outdir}")
    st.stop()

# ── Resolve canonical artifact paths ──────────────────────────────────────
_paths = ArtifactPaths(outdir)

# ── Raw (immutable) pipeline output directories ────────────────────────────
raw_tree_dir    = _paths.trees_dir
raw_fasta_dir   = _paths.fasta_dir
raw_clipkit_dir = _paths.alignments_clipkit_dir

# ── Curated overlay directories ────────────────────────────────────────────
curated_dir         = _paths.curated_dir
cur_tree_dir        = _paths.curated_trees_dir
cur_fasta_dir       = _paths.curated_fasta_dir
cur_clipkit_dir     = _paths.curated_alignments_dir
cur_manifest_dir    = _paths.curated_manifests_dir
cur_audit_dir       = _paths.curated_audit_dir

for _d in [cur_tree_dir, cur_fasta_dir, cur_clipkit_dir,
           cur_manifest_dir, cur_audit_dir]:
    os.makedirs(_d, exist_ok=True)

# ── Discover trees (curated overlay first, fall back to raw) ───────────────
# A tree is listed if it exists in the raw directory.  Curated presence is
# shown as a status indicator per HMM but does not filter the list.
raw_trees = sorted(glob.glob(os.path.join(raw_tree_dir, "*.treefile")))
hmms = [os.path.basename(t).replace(".treefile", "") for t in raw_trees]

if not hmms:
    st.warning("No .treefile files found in outdir/trees_iqtree")
    st.stop()

st.sidebar.markdown("---")
selected_hmm = st.sidebar.selectbox("Select Gene Family (HMM):", hmms)

raw_tree_fp      = os.path.join(raw_tree_dir,    f"{selected_hmm}.treefile")
curated_tree_fp  = os.path.join(cur_tree_dir,    f"{selected_hmm}.treefile")

# Prefer curated tree for display if it exists; otherwise show raw tree.
active_tree_fp = curated_tree_fp if os.path.exists(curated_tree_fp) else raw_tree_fp

# Read current tree (curated or raw)
try:
    tree = Phylo.read(active_tree_fp, "newick")
except Exception as e:
    st.error(f"Error reading tree: {e}")
    st.stop()

# --- Load Taxonomy if available ---
tax_map = {}
tax_file = cfg.get("inputs", {}).get("taxonomy_file", "")
gtdb_dir = cfg.get("inputs", {}).get("gtdb_dir", "")

def load_tax(gtdb_dir, tax_file):
    mapping = {}
    if tax_file and os.path.exists(tax_file):
        try:
            df = pd.read_csv(tax_file, sep="\t", header=None)
            for _, r in df.iterrows():
                mapping[r[0]] = r[1]
        except: pass
    elif gtdb_dir and os.path.exists(gtdb_dir):
        for fp in glob.glob(os.path.join(gtdb_dir, "*.tsv")):
            try:
                df = pd.read_csv(fp, sep="\t", low_memory=False)
                if "user_genome" in df.columns and "classification" in df.columns:
                    for _, r in df.iterrows():
                        mapping[r["user_genome"]] = r["classification"]
            except: pass
    return mapping

tax_map = load_tax(gtdb_dir, tax_file)

# Append taxonomy to tree terminals for display
display_tree = Phylo.read(active_tree_fp, "newick")
for term in display_tree.get_terminals():
    genome = term.name.split("|")[0] if "|" in term.name else term.name
    genome = genome.replace(".faa", "").replace(".fna", "")
    tax = tax_map.get(genome, "")
    if tax:
        # Just grab the lowest reliable taxonomic level for brevity in the plot
        short_tax = tax.split(";")[-1]
        term.name = f"{term.name} [{short_tax}]"

terminals = [t.name for t in tree.get_terminals()]  # keep original names for the logic
display_terminals = [t.name for t in display_tree.get_terminals()]

st.subheader(f"Curation: {selected_hmm}")
if os.path.exists(curated_tree_fp):
    st.info("This tree has a curated overlay. Raw source is preserved in trees_iqtree/.")
else:
    st.info("No curated overlay yet. Displaying raw tree from trees_iqtree/.")

col1, col2 = st.columns([2, 1])

with col1:
    fig, ax = plt.subplots(figsize=(12, max(5, len(terminals) * 0.25)))
    Phylo.draw(display_tree, axes=ax, do_show=False)
    st.pyplot(fig)

with col2:
    st.markdown("### Action Panel")
    st.write(f"Total taxa: **{len(terminals)}**")

    to_drop = st.multiselect("Select Taxa to PRUNE (Drop from tree & downstream data):", options=terminals)
    to_test = st.multiselect("Select Taxa to mark as {Test} for HyPhy RELAX:", options=[t for t in terminals if t not in to_drop])

    if st.button("Save Curated Files", type="primary"):
        # Always read from the raw source to ensure idempotency.
        curated_tree = Phylo.read(raw_tree_fp, "newick")

        for drop_taxa in to_drop:
            try:
                curated_tree.prune(drop_taxa)
            except ValueError:
                pass

        # Add {Test} labels
        for term in curated_tree.get_terminals():
            term.name = term.name.replace("{Test}", "")
            if term.name in to_test:
                term.name = term.name + "{Test}"

        # Write curated tree to overlay directory (raw file untouched)
        Phylo.write(curated_tree, curated_tree_fp, "newick")

        # Pruned FASTA — read from raw, write to curated overlay
        raw_fasta_fp = os.path.join(raw_fasta_dir, f"{selected_hmm}.fasta")
        out_fasta_fp = os.path.join(cur_fasta_dir,  f"{selected_hmm}.fasta")
        if os.path.exists(raw_fasta_fp):
            recs = [r for r in SeqIO.parse(raw_fasta_fp, "fasta") if r.id not in to_drop]
            SeqIO.write(recs, out_fasta_fp, "fasta")

        # Pruned alignment — read from raw, write to curated overlay
        raw_aln_fp = os.path.join(raw_clipkit_dir, f"{selected_hmm}.clipkit.faa")
        out_aln_fp = os.path.join(cur_clipkit_dir,  f"{selected_hmm}.clipkit.faa")
        if os.path.exists(raw_aln_fp):
            recs = [r for r in SeqIO.parse(raw_aln_fp, "fasta") if r.id not in to_drop]
            SeqIO.write(recs, out_aln_fp, "fasta")

        # Write per-HMM manifest
        retained = {t.name.replace("{Test}", "") for t in curated_tree.get_terminals()}
        timestamp = datetime.datetime.now(datetime.timezone.utc).isoformat()
        manifest = {
            "hmm": selected_hmm,
            "timestamp": timestamp,
            "source": "manual",
            "raw_tree": raw_tree_fp,
            "curated_tree": curated_tree_fp,
            "filters_applied": ["manual_prune"],
            "dropped_taxa": sorted(to_drop),
            "test_taxa": sorted(to_test),
            "retained_taxa_count": len(retained),
            "was_modified": bool(to_drop or to_test),
        }
        manifest_fp = os.path.join(cur_manifest_dir, f"{selected_hmm}.json")
        with open(manifest_fp, "w") as _mf:
            json.dump(manifest, _mf, indent=2)

        # Append to audit log
        audit_fp = os.path.join(cur_audit_dir, "manual_curation.jsonl")
        with open(audit_fp, "a") as _af:
            _af.write(json.dumps(manifest) + "\n")

        st.success(
            f"Curated files saved to {curated_dir}. "
            "Raw pipeline outputs are preserved. Ready for downstream steps!"
        )

        try:
            st.rerun()
        except AttributeError:
            st.experimental_rerun()

st.markdown("---")
st.subheader("Alignment Preview")
# Prefer curated alignment if available, otherwise raw
active_aln_fp = (
    os.path.join(cur_clipkit_dir, f"{selected_hmm}.clipkit.faa")
    if os.path.exists(os.path.join(cur_clipkit_dir, f"{selected_hmm}.clipkit.faa"))
    else os.path.join(raw_clipkit_dir, f"{selected_hmm}.clipkit.faa")
)
if os.path.exists(active_aln_fp):
    seqs = list(SeqIO.parse(active_aln_fp, "fasta"))
    if seqs:
        df = pd.DataFrame([{"ID": s.id, "Seq (1-100)": str(s.seq)[:100] + "..."} for s in seqs])
        st.dataframe(df, use_container_width=True)

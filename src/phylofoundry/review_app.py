import os
import glob
import json
import shutil
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt
from Bio import Phylo, SeqIO

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

tree_dir = os.path.join(outdir, "trees_iqtree")
fasta_dir = os.path.join(outdir, "fasta_per_hmm")
clipkit_dir = os.path.join(outdir, "alignments_clipkit")

trees = sorted(glob.glob(os.path.join(tree_dir, "*.treefile")))
# Filter out backups
trees = [t for t in trees if not t.endswith(".raw.treefile")]

hmms = [os.path.basename(t).replace(".treefile", "") for t in trees]

if not hmms:
    st.warning("No .treefile files found in outdir/trees_iqtree")
    st.stop()

st.sidebar.markdown("---")
selected_hmm = st.sidebar.selectbox("Select Gene Family (HMM):", hmms)

tree_fp = os.path.join(tree_dir, f"{selected_hmm}.treefile")
raw_tree_fp = os.path.join(tree_dir, f"{selected_hmm}.raw.treefile")

# Read current tree (which might be already curated, or the raw one)
try:
    tree = Phylo.read(tree_fp, "newick")
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

terminals = [t.name for t in tree.get_terminals()] # keep original names for the logic
display_terminals = [t.name for t in display_tree.get_terminals()]

st.subheader(f"Curation: {selected_hmm}")
if os.path.exists(raw_tree_fp):
    st.info("This tree has been previously curated. The original backup exists.")

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
        # Create backups if they don't exist
        for d, ext in [(tree_dir, ".treefile"), (fasta_dir, ".fasta"), (clipkit_dir, ".clipkit.faa")]:
            orig = os.path.join(d, f"{selected_hmm}{ext}")
            backup = os.path.join(d, f"{selected_hmm}.raw{ext}")
            if os.path.exists(orig) and not os.path.exists(backup):
                shutil.copy2(orig, backup)
                
        # 1. Prune Tree
        # Re-read raw tree to ensure idempotency if user drops/undrops taxa repeatedly
        if os.path.exists(raw_tree_fp):
            curated_tree = Phylo.read(raw_tree_fp, "newick")
        else:
            curated_tree = Phylo.read(tree_fp, "newick")
            
        for drop_taxa in to_drop:
            try:
                curated_tree.prune(drop_taxa)
            except ValueError:
                pass
                
        # 2. Add {Test} labels
        for term in curated_tree.get_terminals():
            # First remove any existing {Test} just in case
            term.name = term.name.replace("{Test}", "")
            if term.name in to_test:
                term.name = term.name + "{Test}"
                
        Phylo.write(curated_tree, tree_fp, "newick")
        
        # 3. Prune fasta
        ref_fasta = os.path.join(fasta_dir, f"{selected_hmm}.raw.fasta")
        if not os.path.exists(ref_fasta): ref_fasta = os.path.join(fasta_dir, f"{selected_hmm}.fasta")
        out_fasta = os.path.join(fasta_dir, f"{selected_hmm}.fasta")
        if os.path.exists(ref_fasta):
            recs = [r for r in SeqIO.parse(ref_fasta, "fasta") if r.id not in to_drop]
            SeqIO.write(recs, out_fasta, "fasta")
            
        # 4. Prune alignment
        ref_aln = os.path.join(clipkit_dir, f"{selected_hmm}.raw.clipkit.faa")
        if not os.path.exists(ref_aln): ref_aln = os.path.join(clipkit_dir, f"{selected_hmm}.clipkit.faa")
        out_aln = os.path.join(clipkit_dir, f"{selected_hmm}.clipkit.faa")
        if os.path.exists(ref_aln):
            recs = [r for r in SeqIO.parse(ref_aln, "fasta") if r.id not in to_drop]
            SeqIO.write(recs, out_aln, "fasta")
            
        st.success("Successfully pruned tree and synced sequences. Ready for module B (start_at: codon)!")
        
        try:
            st.rerun()
        except AttributeError:
            st.experimental_rerun()

st.markdown("---")
st.subheader("Alignment Preview")
# Show a snippet of the alignment
ref_aln = os.path.join(clipkit_dir, f"{selected_hmm}.clipkit.faa")
if os.path.exists(ref_aln):
    seqs = list(SeqIO.parse(ref_aln, "fasta"))
    if seqs:
        df = pd.DataFrame([{"ID": s.id, "Seq (1-100)": str(s.seq)[:100] + "..."} for s in seqs])
        st.dataframe(df, use_container_width=True)

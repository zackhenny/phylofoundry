# PhyloFoundry

**PhyloFoundry** is a robust, HPC-ready bioinformatics pipeline for competitive HMM analysis, phylogenetics, and optional protein language model embeddings. It automates the journey from raw proteomes and HMM profiles to publication-ready phylogenetic trees and functional landscape metrics.

---

## 🚀 Features

-   **Competitive HMM Hits**: Uses both `hmmscan` and `hmmsearch` to identify the best functional assignments for proteins, resolving overlapping hits competitively by bitscore.
-   **Automated Phylogeny**: Per-HMM alignment (MAFFT/HMMER), trimming (ClipKit), and tree inference (IQ-TREE).
-   **Protein Embeddings** (Optional): Generates per-HMM embeddings (ESM) and dimensionality reduction (PCA/UMAP), with HDBSCAN clustering and UMAP scatter plots.
-   **Ancestral Sequence Reconstruction**: Parses IQ-TREE `.state` files to reconstruct ancestral protein sequences, embeds them alongside modern sequences, and visualizes evolutionary trajectories in UMAP space.
-   **Combined Tree Mode**: `--combined` flag to build a single tree from all HMM hits, with combined embeddings and clustering.
-   **Motif Scoring** (Optional): Uses ESM-2 attention weights to score structurally important motifs (e.g., `--motifs HPEVY,HPEVF`).
-   **Motif Discovery** (Optional): Compares attention profiles across all HDBSCAN clades to discover novel structural hubs natively in a 1-vs-All manner.
-   **Synteny Analysis** (Optional): Extracts gene neighborhoods (configurable window), computes similarity (DIAMOND/MMseqs2), and plots synteny tracks ordered by phylogeny.
-   **HDBSCAN Clustering** (Optional): Clusters protein embeddings and outputs `clade_assignment.tsv` with taxonomy.
-   **GTDB Taxonomy Integration**: Merges GTDB-Tk taxonomy into summary tables and cluster assignments.
-   **Resumable**: Smart checkpointing skips already completed steps.
-   **HPC Ready**: Auto-detects Slurm CPU allocations.

---

## 🛠️ Installation

### Option A: Conda (Recommended)

```bash
# 1. Clone the repository
git clone https://github.com/yourusername/phylofoundry.git
cd phylofoundry

# 2. Create the environment
conda env create -f environment.yml
conda activate phylofoundry

# 3. Install the package
pip install -e .
```

### Option B: Docker

Build the image locally:

```bash
docker build -t phylofoundry:latest .
```

### Option C: Apptainer / Singularity (HPC)

Convert the Docker image to an Apptainer (Singularity) image file (`.sif`) for use on HPC systems.

**Method 1: Pull from Docker Daemon (if you built it locally)**
```bash
# Save docker image to tarball
docker save phylofoundry:latest -o phylofoundry.tar
# Build SIF from tarball
apptainer build phylofoundry.sif docker-archive://phylofoundry.tar
```

**Method 2: Build from Recipe (Definition File)**
*Create a `PhyloFoundry.def` file based on the Dockerfile if needed, but converting from Docker is usually easier.*

---

## 📦 Dependencies & External Tools

While the Conda environment installs all necessary software, you can also install these tools manually if needed.

### Core Dependencies (Required)

| Tool | Purpose | Conda Package | Manual Installation / Source |
| :--- | :--- | :--- | :--- |
| **HMMER** | Searching profiles against proteomes. | `hmmer` | [hmmer.org](http://hmmer.org/) |
| **MAFFT** | Multiple sequence alignment. | `mafft` | [mafft.cbrc.jp](https://mafft.cbrc.jp/alignment/software/) |
| **ClipKit** | Alignment trimming. | `clipkit` | `pip install clipkit` |
| **IQ-TREE** | Phylogenetic tree inference. | `iqtree` | [iqtree.org](http://www.iqtree.org/) |

### Optional Dependencies (Advanced Steps)

These are only required if you run the `codon` or `hyphy` steps.

| Tool | Purpose | Conda Package | Manual Installation / Source |
| :--- | :--- | :--- | :--- |
| **PAL2NAL** | Converts protein alignments to codon alignments. | `pal2nal` | [bork.embl.de/pal2nal](http://www.bork.embl.de/pal2nal/) (Project is Perl script, download `pal2nal.pl` and add to `$PATH`) |
| **HyPhy** | Selection pressure analysis. | `hyphy` | [hyphy.org](http://www.hyphy.org/) |

### Synteny Dependencies (Optional)

Required if running `synteny` step.

| Tool | Purpose | Conda Package |
| :--- | :--- | :--- |
| **DIAMOND** | Fast protein alignment. | `diamond` |
| **MMseqs2** | Alternative protein alignment. | `mmseqs2` |
| **pyGenomeViz** | Genome visualization (Python). | `pip install pygenomeviz` |

---

## 📂 Inputs & Outputs

### Inputs

| Argument | Description | Required | format |
| :--- | :--- | :--- | :--- |
| `inputs.faa_dir` | Directory containing proteome files (one per genome) OR a single `.faa` file. | **Yes** | FASTA (`.faa`) |
| `inputs.hmm_input` | Directory of HMM profiles OR a single `.hmm` file to search against. | **Yes** | HMMER3 (`.hmm`) |
| `inputs.cds_dir` | (Optional) Directory of nucleotide coding sequences (CDS). Required only if running `codon` or `hyphy` steps. | No | FASTA (`.fna` / `.ffn`) |
| `inputs.gtdb_dir` | (Optional) Directory of GTDB-Tk output (e.g., `gtdbtk.bac120.summary.tsv`). | No | Directory |
| `inputs.taxonomy_file` | (Optional) Custom TSV mapping `genome` -> `lineage`. | No | TSV |
| `synteny.gbk_dir` | (Optional) Directory of GenBank files for neighborhood extraction. | No | `.gbk` / `.gbff` |
| `synteny.gff_dir` | (Optional) Directory of GFF3 files (requires matching fasta in `inputs.faa_dir` or similar). | No | `.gff` |
| `post.clades_tsv` | (Optional) TSV mapping tip names to clades for dispersion metrics and KL divergence. | No | TSV: `clade_name` `tip` |

### Outputs

The pipeline creates a structured `results` directory:

| Path | Description | Format |
| :--- | :--- | :--- |
| `summary/best_hits.competitive.tsv` | **Key Result**. Table of the best HMM hit for each protein (resolved by bitscore). | TSV |
| `summary/best_hits.with_taxonomy.tsv` | **Key Result + Tax**. Same as above, but with a `taxonomy` column merged from GTDB. | TSV |
| `summary/genome_taxonomy.tsv` | Helper table mapping `genome` -> `classification`. | TSV |
| `summary/detected_clades.tsv` | Auto-detected clades (taxonomy rank or TreeCluster), columns `clade_name`, `tip`. | TSV |
| `summary/resolved_config.json` | The exact configuration used for the run (provenance). | JSON |
| `trees_iqtree/<HMM>.treefile` | The final Maximum Likelihood phylogenetic tree. | Newick |
| `fasta_per_hmm/<HMM>.faa` | Unaligned protein sequences extracted for that HMM. | FASTA |
| `alignments_clipkit/<HMM>.clipkit.faa` | The final trimmed alignment used for tree building. | FASTA |
| `embeddings/<HMM>.pca.tsv` | PCA coordinates of protein embeddings (if enabled). | TSV |
| `embeddings/<HMM>.umap.tsv` | UMAP coordinates of protein embeddings. | TSV |
| `embeddings/<HMM>.umap.png` | UMAP scatter plot (colored by clades if provided). | PNG |
| `embeddings/<HMM>.umap.clustered.png` | UMAP scatter plot colored by HDBSCAN cluster. | PNG |
| `embeddings/<HMM>.dispersion.tsv` | Quantified "functional tightness" of clades in embedding space. | TSV |
| `summary/clade_assignment.tsv` | HDBSCAN cluster assignments with protein, genome, cluster ID, and taxonomy. | TSV |
| `synteny/<HMM>/synteny.<HMM>.pdf` | Synteny plot of gene neighborhoods. | PDF |
| `synteny/<HMM>/neighborhood_proteins.faa` | Sequences of all genes in the extracted neighborhoods. | FASTA |
| `codon_alignments/<HMM>.codon.fasta` | Codon-aware alignment (if enabled). | FASTA |
| `summary/hyphy/<HMM>.<test>.json` | Legacy single-run HyPhy output when clade-aware mode is disabled or no clades are detected. | JSON |
| `summary/hyphy/<HMM>/<TEST>/<CLADE>.json` | Clade-aware HyPhy output (auto-generated from `summary/detected_clades.tsv`). | JSON |
| `summary/hyphy/<HMM>/<TEST>/<CLADE>.log` | Captured HyPhy run log for each clade-aware run. | LOG |
| `summary/hyphy/clade_runs.tsv` | QC manifest of all clade-aware HyPhy attempts and output paths. | TSV |
| `trees_labeled/<HMM>/<CLADE>.<TEST>.nwk` | Auto-generated labeled trees consumed by HyPhy (e.g., `{FG}`, `{test}`, `{reference}`). | Newick |

---

## 🏃 Usage

### 1. Basic Execution (Conda)

```bash
phylofoundry \
  --faa_dir ./data/proteomes \
  --hmm_dir ./data/markers \
  --outdir ./results_run1
```

Or using a config file:

```bash
phylofoundry --config config.json
```

### 2. Running with Docker

Mount your data directories so the container can see them.

```bash
docker run --rm -v $(pwd)/data:/data -v $(pwd)/results:/results \
  phylofoundry:latest \
  --faa_dir /data/proteomes \
  --hmm_dir /data/markers \
  --outdir /results
```

### 3. Running with Apptainer (HPC)

On HPC systems (Slurm/PBS), you typically cannot run Docker directly. Use Apptainer (formerly Singularity).

#### Step A: Build the Image
```bash
# Build SIF image from the Docker directory
apptainer build phylofoundry.sif docker://phylofoundry:latest
# OR build from a definition file if you have one
# apptainer build phylofoundry.sif PhyloFoundry.def
```

#### Step B: Run the Pipeline
You must `bind` (mount) your data directories so the container can access them.

```bash
# Example Slurm script snippet
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G

# Auto-detect CPUs (PhyloFoundry handles this internally if passed, but good to know)
export SLURM_CPUS_PER_TASK

apptainer run \
  --bind /path/to/my/data:/data \
  --bind /path/to/my/results:/results \
  phylofoundry.sif \
  --faa_dir /data/proteomes \
  --hmm_dir /data/markers \
  --outdir /results
```

### 4. CLI Reference

| Flag | Description |
| :--- | :--- |
| `--config <path>` | JSON config file (merged with defaults and CLI overrides). |
| `--faa_dir <path>` | Override `inputs.faa_dir`. |
| `--hmm_dir <path>` | Override `inputs.hmm_input`. |
| `--diamond_query <path>` | Override `inputs.diamond_query` (FASTA file or directory for DIAMOND mode). |
| `--diamond_mode` | Enable DIAMOND search mode (use protein FASTA queries instead of HMMs). |
| `--outdir <path>` | Override `output.outdir`. |
| `--cpu <N>` | Override `resources.cpu`. |
| `--start_at <step>` | Override `workflow.start_at`. |
| `--stop_after <step>` | Override `workflow.stop_after`. |
| `--force` | Override `workflow.force=true` (re-run existing steps). |
| `--combined` | Enable combined tree from all HMMs (`phylo.combined_tree`). |
| `--motifs <list>` | Comma-separated motif list for attention scoring (e.g., `HPEVY,HPEVF`). |
| `--dump_default_config` | Print the default config JSON and exit. |
| `--list-steps` | List all known workflow steps and exit. |
| `--plan` | Show the execution plan for the given config and exit (no steps run). |
| `--validate-config` | Validate the config without running the pipeline and exit. |
| `--doctor` | Check tool availability and environment health, then exit. |

*Note*: Paths inside the container (`/data`) must match where you mounted them, or just map them 1:1 (e.g., `--bind /scratch/user/project:/scratch/user/project`).

---

## 🔄 Workflow Logic

The pipeline runs as a series of sequential **Steps**. You can control execution using `--start_at <STEP>` and `--stop_after <STEP>`.

### Step 1: `prep`
-   **Input**: Directory of `.faa` files (genomes) and `.hmm` files.
-   **Action**:
    -   Concatenates all proteomes into `combined_proteomes.faa`.
    -   Concatenates all HMMs into `combined.hmm` and runs `hmmpress`.
-   **Output**: `combined_proteomes.faa`, `combined.hmm` indices.

### Step 2: `hmmer`
-   **Action**: 
    -   Runs `hmmscan` (Proteins query vs HMM db) for each genome.
    -   Runs `hmmsearch` (HMM query vs Protein db) for each HMM.
    -    parses outputs, filters by bitscore/coverage, and resolves "Best Hits".
-   **Resolution**: If a protein hits multiple HMMs, the specific HMM with the highest bitscore wins.
-   **Output**: `summary/best_hits.competitive.tsv`.

### Step 3: `extract`
-   **Action**: extracting sequences for the "best hits" identified in the previous step.
-   **Output**: `fasta_per_hmm/<hmm_name>.faa` (unaligned sequences).

### Step 4: `embed` (Optional)
-   **Action**: Uses Protein Language Models (ESM-2, etc.) to embed sequences.
-   **Analysis**: Performs PCA and UMAP on the embeddings to visualize sequence space.
-   **Clustering**: Runs HDBSCAN on raw embeddings to auto-discover functional clusters.
-   **Output**: `embeddings/<hmm_name>.pca.tsv`, `.umap.tsv`, `.umap.png`, `.umap.clustered.png`, `summary/clade_assignment.tsv`.

### Step 5: `phylo`
-   **Action**:
    1.  **Align**: Runs `mafft` (or `hmmalign`) on per-HMM FASTAs.
    2.  **Trim**: Runs `clipkit` to remove poor alignment sites.
    3.  **Tree**: Runs `iqtree` (ModelFinder + Tree search + Bootstrap).
-   **Output**: `trees_iqtree/<hmm_name>.treefile`.

### Step 6: `curate` (Optional)
-   **Action**: Prunes outlier branches using TreeShrink and/or ESM-based sequence filtering. Writes curated artifacts to a `curated/` overlay directory without overwriting raw pipeline outputs.
-   **Output**: `curated/trees/`, `curated/fasta_per_hmm/`, `curated/alignments_clipkit/`.

### Step 7: `taxonomy_integrate` (Optional)
-   **Action**: Loads taxonomy from a GTDB-Tk summary directory (`inputs.gtdb_dir`) or a custom genome→lineage TSV (`inputs.taxonomy_file`) and annotates best-hit results.
-   **Output**: `summary/genome_taxonomy.tsv`, `summary/best_hits.with_taxonomy.tsv`.

### Step 8: `conservation_metrics` (Optional)
-   **Action**: Calculates per-site conservation scores and KL divergence using scikit-bio.
-   **Output**: `summary/post_scikitbio/`.

### Step 9: `detect_clades` (Optional)
-   **Action**: Detects clades via taxonomy rank, TreeCluster, or tree+embedding strategies. Writes the clade table consumed by `hyphy` and `discover_motifs`.
-   **Output**: `summary/detected_clades.tsv`, `clade_assignments/`.

### Step 10: `post` (Optional, Legacy)
-   **Action**: Backward-compatibility shim combining conservation metrics, KL divergence, and clade detection in a single step. New workflows should prefer the dedicated `taxonomy_integrate`, `conservation_metrics`, and `detect_clades` steps.
-   **Output**: `summary/post_scikitbio/`, `clade_assignments/`.

### Step 11: `synteny` (Optional)
-   **Action**: Extracts gene neighborhoods (configurable window), computes similarity (DIAMOND/MMseqs2), and plots synteny tracks ordered by phylogeny.
-   **Output**: `synteny/<hmm_name>/synteny.<hmm_name>.pdf`.

### Step 12: `codon` (Optional)
-   **Action**: 
    -   Matches protein sequences to their CDS (nucleotide) sequences.
    -   **Strips terminal stop codons** (`*` from AA, TAA/TAG/TGA from CDS) before running `pal2nal.pl`.
    -   Uses `pal2nal.pl` with `-nogap -nomismatch` flags for robust codon alignment.
-   **Output**: `codon_alignments/<hmm_name>.codon.fasta`.

### Step 13: `hyphy` (Optional)
-   **Action**: Runs selection tests (e.g., RELAX, aBSREL, MEME) on codon alignments and trees. If `summary/detected_clades.tsv` exists and `hyphy.use_detected_clades=true`, HyPhy automatically builds labeled trees and runs per-clade analyses (no manual RELAX labeling required).
-   **Output**: Legacy `summary/hyphy/<hmm_name>.<test>.json` plus clade-aware outputs under `summary/hyphy/<hmm_name>/<TEST>/<clade_name>.json` and labeled trees under `trees_labeled/<hmm_name>/`.

### Step 14: `score_motifs` (Optional)
-   **Action**: Passes sequences through ESM-2 with `output_attentions=True`, extracts attention weights at user-specified motif positions.
-   **CLI**: `--motifs HPEVY,HPEVF`
-   **Output**: `summary/motif_attention_scores.tsv` — columns: `seq_id`, `motif`, `start_pos`, `end_pos`, `attention_score`, `clade_id`, `type`.
-   **HA Output (optional)**: `attention/<HMM>.ha_sites.tsv` + `summary/ha_summary.tsv` when `ha.enabled=true` and `motifs.use_ha=true`.

### Step 15: `discover_motifs` (Optional)
-   **Action**: Iterates over all HDBSCAN clades, comparing the 1D attention profiles of each clade against the combined average of all others. Finds peaks in the attention delta and extracts k-mers as candidate novel structural hubs for that specific clade.
-   **CLI**: N/A, runs automatically if `discover.enabled` is `true`.
-   **Output**: `summary/discovered_motifs.tsv` — columns: `kmer`, `n_sequences`, `mean_attention_delta`, `source_clade`, `reference_clade`.
-   **HA Outputs (optional)**: `discover/<HMM>.ha_enrichment.tsv` and `discover/<HMM>.ha_hubs.tsv` when `ha.enabled=true` and `discover.use_ha=true`.
-   **Candidate residue outputs (optional)**: `discover/<HMM>.candidate_residues.tsv` and `discover/<HMM>.candidate_regions.tsv` when `discover.candidates.enabled=true`.

---

## ⏯️ Resuming & Checkpoints

PhyloFoundry is designed to be highly resilient. It uses a file-existence check to determine if a step needs to be run.

### Scenarios

#### Scenario A: Pipeline Interrupted
If the pipeline crashes or is cancelled (e.g., walltime limit reached on HPC):
1.  Simply run the **exact same command** again.
2.  It will detect existing output files (e.g., `hmmer` tables, `extract` FASTAs) and skip those steps.
3.  It will pick up exactly where it left off (e.g., processing the remaining HMMs in `phylo`).

#### Scenario B: Adding Analysis (e.g., Embeddings)
You ran the pipeline without embeddings, but now want to add them:
1.  Enable embeddings in your config (`"embeddings": { "enabled": true }`).
2.  Run with `--start_at embed`.
    ```bash
    phylofoundry --config config.json --start_at embed
    ```
3.  This skips `prep`, `hmmer`, and `extract`, loading the necessary data to run `embed`.

#### Scenario C: Force Re-run
To overwrite existing results (e.g., if you changed parameters like `mafft_mode`):
```bash
phylofoundry --config config.json --force
```
*Note*: This forces **all** steps in the workflow range. To force only one step, use:
```bash
phylofoundry --config config.json --start_at phylo --stop_after phylo --force
```

---

## ⚙️ Configuration

Generate a template config:
```bash
phylofoundry --dump_default_config > config.json
```

### Key Options

#### `phylo`
-   `mafft_mode`: Alignment strategy.
    -   `auto` (default): Fast/Standard.
    -   `linsi`: High accuracy (L-INS-i).
    -   `ginsi`: High accuracy (G-INS-i).
-   `iq_boot`: Number of bootstrap replicates (default 1000).

#### `embeddings`
-   `enabled`: Set to `true` to run embedding step.
-   `model`: `esm2_t33_650M_UR50D` (default) or other HuggingFace models.
-   `device`: `cuda` (GPU) or `cpu`.

#### `job` / `resources`
-   `cpu`: Number of threads. If running on Slurm, this is auto-detected from `$SLURM_CPUS_PER_TASK`.

---

## 💎 DIAMOND BLAST Search Mode

PhyloFoundry supports an alternative search mode using [DIAMOND](https://github.com/bbuchfink/diamond) blastp. Instead of HMM profiles, you supply protein FASTA query files and DIAMOND searches them against the combined genome proteome database.

### When to Use DIAMOND Mode

-   You have representative protein sequences but **not HMM profiles**.
-   You want fast BLAST-based homology searches with identity/coverage thresholds.
-   You are doing exploratory analysis and haven't built curated HMMs yet.

> **Note**: In DIAMOND mode, the `phylo` step automatically uses **MAFFT** for multiple sequence alignment (since there are no HMM profiles available for `hmmalign`).

### CLI Usage

```bash
phylofoundry \
  --faa_dir /path/to/genomes/ \
  --diamond_query /path/to/queries/ \
  --diamond_mode \
  --outdir /path/to/output/ \
  --cpu 16
```

The `--diamond_query` argument accepts:
-   A **single FASTA file** (`.faa`, `.fasta`, `.fa`, `.fas`) — one query group.
-   A **directory of FASTA files** — one query group per file (basename without extension becomes the group name).

### Config-Based Usage

Add this to your `config.json`:

```json
{
    "inputs": {
        "faa_dir": "/path/to/genomes/",
        "diamond_query": "/path/to/queries/",
        "hmm_input": null
    },
    "diamond": {
        "enabled": true,
        "sensitivity": "sensitive",
        "max_evalue": 1e-5,
        "max_target_seqs": 500,
        "min_identity": 30.0,
        "min_coverage": 0.5,
        "block_size": 2.0,
        "index_chunks": 4
    }
}
```

Then run:

```bash
phylofoundry --config config.json --outdir /path/to/output/
```

### DIAMOND Configuration Options

| Key | Default | Description |
| :--- | :--- | :--- |
| `diamond.enabled` | `false` | Enable DIAMOND mode (replaces HMMER). |
| `diamond.sensitivity` | `"sensitive"` | Sensitivity flag: `fast`, `mid-sensitive`, `sensitive`, `more-sensitive`, `very-sensitive`, `ultra-sensitive`. |
| `diamond.max_evalue` | `1e-5` | Maximum e-value threshold for hits. |
| `diamond.max_target_seqs` | `500` | Maximum number of target sequences per query. |
| `diamond.min_identity` | `30.0` | Minimum percent identity (post-search filter). |
| `diamond.min_coverage` | `0.5` | Minimum query coverage (post-search filter). |
| `diamond.block_size` | `2.0` | DIAMOND `-b` parameter (memory scaling). |
| `diamond.index_chunks` | `4` | DIAMOND `-c` parameter. |

### Pipeline Differences in DIAMOND Mode

| Aspect | HMM Mode | DIAMOND Mode |
| :--- | :--- | :--- |
| Input | `.hmm` profiles | Protein FASTA(s) |
| `prep` step | Builds `combined.hmm` + `combined_proteomes.faa` | Builds `combined_proteomes.faa` only |
| `hmmer` step | Runs `hmmscan` + `hmmsearch` | Runs `diamond blastp` |
| `phylo` alignment | `hmmalign` (default) or MAFFT | MAFFT (auto-enabled) |
| `extract` and downstream | Unchanged | Unchanged |

---


```text
results/
├── summary/
│   ├── best_hits.competitive.tsv  # Main results table
│   ├── best_hits.with_taxonomy.tsv # With GTDB taxonomy
│   ├── clade_assignment.tsv       # HDBSCAN cluster assignments
│   ├── resolved_config.json       # Provenance
│   └── post_scikitbio/            # Conservation metrics
├── fasta_per_hmm/                 # Extracted sequences
├── alignments_hmm/                # Raw alignments
├── alignments_clipkit/            # Trimmed alignments
├── trees_iqtree/                  # Final Newick trees
├── embeddings/                    # PCA/UMAP data + plots
│   ├── <HMM>.umap.png             # UMAP scatter
│   └── <HMM>.umap.clustered.png   # UMAP colored by HDBSCAN
├── synteny/                       # Gene neighborhood plots
└── codon_alignments/              # PAL2NAL codon alignments
```

---

## 📘 Detailed Configuration Guide

This section explains every key option in `config.json`.

### `inputs`
Defines your raw data.
-   `faa_dir`: (Required) Path to directory containing protein FASTA files (`.faa`), or a single merged `.faa` file.
-   `hmm_input`: (Required) Path to directory containing HMM profiles (`.hmm`), or a single `.hmm` file.
-   `cds_dir`: (Optional) Directory of nucleotide coding sequences (`.fna`). Only needed for `codon` / `hyphy` steps.
-   `gtdb_dir`: (Optional) Directory containing GTDB-Tk summary files (e.g. `gtdbtk.bac120.summary.tsv`). Used to add taxonomy to summary tables.
-   `taxonomy_file`: (Optional) Custom TSV (columns: `genome`, `lineage`) if not using GTDB.

### `output`
-   `outdir`: (Required) Where all results go.

### `workflow`
Controls execution flow.
-   `start_at`: Start pipeline at a specific step (e.g., `"phylo"`).
-   `stop_after`: Stop after a specific step (e.g., `"hmmer"`).
-   `force`: (Default: `false`) Overwrite existing output files for the active steps.
-   `hmm_manifest`: (Default: `null`) Path to a text file listing specific HMM names to process (one per line).

### `prep`
Input preparation.
-   `cleanup_combined_faa`: (Default: `false`) Delete the merged `combined_proteomes.faa` after `hmmer` completes to save disk space.

### `filtering`
-   `global_min_score`: (Default: `25.0`) Minimum bitscore for a hit.
-   `min_coverage`: (Default: `0.5`) Minimum query coverage (0.0-1.0).
-   `keep_tbl`: (Default: `false`) Keep raw `hmmscan`/`hmmsearch` output tables (can be large).

### `phylo`
Phylogenetic inference.
-   `mafft_mode`: (Default: `"auto"`) alignment strategy. `auto`, `linsi` (accurate), `ginsi`, `fftnsi` (fast).
-   `iqtree_bin`: (Default: `"iqtree"`) Name/path of IQ-TREE executable.
-   `iq_boot`: (Default: `1000`) Bootstrap replicates.
-   `no_asr`: (Default: `false`) Skip Ancestral Sequence Reconstruction (saves memory).
-   `skip_clipkit`: (Default: `false`) Skip alignment trimming.
-   `mafft`: (Default: `false`) Force MAFFT alignment even if HMM alignment is available.

### `synteny`
Gene neighborhood analysis.
-   `enabled`: Set to `true` to run.
-   `gbk_dir`: (Required if enabled) Directory of GenBank files (`.gbk`) with genomic context.
-   `gff_dir`: (Optional) Directory of GFF3 files. If both GBK and GFF are provided, GBK is tried first, then GFF as fallback.
-   `window_genes`: (Default: `10`) Number of genes to extract upstream and downstream of the hit. May be fewer at contig boundaries.
-   `similarity`: homology search settings.
    -   `method`: `"diamond"` (default) or `"mmseqs"`.
    -   `min_identity`: (Default: `30`) % identity cutoff.
-   `include_tree`: (Default: `true`) Plot tree alongside synteny tracks.
-   `output_format`: (Default: `"pdf"`) `pdf` or `png` or `html`.

### `embeddings`
Protein Language Model analysis.
-   `enabled`: Set to `true` to run.
-   `model`: (Default: `"esm2_t33_650M_UR50D"`) ESM2 model.
-   `device`: `"cuda"` (GPU) or `"cpu"`.
-   `write_full_vectors`: (Default: `false`) Must be `true` for HA attention-based methods.
    - Attention-driven analyses include HA/LoC calling, motif HA scoring, HA enrichment/hubs, and candidate residues.
-   `cluster_embeddings`: (Default: `true`) Run HDBSCAN clustering on embeddings.
-   `hdbscan_min_cluster_size`: (Default: `5`) Minimum cluster size for HDBSCAN.

### `ha`
High-Attention (HA) site calling.
-   `enabled`: (Default: `false`) Enable HA site extraction.
-   `method`: (Default: `"middle"`) `"middle"` (legacy layer-range aggregation) or `"loc"` (paper-style convergence layer).
-   `layer_mode`: `"middle"` (default) or `"range"`.
-   `layer_start`, `layer_end`: Explicit layer interval when `layer_mode="range"`.
-   `agg`: (Default: `"median"`) Layer aggregation statistic (`"median"` or `"mean"`).
-   `call_mode`: `"percentile"` (default) or `"topk"`.
-   `percentile`: (Default: `0.05`) Top fraction when percentile mode is used.
-   `topk`: (Default: `20`) Number of HA sites in top-k mode.
-   `min_sites`, `max_sites`: (Defaults: `8`, `60`) Clamp called HA-site counts.
-   `loc_norm_mode`: (Default: `"max"`) Per-layer normalization mode used in `method="loc"`.
-   `loc_theta_target_deg`: (Default: `90`) Target angle for convergence-layer selection.
-   `loc_break_adjust`: (Default: `-1`) Optional adjustment applied to the PWLF breakpoint-derived site count.

### `taxonomy_integrate`
GTDB / custom taxonomy annotation.
-   `enabled`: (Default: `false`) Set to `true` to run. Reads `inputs.gtdb_dir` or `inputs.taxonomy_file` and writes `summary/genome_taxonomy.tsv` and `summary/best_hits.with_taxonomy.tsv`.

### `conservation_metrics`
Per-site conservation and KL divergence.
-   `enabled`: (Default: `false`) Set to `true` to run.
-   `compute_conservation`: (Default: `false`) Calculate per-column conservation scores.
-   `conservation_metric`: (Default: `"inverse_shannon_uncertainty"`) scikit-bio metric name.
-   `compute_kl`: (Default: `false`) Compute KL divergence between clade pairs.
-   `kl_pairs`: (Default: `null`) Explicit clade pairs (e.g. `"A:B,A:background"`). If `null` and `compute_kl=true`, computes each detected clade vs. all other tips.
-   `clades_tsv`: (Optional) TSV mapping tips to groups; overrides auto-detected clades for dispersion / KL calculations.

### `detect_clades`
Automated clade detection.
-   `enabled`: (Default: `false`) Set to `true` to run.
-   `clades_tsv`: (Optional) Pre-existing clade TSV (columns: `clade_name`, `tip`) to load instead of auto-detecting.
-   `detect_method`: (Default: `null`) `"taxonomy"`, `"treecluster"`, or `"tree_embed"` to auto-generate `summary/detected_clades.tsv`.
-   `taxonomy_clade_level`: (Default: `"genus"`) Taxonomic rank used when `detect_method="taxonomy"`.
-   `treecluster_threshold`: (Default: `0.045`) Distance threshold passed to TreeCluster.
-   `treecluster_method`: (Default: `"max_clade"`) TreeCluster clustering method.
-   `embedtree_support_min`: (Default: `80`) Minimum internal-node support for candidate splits when `detect_method="tree_embed"`.
-   `embedtree_min_size`: (Default: `5`) Minimum tips in a candidate clade.
-   `embedtree_max_size`: (Default: `5000`) Maximum tips in a candidate clade (`null` to disable).
-   `embedtree_top_k`: (Default: `10`) Maximum non-overlapping embedding-shift clades emitted per HMM.
-   `embedtree_pcs`: (Default: `10`) Number of embedding PCs used for split scoring.
-   `embedtree_distance`: (Default: `"euclidean"`) Distance metric for centroid and dispersion calculations (`"euclidean"` or `"cosine"`).
-   `embedtree_allow_nested`: (Default: `false`) If `true`, descendants of selected splits may also be emitted.
-   `embedtree_require_monophyly`: (Default: `true`) Enforces clades to be internal tree nodes (monophyletic by construction).
-   `embedtree_emit_all`: (Default: `false`) If `true`, emit all accepted nodes instead of truncating at `embedtree_top_k`.

### `post`
> ⚠️ **Legacy / backward-compatibility shim.** New workflows should use the dedicated `taxonomy_integrate`, `conservation_metrics`, and `detect_clades` steps instead.

Post-processing metrics (combined legacy step).
-   `enabled`: Set to `true` to run.
-   `compute_conservation`: (Default: `false`) Calculate conservation scores.
-   `clades_tsv`: (Optional) TSV mapping tips to groups for dispersion analysis.
-   `detect_clades_method`: (Optional) `taxonomy`, `treecluster`, or `tree_embed` to auto-generate `summary/detected_clades.tsv`.
-   `taxonomy_clade_level`: (Default: `"genus"`) Taxonomic rank used when `detect_clades_method=taxonomy`.
-   `treecluster_threshold`: (Default: `0.045`) Distance threshold passed to TreeCluster.
-   `treecluster_method`: (Default: `"max_clade"`) TreeCluster clustering method.
-   `embedtree_support_min`: (Default: `80`) Minimum internal-node support for candidate splits when `detect_clades_method=tree_embed`.
-   `embedtree_min_size`: (Default: `5`) Minimum tips in a candidate clade.
-   `embedtree_max_size`: (Default: `5000`) Maximum tips in a candidate clade (`null` to disable).
-   `embedtree_top_k`: (Default: `10`) Maximum non-overlapping embedding-shift clades emitted per HMM.
-   `embedtree_pcs`: (Default: `10`) Number of embedding PCs used for split scoring.
-   `embedtree_distance`: (Default: `"euclidean"`) Distance metric for centroid and dispersion calculations (`"euclidean"` or `"cosine"`).
-   `embedtree_allow_nested`: (Default: `false`) If `true`, descendants of selected splits may also be emitted.
-   `embedtree_require_monophyly`: (Default: `true`) Enforces clades to be internal tree nodes (monophyletic by construction).
-   `embedtree_emit_all`: (Default: `false`) If `true`, emit all accepted nodes instead of truncating at `embedtree_top_k`.
-   `summary/node_scores.embedtree.tsv`: Per-node QC table containing support, dispersion, separation, effect size, and tree diameter for tree-embedding scoring.
-   `compute_kl`: If enabled and no explicit `kl_pairs`, computes `clade vs all other tips` for each detected clade.


### HA Sites (High-Attention residues, paper-inspired)
- HA sites summarize ESM2 attention into per-residue "attention received" scores.
- Two modes are supported:
  - `ha.method="middle"`: aggregate normalized received-attention vectors across middle/configured layers, then call HA by percentile/top-k.
  - `ha.method="loc"`: fit a piecewise linear model per layer, select a convergence layer (LoC), and call HA above the layer-specific breakpoint.
- HA outputs:
  - `attention/<HMM>.ha_sites.tsv` (`pos_ungapped` is 1-based; includes `loc_layer` when LoC mode is active)
  - `summary/ha_summary.tsv` (includes `method`, and `loc_layer`/`norm_mode` for LoC)

> ⚠️ Attention-driven analyses (HA / convergence layer / motif scoring / HA discovery / candidate residues) require `embeddings.write_full_vectors: true`.
> If HA is enabled while `write_full_vectors` is false, the pipeline exits with a clear error.

Example HA config:
```json
{
  "embeddings": { "write_full_vectors": true },
  "ha": {
    "enabled": true,
    "method": "loc",
    "loc_norm_mode": "max"
  },
  "motifs": { "use_ha": true },
  "discover": { "enabled": true, "use_ha": true, "candidates": { "enabled": true } }
}
```

### `codon`
Codon alignments.
-   `enabled`: Set to `true` to run.
-   `pal2nal_cmd`: (Default: `"pal2nal.pl"`) Path to PAL2NAL script.

### `motifs`
Targeted motif scoring.
-   `use_ha`: (Default: `false`) Add HA-overlap and HA-score metrics to `summary/motif_ha_scores.tsv`.

### `discover`
Unsupervised motif discovery.
-   `standard_clade`: (Default: `null`) Manually specify a reference clade ID/name for 1-vs-1 comparisons.
-   `novel_clade`: (Default: `null`) Manually specify the focal clade ID/name for 1-vs-1 comparisons.
-   `use_ha`: (Default: `false`) Enable alignment-aware HA enrichment/hub calling outputs.
-   `ha_window`: (Default: `9`) Smoothing window for HA delta profiles.
-   `ha_delta_min`: (Default: `0.2`) Minimum smoothed delta for HA hub calls.
-   `ha_gap_frac_max`: (Default: `0.6`) Maximum mean gap fraction allowed in HA hubs.
-   `candidates.enabled`: (Default: `true`) Emit clade-aware functional candidate residues/regions.
-   `candidates.min_delta_ha`: (Default: `0.1`) Minimum clade-vs-rest HA enrichment.
-   `candidates.min_cons_frac`: (Default: `0.6`) Minimum consensus frequency gate.
-   `candidates.max_gap_frac`: (Default: `0.6`) Maximum gap fraction at candidate columns.
-   `candidates.w_delta_ha`, `w_aa_shift`, `w_js`: Score weights for HA enrichment, consensus shift, and JS divergence.

### `hyphy`
Selection tests.
-   `enabled`: Set to `true` to run.
-   `hyphy_tests`: (Default: `"RELAX,aBSREL,MEME"`) List of tests to run.
-   `use_detected_clades`: (Default: `true`) Auto-load `summary/detected_clades.tsv` for clade-aware HyPhy runs.
-   `min_clade_size`: (Default: `4`) Skip per-HMM clades smaller than this number of tips.
-   `label_mode`: (Default: `"crown"`) Branch-label strategy for clade foreground (`"crown"` or `"stem"`).
-   `relax_label_reference`: (Default: `true`) Add `{reference}` labels to non-foreground branches for RELAX.
-   `hyphy_args`: Existing per-test args are still supported; in clade-aware mode, `aBSREL`/`BUSTED` branch labels are forced to `FG` and RELAX labels are synchronized to `test`/`reference`.

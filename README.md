# PhyloFoundry

**PhyloFoundry** is a robust, HPC-ready bioinformatics pipeline for competitive HMM analysis, phylogenetics, and optional protein language model embeddings. It automates the journey from raw proteomes and HMM profiles to publication-ready phylogenetic trees and functional landscape metrics.

---

## 🚀 Features

-   **Competitive HMM Hits**: Uses both `hmmscan` and `hmmsearch` to identify the best functional assignments for proteins, resolving overlapping hits competitively by bitscore.
-   **Automated Phylogeny**: Per-HMM alignment (MAFFT/HMMER), trimming (ClipKit), and tree inference (IQ-TREE).
-   **Protein Embeddings** (Optional): Generates per-HMM embeddings (ESM) and dimensionality reduction (PCA/UMAP), with HDBSCAN clustering and UMAP scatter plots.
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
| `summary/hyphy/<HMM>.<test>.json` | Selection test results (e.g., RELAX, MEME). | JSON |

---

## 🏃 Usage

### 1. Basic Execution (Conda)

```bash
phylofoundry \
  --inputs.faa_dir ./data/proteomes \
  --inputs.hmm_input ./data/markers \
  --output.outdir ./results_run1
```

### 2. Running with Docker

Mount your data directories so the container can see them.

```bash
docker run --rm -v $(pwd)/data:/data -v $(pwd)/results:/results \
  phylofoundry:latest \
  --inputs.faa_dir /data/proteomes \
  --inputs.hmm_input /data/markers \
  --output.outdir /results
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
  --inputs.faa_dir /data/proteomes \
  --inputs.hmm_input /data/markers \
  --output.outdir /results
```

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

### Step 6: `post` (Optional)
-   **Action**: Calculates conservation scores and Scikit-bio metrics.
-   **Output**: `summary/post_scikitbio/`.

### Step 7: `codon` (Optional)
-   **Action**: 
    -   Matches protein sequences to their CDS (nucleotide) sequences.
    -   Uses `pal2nal.pl` to generate codon alignments from the protein alignments.
-   **Output**: `codon_alignments/<hmm_name>.codon.fasta`.

### Step 8: `hyphy` (Optional)
-   **Action**: Runs selection tests (e.g., RELAX, aBSREL, MEME) on the codon alignments and trees.
-   **Output**: `summary/hyphy/<hmm_name>.<test>.json`.

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

## 📂 Output Structure

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
-   `cluster_embeddings`: (Default: `true`) Run HDBSCAN clustering on embeddings.
-   `hdbscan_min_cluster_size`: (Default: `5`) Minimum cluster size for HDBSCAN.

### `post`
Post-processing metrics.
-   `enabled`: Set to `true` to run.
-   `compute_conservation`: (Default: `false`) Calculate conservation scores.
-   `clades_tsv`: (Optional) TSV mapping tips to groups for dispersion analysis.

### `codon`
Codon alignments.
-   `enabled`: Set to `true` to run.
-   `pal2nal_cmd`: (Default: `"pal2nal.pl"`) Path to PAL2NAL script.

### `hyphy`
Selection tests.
-   `enabled`: Set to `true` to run.
-   `hyphy_tests`: (Default: `"RELAX,aBSREL,MEME"`) List of tests to run.


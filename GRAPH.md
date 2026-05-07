# PhyloFoundry — Codebase Graph

Generated from `/graphify .` on the `src/` tree.

---

## Pipeline Step Dependency Graph

Each node is a workflow step. Solid arrows (`-->`) represent hard dependencies;
dashed arrows (`-.->`) represent optional steps that depend on an upstream step.
Steps marked `[opt]` are skipped unless their `enabled` config flag is `true`.

```mermaid
flowchart TD
    prep["⚙️ prep\nCombine FAA + HMM inputs"]
    hmmer["🔍 hmmer\nhmmscan / hmmsearch / DIAMOND"]
    extract["✂️ extract\nExtract hit sequences"]
    embed["🧠 embed [opt]\nESM-2 / HuggingFace embeddings"]
    maape["🔀 maape [opt]\nMAAPE embedding network graphs"]
    phylo["🌳 phylo\nMSA → ClipKit → IQ-TREE"]
    curate["🔧 curate [opt]\nTreeShrink pruning / ESM filter"]
    taxonomy_integrate["🗂️ taxonomy_integrate [opt]\nGTDB / custom taxonomy"]
    conservation_metrics["📊 conservation_metrics [opt]\nscikit-bio conservation + KL div"]
    detect_clades["🔎 detect_clades [opt]\nTaxonomy / TreeCluster / embedding"]
    aa_composition["🧪 aa_composition [opt]\nAA composition + biochemical metrics"]
    post["🔄 post [opt]\nLegacy shim — prefers dedicated steps"]
    tree_viz["🎨 tree_viz\nggtree annotated tree plots"]
    synteny["🗺️ synteny [opt]\nclinker / pygenomeviz"]
    codon["🧬 codon [opt]\npal2nal codon alignments"]
    hyphy["⚡ hyphy [opt]\nRELAX / aBSREL / MEME"]
    score_motifs["🎯 score_motifs [opt]\nESM-2 attention motif scoring"]
    discover_motifs["🔬 discover_motifs [opt]\nESM-2 attention motif discovery"]

    prep --> hmmer
    hmmer --> extract
    extract --> embed
    embed -.-> maape
    extract --> phylo
    extract -.-> score_motifs
    extract -.-> discover_motifs
    hmmer -.-> synteny
    phylo --> curate
    phylo --> taxonomy_integrate
    phylo --> conservation_metrics
    phylo --> detect_clades
    phylo --> aa_composition
    phylo -.-> post
    phylo --> tree_viz
    taxonomy_integrate -.-> tree_viz
    detect_clades -.-> tree_viz
    phylo --> codon
    codon --> hyphy
    detect_clades --> hyphy
    detect_clades -.-> discover_motifs

    style embed        fill:#d0e8ff,stroke:#5599cc
    style maape        fill:#d0e8ff,stroke:#5599cc
    style curate       fill:#d0e8ff,stroke:#5599cc
    style taxonomy_integrate fill:#d0e8ff,stroke:#5599cc
    style conservation_metrics fill:#d0e8ff,stroke:#5599cc
    style detect_clades fill:#d0e8ff,stroke:#5599cc
    style aa_composition fill:#d0e8ff,stroke:#5599cc
    style post         fill:#ffe4b5,stroke:#cc8800
    style tree_viz     fill:#d4edda,stroke:#28a745
    style synteny      fill:#d0e8ff,stroke:#5599cc
    style codon        fill:#d0e8ff,stroke:#5599cc
    style hyphy        fill:#d0e8ff,stroke:#5599cc
    style score_motifs fill:#d0e8ff,stroke:#5599cc
    style discover_motifs fill:#d0e8ff,stroke:#5599cc
```

---

## Source Module Structure

```mermaid
graph TD
    subgraph entry["Entry Points"]
        main["main.py\nCLI (argparse)"]
        review_app["review_app.py\nStreamlit curation UI"]
    end

    subgraph core["Core Infrastructure"]
        config["config.py\nConfig resolution + validation"]
        constants["constants.py\nSTEPS list + DEFAULT_CONFIG"]
        artifact_paths["artifact_paths.py\nArtifact path helpers"]
        pipeline["pipeline.py\nPipeline runner (run_pipeline)"]
        execution_planner["execution_planner.py\nBuild ExecutionPlan from config"]
        execution_schema["execution_schema.py\nExecutionPlan / PlannedStep / StepState"]
        workflow_registry["workflow_registry.py\nRegistered StepDefinitions"]
        workflow_schema["workflow_schema.py\nStepDefinition / ArtifactSpec / ToolRequirement"]
        failure_policy["failure_policy.py\nStep failure + dependency propagation"]
        checkpoint["checkpoint.py\nTwo-tier NDJSON + SQLite checkpointing"]
        logging_utils["logging_utils.py\nStructured pipeline + step logs"]
        preflight["preflight.py\nTool + config preflight checks"]
        provenance["provenance.py\nInput provenance recording"]
    end

    subgraph tasks["tasks/ (one module per step)"]
        t_prep["prep.py"]
        t_hmmer["hmmer.py"]
        t_diamond["diamond.py"]
        t_extract["extract.py"]
        t_embed["embed.py"]
        t_maape["maape.py"]
        t_phylo["phylo.py"]
        t_curate["curate.py"]
        t_post["post.py"]
        t_tree_viz["tree_viz.py"]
        t_synteny["synteny.py"]
        t_codon["codon.py"]
        t_hyphy["hyphy.py"]
        t_asr["asr.py"]
        t_tax["taxonomy_integrate.py"]
        t_cons["conservation_metrics.py"]
        t_clades["detect_clades.py"]
        t_aacomp["aa_composition.py"]
        t_score["score_motifs.py"]
        t_discover["discover_motifs.py"]
    end

    subgraph utils["utils/ (shared helpers)"]
        u_bio["bio.py\nFASTA I/O + sequence utils"]
        u_ha["ha.py\nHA-specific utilities"]
        u_helpers["helpers.py\nsafe_mkdir, write_json, …"]
        u_pfam["pfam_scan.py\nPfam scan wrappers"]
        u_test["test.py\nTest helpers"]
    end

    main --> config
    main --> pipeline
    main --> execution_planner
    main --> preflight
    main --> provenance

    pipeline --> artifact_paths
    pipeline --> constants
    pipeline --> execution_planner
    pipeline --> execution_schema
    pipeline --> failure_policy
    pipeline --> checkpoint
    pipeline --> logging_utils
    pipeline --> tasks

    execution_planner --> workflow_registry
    execution_planner --> execution_schema
    execution_planner --> constants
    workflow_registry --> workflow_schema

    tasks --> u_bio
    tasks --> u_helpers
```

---

## Artifact Flow (Data Lineage)

```mermaid
flowchart LR
    faa_dir["📁 faa_dir\n(input proteomes)"]
    hmm_input["📄 hmm_input\n(HMM profiles)"]

    combined_faa["combined_proteomes.faa"]
    combined_hmm["combined.hmm"]
    hits_tsv["hmmscan/hmmsearch\n_hits.filtered.tsv"]
    fasta_per_hmm["fasta_per_hmm/\n*.faa"]
    alignments_hmm["alignments_hmm/\n*.hmm.faa"]
    alignments_clipkit["alignments_clipkit/\n*.clipkit.faa"]
    trees["trees_iqtree/\n*.treefile"]
    embeddings["embeddings/\n*.pca.tsv\n*.umap.tsv\n*.png"]
    maape_out["maape/\n*.pkl *.txt *.png"]
    curated["curated/\ntrees/ fasta/ alignments/"]
    tree_viz_out["tree_viz/\n*.tree.png/pdf"]
    codon_aln["codon_alignments/\n*.codon.fasta"]
    hyphy_out["summary/hyphy/\n*.json"]
    synteny_out["synteny/\n*.pdf"]
    post_bio["summary/post_scikitbio/"]
    clades_tsv["summary/detected_clades.tsv"]
    tax_tsv["summary/genome_taxonomy.tsv"]
    best_hits_tax["summary/best_hits.with_taxonomy.tsv"]
    aa_comp["summary/aa_composition/"]
    motif_scores["motifs/motif_attention_scores.tsv"]
    disc_motifs["discover/discovered_motifs.tsv"]

    faa_dir --> combined_faa
    hmm_input --> combined_hmm
    combined_faa --> hits_tsv
    combined_hmm --> hits_tsv
    hits_tsv --> fasta_per_hmm
    fasta_per_hmm --> alignments_hmm
    alignments_hmm --> alignments_clipkit
    alignments_clipkit --> trees
    fasta_per_hmm --> embeddings
    embeddings --> maape_out
    trees --> curated
    trees --> tree_viz_out
    clades_tsv --> tree_viz_out
    tax_tsv --> tree_viz_out
    trees --> codon_aln
    trees --> post_bio
    trees --> clades_tsv
    trees --> tax_tsv
    tax_tsv --> best_hits_tax
    codon_aln --> hyphy_out
    clades_tsv --> hyphy_out
    hits_tsv --> synteny_out
    trees --> aa_comp
    fasta_per_hmm --> motif_scores
    fasta_per_hmm --> disc_motifs
    clades_tsv --> disc_motifs
```

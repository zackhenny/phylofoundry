#!/usr/bin/env Rscript
# plot_synteny_tree.R
#
# Produce a combined tree + synteny neighborhood figure using ggtree and
# ggtreeExtra.  The tree occupies the left panel; gene neighborhood arrows
# are drawn as external panels on the right using geom_fruit.
#
# The synteny data TSV written by phylofoundry's Python synteny step has
# the following required columns:
#   genome, focal_protein, position_relative, strand, gene_label,
#   hmm_annotation, is_focal
#
# Optional columns:
#   clade_name  (if clade_assign_dir output is joined in)
#   tax_label   (if taxonomy_integrate output is joined in)
#
# Usage:
#   Rscript plot_synteny_tree.R \
#     --treefile      /path/to/hmm.treefile \
#     --synteny_tsv   /path/to/hmm_gene_clusters.tsv \
#     --outdir        /path/to/synteny/hmm/tree_panel/ \
#     --hmm_name      MyHMM \
#     [--clades_tsv   /path/to/hmm.clades.tsv] \
#     [--taxonomy_tsv /path/to/genome_taxonomy.tsv] \
#     [--tax_level    genus] \
#     [--format       png,pdf] \
#     [--width        20] \
#     [--height       0] \
#     [--bootstrap_min 80] \
#     [--show_tip_labels auto] \
#     [--color_palette Set3]

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(ape)
  library(treeio)
  library(ggtree)
  library(ggtreeExtra)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(RColorBrewer)
})

# ---------------------------------------------------------------------------
# CLI argument parsing
# ---------------------------------------------------------------------------
option_list <- list(
  make_option("--treefile",        type = "character", default = NULL,
              help = "Path to IQ-TREE .treefile [required]"),
  make_option("--synteny_tsv",     type = "character", default = NULL,
              help = "Synteny gene-cluster TSV from phylofoundry [required]"),
  make_option("--clades_tsv",      type = "character", default = NULL,
              help = "Clade assignment TSV (clade_name, tip)"),
  make_option("--taxonomy_tsv",    type = "character", default = NULL,
              help = "Taxonomy TSV (genome, taxonomy)"),
  make_option("--tax_level",       type = "character", default = "genus",
              help = "Taxonomy rank to display [default: genus]"),
  make_option("--outdir",          type = "character", default = ".",
              help = "Output directory [default: .]"),
  make_option("--hmm_name",        type = "character", default = "tree",
              help = "Label for output filenames [default: tree]"),
  make_option("--format",          type = "character", default = "png",
              help = "Comma-separated formats: png,pdf,svg [default: png]"),
  make_option("--width",           type = "double",    default = 20,
              help = "Figure width in inches [default: 20]"),
  make_option("--height",          type = "double",    default = 0,
              help = "Figure height in inches; 0 = auto [default: 0]"),
  make_option("--bootstrap_min",   type = "double",    default = 80,
              help = "Min bootstrap to annotate [default: 80]"),
  make_option("--show_tip_labels", type = "character", default = "auto",
              help = "Show tip labels: auto/true/false [default: auto]"),
  make_option("--color_palette",   type = "character", default = "Set3",
              help = "RColorBrewer palette for clade colors [default: Set3]")
)

parser <- OptionParser(option_list = option_list)
opt    <- parse_args(parser)

if (is.null(opt$treefile) || is.null(opt$synteny_tsv)) {
  stop("--treefile and --synteny_tsv are required.")
}
if (!file.exists(opt$treefile)) stop(paste("treefile not found:", opt$treefile))
if (!file.exists(opt$synteny_tsv)) stop(paste("synteny_tsv not found:", opt$synteny_tsv))

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
formats <- trimws(strsplit(opt$format, ",")[[1]])

# ---------------------------------------------------------------------------
# Helper: extract taxonomy rank label
# ---------------------------------------------------------------------------
extract_tax_label <- function(lineage, level) {
  rank_prefix <- switch(tolower(level),
    domain  = "d__", phylum = "p__", class  = "c__", order = "o__",
    family  = "f__", genus  = "g__", species = "s__",
    paste0(level, "__"))
  tokens <- trimws(strsplit(as.character(lineage), ";")[[1]])
  hit <- tokens[startsWith(tolower(tokens), tolower(rank_prefix))]
  if (length(hit) == 0) return(NA_character_)
  sub(paste0("(?i)", rank_prefix), "", hit[1])
}

# ---------------------------------------------------------------------------
# Helper: build clade color palette
# ---------------------------------------------------------------------------
get_clade_colors <- function(names_vec, palette_name) {
  n <- length(names_vec)
  if (n == 0) return(setNames(character(0), character(0)))
  max_b <- brewer.pal.info[palette_name, "maxcolors"]
  if (n <= max_b) {
    cols <- brewer.pal(max(3, n), palette_name)[seq_len(n)]
  } else {
    cols <- colorRampPalette(brewer.pal(max_b, palette_name))(n)
  }
  setNames(cols, names_vec)
}

# ---------------------------------------------------------------------------
# Load tree
# ---------------------------------------------------------------------------
message("[synteny_tree] Reading tree: ", opt$treefile)
tree     <- treeio::read.newick(opt$treefile)
n_tips   <- length(tree@phylo$tip.label)

show_tips <- switch(tolower(opt$show_tip_labels),
  "true"  = TRUE, "false" = FALSE, "auto" = (n_tips <= 100), TRUE)
fig_height <- if (opt$height > 0) opt$height else max(6, n_tips * 0.15)

# ---------------------------------------------------------------------------
# Load synteny TSV
# ---------------------------------------------------------------------------
message("[synteny_tree] Reading synteny data: ", opt$synteny_tsv)
syn <- read.delim(opt$synteny_tsv, stringsAsFactors = FALSE)

required_syn_cols <- c("genome", "focal_protein", "position_relative",
                       "strand", "gene_label", "is_focal")
if (!all(required_syn_cols %in% colnames(syn))) {
  stop(paste("synteny_tsv missing columns:",
             paste(setdiff(required_syn_cols, colnames(syn)), collapse = ", ")))
}

# Build a tip-key matching genome|protein → tip label in the tree
# Tips are typically "genome|protein" form
syn$tip_key <- paste0(syn$genome, "|", syn$focal_protein)

# Also try just genome part if full key not in tree
tree_tips_set <- tree@phylo$tip.label

syn$matched_tip <- vapply(seq_len(nrow(syn)), function(i) {
  full_key <- syn$tip_key[i]
  if (full_key %in% tree_tips_set) return(full_key)
  # Try just the genome part
  genome_key <- syn$genome[i]
  hits <- tree_tips_set[startsWith(tree_tips_set, paste0(genome_key, "|"))]
  if (length(hits) > 0) return(hits[1])
  return(NA_character_)
}, character(1))

# Only keep tracks whose focal protein is in the tree
syn_matched <- syn %>%
  filter(!is.na(matched_tip)) %>%
  # Deduplicate: keep each (matched_tip, position_relative) once
  distinct(matched_tip, position_relative, .keep_all = TRUE)

if (nrow(syn_matched) == 0) {
  message("[synteny_tree] WARNING: no synteny rows matched tree tips — ",
          "synteny panel will be empty.")
}

# ---------------------------------------------------------------------------
# Load clades TSV
# ---------------------------------------------------------------------------
clades_df <- NULL
if (!is.null(opt$clades_tsv) && file.exists(opt$clades_tsv)) {
  clades_df <- read.delim(opt$clades_tsv, stringsAsFactors = FALSE)
  if (all(c("clade_name", "tip") %in% colnames(clades_df))) {
    clades_df$clade_name <- sub("^[^|]+\\|", "", clades_df$clade_name)
    tip_clades <- clades_df %>% group_by(tip) %>% slice(1) %>% ungroup()
  } else {
    clades_df <- NULL
  }
}

# ---------------------------------------------------------------------------
# Load taxonomy TSV
# ---------------------------------------------------------------------------
tax_df <- NULL
if (!is.null(opt$taxonomy_tsv) && file.exists(opt$taxonomy_tsv)) {
  raw_tax <- read.delim(opt$taxonomy_tsv, stringsAsFactors = FALSE)
  if (all(c("genome", "taxonomy") %in% colnames(raw_tax))) {
    tax_lookup <- setNames(raw_tax$taxonomy, raw_tax$genome)
    tips_vec   <- tree@phylo$tip.label
    tax_labels <- vapply(tips_vec, function(tip) {
      gpart  <- sub("\\|.*$", "", tip)
      lin    <- tax_lookup[gpart]
      if (is.na(lin)) lin <- tax_lookup[tip]
      if (is.na(lin)) return(NA_character_)
      extract_tax_label(lin, opt$tax_level)
    }, character(1))
    tax_df <- data.frame(tip = tips_vec, tax_label = tax_labels,
                         stringsAsFactors = FALSE)
  }
}

# ---------------------------------------------------------------------------
# Base tree plot
# ---------------------------------------------------------------------------
message("[synteny_tree] Building tree panel...")

p <- ggtree(tree, layout = "rectangular") + theme_tree2()

# Bootstrap annotations
tree_data <- as_tibble(tree)
boot_data <- tree_data %>%
  filter(!isTip) %>%
  mutate(support = suppressWarnings(as.numeric(label))) %>%
  filter(!is.na(support) & support >= opt$bootstrap_min)

if (nrow(boot_data) > 0) {
  p <- p + geom_nodepoint(
    data = boot_data,
    aes(subset = !isTip & !is.na(suppressWarnings(as.numeric(label))) &
          suppressWarnings(as.numeric(label)) >= opt$bootstrap_min),
    color = "steelblue", alpha = 0.7, size = 1.2
  )
}

if (show_tips) {
  p <- p + geom_tiplab(size = 1.8, align = FALSE)
}

# ---------------------------------------------------------------------------
# Clade highlight blocks
# ---------------------------------------------------------------------------
if (!is.null(clades_df)) {
  cn_sorted    <- sort(unique(clades_df$clade_name))
  clade_colors <- get_clade_colors(cn_sorted, opt$color_palette)

  for (cn in cn_sorted) {
    tips_in <- intersect(
      clades_df %>% filter(clade_name == cn) %>% pull(tip),
      tree@phylo$tip.label
    )
    if (length(tips_in) == 0) next
    nid <- if (length(tips_in) == 1) {
      which(tree@phylo$tip.label == tips_in[1])
    } else {
      tryCatch(getMRCA(tree@phylo, tips_in), error = function(e) NULL)
    }
    if (is.null(nid) || is.na(nid)) next
    p <- p +
      geom_highlight(node = nid, fill = clade_colors[cn], alpha = 0.18, extend = 0.3) +
      geom_cladelab(node = nid, label = cn, fontsize = 2.5,
                    offset = 0.01, barsize = 0.6,
                    barcolor = clade_colors[cn], color = clade_colors[cn])
  }
}

# ---------------------------------------------------------------------------
# Taxonomy tile strip
# ---------------------------------------------------------------------------
if (!is.null(tax_df)) {
  tax_ann <- tax_df %>% filter(!is.na(tax_label))
  if (nrow(tax_ann) > 0) {
    n_taxa   <- length(unique(tax_ann$tax_label))
    tax_cols <- get_clade_colors(sort(unique(tax_ann$tax_label)), "Paired")
    p <- p + new_scale_fill() +
      geom_fruit(
        data    = tax_ann,
        geom    = geom_tile,
        mapping = aes(y = tip, fill = tax_label),
        offset  = 0.04, pwidth = 0.05
      ) +
      scale_fill_manual(name = opt$tax_level, values = tax_cols,
                        na.value = "grey90",
                        guide = guide_legend(ncol = max(1, ceiling(n_taxa / 20))))
  }
}

# ---------------------------------------------------------------------------
# Synteny annotation panels using geom_fruit
# ---------------------------------------------------------------------------
# Each distinct position_relative gets its own geom_fruit tile/segment layer.
# Color focal gene red; hmm-annotated genes in teal; others in grey.

if (nrow(syn_matched) > 0) {
  # Determine color for each gene
  syn_matched <- syn_matched %>%
    mutate(gene_color = case_when(
      is_focal == TRUE | is_focal == "True" | is_focal == 1 ~ "focal",
      hmm_annotation != "" & !is.na(hmm_annotation) ~ "hmm_hit",
      TRUE ~ "other"
    ))

  gene_fill_scale <- c(
    focal   = "#e63946",   # red — focal hit gene
    hmm_hit = "#2a9d8f",   # teal — other HMM-annotated neighbor
    other   = "#adb5bd"    # grey — unannotated neighbor
  )

  # Add a tile panel per distinct relative position
  positions <- sort(unique(syn_matched$position_relative))

  for (pos in positions) {
    pos_data <- syn_matched %>%
      filter(position_relative == pos) %>%
      select(matched_tip, gene_color, gene_label) %>%
      rename(tip = matched_tip)

    p <- p + new_scale_fill() +
      geom_fruit(
        data    = pos_data,
        geom    = geom_tile,
        mapping = aes(y = tip, fill = gene_color),
        offset  = if (pos == positions[1]) 0.06 else 0.005,
        pwidth  = 0.04
      ) +
      scale_fill_manual(
        values   = gene_fill_scale,
        na.value = "white",
        guide    = if (pos == 0) guide_legend(title = "Gene type") else "none"
      )
  }
}

# ---------------------------------------------------------------------------
# Final theme
# ---------------------------------------------------------------------------
p <- p +
  theme(
    legend.position = "right",
    legend.text     = element_text(size = 6),
    legend.title    = element_text(size = 7, face = "bold"),
    legend.key.size = unit(0.35, "cm"),
    plot.title      = element_text(size = 9, face = "bold"),
    plot.margin     = margin(5, 40, 5, 5)
  ) +
  ggtitle(paste("Synteny + Tree:", opt$hmm_name))

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
for (fmt in formats) {
  out_fp <- file.path(opt$outdir, paste0("synteny_tree.", opt$hmm_name, ".", fmt))
  message("[synteny_tree] Saving: ", out_fp)
  tryCatch(
    ggsave(out_fp, plot = p, width = opt$width, height = fig_height,
           limitsize = FALSE, dpi = 150),
    error = function(e) warning("[synteny_tree] ggsave failed: ", conditionMessage(e))
  )
}

message("[synteny_tree] Done.")

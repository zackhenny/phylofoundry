#!/usr/bin/env Rscript
# plot_iqtree_ggtree.R
#
# Produce annotated phylogenetic tree plots from IQ-TREE output using ggtree.
# Supports clade highlights (from detect_clades TSV) and taxonomy annotation
# panels (from taxonomy_integrate TSV).
#
# Usage:
#   Rscript plot_iqtree_ggtree.R \
#     --treefile /path/to/hmm.treefile \
#     --outdir   /path/to/tree_viz/ \
#     --hmm_name MyHMM \
#     [--clades_tsv  /path/to/hmm.clades.tsv] \
#     [--taxonomy_tsv /path/to/genome_taxonomy.tsv] \
#     [--tax_level genus] \
#     [--format png,pdf] \
#     [--width 10] \
#     [--height 0] \
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
  make_option("--treefile",       type = "character", default = NULL,
              help = "Path to IQ-TREE .treefile (Newick format) [required]"),
  make_option("--clades_tsv",     type = "character", default = NULL,
              help = "TSV with columns clade_name,tip (from detect_clades step)"),
  make_option("--taxonomy_tsv",   type = "character", default = NULL,
              help = "TSV with columns genome,taxonomy (from taxonomy_integrate step)"),
  make_option("--tax_level",      type = "character", default = "genus",
              help = "Taxonomy rank to display [default: genus]"),
  make_option("--outdir",         type = "character", default = ".",
              help = "Output directory [default: current directory]"),
  make_option("--hmm_name",       type = "character", default = "tree",
              help = "Label used for output filenames [default: tree]"),
  make_option("--format",         type = "character", default = "png",
              help = "Comma-separated output formats: png,pdf,svg [default: png]"),
  make_option("--width",          type = "double",    default = 10,
              help = "Figure width in inches [default: 10]"),
  make_option("--height",         type = "double",    default = 0,
              help = "Figure height in inches; 0 = auto-scale with tip count [default: 0]"),
  make_option("--bootstrap_min",  type = "double",    default = 80,
              help = "Min bootstrap value to annotate on nodes [default: 80]"),
  make_option("--show_tip_labels",type = "character", default = "auto",
              help = "Show tip labels: auto (hide if >100 tips), true, false [default: auto]"),
  make_option("--color_palette",  type = "character", default = "Set3",
              help = "RColorBrewer palette for clade highlighting [default: Set3]")
)

parser <- OptionParser(option_list = option_list,
                       description = "Annotated ggtree plot from IQ-TREE output")
opt    <- parse_args(parser)

if (is.null(opt$treefile)) {
  stop("--treefile is required.")
}
if (!file.exists(opt$treefile)) {
  stop(paste("treefile not found:", opt$treefile))
}

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

formats <- trimws(strsplit(opt$format, ",")[[1]])

# ---------------------------------------------------------------------------
# Helper: extract one taxonomy rank token from a GTDB lineage string
# e.g. "d__Bacteria;p__Proteobacteria;...;g__Escherichia;s__..." -> "g__Escherichia"
# ---------------------------------------------------------------------------
extract_tax_label <- function(lineage, level) {
  rank_prefix <- switch(tolower(level),
    domain  = "d__",
    phylum  = "p__",
    class   = "c__",
    order   = "o__",
    family  = "f__",
    genus   = "g__",
    species = "s__",
    paste0(level, "__")
  )
  tokens <- trimws(strsplit(as.character(lineage), ";")[[1]])
  hit <- tokens[startsWith(tolower(tokens), tolower(rank_prefix))]
  if (length(hit) == 0) return(NA_character_)
  # Strip prefix for a cleaner label (e.g. "Escherichia" from "g__Escherichia")
  sub(paste0("(?i)", rank_prefix), "", hit[1])
}

# ---------------------------------------------------------------------------
# Load tree
# ---------------------------------------------------------------------------
message("[ggtree] Reading tree: ", opt$treefile)
tree <- treeio::read.newick(opt$treefile)
n_tips <- length(tree@phylo$tip.label)
message("[ggtree] Tips: ", n_tips)

# ---------------------------------------------------------------------------
# Load optional clades TSV
# ---------------------------------------------------------------------------
clades_df <- NULL
if (!is.null(opt$clades_tsv) && file.exists(opt$clades_tsv)) {
  clades_df <- read.delim(opt$clades_tsv, stringsAsFactors = FALSE)
  # Expected columns: clade_name, tip
  if (!all(c("clade_name", "tip") %in% colnames(clades_df))) {
    warning("[ggtree] clades_tsv missing 'clade_name' or 'tip' columns — ignoring.")
    clades_df <- NULL
  } else {
    # Remove clade prefix of the form "HMM|" so the clade name is cleaner
    clades_df$clade_name <- sub("^[^|]+\\|", "", clades_df$clade_name)
    # Build tip-level clade assignment (first clade wins for each tip)
    tip_clades <- clades_df %>%
      group_by(tip) %>%
      slice(1) %>%
      ungroup()
    message("[ggtree] Loaded ", nrow(tip_clades), " clade assignments for ",
            length(unique(tip_clades$clade_name)), " clades.")
  }
}

# ---------------------------------------------------------------------------
# Load optional taxonomy TSV
# ---------------------------------------------------------------------------
tax_df <- NULL
if (!is.null(opt$taxonomy_tsv) && file.exists(opt$taxonomy_tsv)) {
  raw_tax <- read.delim(opt$taxonomy_tsv, stringsAsFactors = FALSE)

  # Support two formats:
  #  1. genome_taxonomy.tsv: columns genome, taxonomy
  #  2. best_hits.with_taxonomy.tsv: columns genome, protein, taxonomy
  if ("genome" %in% colnames(raw_tax) && "taxonomy" %in% colnames(raw_tax)) {
    # Try to match tips: tips are "genome|protein" or just "genome"
    # Build both genome-level and tip-level lookups
    tax_lookup <- setNames(raw_tax$taxonomy, raw_tax$genome)

    tips_vec <- tree@phylo$tip.label
    tax_labels <- vapply(tips_vec, function(tip) {
      # Strip path separators and try genome part before the first "|"
      genome_part <- sub("\\|.*$", "", tip)
      lineage <- tax_lookup[genome_part]
      if (is.na(lineage)) lineage <- tax_lookup[tip]
      if (is.na(lineage)) return(NA_character_)
      extract_tax_label(lineage, opt$tax_level)
    }, character(1))

    tax_df <- data.frame(tip = tips_vec, tax_label = tax_labels,
                         stringsAsFactors = FALSE)
    n_annotated <- sum(!is.na(tax_df$tax_label))
    message("[ggtree] Taxonomy annotations: ", n_annotated, "/", n_tips,
            " tips matched at level '", opt$tax_level, "'.")
  } else {
    warning("[ggtree] taxonomy_tsv missing expected columns — ignoring.")
  }
}

# ---------------------------------------------------------------------------
# Decide whether to show tip labels
# ---------------------------------------------------------------------------
show_tips <- switch(tolower(opt$show_tip_labels),
  "true"  = TRUE,
  "false" = FALSE,
  "auto"  = (n_tips <= 100),
  TRUE
)

# ---------------------------------------------------------------------------
# Figure height: auto-scale to 0.15 in per tip, minimum 6 in
# ---------------------------------------------------------------------------
fig_height <- if (opt$height > 0) opt$height else max(6, n_tips * 0.15)

# ---------------------------------------------------------------------------
# Color palette for clades
# ---------------------------------------------------------------------------
n_clades <- if (!is.null(clades_df)) length(unique(clades_df$clade_name)) else 0

get_clade_colors <- function(clade_names, palette_name) {
  n <- length(clade_names)
  if (n == 0) return(character(0))
  # brewer.pal has limits; fall back to colorRampPalette for large n
  max_brewer <- brewer.pal.info[palette_name, "maxcolors"]
  if (n <= max_brewer) {
    cols <- brewer.pal(max(3, n), palette_name)[seq_len(n)]
  } else {
    base_cols <- brewer.pal(max_brewer, palette_name)
    cols <- colorRampPalette(base_cols)(n)
  }
  setNames(cols, clade_names)
}

# ---------------------------------------------------------------------------
# Build the base ggtree plot
# ---------------------------------------------------------------------------
message("[ggtree] Building tree plot...")

p <- ggtree(tree, layout = "rectangular") +
  theme_tree2()

# Bootstrap support annotations on internal nodes
# treeio stores bootstrap values as node labels; parse them numerically
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
    color = "steelblue", alpha = 0.7, size = 1.5
  ) +
  geom_nodelab(
    aes(subset = !isTip & !is.na(suppressWarnings(as.numeric(label))) &
          suppressWarnings(as.numeric(label)) >= opt$bootstrap_min,
        label = label),
    size = 1.8, hjust = -0.1, color = "steelblue3"
  )
}

# Tip labels
if (show_tips) {
  p <- p + geom_tiplab(size = 2, align = FALSE)
}

# ---------------------------------------------------------------------------
# Clade highlights and labels
# ---------------------------------------------------------------------------
if (!is.null(clades_df) && n_clades > 0) {
  clade_names_sorted <- sort(unique(clades_df$clade_name))
  clade_colors <- get_clade_colors(clade_names_sorted, opt$color_palette)

  # geom_highlight requires MRCA node for each clade
  for (cn in clade_names_sorted) {
    tips_in_clade <- clades_df %>%
      filter(clade_name == cn) %>%
      pull(tip)
    # Keep only tips that are actually in the tree
    tips_in_clade <- intersect(tips_in_clade, tree@phylo$tip.label)
    if (length(tips_in_clade) < 1) next

    if (length(tips_in_clade) == 1) {
      # Single-tip clade: highlight just the tip
      node_id <- which(tree@phylo$tip.label == tips_in_clade[1])
    } else {
      node_id <- tryCatch(
        getMRCA(tree@phylo, tips_in_clade),
        error = function(e) NULL
      )
    }
    if (is.null(node_id) || is.na(node_id)) next

    col <- clade_colors[cn]
    p <- p + geom_highlight(
      node = node_id, fill = col, alpha = 0.2, extend = 0.5
    )
    p <- p + geom_cladelab(
      node = node_id, label = cn,
      fontsize = 2.8, offset = 0.02, barsize = 0.8, barcolor = col, color = col
    )
  }

  # Add a legend for clade colors via a dummy data layer
  dummy_df <- data.frame(
    x = NA_real_, y = NA_real_,
    clade = factor(clade_names_sorted, levels = clade_names_sorted)
  )
  p <- p + geom_point(
    data = dummy_df, aes(x = x, y = y, fill = clade),
    shape = 22, size = 4, inherit.aes = FALSE
  ) +
  scale_fill_manual(
    name  = "Clade",
    values = clade_colors,
    guide  = guide_legend(override.aes = list(alpha = 0.8))
  )
}

# ---------------------------------------------------------------------------
# Taxonomy annotation panel (geom_fruit tile strip)
# ---------------------------------------------------------------------------
if (!is.null(tax_df)) {
  # Only annotate tips that have a taxonomy label
  tax_annotated <- tax_df %>% filter(!is.na(tax_label))
  if (nrow(tax_annotated) > 0) {
    n_taxa <- length(unique(tax_annotated$tax_label))
    tax_colors <- get_clade_colors(sort(unique(tax_annotated$tax_label)), "Paired")

    p <- p + new_scale_fill() +
      geom_fruit(
        data    = tax_annotated,
        geom    = geom_tile,
        mapping = aes(y = tip, fill = tax_label),
        offset  = 0.05, pwidth = 0.06
      ) +
      scale_fill_manual(
        name   = opt$tax_level,
        values = tax_colors,
        guide  = guide_legend(
          ncol           = max(1, ceiling(n_taxa / 20)),
          override.aes   = list(size = 4)
        ),
        na.value = "grey90"
      )
  }
}

# ---------------------------------------------------------------------------
# Final theme tweaks
# ---------------------------------------------------------------------------
p <- p +
  theme(
    legend.position    = "right",
    legend.text        = element_text(size = 7),
    legend.title       = element_text(size = 8, face = "bold"),
    legend.key.size    = unit(0.4, "cm"),
    plot.title         = element_text(size = 10, face = "bold"),
    plot.margin        = margin(5, 40, 5, 5)
  ) +
  ggtitle(paste("IQ-TREE:", opt$hmm_name))

# ---------------------------------------------------------------------------
# Save outputs
# ---------------------------------------------------------------------------
for (fmt in formats) {
  out_fp <- file.path(opt$outdir, paste0(opt$hmm_name, ".tree.", fmt))
  message("[ggtree] Saving: ", out_fp)
  tryCatch(
    ggsave(out_fp, plot = p, width = opt$width, height = fig_height,
           limitsize = FALSE, dpi = 150),
    error = function(e) warning("[ggtree] ggsave failed: ", conditionMessage(e))
  )
}

message("[ggtree] Done.")

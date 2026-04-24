#!/usr/bin/env Rscript
# core_microbiome.R
# Detects core microbiome members and writes CSV outputs (no plotting).

suppressPackageStartupMessages({
  library(optparse)
  library(phyloseq)
  library(microbiome)
})

# Source loader (resolved relative to this script's location)
script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    f    <- sub("--file=", "", args[grep("--file=", args)])
    if (length(f)) dirname(normalizePath(f)) else getwd()
  }
)
source(file.path(script_dir, "load_feature_table.R"))

# ---- CLI options -------------------------------------------------------------

option_list <- list(
  # Input
  make_option(c("--feature_table"),  type = "character", default = NULL,
              help = "Path to feature table (BIOM, TSV, or GTDB)"),
  make_option(c("--biom_file"),      type = "character", default = NULL,
              help = "Alias for --feature_table (legacy)"),
  make_option(c("--taxonomy_table"), type = "character", default = "",
              help = "Path to taxonomy table (required for tsv/gtdb)"),
  make_option(c("--input_format"),   type = "character", default = "biom",
              help = "Input format: biom | tsv | gtdb [default: biom]"),
  make_option(c("--meta_table"),     type = "character", default = NULL,
              help = "Path to metadata CSV"),
  make_option(c("--output_dir"),     type = "character", default = ".",
              help = "Output directory [default: .]"),

  # Filtering
  make_option(c("--which_level"),       type = "character", default = "Phylum",
              help = "Taxonomic level: Kingdom|Phylum|Class|Order|Family|Genus|Otus [default: Phylum]"),
  make_option(c("--label"),             type = "character", default = "analysis",
              help = "Label prefix for output files [default: analysis]"),
  make_option(c("--min_library_size"),  type = "integer",   default = 5000L,
              help = "Minimum library size to retain sample [default: 5000]"),
  make_option(c("--exclude_column"),    type = "character", default = "",
              help = "Metadata column to use for sample exclusion"),
  make_option(c("--exclude_values"),    type = "character", default = "",
              help = "Comma-separated values in exclude_column to drop"),

  # Grouping
  make_option(c("--groups_column"),       type = "character", default = "",
              help = "Metadata column to use as Groups"),
  make_option(c("--groups_paste_columns"), type = "character", default = "",
              help = "Comma-separated metadata columns to paste together as Groups"),
  make_option(c("--type_column"),         type = "character", default = ""),
  make_option(c("--type2_column"),        type = "character", default = ""),
  make_option(c("--type2_levels"),        type = "character", default = ""),

  # Core microbiome options
  make_option(c("--prevalence_min"),      type = "double",    default = 0.85,
              help = "Minimum prevalence for core detection [default: 0.85]"),
  make_option(c("--what_detection"),      type = "character", default = "absolute",
              help = "Detection mode: absolute | relative [default: absolute]"),
  make_option(c("--detection_range_min"), type = "double",    default = -1,
              help = "Min detection threshold (-1 = auto) [default: -1]"),
  make_option(c("--detection_range_max"), type = "double",    default = -1,
              help = "Max detection threshold (-1 = auto) [default: -1]"),
  make_option(c("--detection_steps"),     type = "integer",   default = 20L,
              help = "Number of detection threshold steps [default: 20]")
)

opt_parser <- OptionParser(option_list = option_list)
opt        <- parse_args(opt_parser)

# Resolve --biom_file alias
if (is.null(opt$feature_table) && !is.null(opt$biom_file)) {
  opt$feature_table <- opt$biom_file
}
if (is.null(opt$feature_table)) {
  print_help(opt_parser)
  stop("--feature_table (or --biom_file) is required.", call. = FALSE)
}
if (is.null(opt$meta_table)) {
  print_help(opt_parser)
  stop("--meta_table is required.", call. = FALSE)
}

# ---- Prepare output dir -----------------------------------------------------

dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)
which_level <- opt$which_level

# ---- Load data --------------------------------------------------------------

message("==> Loading feature table ...")
ft_result <- load_feature_table(
  feature_table  = opt$feature_table,
  input_format   = opt$input_format,
  taxonomy_table = if (nchar(opt$taxonomy_table) > 0) opt$taxonomy_table else NULL
)
abund_table  <- ft_result$abund_table
OTU_taxonomy <- ft_result$OTU_taxonomy

message("==> Loading metadata from: ", opt$meta_table)
meta_table <- utils::read.csv(opt$meta_table, row.names = 1,
                              check.names = FALSE, stringsAsFactors = FALSE)

# ---- Sample exclusion -------------------------------------------------------

if (nchar(opt$exclude_column) > 0 && opt$exclude_column %in% colnames(meta_table)) {
  excl_vals <- trimws(unlist(strsplit(opt$exclude_values, ",")))
  keep_rows <- !(meta_table[[opt$exclude_column]] %in% excl_vals)
  meta_table  <- meta_table[keep_rows, , drop = FALSE]
  message("  Excluded ", sum(!keep_rows), " sample(s) via exclude_column '",
          opt$exclude_column, "'")
}

# Align samples
shared_samples <- intersect(rownames(abund_table), rownames(meta_table))
if (length(shared_samples) == 0)
  stop("No shared sample IDs between feature table and metadata.")
abund_table <- abund_table[shared_samples, , drop = FALSE]
meta_table  <- meta_table[shared_samples,  , drop = FALSE]

# ---- Library size filter ----------------------------------------------------

lib_sizes <- rowSums(abund_table)
pass_lib   <- lib_sizes >= opt$min_library_size
if (any(!pass_lib)) {
  message("  Dropping ", sum(!pass_lib), " sample(s) below min_library_size (",
          opt$min_library_size, "): ", paste(names(lib_sizes)[!pass_lib], collapse = ", "))
  abund_table <- abund_table[pass_lib, , drop = FALSE]
  meta_table  <- meta_table[pass_lib,  , drop = FALSE]
}

# Hardening: require >= 2 samples
if (nrow(abund_table) < 2) {
  stop("After filtering, fewer than 2 samples remain. Cannot compute core microbiome.")
}

# ---- Groups -----------------------------------------------------------------

.validate_col <- function(df, col, arg_name) {
  if (nchar(col) > 0 && !(col %in% colnames(df)))
    stop("Column '", col, "' (", arg_name, ") not found in metadata. ",
         "Available columns: ", paste(colnames(df), collapse = ", "))
}

.validate_col(meta_table, opt$groups_column,        "--groups_column")
.validate_col(meta_table, opt$type_column,           "--type_column")
.validate_col(meta_table, opt$type2_column,          "--type2_column")

if (nchar(opt$groups_paste_columns) > 0) {
  paste_cols <- trimws(unlist(strsplit(opt$groups_paste_columns, ",")))
  for (pc in paste_cols) .validate_col(meta_table, pc, "--groups_paste_columns")
}

if (nchar(opt$groups_column) > 0) {
  meta_table$Groups <- meta_table[[opt$groups_column]]
} else if (nchar(opt$groups_paste_columns) > 0) {
  paste_cols <- trimws(unlist(strsplit(opt$groups_paste_columns, ",")))
  meta_table$Groups <- apply(meta_table[, paste_cols, drop = FALSE], 1, paste, collapse = "_")
} else {
  message("  No groups_column or groups_paste_columns defined — assigning all samples to 'All'")
  meta_table$Groups <- "All"
}

# ---- Taxonomic collation ----------------------------------------------------

valid_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Otus")
if (!(which_level %in% valid_levels)) {
  stop("--which_level must be one of: ", paste(valid_levels, collapse = ", "))
}

if (which_level == "Otus") {
  # Keep features as-is; label from taxonomy
  rownames(OTU_taxonomy) <- make.unique(
    ifelse(OTU_taxonomy$Genus != "", OTU_taxonomy$Genus, rownames(OTU_taxonomy))
  )
  message("  Using OTU-level features (", ncol(abund_table), " features)")
} else {
  # Aggregate by taxonomic level
  level_idx <- match(which_level, valid_levels)
  if (level_idx < 1)
    stop("Invalid which_level: ", which_level)

  # Build a grouping vector
  tax_group <- OTU_taxonomy[[which_level]]
  tax_group[is.na(tax_group) | tax_group == ""] <- "Unknown"

  # Aggregate columns of abund_table by tax_group
  groups_unique <- unique(tax_group)
  agg_mat <- sapply(groups_unique, function(g) {
    cols <- colnames(abund_table)[tax_group == g]
    if (length(cols) == 1) abund_table[, cols]
    else rowSums(abund_table[, cols, drop = FALSE])
  })
  rownames(agg_mat) <- rownames(abund_table)
  abund_table <- agg_mat
  message("  Aggregated to ", which_level, " level: ", ncol(abund_table), " taxa")
}

# ---- Build phyloseq ---------------------------------------------------------

OTU_ps  <- otu_table(t(abund_table), taxa_are_rows = TRUE)
SAM_ps  <- sample_data(meta_table)
physeq  <- merge_phyloseq(phyloseq(OTU_ps), SAM_ps)
pseq.2  <- prune_taxa(taxa_sums(physeq) > 0, physeq)

# Relative transform
pseq.rel <- microbiome::transform(pseq.2, "compositional")

# ---- Detection range --------------------------------------------------------

prevalences <- seq(0.05, 1, 0.05)
steps       <- max(2L, opt$detection_steps)

if (opt$what_detection == "relative") {
  d_min_default <- 1e-3
  d_max_default <- 0.2

  d_min <- if (opt$detection_range_min > 0) opt$detection_range_min else d_min_default
  d_max <- if (opt$detection_range_max > 0) opt$detection_range_max else d_max_default

  detections   <- 10^seq(log10(d_min), log10(d_max), length.out = steps)
  pseq_to_plot <- pseq.rel

} else {
  # absolute
  max_abund <- max(microbiome::abundances(pseq.2), na.rm = TRUE)
  if (is.na(max_abund) || max_abund <= 0) {
    warning("max(abundances) is 0 or NA; using fallback detection range 1–100.")
    max_abund <- 100
  }

  d_min_default <- 1
  d_max_default <- max(max_abund / 10, d_min_default * 2)

  d_min <- if (opt$detection_range_min > 0) opt$detection_range_min else d_min_default
  d_max <- if (opt$detection_range_max > 0) opt$detection_range_max else d_max_default

  if (d_min >= d_max) {
    warning("detection_range_min >= detection_range_max; using fallback 1–", d_max_default)
    d_min <- d_min_default
    d_max <- d_max_default
  }

  detections   <- round(10^seq(log10(d_min), log10(d_max), length.out = steps))
  detections   <- unique(detections)
  pseq_to_plot <- pseq.2
}

message("  Detection range: ", round(min(detections), 5), " to ", round(max(detections), 5),
        " (", length(detections), " steps), mode = ", opt$what_detection)

# ---- Run plot_core ----------------------------------------------------------

core_df <- tryCatch({
  datacore <- plot_core(
    pseq_to_plot,
    plot.type       = "heatmap",
    prevalences     = prevalences,
    detections      = detections,
    min.prevalence  = opt$prevalence_min,
    horizontal      = TRUE
  )
  df      <- datacore$data
  df$Taxa              <- as.character(df$Taxa)
  df$DetectionThreshold <- as.numeric(as.character(df$DetectionThreshold))
  df$Prevalence        <- as.numeric(as.character(df$Prevalence))
  df
}, error = function(e) {
  warning("plot_core() failed: ", conditionMessage(e),
          ". Writing empty heatmap CSV with warning column.")
  data.frame(Taxa = character(0), DetectionThreshold = numeric(0),
             Prevalence = numeric(0), Warning = character(0))
})

# ---- Output 1: Core heatmap CSV ---------------------------------------------

heatmap_file <- file.path(opt$output_dir,
  paste0("Core_heatmap_", opt$label, "_", which_level, "_", opt$what_detection, ".csv"))
utils::write.csv(core_df, heatmap_file, row.names = FALSE)
message("Wrote heatmap data: ", heatmap_file)

# ---- Output 2: Core members CSV ---------------------------------------------

if (nrow(core_df) > 0 && !("Warning" %in% colnames(core_df))) {
  # Core members: taxa that meet min.prevalence at ALL detection thresholds
  # i.e., prevalence >= prevalence_min at every detection threshold where they appear.
  taxa_by_detection <- split(core_df, core_df$DetectionThreshold)

  core_taxa_per_step <- lapply(taxa_by_detection, function(sub) {
    sub$Taxa[sub$Prevalence >= opt$prevalence_min]
  })

  if (length(core_taxa_per_step) > 0) {
    # Taxa that pass at every step
    core_members_all <- Reduce(intersect, core_taxa_per_step)
  } else {
    core_members_all <- character(0)
  }

  if (length(core_members_all) > 0) {
    # For each core member, find the minimum detection threshold at which it is core
    members_df <- do.call(rbind, lapply(core_members_all, function(taxon) {
      sub <- core_df[core_df$Taxa == taxon & core_df$Prevalence >= opt$prevalence_min, ]
      data.frame(
        Taxa                     = taxon,
        min_detection_in_core    = min(sub$DetectionThreshold, na.rm = TRUE),
        prevalence_at_min_detection = sub$Prevalence[which.min(sub$DetectionThreshold)],
        stringsAsFactors         = FALSE
      )
    }))
    members_df <- members_df[order(members_df$min_detection_in_core), ]
  } else {
    message("  No core members found at min.prevalence = ", opt$prevalence_min,
            " across all detection thresholds.")
    members_df <- data.frame(
      Taxa                     = character(0),
      min_detection_in_core    = numeric(0),
      prevalence_at_min_detection = numeric(0),
      stringsAsFactors         = FALSE
    )
  }
} else {
  members_df <- data.frame(
    Taxa                     = character(0),
    min_detection_in_core    = numeric(0),
    prevalence_at_min_detection = numeric(0),
    stringsAsFactors         = FALSE
  )
}

members_file <- file.path(opt$output_dir,
  paste0("Core_members_", opt$label, "_", which_level, "_", opt$what_detection, ".csv"))
utils::write.csv(members_df, members_file, row.names = FALSE)
message("Wrote core members: ", members_file, " (", nrow(members_df), " taxa)")

message("==> core_microbiome.R complete.")

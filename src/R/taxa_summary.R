#!/usr/bin/env Rscript
# taxa_summary.R
# Produces a tidy long-format taxa summary CSV from a feature table.

suppressPackageStartupMessages({
  library(optparse)
  library(phyloseq)
  library(vegan)
})

# Source loader
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
  make_option(c("--groups_column"),         type = "character", default = "",
              help = "Metadata column to use as Groups"),
  make_option(c("--groups_paste_columns"),  type = "character", default = "",
              help = "Comma-separated metadata columns to paste together as Groups"),
  make_option(c("--taxa_groups_paste_columns"), type = "character", default = "",
              help = "Alternative groups_paste_columns for taxa summary"),

  # Taxa summary options
  make_option(c("--top_n"),      type = "integer",   default = 25L,
              help = "Number of top taxa to retain [default: 25]"),
  make_option(c("--normalise"),  type = "character", default = "relative",
              help = "Normalisation: relative | rarefied | clr [default: relative]")
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

normalise <- tolower(trimws(opt$normalise))
if (!(normalise %in% c("relative", "rarefied", "clr"))) {
  stop("--normalise must be one of: relative, rarefied, clr")
}

# Merge taxa_groups_paste_columns into groups_paste_columns if provided
if (nchar(opt$taxa_groups_paste_columns) > 0 && nchar(opt$groups_paste_columns) == 0) {
  opt$groups_paste_columns <- opt$taxa_groups_paste_columns
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
  excl_vals  <- trimws(unlist(strsplit(opt$exclude_values, ",")))
  keep_rows  <- !(meta_table[[opt$exclude_column]] %in% excl_vals)
  meta_table <- meta_table[keep_rows, , drop = FALSE]
  message("  Excluded ", sum(!keep_rows), " sample(s) via exclude_column")
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
  stop("After filtering, fewer than 2 samples remain. Cannot compute taxa summary.")
}

# ---- Groups -----------------------------------------------------------------

.validate_col <- function(df, col, arg_name) {
  if (nchar(col) > 0 && !(col %in% colnames(df)))
    stop("Column '", col, "' (", arg_name, ") not found in metadata. ",
         "Available columns: ", paste(colnames(df), collapse = ", "))
}

.validate_col(meta_table, opt$groups_column, "--groups_column")

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
  message("  Using OTU-level features (", ncol(abund_table), " features)")
} else {
  tax_group     <- OTU_taxonomy[[which_level]]
  tax_group[is.na(tax_group) | tax_group == ""] <- "Unknown"
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

# ---- Normalisation ----------------------------------------------------------

message("  Applying normalisation: ", normalise)

if (normalise == "relative") {
  rs <- rowSums(abund_table)
  rs[rs == 0] <- 1
  x <- sweep(abund_table, 1, rs, "/")

} else if (normalise == "rarefied") {
  min_depth <- min(rowSums(abund_table))
  if (min_depth < 1) stop("Cannot rarefy: min library size after filtering is 0.")
  message("  Rarefying to depth: ", min_depth)
  set.seed(42L)
  rarefied <- vegan::rrarefy(abund_table, min_depth)
  rs <- rowSums(rarefied)
  rs[rs == 0] <- 1
  x <- sweep(rarefied, 1, rs, "/")

} else {
  # clr
  # Add pseudocount of 0.5 to avoid log(0)
  pseudo   <- abund_table + 0.5
  log_x    <- log(pseudo)
  geo_mean <- rowMeans(log_x)
  x        <- sweep(log_x, 1, geo_mean, "-")
  # CLR values can be negative; keep as-is
}

# ---- Top-N selection --------------------------------------------------------

N <- opt$top_n

# Sort taxa by total abundance (column sums)
col_totals <- colSums(x)
x <- x[, order(col_totals, decreasing = TRUE), drop = FALSE]

taxa_list <- colnames(x)[seq_len(min(ncol(x), N))]

# Remove __Unknowns__ from top-N (push it to Others)
unknown_pat <- "__Unknowns__|^Unknown$"
if (any(grepl(unknown_pat, taxa_list))) {
  extended  <- colnames(x)[seq_len(min(ncol(x), N + sum(grepl(unknown_pat, taxa_list))))]
  taxa_list <- extended[!grepl(unknown_pat, extended)]
  taxa_list <- taxa_list[seq_len(min(length(taxa_list), N))]
}
N <- length(taxa_list)
message("  Top taxa selected: ", N)

# Build the final matrix with an 'Others' column if needed
in_top  <- colnames(x) %in% taxa_list
out_top <- !in_top

if (all(in_top)) {
  new_x <- as.data.frame(x[, taxa_list, drop = FALSE])
} else {
  others <- rowSums(x[, out_top, drop = FALSE])
  new_x  <- as.data.frame(cbind(x[, taxa_list, drop = FALSE], Others = others))
}

# ---- Build tidy output dataframe --------------------------------------------

df <- do.call(rbind, lapply(seq_len(ncol(new_x)), function(i) {
  data.frame(
    Sample = rownames(new_x),
    Taxa   = colnames(new_x)[i],
    Value  = new_x[, i],
    Groups = meta_table$Groups,
    stringsAsFactors = FALSE
  )
}))
rownames(df) <- NULL

# ---- Write output -----------------------------------------------------------

out_file <- file.path(opt$output_dir,
  paste0("Taxa_summary_", which_level, "_", opt$label, ".csv"))
utils::write.csv(df, out_file, row.names = FALSE)
message("Wrote taxa summary: ", out_file)
message("==> taxa_summary.R complete.")

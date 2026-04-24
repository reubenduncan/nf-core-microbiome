suppressPackageStartupMessages({
  library(optparse)
  library(arrow)
})

option_list <- list(
  make_option("--label",      type = "character", default = "analysis",
              help = "Label used in output filename [default: analysis]"),
  make_option("--output_dir", type = "character", default = ".",
              help = "Directory for the output Parquet file [default: .]")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (!dir.exists(opt$output_dir))
  dir.create(opt$output_dir, recursive = TRUE)

classify_csv <- function(fname) {
  b <- basename(fname)
  if      (grepl("^Core_heatmap_", b)) list(analysis = "core", table = "core_heatmap")
  else if (grepl("^Core_members_", b)) list(analysis = "core", table = "core_members")
  else                                  list(analysis = "unknown", table = "unknown")
}

csv_files <- list.files(".", pattern = "\\.csv$", full.names = TRUE, recursive = FALSE)
if (length(csv_files) == 0)
  stop("No CSV files found in working directory.")

message("Found ", length(csv_files), " CSV file(s) to merge.")

frames <- lapply(csv_files, function(f) {
  cls <- classify_csv(f)
  df  <- tryCatch(
    read.csv(f, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) { warning("Failed to read ", basename(f), ": ", conditionMessage(e)); NULL }
  )
  if (is.null(df)) return(NULL)

  if (ncol(df) > 0 && colnames(df)[1] == "") df <- df[, -1, drop = FALSE]
  idx_cols <- grepl("^X$|^X\\.\\d+$", colnames(df))
  if (any(idx_cols)) df <- df[, !idx_cols, drop = FALSE]

  df$analysis    <- cls$analysis
  df$table       <- cls$table
  df$source_file <- basename(f)

  message("  ", basename(f), " → analysis=", cls$analysis, ", table=", cls$table,
          " (", nrow(df), " rows)")
  df
})

frames <- Filter(Negate(is.null), frames)
if (length(frames) == 0)
  stop("All CSV files failed to read — cannot write Parquet.")

all_cols       <- unique(unlist(lapply(frames, colnames)))
frames_aligned <- lapply(frames, function(df) {
  df[setdiff(all_cols, colnames(df))] <- NA
  df[, all_cols, drop = FALSE]
})

merged <- do.call(rbind, frames_aligned)
rownames(merged) <- NULL

out_path <- file.path(opt$output_dir, paste0("core_microbiome_", opt$label, ".parquet"))
arrow::write_parquet(merged, out_path)
message("Written: ", out_path, " (", nrow(merged), " rows x ", ncol(merged), " columns)")

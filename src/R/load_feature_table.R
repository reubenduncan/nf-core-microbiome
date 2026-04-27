# load_feature_table.R
# Loads a feature table from BIOM, TSV, or GTDB format and returns a
# standardised list: abund_table (samples x features) + feature_taxonomy (data.frame).

load_feature_table <- function(feature_table, input_format, taxonomy_table = NULL) {

  input_format <- tolower(trimws(input_format))
  stopifnot(input_format %in% c("biom", "tsv", "gtdb"))

  # ---- helpers ----------------------------------------------------------------

  .strip_prefixes <- function(tax_string, format = "qiime") {
    # Removes rank prefixes from a semicolon-separated taxonomy string.
    # format = "qiime" handles both D_0__ and d__ styles
    # format = "gtdb"  handles d__, p__, c__, o__, f__, g__, s__
    prefixes <- c(
      "D_0__", "D_1__", "D_2__", "D_3__", "D_4__", "D_5__", "D_6__",
      "d__",   "p__",   "c__",   "o__",   "f__",   "g__",   "s__",
      "k__"
    )
    parts <- trimws(unlist(strsplit(tax_string, ";")))
    for (p in prefixes) {
      parts <- sub(paste0("^", p), "", parts)
    }
    parts
  }

  .build_taxonomy_df <- function(feature_ids, tax_strings, format = "qiime") {
    ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Feature")
    rows <- lapply(tax_strings, function(ts) {
      parts <- .strip_prefixes(ts, format)
      # Pad or truncate to exactly 7 ranks
      length(parts) <- 7
      as.list(parts)
    })
    df <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
    colnames(df) <- ranks
    rownames(df) <- feature_ids
    df[is.na(df)] <- ""
    df
  }

  .prune_taxonomy <- function(abund_table, tax_df) {
    # Remove features with Unassigned Kingdom, blank Phylum,
    # Chloroplast Order, or Mitochondria Family.
    keep <- rep(TRUE, nrow(tax_df))

    if ("Kingdom" %in% colnames(tax_df)) {
      unassigned <- grepl("^[Uu]nassigned$|^[Uu]nclassified$|^$", tax_df$Kingdom)
      if (any(unassigned)) {
        message("  Pruning ", sum(unassigned), " feature(s) with Unassigned/blank Kingdom")
        keep <- keep & !unassigned
      }
    }

    if ("Phylum" %in% colnames(tax_df)) {
      blank_phylum <- tax_df$Phylum == "" | is.na(tax_df$Phylum)
      if (any(blank_phylum)) {
        message("  Pruning ", sum(blank_phylum), " feature(s) with blank Phylum")
        keep <- keep & !blank_phylum
      }
    }

    if ("Order" %in% colnames(tax_df)) {
      chloroplast <- grepl("Chloroplast", tax_df$Order, ignore.case = TRUE)
      if (any(chloroplast)) {
        message("  Pruning ", sum(chloroplast), " Chloroplast feature(s)")
        keep <- keep & !chloroplast
      }
    }

    if ("Family" %in% colnames(tax_df)) {
      mitochondria <- grepl("Mitochondria|Mitochondri[a-z]", tax_df$Family, ignore.case = TRUE)
      if (any(mitochondria)) {
        message("  Pruning ", sum(mitochondria), " Mitochondria feature(s)")
        keep <- keep & !mitochondria
      }
    }

    feat_keep <- rownames(tax_df)[keep]
    # abund_table is samples x features; subset columns
    shared <- intersect(feat_keep, colnames(abund_table))
    list(
      abund_table  = abund_table[, shared, drop = FALSE],
      feature_taxonomy = tax_df[shared, , drop = FALSE]
    )
  }

  # ---- BIOM -------------------------------------------------------------------

  if (input_format == "biom") {
    if (!requireNamespace("phyloseq", quietly = TRUE))
      stop("Package 'phyloseq' is required to import BIOM files.")

    message("Loading BIOM file: ", feature_table)
    ps <- phyloseq::import_biom(feature_table)

    # abund_table: samples x features
    abund_mat <- as.matrix(t(phyloseq::otu_table(ps)))
    if (!phyloseq::taxa_are_rows(ps)) {
      abund_mat <- as.matrix(phyloseq::otu_table(ps))
    }

    ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Feature")

    if (inherits(ps, "phyloseq")) {
      tax_mat <- as.matrix(phyloseq::tax_table(ps))
      colnames(tax_mat)[seq_len(min(ncol(tax_mat), 7))] <- ranks[seq_len(min(ncol(tax_mat), 7))]
      tax_df <- as.data.frame(tax_mat, stringsAsFactors = FALSE)
      for (r in colnames(tax_df)) {
        tax_df[[r]] <- sub("^D_[0-9]+__", "", tax_df[[r]])
        tax_df[[r]] <- sub("^[dpcofgs]__", "", tax_df[[r]])
        tax_df[[r]][is.na(tax_df[[r]])] <- ""
      }
      shared_feats <- intersect(rownames(tax_df), colnames(abund_mat))
      abund_mat <- abund_mat[, shared_feats, drop = FALSE]
      tax_df    <- tax_df[shared_feats, , drop = FALSE]
      result <- .prune_taxonomy(abund_mat, tax_df)
      message("BIOM loaded: ", nrow(result$abund_table), " samples x ",
              ncol(result$abund_table), " features (after pruning)")
      return(result)
    } else if (!is.null(taxonomy_table) && taxonomy_table != "" && file.exists(taxonomy_table)) {
      message("BIOM has no embedded taxonomy — loading separate taxonomy table: ", taxonomy_table)
      tax_raw <- tryCatch(
        utils::read.table(taxonomy_table, header = FALSE, sep = "\t",
                          check.names = FALSE, stringsAsFactors = FALSE, fill = TRUE),
        error = function(e) stop("Failed to read taxonomy table: ", conditionMessage(e))
      )
      feat_ids <- as.character(tax_raw[[1]])
      tax_strs <- as.character(tax_raw[[2]])
      tax_df   <- .build_taxonomy_df(feat_ids, tax_strs, "qiime")
      common   <- intersect(colnames(abund_mat), rownames(tax_df))
      if (length(common) == 0)
        stop("No feature IDs overlap between BIOM file and taxonomy table.")
      abund_mat <- abund_mat[, common, drop = FALSE]
      tax_df    <- tax_df[common, , drop = FALSE]
      return(.prune_taxonomy(abund_mat, tax_df))
    } else {
      message(
        "BIOM file has no embedded taxonomy. All features treated as features. ",
        "Use --taxon_rank Feature, or supply a separate --taxonomy_table."
      )
      tax_df <- as.data.frame(
        matrix("", nrow = ncol(abund_mat), ncol = 7,
               dimnames = list(colnames(abund_mat), ranks)),
        stringsAsFactors = FALSE
      )
      tax_df$Feature <- colnames(abund_mat)
      return(list(abund_table = abund_mat, feature_taxonomy = tax_df))
    }
  }

  # ---- TSV / GTDB -------------------------------------------------------------

  if (input_format %in% c("tsv", "gtdb")) {
    message("Loading TSV feature table: ", feature_table)

    ft <- utils::read.table(feature_table, header = TRUE, sep = "\t",
                            row.names = 1, check.names = FALSE,
                            comment.char = "", stringsAsFactors = FALSE)

    # Auto-detect orientation: if features > samples → transpose
    # Heuristic: a feature table typically has more features than samples.
    # If nrow > ncol we assume rows = features → transpose to get samples x features.
    if (nrow(ft) > ncol(ft)) {
      message("  Auto-detecting orientation: transposing (features were rows)")
      ft <- t(ft)
    }
    abund_mat <- as.matrix(ft)
    storage.mode(abund_mat) <- "numeric"

    # Taxonomy
    if (is.null(taxonomy_table) || taxonomy_table == "") {
      stop("--taxonomy_table is required for tsv/gtdb input format.")
    }
    message("Loading taxonomy table: ", taxonomy_table)
    tax_raw <- utils::read.table(taxonomy_table, header = TRUE, sep = "\t",
                                 check.names = FALSE, stringsAsFactors = FALSE)

    # col1 = featureID, col2 = taxonomy string
    feat_col <- colnames(tax_raw)[1]
    tax_col  <- colnames(tax_raw)[2]

    fmt <- if (input_format == "gtdb") "gtdb" else "qiime"
    tax_df <- .build_taxonomy_df(tax_raw[[feat_col]], tax_raw[[tax_col]], fmt)

    # Align feature IDs between abund_mat and tax_df
    shared <- intersect(colnames(abund_mat), rownames(tax_df))
    if (length(shared) == 0) {
      stop("No shared feature IDs between the feature table and taxonomy table.")
    }
    abund_mat <- abund_mat[, shared, drop = FALSE]
    tax_df    <- tax_df[shared, , drop = FALSE]

    result <- .prune_taxonomy(abund_mat, tax_df)
    message(toupper(input_format), " loaded: ", nrow(result$abund_table), " samples x ",
            ncol(result$abund_table), " features (after pruning)")
    return(result)
  }
}

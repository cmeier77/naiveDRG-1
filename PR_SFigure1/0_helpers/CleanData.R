#!/usr/bin/env Rscript
#-------------------------------------------------
# Title: Make DRG dataframe cleaner helper
# Author: Christina Meier
# Date: 2025-05-22
# Email: 19cm51@queensu.ca
#-------------------------------------------------
# Notes ------------------------------------------
# module load StdEnv/2023 r/4.5.0
#
# Options ----------------------------------------

# Packages ---------------------------------------
library(data.table) # 1.18.2.1
library(dplyr) # 1.1.4

# Function ---------------------------------------
cleanData <- function(df_id, input_dir, output_dir) {
  #' Custom function to only keep naive samples and remove other timepoints from DRG single cell datasets.
  #' @param df_id DRG dataset number (integer); matches in file names
  #' @param input_dir path to input folder (string)
  #' @param output_dir path to output folder (string)
  #' @return saves cleaned count and metadata RDS files to output_dir
  #' @example cleanData(df_id = 1, input_dir = "data/raw/GSE155622", output_dir = "data/sc_cleaned")
  #'
  counts <- fread(
    file.path(input_dir, sprintf("GSE155622_raw_UMI_counts_%d.txt", df_id)),
    sep = "\t"
  )

  metadata <- fread(
    file.path(input_dir, sprintf("GSE155622_raw_UMI_counts_%d_metadata.txt", df_id)),
    sep = "\t",
  )

  print(names(metadata))
  print(head(metadata))

  # make all Condition column names the same for filtering
  condition_col <- intersect(
    c("Conditions", "condition", "modify.ident"),
    names(metadata)
  )

  if (length(condition_col) != 1) {
    stop(
      "Expected exactly one condition column, found: ",
      paste(condition_col, collapse = ", "),
      "\nMetadata columns are: ",
      paste(names(metadata), collapse = ", ")
    )
  }

  metadata <- metadata |>
    rename(Conditions = all_of(condition_col))

  conditions_to_remove <- c("SNI 6h", "SNI 24h", "SNI 14d", "SNI 2d", "SNI 7d", "SNI 28d", "SNI 6h_2", "SNI 24h_2")

  metadata_filtered <- metadata |>
    filter(!(Conditions %in% conditions_to_remove))

  meta_ids <- metadata_filtered$V1
  count_cols <- colnames(counts)[-1] # first column is gene names

  missing <- !(meta_ids %in% count_cols)

  ids_fixed <- meta_ids

  ids_fixed[missing] <- ifelse(
    sub("-\\d+$", "-1", meta_ids[missing]) %in% count_cols,
    sub("-\\d+$", "-1", meta_ids[missing]),
    meta_ids[missing]
  )

  metadata_filtered$V1 <- ids_fixed

  still_missing <- setdiff(metadata_filtered$V1, colnames(counts))

  if (length(still_missing) > 0) {
    warning(
      length(still_missing), " metadata cells still not found in counts: ",
      paste(head(still_missing), collapse = ", ")
    )
  }

  counts_filtered <- counts |>
    select(V1, all_of(metadata_filtered$V1))

  metadata_filtered <- as.data.frame(metadata_filtered)
  rownames(metadata_filtered) <- metadata_filtered$V1

  # 1. Check sizes match
  if (ncol(counts_filtered) - 1 != nrow(metadata_filtered)) {
    stop("Mismatch: number of cells in counts and metadata differ")
  } else {
    message("DRG ", df_id, ": counts and metadata have matching number of cells (", ncol(counts_filtered) - 1, ")")
  }

  # 2. Check all IDs match (ignoring order)
  if (!setequal(colnames(counts_filtered)[-1], rownames(metadata_filtered))) {
    stop("Mismatch: counts and metadata contain different cell IDs")
  } else {
    message("DRG ", df_id, ": counts and metadata contain the same cell IDs")
  }

  # 3. Force correct ordering
  metadata_filtered <- metadata_filtered[
    colnames(counts_filtered)[-1], ,
    drop = FALSE
  ]

  # 4. Final strict check (order + identity)
  if (!identical(colnames(counts_filtered)[-1], rownames(metadata_filtered))) {
    stop("Ordering failed: metadata does not match counts")
  } else {
    message("DRG ", df_id, ": counts and metadata are in the same order")
  }

  message("DRG ", df_id, ": kept ", ncol(counts_filtered) - 1, " cell columns")
  message("DRG ", df_id, ": kept ", nrow(metadata_filtered), " metadata rows")

  saveRDS(
    counts_filtered,
    file.path(output_dir, sprintf("naive_drg%d_good.rds", df_id))
  )

  saveRDS(
    metadata_filtered,
    file.path(output_dir, sprintf("naive_drg%d_mt_good.rds", df_id))
  )
}

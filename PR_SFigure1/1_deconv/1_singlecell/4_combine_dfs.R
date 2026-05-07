#!/usr/bin/env Rscript
#-------------------------------------------------
# Title: combine all the single celldata frames
# Author: Christina Meier
# Date: 2025-05-27
# Email: 19cm51@queensu.ca
#-------------------------------------------------
# Notes ------------------------------------------
# module load StdEnv/2023 r/4.5.0
#
# Options-----------------------------------------
message("Script started at: ", Sys.time())

# Packages ---------------------------------------
library(dplyr) # 1.1.4


# Input paths--------------------------------------
setwd("/global/project/hpcg1604/Christina_Meier/DRG_DEA/1_deconv")
rdsPath <- "1_singlecell/data/sc_cleaned"


# Load in data-------------------------------------
df_ids <- 1:4
counts_list <- list()

for (id in df_ids) {
  counts <- readRDS(sprintf(file.path(rdsPath, "naive_drg%d_good.rds"), id))
  print(paste("Loading in DRG", id, "counts"))
  counts <- as.data.frame(counts)
  rownames(counts) <- counts$V1
  counts <- counts |>
    select(-V1)
  colnames(counts) <- paste0("drg", id, "_", colnames(counts))
  counts_list[[id]] <- counts
  print(paste("DRG", id, "counts loaded and processed."))
}

# List of unique genes across all datasets
all_genes <- unique(unlist(lapply(counts_list, rownames)))

# Function to add missing genes and fill with zeros
fill_missing_genes <- function(df, all_genes) {
  #' @param df data frame of counts for one DRG dataset
  #' @param all_genes vector of all unique gene names across datasets
  #' @return data frame with all genes as rows, missing genes filled with zeros
  missing_genes <- setdiff(all_genes, rownames(df))
  if (length(missing_genes) > 0) {
    missing_df <- matrix(
      0,
      nrow = length(missing_genes),
      ncol = ncol(df),
      dimnames = list(missing_genes, colnames(df))
    )
    df <- rbind(df, missing_df)
  }

  df[all_genes, , drop = FALSE]
}

counts_list_filled <- lapply(counts_list, fill_missing_genes, all_genes = all_genes)
allscDRGcounts <- do.call(cbind, counts_list_filled)

# Ensure no duplicates
sum(duplicated(colnames(allscDRGcounts)))
colnames(allscDRGcounts)[duplicated(colnames(allscDRGcounts))]

# Save RDS ----------------------------------------
saveRDS(allscDRGcounts, file.path(rdsPath, "naive_allscDRGcounts.rds"))

message("Script finished at: ", Sys.time())

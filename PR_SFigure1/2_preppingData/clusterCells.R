#!/usr/bin/env Rscript
#-------------------------------------------------
# Title: Cluster cell states for DEA
# Author: Christina Meier
# Date: 2025-07-02-17
# Email: 19cm51@queensu.ca
#-------------------------------------------------
# Notes ------------------------------------------
# Prep for DEA by making cell type data frames
# module load StdEnv/2023 r/4.5.0
#
# Options ----------------------------------------
message("Script started at ", Sys.time())

# Packages ---------------------------------------
library(dplyr) # 1.3.1
library(reshape2) # 1.4.4

setwd("/global/project/hpcg1604/Christina_Meier/DRG_DEA")
source("0_helpers/CombineStates.R")

# Pathways ---------------------------------------
data_path <- "1_deconv/2_instaprism/data/naive_celltype_geneExpressionsDRG.rds"
output_path <- "2_preppingData/clusters"
list_path <- "2_preppingData"
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

# Load data ---------------------------------------
data <- readRDS(file.path(data_path))
drg_col_data <- read.csv(
  "/global/project/hpcg1604/Amanda_Zacharias/mouNaiveSNI/3_naiveSNI/2_dataPrep/gene/cleanData/d.coldata.csv"
)

# Clean coldata ------------------------------------
samples_to_keep <- c("X11", "X12", "X13", "X3", "X4", "X5")

coldat <- drg_col_data |>
  filter(naiveVsSni == "naive") |>
  filter(ztTime %in% c(2, 14)) |>
  dplyr::select(-any_of(c(
    "X.1", "X", "sampleGrps", "day", "tissue", "sex", "readType"
  )))

rownames(coldat) <- paste0("X", coldat$sampleNum)
stopifnot(identical(rownames(coldat), samples_to_keep))

# Define merged cell type groups ------------------
samples <- dimnames(data)[[1]]
genes <- dimnames(data)[[2]]
cell_types <- dimnames(data)[[3]]

merged_map <- list(
  Mrgprd = c("Mrgprd/Gm7271", "Mrgprd/Lpar3"),
  Mrgpra = c("Mrgpra3/Mrgprb4", "Mrgpra3"),
  S100b = c("S100b/Smr2", "S100b/Ntrk3/Gfra1", "S100b/Prokr2", "S100b/Baiap2l1", "S100b/Wnt7a"),
  Atf3 = c("Atf3/Mrgprd", "Atf3/Gfra3/Gal"),
  Zcchc12 = c("Zcchc12/Sstr2", "Zcchc12/Dcn", "Zcchc12/Trpm8", "Zcchc12/Rxfp1")
)

# Keep merged clusters + all remaining cell types--
merged_subtypes <- unlist(merged_map)
remaining_cell_types <- setdiff(cell_types, merged_subtypes)

single_map <- as.list(remaining_cell_types)
names(single_map) <- remaining_cell_types

final_map <- c(merged_map, single_map)

message("Groups to create:")
print(names(final_map))

# Create all matrices -------------------------
de_dfs <- lapply(names(final_map), function(group_name) {
  CombineStates(
    group_name = group_name,
    subtypes = final_map[[group_name]],
    data = data,
    samples_to_keep = samples_to_keep
  )
})

names(de_dfs) <- names(final_map)

# Remove any null results
de_dfs <- de_dfs[!sapply(de_dfs, is.null)]

# Save each matrix -------------------------------
for (group_name in names(de_dfs)) {
  safe_name <- gsub("[^A-Za-z0-9_]", "_", group_name)
  out_dir <- file.path(output_path, safe_name, "cleanData")

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  raw_counts <- de_dfs[[group_name]]
  # coldat_group <- coldat[colnames(raw_counts), , drop = FALSE]

  stopifnot(identical(rownames(coldat), colnames(raw_counts)))

  write.csv(
    de_dfs[[group_name]],
    file = file.path(out_dir, "rawCounts.csv"),
    row.names = TRUE
  )

  write.csv(
    coldat,
    file = file.path(out_dir, "coldata.csv"),
    row.names = TRUE
  )

  message("Saved raw counts and coldata for ", group_name, ": ", file.path(out_dir, "rawCounts.csv"))
}

# Save the list -------------------------------
saveRDS(de_dfs, file.path(list_path, "all_DE_ready_matrices.rds"))

message("Script ended at: ", Sys.time())

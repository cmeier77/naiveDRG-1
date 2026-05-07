#!/usr/bin/env Rscript
#-------------------------------------------------
# Title: Prep bulk data for deconvolution
# Author: Christina Meier
# Date: 2025-06-23
# Email: 19cm51@queensu.ca
#-------------------------------------------------
# Notes ------------------------------------------
#
# module load StdEnv/2023 r/4.5.0
# Options ----------------------------------------
message("Script started at ", Sys.time())

# Packages ---------------------------------------
library(dplyr) # 1.1.4

# Paths-------------------------------------------
setwd("/global/project/hpcg1604/Christina_Meier/DRG_DEA/1_deconv")
output_dir <- file.path("0_bulk", "data")
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# 1. Load raw DRG counts---------------------------
drg_rawCounts <- read.csv("/global/project/hpcg1604/Amanda_Zacharias/mouNaiveSNI/3_naiveDrg/2_dataPrep/gene/cleanData/rawCounts.csv")

# 2. Remove novel genes, 'MSTRG'--------------------
bulk_wo_mstrg_genes <- drg_rawCounts[grepl("^ENSMUS", drg_rawCounts$X), ]

# 3. Replace ENSMUSC with gene names, use id2name---
# load bulk id2name info
drg_id2name <- read.csv("/global/project/hpcg1604/Amanda_Zacharias/mouNaiveSNI/3_naiveDrg/2_dataPrep/id2name.csv")

# 4. Remove Mstrg genes-----------------------------
drg_gene_names_wo_mstrg_genes <- drg_id2name[grepl("^ENSMUS", drg_id2name$isoform_id), ]

# 5. Join by the X column (Ensembl ID's) and remove X column now with ensemble ID's
drg_wGeneNames <- bulk_wo_mstrg_genes |>
    left_join(drg_gene_names_wo_mstrg_genes, by = c("X" = "gene_id")) |>
    dplyr::select(gene_name, everything()) |>
    dplyr::select(-isoform_id, -X, -X.y)

# 6. Keep the first occurrence of each gene name-----
drg_wGeneNames_asRowNames <- drg_wGeneNames |>
    distinct(gene_name, .keep_all = TRUE)

rownames(drg_wGeneNames_asRowNames) <- drg_wGeneNames_asRowNames$gene_name

drg_wGeneNames_asRowNames <- drg_wGeneNames_asRowNames |>
    select(-gene_name)

# 7. Save rds-----------------------------------------
saveRDS(drg_wGeneNames_asRowNames, file.path(output_dir, "naiveDrg_wGeneNames_wRowNames.rds"))

message("Script ended at ", Sys.time())

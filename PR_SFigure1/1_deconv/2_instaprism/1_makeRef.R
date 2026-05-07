#!/usr/bin/env Rscript
#-------------------------------------------------
# Title: prepare reference sequence for DRG deconvolution
# Author: Christina Meier
# Date: 2025-06-05
# Email: 19cm51@queensu.ca
#-------------------------------------------------
# Notes ------------------------------------------
# module load StdEnv/2023 r/4.5.0
#
# Options ----------------------------------------
message("Script started at ", Sys.time())

# Packages ---------------------------------------
library(InstaPrism) # 0.1.6
library(Biobase) # 2.64.0
library(dplyr) # 1.1.4
library(ggplot2) # 3.5.2

# Input Paths-------------------------------------
setwd("/global/project/hpcg1604/Christina_Meier/DRG_DEA/1_deconv")

output_path <- "2_instaprism"
if (!dir.exists(output_path)) {
    dir.create(output_path)
}

data_path <- file.path(output_path, "data")
if (!dir.exists(data_path)) {
    dir.create(data_path)
}

# 1. Load in data---------------------------------
bulkDRG <- readRDS("./0_bulk/data/naiveDrg_wGeneNames_wRowNames.rds")
scDRG <- readRDS("./1_singlecell/data/sc_cleaned/naive_allscDRGcounts.rds")
scDRGmeta <- readRDS("./1_singlecell/data/sc_cleaned/allscDRGmeta.rds")

# 2. Extract info for reference------------------
cell_type_labels <- scDRGmeta$cell_type
cell_state_labels <- scDRGmeta$cell_state

ncol(scDRG) == nrow(scDRGmeta)
length(cell_type_labels) == ncol(scDRG) 
length(cell_state_labels) == ncol(scDRG)

cell.state.labels <- scDRGmeta[colnames(scDRG), "cell_state"]
cell.type.labels <- scDRGmeta[colnames(scDRG), "cell_type"]

# 3. Prepare reference with refprepare function----
# sc_Expr has to have only numeric arguments, no gene name columns
refPhi_obj <- refPrepare(sc_Expr = scDRG, cell.type.labels = cell.type.labels, cell.state.labels = cell.state.labels)
saveRDS(refPhi_obj, file.path(data_path, "refPhi_objDRG.rds"))

# 4. Plot gene overlap----------------------------
bulk <- as.data.frame(bulkDRG)
bulk_genes <- rownames(bulk)

ref_genes <- rownames(refPhi_obj@phi.cs)
common_genes_clean <- length(intersect(ref_genes, bulk_genes))
unique_sc_clean <- length(setdiff(ref_genes, bulk_genes))
unique_bulk_clean <- length(setdiff(bulk_genes, ref_genes))

gene_summary_clean <- data.frame(
  Category = c("Shared", "Single-Cell Only", "Bulk Only"),
  Count = c(common_genes_clean, unique_sc_clean, unique_bulk_clean)
)

plot_clean <- ggplot(gene_summary_clean, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 0.5) +
  labs(
    title = "Gene Distribution Between SC and Bulk DRG",
    y = "Number of Genes",
    x = "Category"
  ) +
  scale_fill_manual(values = c("Shared" = "blue", "Single-Cell Only" = "pink", "Bulk Only" = "forestgreen")) +
  ylim(0, 30000) +
  theme_minimal()

# Save plot-----------------------------------------
ggsave(file.path(data_path, "geneOverlapDRG.pdf"), plot = plot_clean, width = 8, height = 6)

message("Script ended at: ", Sys.time())

#!/usr/bin/env Rscript
#-------------------------------------------------
# Title: run DRG deconvolution with Instaprism 
# Author: Christina Meier
# Date: 2025-06-10
# Email: 19cm51@queensu.ca
#-------------------------------------------------
# Notes ------------------------------------------
# module load StdEnv/2023 r/4.5.0
#
# Options ----------------------------------------
message("Script started at ", Sys.time())

# Packages ----------------------------------------
library(InstaPrism) # 0.1.6
library(Biobase) # 2.64.0
library(dplyr) # 1.3.1

# Pathways-----------------------------------------
setwd("/global/project/hpcg1604/Christina_Meier/DRG_DEA/1_deconv")
bulk_path <- "0_bulk/data/naiveDrg_wGeneNames_wRowNames.rds"
sc_path <- "1_singlecell/data/sc_cleaned/naive_allscDRGcounts.rds"
ref_path <- "2_instaprism/data/refPhi_objDRG.rds"
output_path <- "2_instaprism/data"

# Load input data----------------------------------
bulkDRG <- readRDS(file.path(bulk_path))
scDRG <- readRDS(file.path(sc_path))
refDRG <- readRDS(file.path(ref_path))

# 1. Run InstaPrism, visualize convergence
InstaPrism.res <- InstaPrism(
  bulk_Expr = bulkDRG, refPhi_cs = refDRG,
  verbose = T, convergence.plot = T, max_n_per_plot = 100
)

# The deconvolved cell type fraction is accessible with:
estimated_frac <- t(InstaPrism.res@Post.ini.ct@theta)
saveRDS(estimated_frac, file.path(output_path, "naive_deconv_results_DRG.rds"))

# The sample by gene by cell type array contains the cell-type specific gene expression deconvolution results and can be accessed with:
Z <- get_Z_array(InstaPrism.res)
saveRDS(Z, file.path(output_path, "naive_celltype_geneExpressionsDRG.rds"))

message("Script ended at ", Sys.time())

#!/usr/bin/env Rscript
#-------------------------------------------------
# Title: combine all the single data metadatas
# Author: Christina Meier
# Date: 2025-06-02
# Email: 19cm51@queensu.ca
#-------------------------------------------------
# Notes ------------------------------------------
# module load StdEnv/2023 r/4.5.0
#
# Options-----------------------------------------
message("Script started at: ", Sys.time())


# Packages ----------------------------------------
library(dplyr) # 1.1.4

# Input paths--------------------------------------
setwd("/global/project/hpcg1604/Christina_Meier/DRG_DEA/1_deconv")
rdsPath <- "1_singlecell/data/sc_cleaned/"

df_ids <- 1:4
meta_list <- list()

# Load in data-------------------------------------
for (id in df_ids) {
    meta <- readRDS(file.path(rdsPath, sprintf("naive_drg%d_mt_good.rds", id)))
    print(paste("Loading in DRG", id, "metadata"))
    meta <- as.data.frame(meta)
    meta$V1 <- paste0("drg", id, "_", meta$V1)
    rownames(meta) <- meta$V1
    meta_list[[id]] <- meta
    print(paste("DRG", id, "metadata loaded."))
}

sc_counts <- readRDS(file.path(rdsPath, "naive_allscDRGcounts.rds"))

unique(meta_list[[1]]$Celltype)
unique(meta_list[[2]]$celltype)
unique(meta_list[[3]]$celltype)
unique(meta_list[[4]]$celltype)

# Clean up drg1------------------------------------
meta_list[[1]]$Celltype[meta_list[[1]]$Celltype == "Immune"] <- "Immune_Cell"
meta_list[[1]]$Celltype[meta_list[[1]]$Celltype == "Schwann"] <- "Schwann_Cell"
meta_list[[1]]$Celltype[meta_list[[1]]$Celltype == "Satellite"] <- "Satellite_Cell"
meta_list[[1]]$Celltype[meta_list[[1]]$Celltype %in% c("VEC", "VECC")] <- "VECC"

# Clean up drg2------------------------------------
meta_list[[2]]$celltype[meta_list[[2]]$celltype == "Immune Cell"] <- "Immune_Cell"
meta_list[[2]]$celltype[meta_list[[2]]$celltype == "Schwann Cell"] <- "Schwann_Cell"
meta_list[[2]]$celltype[meta_list[[2]]$celltype == "Satellite Cell"] <- "Satellite_Cell"
meta_list[[2]]$celltype[meta_list[[2]]$celltype %in% c("VEC", "VECC")] <- "VECC"
meta_list[[2]]$celltype[meta_list[[2]]$celltype == "Red Blood Cell"] <- "Red_blood_cell"

# Clean up drg3------------------------------------
meta_list[[3]]$celltype[meta_list[[3]]$celltype == "Immune"] <- "Immune_Cell"
meta_list[[3]]$celltype[meta_list[[3]]$celltype == "Schwann"] <- "Schwann_Cell"
meta_list[[3]]$celltype[meta_list[[3]]$celltype == "Satellite"] <- "Satellite_Cell"
meta_list[[3]]$celltype[meta_list[[3]]$celltype %in% c("VEC", "VECC")] <- "VECC"
meta_list[[3]]$celltype[meta_list[[3]]$celltype == "RBC"] <- "Red_blood_cell"

# Clean up drg4------------------------------------
meta_list[[4]]$celltype[meta_list[[4]]$celltype == "Immune"] <- "Immune_Cell"
meta_list[[4]]$celltype[meta_list[[4]]$celltype == "Schwann"] <- "Schwann_Cell"
meta_list[[4]]$celltype[meta_list[[4]]$celltype == "Satellite"] <- "Satellite_Cell"
meta_list[[4]]$celltype[meta_list[[4]]$celltype %in% c("VEC", "VECC")] <- "VECC"
meta_list[[4]]$celltype[meta_list[[4]]$celltype == "RBC"] <- "Red_blood_cell"

# cell type information-----------------------------
celltype_1 <- data.frame(cell_type = meta_list[[1]]$Celltype, row.names = meta_list[[1]]$V1)
celltype_2 <- data.frame(cell_type = meta_list[[2]]$celltype, row.names = meta_list[[2]]$V1)
celltype_3 <- data.frame(cell_type = meta_list[[3]]$celltype, row.names = meta_list[[3]]$V1)
celltype_4 <- data.frame(cell_type = meta_list[[4]]$celltype, row.names = meta_list[[4]]$V1)

unique(meta_list[[1]]$Celltype)
unique(meta_list[[2]]$celltype)
unique(meta_list[[3]]$celltype)
unique(meta_list[[4]]$celltype)

# cell state information----------------------------
cellstate_1 <- data.frame(cell_state = meta_list[[1]]$Celltype, row.names = meta_list[[1]]$V1)
cellstate_2 <- data.frame(cell_state = meta_list[[2]]$celltype, row.names = meta_list[[2]]$V1)
cellstate_3 <- data.frame(cell_state = meta_list[[3]]$celltype_2, row.names = meta_list[[3]]$V1)
cellstate_4 <- data.frame(cell_state = meta_list[[4]]$celltype, row.names = meta_list[[4]]$V1)

# Combine-------------------------------------------
scDRGcelltypes <- do.call(rbind, list(celltype_1, celltype_2, celltype_3, celltype_4))
scDRGcellstates <- do.call(rbind, list(cellstate_1, cellstate_2, cellstate_3, cellstate_4))

finalDRGmeta <- cbind(scDRGcelltypes, scDRGcellstates)
unique(finalDRGmeta$cell_type)
unique(finalDRGmeta$cell_state)

# Make cell states consistent too--------------------
finalDRGmeta$cell_state[finalDRGmeta$cell_state == "Immune"] <- "Immune_Cell"
finalDRGmeta$cell_state[finalDRGmeta$cell_state == "VEC"] <- "VECC"
finalDRGmeta$cell_state[finalDRGmeta$cell_state == "Schwann"] <- "Schwann_Cell"
finalDRGmeta$cell_state[finalDRGmeta$cell_state == "Satellite"] <- "Satellite_Cell"
finalDRGmeta$cell_state[finalDRGmeta$cell_state == "RBC"] <- "Red_blood_cell"

# Check matches with counts--------------------------
rownames(finalDRGmeta) <- rownames(scDRGcelltypes)
finalscDRGmeta <- finalDRGmeta[
    colnames(sc_counts), ,
    drop = FALSE
]
stopifnot(identical(colnames(sc_counts), rownames(finalscDRGmeta)))

# Save RDS -------------------------------------------
saveRDS(finalscDRGmeta, file.path(rdsPath, "allscDRGmeta.rds"))

message("Script ended at: ", Sys.time())

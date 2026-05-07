#!/usr/bin/env Rscript
#-------------------------------------------------
# Title: Download single cell DRG data
# Author: Christina Meier
# Date: 2025-05-22
# Email: 19cm51@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------
# Load in data downloaded from GEO Series GSE155622.
# This contains single cell RNA-seq data from mouse dorsal root ganglia (DRG) after sciatic nerve injury (SNI) and sham surgery.
# The data includes raw UMI counts for 4 samples (2 SNI and 2 sham) and associated metadata.
#
# module load StdEnv/2023 r/4.5.0
#

# Options -----------------------------------------
message("Script started at ", Sys.time())

# Packages -----------------------------------------
library(data.table) # 1.18.2.1

gse155622_dir <- file.path("data", "raw", "GSE155622")

read_drg_counts <- function(i, data_dir = gse155622_dir) {
    stopifnot(i %in% 1:4)

    fread(
        file.path(data_dir, sprintf("GSE155622_raw_UMI_counts_%d.txt.gz", i)),
        sep = "\t"
    )
}

read_drg_metadata <- function(i, data_dir = gse155622_dir) {
    stopifnot(i %in% 1:4)

    fread(
        file.path(data_dir, sprintf("GSE155622_raw_UMI_counts_%d_metadata.txt.gz", i)),
        sep = "\t"
    )
}

message("Script ended at ", Sys.time())

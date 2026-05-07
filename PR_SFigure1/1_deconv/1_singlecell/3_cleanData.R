#!/usr/bin/env Rscript
#-------------------------------------------------
# Title: Clean sc DRG data; only keep naive samples
# Author: Christina Meier
# Date: 2025-05-25
# Email: 19cm51@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------
# module load StdEnv/2023 r/4.5.0
#
# Options -----------------------------------------
message("Script started at ", Sys.time())

# Packages -----------------------------------------
library(dplyr) # 1.1.4

# Load helper function -----------------------------
source("../../0_helpers/CleanData.R")

df_ids <- 1:4

for (id in df_ids) {
    message("cleaning df ", id, " now")
    cleanData(
        df_id = id,
        input_dir = "data/raw/GSE155622",
        output_dir = "data/sc_cleaned"
    )
    message("done cleaning df ", id)
}

message("Script ended at ", Sys.time())

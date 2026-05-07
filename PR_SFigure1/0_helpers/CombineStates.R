#!/usr/bin/env Rscript
#-------------------------------------------------
# Title: Combine cell type states
# Author: Christina Meier
# Date: 2025-07-15
# Email: 19cm51@queensu.ca
#-------------------------------------------------
# Notes ------------------------------------------
# module load StdEnv/2023 r/4.5.0
#
#
#
# Options ----------------------------------------

# Packages ---------------------------------------
library(dplyr) # 1.1.4
library(reshape2) # 1.4.4

# Function ---------------------------------------
CombineStates <- function(group_name, subtypes, data, samples_to_keep) {
    #' @param group_name: Name of the group (e.g., "Mrgprd")
    #' @param subtypes: Vector of subtypes to merge (e.g., c("Mrgprd/Gm7271", "Mrgprd/Lpar3"))
    #' @param data: The original 3D array of gene expression data
    #' @param samples_to_keep: Vector of sample names to keep (e.g., c("X11", "X12", "X13", "X3", "X4", "X5"))
    #' @return: A data frame with genes as rows and samples as columns, with subtypes merged

    subtypes_present <- subtypes[subtypes %in% dimnames(data)[[3]]]

    if (length(subtypes_present) == 0) {
        warning("No matching subtypes found for: ", group_name)
        return(NULL)
    }

    message("Processing ", group_name, ": ", paste(subtypes_present, collapse = ", "))

    subtype_matrices <- lapply(subtypes_present, function(subtype) {
        mat <- t(data[samples_to_keep, , subtype])
        colnames(mat) <- samples_to_keep
        mat
    })

    combined_mat <- do.call(cbind, subtype_matrices)

    # Sum duplicate sample columns after combining subtypes
    combined_t <- t(combined_mat)
    combined_t <- rowsum(combined_t, group = rownames(combined_t))

    final_mat <- as.data.frame(t(combined_t))
    rownames(final_mat) <- dimnames(data)[[2]]

    return(final_mat)
}

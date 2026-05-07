#!/usr/bin/env Rscript
#-------------------------------------------------
# Title: Write individual bash scripts
# Author: Amanda Zacharias & Christina Meier
# Date: 2023-08-11 & 2025-08-11
# Email: 19cm51@queensu.ca
#-------------------------------------------------
# Notes ------------------------------------------
# Run this script from working directory
# module load StdEnv/2023 gcc/12.3 r/4.5.0
#
#
#
# Options ----------------------------------------
message("Script started at ", Sys.time())

# Packages ---------------------------------------
library(dplyr) # 1.1.4

# Pathways ----------------------------------------
# Input ===========
projectName <- "3_naiveDrg"

cellDir <- "2_preppingData/clusters"
cellPaths <- list.files(cellDir, full.names = TRUE)
cellType <- gsub(".csv", "", basename(cellPaths))

rscriptNames <- list("2vs14" = "0_2vs14.R")

baseDir <- file.path(getwd(), projectName, "3_deseqCandidate")
candidatesDir <- file.path(baseDir, "candidates", "cleanCandidates")
candPaths <- list.files(candidatesDir, full.names = TRUE)
names(candPaths) <- gsub(".csv", "", basename(candPaths))

# Output ===========
for (cell in cellType) {
  bashDir <- file.path(baseDir, "bash")
  if (!dir.exists(bashDir)) {
    dir.create(bashDir, recursive = TRUE)
  }
}


# Load data -----------------------------------------
baseScript <- as.character(
  read.table(file.path(baseDir, "1_baseScript.sh"),
    sep = "\n", blank.lines.skip = FALSE,
    comment.char = "", quote = "\'",
    stringsAsFactors = FALSE
  )$V1
)

# Modify lines -----------------------------------------
# Functions ============
ModifyScript <- function(name, rscriptPath, candPath, projName, cellType, thresh, script) {
  #' Modify the base script
  #'
  #' @param name string; what to call the SLURM job & files
  #' @param rscriptPath Path to Rscript that will be run (string)
  #' @param candPath Path to candidate genes txt file (string)
  #' @param projName string; name of project folder
  #' @param cellType string; name of project sub-folder
  #' @param thresh double; threshold for non-specific filtering
  #' @param script unmodified script (vector)
  #' @return a modified base script (vector)
  #' @example
  #'
  # Header
  script[2] <- paste(script[2], paste0(name), sep = "")
  script[9] <- paste(script[9], paste0(name, ".out"), sep = "")
  script[10] <- paste(script[10], paste0(name, ".err"), sep = "")
  # Content
  script[25] <- paste(script[25], rscriptPath, sep = "")
  script[26] <- paste(script[26], name, sep = "")
  script[27] <- paste(script[27], projName, sep = "")
  script[28] <- paste(script[28], cellType, sep = "")
  script[29] <- paste(script[29], thresh, sep = "")
  script[30] <- paste(script[30], candPath, sep = "")
  script[31] <- paste(script[31], getwd(), sep = "")
  # Return
  return(script)
}

ProcessInfo <- function(cellPaths, rscriptNames) {
  #' Prepare information for modifying files,
  #'  so don't have to copy and paste code
  #' @param coldat Dataframe with sample metadata
  #' @return NA, writes bash scripts and messages to console
  #' @example
  #'
  # Loop through samples
  cellNames <- gsub(".csv", "", basename(cellPaths))

  for (i in seq_along(cellPaths)) {
    cellName <- cellNames[i]
    candPath <- candPaths
    cat("\n", cellName)

    for (rName in names(rscriptNames)) {
      cat("\n\t", rName)
      sampleLines <- ModifyScript(
        rscriptPath = rscriptNames[[rName]],
        candPath = candPath,
        name = paste(cellName, rName, sep = "."),
        projName = projectName,
        cellType = cellName,
        script = baseScript,
        thresh = 0.015571
      )
      # Save
      outFilename <- paste(cellName, rName, "sh", sep = ".")
      fileConn <- file(file.path(baseDir, "bash", outFilename), open = "w")
      writeLines(sampleLines, fileConn)
      close(fileConn)
    } # end loop through rNames
  } # end loop through candNames
} # end function

# Execute -----------------------------------------
ProcessInfo(cellPaths, rscriptNames)

# To Run
toRun <- tidyr::crossing(cellType, names(rscriptNames))
toRunStr <- paste(toRun$cellType, toRun$`names(rscriptNames)`, "sh", sep = ".")
toRunLines <- paste("sbatch", toRunStr)

fileConn <- file(file.path(baseDir, "jobsToRun.sh"))
writeLines(toRunLines, fileConn)
close(fileConn)

message("Script ended at ", Sys.time())

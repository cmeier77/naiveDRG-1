#!/usr/bin/env Rscript
#-------------------------------------------------
# Title: edgeR
# Author: Amanda Zacharias & CM
# Date: 2025-07-17
# Email: 19cm51@queensu.ca
#-------------------------------------------------
# Notes ------------------------------------------
# module load StdEnv/2023 r/4.5.0
#
#
#
# Options -----------------------------------------
message("Script started at ", Sys.time())
options(ggrepel.max.overlaps = Inf, box.padding = 0)

# Packages ----------------------------------------
library(optparse) # 1.7.5
library(dplyr) # 1.1.4
library(tibble) # 3.3.0
library(DESeq2) # 1.46.0
library(ggplot2) # 3.5.2
library(ggrepel) # 0.9.6

# Run optparser ------------------------------------
# Define arguments ==================
optionList <- list(
  make_option(c("-w", "--workingDir"),
    type = "character", default = getwd(),
    help = "the project working directory path", metavar = "character"
  ),
  make_option(c("-p", "--projectName"),
    type = "character", default = NA,
    help = "name of folder for analysis", metavar = "character"
  ),
  make_option(c("-x", "--cellType"),
    type = "character", default = NA,
    help = "name of the cell type for analysis", metavar = "character"
  ),
  make_option(c("-c", "--countPath"),
    type = "character", default = NA,
    help = "path to your count matrix", metavar = "character"
  ),
  make_option(c("-n", "--basename"),
    type = "character", default = NA,
    help = "the base filename for this analysis' outputs", metavar = "character"
  ),
  make_option(c("-t", "--threshold"),
    type = "double", default = NA,
    help = "non-specific filtering MAD threshold", metavar = "double"
  ),
  make_option(c("-a", "--candPath"),
    type = "character", default = NA,
    help = "path to your candidates list", metavar = "character"
  )
)
# Get parameters ==================
optParser <- OptionParser(option_list = optionList)
opt <- parse_args(optParser)

# There is no parameter checking!!!!!
cat(
  "\nBase filename:", opt$basename,
  "\nProject name:", opt$projectName,
  "\nPrefix:", opt$cellType,
  "\nThreshold for filtering:", opt$threshold,
  "\n"
)

# Set the working directory -----------------------------------------
setwd(opt$workingDir)

# Helper functions -----------------------------------------
source("0_helpers/MakeVolcano.R")

# Pathways -----------------------------------------
# Input ===========
coldataPath <- file.path("2_preppingData", "clusters", opt$cellType, "cleanData", "coldata.csv")

# Output ===========
baseDir <- file.path(opt$projectName, "3_deseqCandidate", "gene")
dfsDir <- file.path(baseDir, "dataframes")
baseDfsDir <- file.path(dfsDir, opt$basename)
rDataDir <- file.path(baseDir, "rData")
plotsDir <- file.path(baseDir, "plots")
basePlotsDir <- file.path(plotsDir, opt$basename)

system(paste(
  "mkdir", baseDir, dfsDir, baseDfsDir, rDataDir,
  plotsDir, basePlotsDir
))

# Load data -----------------------------------------
coldata <- read.csv(coldataPath, row.names = 1, stringsAsFactors = FALSE)
counts <- read.csv(opt$countPath, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
candidates <- read.csv(opt$candPath, row.names = 1, stringsAsFactors = FALSE)

# Execute DeSeq2 ########################################
# Normalize -----------------------------------------
cat("\n", rep("-", 10), "Normalizing", rep("-", 10), "\n")
coldata$sampleGrps <- paste0("naivezt", coldata$ztTime, "d")
dds <- DESeqDataSetFromMatrix(round(data.matrix(counts)), coldata, ~sampleGrps)
ddsSE <- estimateSizeFactors(dds)

# Dispersion estimation -----------------------------------------
cat("\n", rep("-", 10), "Dispersion Estimation", rep("-", 10), "\n")
ddsDisp <- estimateDispersions(ddsSE)

pdf(file.path(plotsDir, "dispersion.pdf"),
  width = 6, height = 6
)
plotDispEsts(ddsDisp,
  ylab = "Dispersion", xlab = "Mean of normalized counts",
  main = opt$basename
)
dev.off()

# Extract only candidates -----------------------------------------
cat("\n", rep("-", 10), "Extracting Candidates", rep("-", 10), "\n")
idCol <- "isoform_id"
if (opt$cellType == "gene") {
  idCol <- "gene_id"
}

ddsCand <- ddsDisp[rownames(ddsDisp) %in% rownames(candidates), ]

candGenesIncluded <- rownames(ddsCand)
candGenesNotIncluded <- setdiff(rownames(candidates), candGenesIncluded)

cat(
  "\nNumber of candidates:", nrow(candidates),
  # "\nNumber of features after subset:", length(rownames(ddsCand)),
  "\nNumber of features after subset: ", length(candGenesIncluded),
  # "\nNumber of unique features after subset:", length(unique(rownames(ddsCand))),
  # "\nCandidate genes included:", candGenesIncluded,
  "\nCandidate genes not included:", candGenesNotIncluded,
  # "\nCandidate genes not included at all:", setdiff(candGenesNotIncluded, candGenesIncluded),
  "\n"
)

# Differential expression analysis -----------------------------------------
cat("\n", rep("-", 10), "DE Test", rep("-", 10), "\n")
ddsWald <- nbinomWaldTest(ddsCand)
cont <- c("sampleGrps", "naivezt14d", "naivezt2d")
coefs <- "sampleGrps_naivezt14d_vs_naivezt2d"
res <- results(
  ddsWald,
  contrast = cont,
  independentFilter = FALSE,
  pAdjustMethod = "bonferroni",
  alpha = 0.05
)
summary(res)

# Format results, no id2name -----------------------------------------
cat("\n", rep("-", 10), "Formatting Results", rep("-", 10), "\n")
resDf <- data.frame(res@listData, row.names = res@rownames) %>%
  arrange(padj) %>%
  tibble::rownames_to_column("gene_id")
sigResDf <- resDf %>% filter(padj < 0.05)

cat("\nNumber of Sig:", nrow(sigResDf), "Out of:", nrow(resDf), "\n")

# Save results -----------------------------------------
cat("\n", rep("-", 10), "Saving Dataframes", rep("-", 10), "\n")
write.csv(resDf, file.path(baseDfsDir, "resDf.csv"))
write.csv(sigResDf, file.path(baseDfsDir, "sigResDf.csv"))

# Plot results -----------------------------------------
cat("\n", rep("-", 10), "Plotting", rep("-", 10), "\n")
# Volcano plot ========
MakeVolcano(
  df = resDf,
  showName = TRUE,
  newTitle = opt$basename,
  newPath = basePlotsDir,
  idColumn = "gene_id"
)
MakeVolcano(
  df = resDf,
  showName = FALSE,
  newTitle = opt$basename,
  newPath = basePlotsDir,
  idColumn = "gene_id"
)

# Save image -----------------------------------------
cat("\n", rep("-", 10), "Saving RData", rep("-", 10), "\n")
save.image(file.path(rDataDir, paste(opt$basename, "RData", sep = ".")))

message("Script ended at ", Sys.time())

#!/bin/bash
#SBATCH --job-name=cleanData
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=19cm51@queensu.ca
#SBATCH --qos=privileged # or SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=10GB  # Job memory request
#SBATCH --time=0-1:00:00  # Day-Hours-Minutes-Seconds
#SBATCH --output=cleanData.out
#SBATCH --error=cleanData.err
# Title: Clean scRNA-seq reads from GEO Series GSE155622.
# Author: Christina Meier
# Date: 2025-05-25
# Email: 19cm51@queensu.ca
#-------------------------------------------------
# Notes ------------------------------------------
#
echo Job started at $(date +%T)

# Dependencies------------------------------------
module load StdEnv/2023 r/4.5.0

mkdir -p data/sc_cleaned

Rscript 3_cleanData.R

echo Job finished at $(date +%T)
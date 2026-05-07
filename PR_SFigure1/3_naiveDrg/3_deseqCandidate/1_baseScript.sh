#!/bin/bash
#SBATCH --job-name=dea
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=19cm51@queensu.ca
#SBATCH --qos=privileged # or SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=5GB  # Job memory request
#SBATCH --time=0-5:00:00  # Day-Hours-Minutes-Seconds
#SBATCH --output=
#SBATCH --error=

# Title: Run edgeR
# Author: Amanda Zacharias
# Date: 2025-07-08
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------

echo Job started at $(date +'%T')

# Load dependencies
module load StdEnv/2023 gcc/12.3 r/4.5.0 

# Variables
RSCRIPTNAME=
BASENAME=
PROJNAME=
PREFIX=
THRESH=
CANDPATH=
CWD=

# Execute R script
Rscript ../${RSCRIPTNAME} \
--basename $BASENAME \
--projectName $PROJNAME \
--cellType $PREFIX \
--threshold $THRESH \
--candPath $CANDPATH \
--countPath ${CWD}/2_preppingData/clusters/${PREFIX}/cleanData/rawCounts.csv \
--workingDir $CWD
  
echo Job ended at $(date +'%T')

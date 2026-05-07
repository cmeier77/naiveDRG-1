#!/bin/bash
#SBATCH --job-name=deaSchwann_Cell.2vs14
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=19cm51@queensu.ca
#SBATCH --qos=privileged # or SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=5GB  # Job memory request
#SBATCH --time=0-5:00:00  # Day-Hours-Minutes-Seconds
#SBATCH --output=Schwann_Cell.2vs14.out
#SBATCH --error=Schwann_Cell.2vs14.err

# Title: Run edgeR
# Author: Amanda Zacharias
# Date: 2025-07-08
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------

echo Job started at $(date +%T)

# Load dependencies
module load StdEnv/2023 gcc/12.3 r/4.5.0 

# Variables
RSCRIPTNAME=0_2vs14.R
BASENAME=Schwann_Cell.2vs14
PROJNAME=3_naiveDrg
PREFIX=Schwann_Cell
THRESH=0.015571
CANDPATH=/global/project/hpcg1604/Christina_Meier/DRG_DEA/3_naiveDrg/3_deseqCandidate/candidates/cleanCandidates/brmmu04040.csv
CWD=/global/project/hpcg1604/Christina_Meier/DRG_DEA

# Execute R script
Rscript ../${RSCRIPTNAME} \
--basename $BASENAME \
--projectName $PROJNAME \
--cellType $PREFIX \
--threshold $THRESH \
--candPath $CANDPATH \
--countPath ${CWD}/2_preppingData/clusters/${PREFIX}/cleanData/rawCounts.csv \
--workingDir $CWD
  
echo Job ended at $(date +%T)

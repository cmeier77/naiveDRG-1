#!/bin/bash
#SBATCH --job-name=downloadReads
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=19cm51@queensu.ca
#SBATCH --qos=privileged # or SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=2GB  # Job memory request
#SBATCH --time=0-1:00:00  # Day-Hours-Minutes-Seconds
#SBATCH --output=downloadReads.out
#SBATCH --error=downloadReads.err
# Title: Download scRNA-seq reads from GEO Series GSE155622.
# Author: Christina Meier
# Date: 2025-05-22
# Email: 19cm51@queensu.ca
#-------------------------------------------------
# Notes ------------------------------------------
#
echo Job started at $(date +%T)

# Dependencies------------------------------------
module load StdEnv/2023

# Code ------------------------------------------
set -euo pipefail

mkdir -p data/raw/GSE155622

while read -r url; do
  file="data/raw/GSE155622/$(basename "$url")"

  if [[ ! -f "$file" ]]; then
    echo "Downloading: $url"
    wget -O "$file" "$url"
  else
    echo "Already exists: $file"
  fi

  echo "Unzipping: $file"
    gunzip "$file"


done < metadata/urls.txt

echo Job finished at $(date +%T)
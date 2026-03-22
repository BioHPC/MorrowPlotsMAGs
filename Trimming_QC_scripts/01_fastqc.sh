#!/usr/bin/env bash
set -euo pipefail

# Paths
RAW="/lab/biohpc/Metagenomics_Huang/new_data/Project_Huang_33_gDNA_merged"
WORK="/lab/biohpc/Metagenomics_Huang/new_data/work_new"
THREADS=32
LOG="$WORK/logs"

# Create directories
mkdir -p "$LOG" "$WORK/qc" "$WORK/trim"

fastqc -t "$THREADS" \
    -o "$WORK/qc" \
    "$RAW"/*_R1_*.fastq.gz "$RAW"/*_R2_*.fastq.gz

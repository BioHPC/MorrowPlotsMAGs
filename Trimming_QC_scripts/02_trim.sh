#!/usr/bin/env bash
set -euo pipefail

RAW="/lab/biohpc/Metagenomics_Huang/new_data/Project_Huang_33_gDNA_merged"
WORK="/lab/biohpc/Metagenomics_Huang/new_data/work_new"
TRIM="$WORK/trim"
QC="$WORK/qc"
THREADS=32

# Adapter file — set TRIMMOMATIC_ADAPT or edit this path to match your Trimmomatic installation
ADAPT="${TRIMMOMATIC_ADAPT:-$(dirname "$(command -v trimmomatic)" 2>/dev/null)/../share/trimmomatic/adapters/TruSeq3-PE-2.fa}"
if [[ ! -f "$ADAPT" ]]; then
  echo "ERROR: Adapter file not found at '$ADAPT'." >&2
  echo "Set TRIMMOMATIC_ADAPT to the full path of TruSeq3-PE-2.fa" >&2
  exit 1
fi

mkdir -p "$TRIM" "$QC"

# Trim reads
for R1 in "$RAW"/*_R1_001.fastq.gz; do
    base=$(basename "$R1" | sed 's/_R1_001\.fastq\.gz//')
    R2="$RAW/${base}_R2_001.fastq.gz"

    out1P="$TRIM/${base}_1P.fq.gz"
    out1U="$TRIM/${base}_1U.fq.gz"
    out2P="$TRIM/${base}_2P.fq.gz"
    out2U="$TRIM/${base}_2U.fq.gz"

    trimmomatic PE -threads "$THREADS" -phred33 \
        "$R1" "$R2" \
        "$out1P" "$out1U" \
        "$out2P" "$out2U" \
        ILLUMINACLIP:"$ADAPT":2:30:10 SLIDINGWINDOW:4:30 MINLEN:70
done

# FastQC
fastqc -t "$THREADS" -o "$QC" "$TRIM"/*_1P.fq.gz "$TRIM"/*_2P.fq.gz

# MultiQC
multiqc -q -o "$QC" "$QC"

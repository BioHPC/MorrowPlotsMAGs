#!/usr/bin/env bash
set -euo pipefail

BASE="/lab/biohpc/Metagenomics_Huang/new_data_work"
MASTER="$BASE/group_samples_master_newdata.tsv"

DB="/lab/biohpc/Metagenomics_Huang/db/k2_standard_20230605_db"
THREADS=16
CONF=0.05
MINHIT=2

OUTSUB="kraken2_reads_group"

FS_RE='[[:space:]]+'

mapfile -t GRP_LIST < <(
  awk -v FS="$FS_RE" 'NR>1{print $1}' "$MASTER" | sort -u
)

for grp in "${GRP_LIST[@]}"; do
  OUTDIR="$BASE/$grp/$OUTSUB"
  mkdir -p "$OUTDIR"

  report="$OUTDIR/${grp}_group_reads.k2.report"
  outcls="$OUTDIR/${grp}_group_reads.k2.classified"

  ARGS=()

  while read -r r1 r2; do
    ARGS+=("$r1" "$r2")
  done < <(
    awk -v FS="$FS_RE" -v G="$grp" 'NR>1 && $1==G {print $3"\t"$4}' "$MASTER"
  )

  kraken2 \
    --db "$DB" \
    --threads "$THREADS" \
    --use-names \
    --gzip-compressed \
    --paired \
    --minimum-hit-groups "$MINHIT" \
    --confidence "$CONF" \
    --report "$report" \
    --output "$outcls" \
    "${ARGS[@]}"
done
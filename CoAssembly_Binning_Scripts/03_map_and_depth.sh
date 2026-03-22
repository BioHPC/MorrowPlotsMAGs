#!/usr/bin/env bash
set -euo pipefail

BASE="/lab/biohpc/Metagenomics_Huang/CoAssemble"
LIST_TSV="${1:-$BASE/group_samples_master.tsv}"

THREADS_TOTAL=48
GROUP_JOBS=6
MAP_THREADS=$(( THREADS_TOTAL / GROUP_JOBS ))
HEAP="16g"

LOG="$BASE/logs"
mkdir -p "$LOG"

# Groups
mapfile -t GROUPS_LIST < <(
  awk -F'\t' 'NR>1{print $1}' "$LIST_TSV" | sort -u
)

process_group() {
  local group="$1"

  local gdir="$BASE/$group"
  local asm="$gdir/asm/final.contigs.fa"
  local clean="$gdir/asm/final.contigs.clean.fa"
  local idx="$gdir/bbmap_index"
  local mapdir="$gdir/map"
  local covdir="$gdir/coverage"

  mkdir -p "$mapdir" "$covdir"

  # Clean FASTA headers
  awk '{
    if ($0 ~ /^>/) {
      sub(/^>/,""); split($0,a,/[\t ]+/); print ">" a[1]
    } else { print }
  }' "$asm" > "$clean"

  # Build index
  rm -rf "$idx"
  bbmap.sh ref="$clean" path="$idx" -Xmx"$HEAP" >"$LOG/bbmap_ref_${group}.log" 2>&1

  # Samples
  mapfile -t ROWS < <(
    awk -F'\t' -v g="$group" 'NR>1 && $1==g {print $2"\t"$3"\t"$4}' "$LIST_TSV"
  )

  for row in "${ROWS[@]}"; do
    sample=$(echo "$row" | cut -f1)
    r1=$(echo "$row" | cut -f2)
    r2=$(echo "$row" | cut -f3)

    bam="$mapdir/${sample}.bam"
    sbam="$mapdir/${sample}.sorted.bam"

    bbmap.sh in1="$r1" in2="$r2" out="$bam" t="$MAP_THREADS" path="$idx" -Xmx"$HEAP" \
      >"$LOG/bbmap_${group}_${sample}.log" 2>&1

    samtools sort -@ "$MAP_THREADS" -o "$sbam" "$bam"
    samtools index "$sbam"
    rm -f "$bam"
  done

  jgi_summarize_bam_contig_depths \
    --outputDepth "$covdir/${group}.depth.txt" \
    "$mapdir"/*.sorted.bam \
    >"$LOG/jgi_${group}.log" 2>&1
}

export -f process_group
export BASE LIST_TSV MAP_THREADS HEAP LOG

parallel --jobs "$GROUP_JOBS" \
  process_group {} ::: "${GROUPS_LIST[@]}"
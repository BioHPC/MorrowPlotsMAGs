#!/usr/bin/env bash
set -euo pipefail

BASE="/lab/biohpc/Metagenomics_Huang/CoAssemble"
LOG="$BASE/logs"
THREADS_TOTAL=48
GROUPS_JOBS=6
CORES=$(( THREADS_TOTAL / GROUPS_JOBS ))
LIST_TSV="${1:-$BASE/group_samples_master.tsv}"

mkdir -p "$LOG"

# Get unique group names
mapfile -t GROUP_NAMES < <(
  awk -F'\t' 'NR>1{print $1}' "$LIST_TSV" | sort -u
)

coassemble_one() {
  local group="$1"
  local list_tsv="$2"
  local base="$3"
  local logdir="$4"
  local cores="$5"

  local outdir="$base/${group}"
  local asmdir="$outdir/asm"
  local tmpdir
  tmpdir=$(mktemp -d "${asmdir}.tmp.XXXXXX")

  mkdir -p "$outdir" "$logdir"

  local R1S R2S
  R1S=$(awk -F'\t' -v g="$group" 'NR>1 && $1==g {print $3}' "$list_tsv" | paste -sd, -)
  R2S=$(awk -F'\t' -v g="$group" 'NR>1 && $1==g {print $4}' "$list_tsv" | paste -sd, -)

  megahit --presets meta-large \
    -1 "$R1S" -2 "$R2S" \
    --min-contig-len 2000 \
    -t "$cores" \
    -o "$tmpdir" \
    >"$logdir/megahit_${group}.log" 2>&1

  rm -rf "$asmdir"
  mv "$tmpdir" "$asmdir"
}

export -f coassemble_one
export BASE LOG

parallel --jobs "$GROUPS_JOBS" \
  coassemble_one {} "$LIST_TSV" "$BASE" "$LOG" "$CORES" ::: "${GROUP_NAMES[@]}"
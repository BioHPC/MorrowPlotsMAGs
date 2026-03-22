#!/usr/bin/env bash
set -euo pipefail

BASE="/lab/biohpc/Metagenomics_Huang/CoAssemble"
THREADS=8

# Groups
mapfile -t GROUPS_LIST < <(
  find "$BASE" -maxdepth 1 -mindepth 1 -type d -printf "%f\n"
)

process_group() {
  local g="$1"
  local gdir="$BASE/$g"

  local bins="$gdir/refine/${g}_DASTool_DASTool_bins"
  local outdir="$gdir/checkm"

  mkdir -p "$outdir"

  # -------- CheckM lineage --------
  checkm lineage_wf \
    -t "$THREADS" \
    -x fa \
    "$bins" \
    "$outdir"

  # -------- CheckM QA --------
  checkm qa \
    "$outdir/lineage.ms" \
    "$outdir" \
    -o 2 \
    -f "$outdir/qa_summary.tsv"
}

export -f process_group
export BASE THREADS

parallel --jobs 4 process_group {} ::: "${GROUPS_LIST[@]}"
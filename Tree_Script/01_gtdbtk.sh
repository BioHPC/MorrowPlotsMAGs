#!/usr/bin/env bash
set -euo pipefail

BASE="/lab/biohpc/Metagenomics_Huang/CoAssemble"
CPUS=8
EXT="fa"

export GTDBTK_DATA_PATH="/lab/biohpc/Metagenomics_Huang/gtdbtk_db/release226"

for gdir in "$BASE"/*; do
  group=$(basename "$gdir")

  bins="$gdir/refine/${group}_DASTool_DASTool_bins"
  out="$gdir/gtdbtk"
  tmp="$out/tmp"

  mkdir -p "$out" "$tmp"

  gtdbtk classify_wf \
    --genome_dir "$bins" \
    --out_dir "$out" \
    --extension "$EXT" \
    --cpus "$CPUS" \
    --skip_ani_screen \
    --tmpdir "$tmp" \
    --force
done
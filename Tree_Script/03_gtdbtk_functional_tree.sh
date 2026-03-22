#!/usr/bin/env bash
set -euo pipefail

BASE="/lab/biohpc/Metagenomics_Huang/new_data_work/tree_all_groups/combined/itol_export_overwrite/tree_AR_only_plus_NitrospiraC"

export GTDBTK_DATA_PATH="/lab/biohpc/Metagenomics_Huang/gtdbtk_db/release226"

GEN_AR="$BASE/genomes_ar"
GEN_BAC="$BASE/genomes_bac"

OUT_AR="$BASE/ar53/gtdbtk_denovo_archaea"
OUT_BAC="$BASE/bac120/gtdbtk_denovo_bacteria_NitrospiraC"

LOGDIR="$BASE/logs"

mkdir -p "$LOGDIR" "$OUT_AR" "$OUT_BAC"

# -------- GTDB-Tk --------
gtdbtk de_novo_wf \
  --genome_dir "$GEN_AR" \
  --archaea \
  --outgroup_taxon "p__Halobacteriota" \
  --out_dir "$OUT_AR" \
  --extension fa \
  --cpus 24 \
  --force \
  > "$LOGDIR/gtdbtk_ar53.log" 2>&1

gtdbtk de_novo_wf \
  --genome_dir "$GEN_BAC" \
  --bacteria \
  --outgroup_taxon "p__Cyanobacteriota" \
  --out_dir "$OUT_BAC" \
  --extension fa \
  --cpus 24 \
  --force \
  > "$LOGDIR/gtdbtk_bac120.log" 2>&1

# -------- Clean + prune --------
TREE_AR="$OUT_AR/gtdbtk.ar53.decorated.tree"
TREE_BAC="$OUT_BAC/gtdbtk.bac120.decorated.tree"

KEEP_AR="$BASE/itol/ar53/ar53.keep_final.K20plusUser.plusSG.txt"
KEEP_BAC="$BASE/itol/bac120/bac120.keep_final.K4plusUser.plusSG.txt"

CLEAN_AR="$BASE/itol/ar53/gtdbtk.ar53.decorated.clean.tree"
CLEAN_BAC="$BASE/itol/bac120/gtdbtk.bac120.decorated.clean.tree"

PRUNED_AR="$BASE/itol/ar53/gtdbtk.ar53.decorated.pruned.tree"
PRUNED_BAC="$BASE/itol/bac120/gtdbtk.bac120.decorated.pruned.tree"

tr -d "'" < "$TREE_AR" > "$CLEAN_AR"
tr -d "'" < "$TREE_BAC" > "$CLEAN_BAC"

echo ";" >> "$CLEAN_AR"
echo ";" >> "$CLEAN_BAC"

nw_prune -f "$KEEP_AR" "$CLEAN_AR" > "$PRUNED_AR"
nw_prune -f "$KEEP_BAC" "$CLEAN_BAC" > "$PRUNED_BAC"
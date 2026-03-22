#!/usr/bin/env bash
set -euo pipefail

BASE="/lab/biohpc/Metagenomics_Huang/new_data_work/tree_all_groups/combined/itol_export_overwrite/tree_AR_only_plus_NitrospiraC"

PRUNED_AR="$BASE/itol/ar53/gtdbtk.ar53.decorated.pruned.K20plusUser.plusSG.tree"
PRUNED_BAC="$BASE/itol/bac120/gtdbtk.bac120.decorated.pruned.K4plusUser.plusSG.tree"

ALIGN_AR_GZ="$BASE/ar53/gtdbtk_denovo_archaea/align/gtdbtk.ar53.msa.fasta.gz"
ALIGN_BAC_GZ="$BASE/bac120/gtdbtk_denovo_bacteria_NitrospiraC/align/gtdbtk.bac120.msa.fasta.gz"

KEEP_AR_TIPS="$BASE/itol/ar53/ar53.keep_final.K20plusUser.plusSG.TREEIDS.txt"
KEEP_BAC_TIPS="$BASE/itol/bac120/bac120.keep_final.K4plusUser.plusSG.TREEIDS.txt"

SUB_AR="$BASE/itol/ar53/ar53.K20plusUser.plusSG.subset.msa.fasta"
SUB_BAC="$BASE/itol/bac120/bac120.K4plusUser.plusSG.subset.msa.fasta"

# Archaea subset
zcat "$ALIGN_AR_GZ" | seqkit grep -f "$KEEP_AR_TIPS" - > "$SUB_AR"
echo "[AR subset seqs]"; grep -c '^>' "$SUB_AR"

# Bacteria subset
zcat "$ALIGN_BAC_GZ" | seqkit grep -f "$KEEP_BAC_TIPS" - > "$SUB_BAC"
echo "[BAC subset seqs]"; grep -c '^>' "$SUB_BAC"

# Archaea
iqtree -s "$SUB_AR" \
  -m MFP \
  -bb 1000 -alrt 1000 \
  -nt AUTO \
  -pre "$BASE/itol/ar53/iqtree_ar53_K20plusUser.plusSG"

# Bacteria
iqtree -s "$SUB_BAC" \
  -m MFP \
  -bb 1000 -alrt 1000 \
  -nt AUTO \
  -pre "$BASE/itol/bac120/iqtree_bac120_K4plusUser.plusSG"
  
ls -lh "$BASE/itol/ar53/"iqtree_ar53_K20plusUser.plusSG*
ls -lh "$BASE/itol/bac120/"iqtree_bac120_K4plusUser.plusSG*

# show a few internal-node supports (numbers before ':' after a ')')
grep -oE '\)[0-9]+(\.[0-9]+)?:' "$BASE/itol/ar53/iqtree_ar53_K20plusUser.plusSG.treefile" | head
grep -oE '\)[0-9]+(\.[0-9]+)?:' "$BASE/itol/bac120/iqtree_bac120_K4plusUser.plusSG.treefile" | head

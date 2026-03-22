#!/usr/bin/env bash
set -euo pipefail

BASE="/lab/biohpc/Metagenomics_Huang/CoAssemble"
OUT="$BASE/tree_all_groups"

mkdir -p "$OUT"

FT=$(command -v FastTree || command -v fasttree)
ZCAT=$(command -v zcat || command -v gzcat)

# Collect MSA 
build_msa() {
  local domain="$1"
  local out="$OUT/${domain}_user_all.faa"
  : > "$out"

  for g in "$BASE"/*; do
    grp=$(basename "$g")
    ad="$g/gtdbtk/align"

    f=""
    [[ -s "$ad/gtdbtk.${domain}.user_msa.fasta.gz" ]] && f="$ad/gtdbtk.${domain}.user_msa.fasta.gz"
    [[ -z "$f" && -s "$ad/gtdbtk.${domain}.user_msa.fasta" ]] && f="$ad/gtdbtk.${domain}.user_msa.fasta"

    if [[ "$f" == *.gz ]]; then
      $ZCAT "$f" | awk -v G="$grp" '
        /^>/ {sub(/^>/,""); sub(/ .*/,""); print ">"G"__"$0; next}
        {print}
      ' >> "$out"
    else
      awk -v G="$grp" '
        /^>/ {sub(/^>/,""); sub(/ .*/,""); print ">"G"__"$0; next}
        {print}
      ' "$f" >> "$out"
    fi
  done
}

build_msa bac120
build_msa ar53

# Build trees 
for domain in bac120 ar53; do
  aln="$OUT/${domain}_user_all.faa"
  $FT -lg -gamma "$aln" > "$OUT/${domain}_all.newick"
done

#Taxonomy map (phylum) 
TMP="$OUT/_gtdbtk_all.tsv"
: > "$TMP"

for g in "$BASE"/*; do
  grp=$(basename "$g")
  for s in \
    "$g/gtdbtk/classify/summary.tsv" \
    "$g/gtdbtk/gtdbtk.bac120.summary.tsv" \
    "$g/gtdbtk/gtdbtk.ar53.summary.tsv"
  do
    [[ -s "$s" ]] && awk -v G="$grp" 'BEGIN{FS=OFS="\t"} NR>1{print G,$0}' "$s" >> "$TMP"
  done
done

awk -F'\t' '
{
  split($3,a,";");
  ph=a[2]; gsub(/^p__/,"",ph);
  print $2 "\t" $1 "\t" ph
}' "$TMP" | sort -u > "$OUT/mag_taxonomy.tsv"

rm -f "$TMP"
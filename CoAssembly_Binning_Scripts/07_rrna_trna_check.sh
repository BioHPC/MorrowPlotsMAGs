#!/usr/bin/env bash
set -euo pipefail

ROOT="/lab/biohpc/Metagenomics_Huang/new_data_work"

for group_dir in "$ROOT"/*; do
  group_name="$(basename "$group_dir")"

  dastool_bins_dir="$(find "$group_dir" -maxdepth 3 -type d -name "*DASTool*bins*" | head -n 1)"

  mapfile -d '' -t fastas < <(
    find "$dastool_bins_dir" -maxdepth 1 -type f \
    \( -name "*.fa" -o -name "*.fna" -o -name "*.fasta" \) -print0
  )

  cd "$group_dir"
  mkdir -p qc_rrna qc_trna

  # -------- rRNA --------
  for fa in "${fastas[@]}"; do
    base="$(basename "$fa")"
    label="${base%.*}"
    out="qc_rrna/${label}.barrnap.gff"

    barrnap --kingdom bac "$fa" > "$out"
  done

  awk 'BEGIN{OFS="\t"; print "bin_label","rrna_5S","rrna_16S","rrna_23S"}
  FNR==1{
    file=FILENAME
    sub(/^.*\//,"",file)
    sub(/\.barrnap\.gff$/,"",file)
    rr5=rr16=rr23=0
  }
  $0 ~ /Name=5S_rRNA/  {rr5++}
  $0 ~ /Name=16S_rRNA/ {rr16++}
  $0 ~ /Name=23S_rRNA/ {rr23++}
  ENDFILE{ print file, rr5, rr16, rr23 }' qc_rrna/*.barrnap.gff > rrna_counts.tsv

  # -------- tRNA --------
  for fa in "${fastas[@]}"; do
    base="$(basename "$fa")"
    label="${base%.*}"
    out="qc_trna/${label}.trnascan.txt"

    tRNAscan-SE -Q -o "$out" "$fa"
  done

  awk 'BEGIN{OFS="\t"; print "bin_label","trna_count"}
  FNR==1{
    file=FILENAME
    sub(/^.*\//,"",file)
    sub(/\.trnascan\.txt$/,"",file)
    n=0
  }
  $0 !~ /^Sequence/ && $0 !~ /^Name/ && $0 !~ /^-+$/ && $0 !~ /^#/ && $0 !~ /^$/ { n++ }
  ENDFILE{ print file, n }' qc_trna/*.trnascan.txt > trna_counts.tsv

done
#!/usr/bin/env bash
set -euo pipefail

BASE="/lab/biohpc/Metagenomics_Huang/CoAssemble"
THREADS=8
SCORE=0.6

# Groups
mapfile -t GROUPS_LIST < <(
  find "$BASE" -maxdepth 1 -mindepth 1 -type d -printf "%f\n"
)

process_group() {
  local g="$1"
  local gdir="$BASE/$g"

  local maps="$gdir/tmp_maps"
  local outdir="$gdir/refine"
  local contigs="$gdir/asm/final.contigs.fa"

  [[ -s "$gdir/asm/final.contigs.clean.fa" ]] && \
    contigs="$gdir/asm/final.contigs.clean.fa"

  mkdir -p "$maps" "$outdir"

  # -------- MetaBAT2 map --------
  : > "$maps/metabat2_scaffolds2bin.tsv"
  for f in "$gdir"/bins/metabat2/*.fa; do
    b=$(basename "$f" .fa)
    awk -v bin="$b" '/^>/{sub(/^>/,""); print $0"\t"bin}' "$f" \
      >> "$maps/metabat2_scaffolds2bin.tsv"
  done

  # -------- MaxBin2 map --------
  : > "$maps/maxbin2_scaffolds2bin.tsv"
  for f in "$gdir"/bins/maxbin2/*.fasta; do
    b=$(basename "$f" .fasta)
    awk -v bin="$b" '/^>/{sub(/^>/,""); print $0"\t"bin}' "$f" \
      >> "$maps/maxbin2_scaffolds2bin.tsv"
  done

  # -------- CONCOCT map --------
  if [[ -s "$gdir/bins/concoct/scaffolds2bin.tsv" ]]; then
    cp "$gdir/bins/concoct/scaffolds2bin.tsv" \
       "$maps/concoct_scaffolds2bin.tsv"
  else
    awk -F',' 'NR>1{print $1"\tconcoct."$2}' \
      "$gdir/concoct/${g}_clustering_merged.csv" \
      > "$maps/concoct_scaffolds2bin.tsv"
  fi

  # -------- DAS Tool --------
  i_files=()
  l_labels=()

  for tool in metabat2 maxbin2 concoct; do
    f="$maps/${tool}_scaffolds2bin.tsv"
    [[ -s "$f" ]] && {
      i_files+=("$f")
      l_labels+=("$tool")
    }
  done

  joined_i=$(IFS=,; echo "${i_files[*]}")
  joined_l=$(IFS=,; echo "${l_labels[*]}")

  outfile="$outdir/${g}_DASTool"

  R_DOCOPT_USE_S4=FALSE DAS_Tool \
    -i "$joined_i" \
    -l "$joined_l" \
    -c "$contigs" \
    -o "$outfile" \
    --search_engine blastp \
    --score_threshold "$SCORE" \
    --threads "$THREADS" \
    --write_bins \
    --write_bin_evals
}

export -f process_group
export BASE THREADS SCORE

parallel --jobs 4 process_group {} ::: "${GROUPS_LIST[@]}"
#!/usr/bin/env bash
set -euo pipefail

BASE="/lab/biohpc/Metagenomics_Huang/CoAssemble"
THREADS=8
LOG="$BASE/logs"

mkdir -p "$LOG"

# Groups
mapfile -t GROUPS_LIST < <(
  find "$BASE" -maxdepth 1 -mindepth 1 -type d -printf "%f\n"
)

process_group() {
  local g="$1"

  local gdir="$BASE/$g"
  local clean="$gdir/asm/final.contigs.clean.fa"
  local depth="$gdir/coverage/${g}.depth.txt"

  local bindir="$gdir/bins"
  local concoct_dir="$gdir/concoct"
  local mapdir="$gdir/map_concoct"
  local covtsv="$concoct_dir/${g}_coverage.tsv"
  local cutfa="$concoct_dir/${g}_10K.fa"
  local bed="$concoct_dir/${g}_10K.bed"

  mkdir -p "$bindir"/{metabat2,maxbin2,concoct} "$concoct_dir"

  # -------- MetaBAT2 --------
  metabat2 -i "$clean" -a "$depth" \
    -o "$bindir/metabat2/${g}.bin" \
    -m 2000 -t "$THREADS" \
    >"$LOG/metabat2_${g}.log" 2>&1 || true

  # -------- MaxBin2 --------
  run_MaxBin.pl \
    -contig "$clean" \
    -abund "$depth" \
    -out "$bindir/maxbin2/${g}" \
    -min_contig_length 2000 \
    -prob 0.8 \
    -thread "$THREADS" \
    >"$LOG/maxbin2_${g}.log" 2>&1 || true

  # -------- CONCOCT --------
  cut_up_fasta.py "$clean" -c 10000 -o 0 --merge_last -b "$bed" > "$cutfa"

  if ls "$gdir/map_concoct"/*.sorted.bam >/dev/null 2>&1; then
    mapdir="$gdir/map_concoct"
  else
    mapdir="$gdir/map"
  fi

  concoct_coverage_table.py "$bed" "$mapdir"/*.sorted.bam > "$covtsv"

  concoct \
    --composition_file "$cutfa" \
    --coverage_file "$covtsv" \
    -b "$concoct_dir/" \
    -t "$THREADS" \
    >"$LOG/concoct_${g}.log" 2>&1

  # -------- scaffolds2bin --------
  csv="$concoct_dir/clustering_gt1000.csv"
  map="$bindir/concoct/scaffolds2bin.tsv"

  awk -F, 'NR>1{
    cont=$1; sub(/\.concoct_part_[0-9]+$/,"",cont)
    bin=$2; key=cont"\t"bin; cnt[key]++
  } END{
    for(k in cnt){
      split(k,a,"\t"); c=a[1]; b=a[2]
      if(cnt[k]>max[c]){max[c]=cnt[k]; best[c]=b}
    }
    for(c in best){ print c"\tconcoct."best[c] }
  }' "$csv" > "$map"

  # -------- Export bins --------
  awk '{print $2}' "$map" | sort -u | while read -r BIN; do
    awk -v b="$BIN" '$2==b{print $1}' "$map" > "$bindir/concoct/.list.$BIN"
    seqtk subseq "$clean" "$bindir/concoct/.list.$BIN" \
      > "$bindir/concoct/${BIN}.fa"
    rm -f "$bindir/concoct/.list.$BIN"
  done
}

export -f process_group
export BASE THREADS LOG

parallel --jobs 6 process_group {} ::: "${GROUPS_LIST[@]}"

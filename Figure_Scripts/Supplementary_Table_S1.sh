#!/usr/bin/env bash
set -euo pipefail

BASE="/lab/biohpc/Metagenomics_Huang/new_data_work"
OUTDIR="$BASE/summary"
mkdir -p "$OUTDIR"

ASM_TSV="$OUTDIR/assembly_stats.tsv"
MAP_TSV="$OUTDIR/mapping_stats.tsv"
BIN_TSV="$OUTDIR/bin_counts.tsv"

fa_stats_ge() {
  local fa="$1" minlen="$2"
  cat "$fa" |
  awk -v MIN="$minlen" '
    /^>/ {
      if (len>=MIN) { n++; total+=len; lens[n]=len; if (len>maxlen) maxlen=len }
      len=0; next
    }
    { gsub(/[ \t\r]/,"",$0); len+=length($0) }
    END {
      if (len>=MIN) { n++; total+=len; lens[n]=len; if (len>maxlen) maxlen=len }
      n50=0;
      if (n>0) {
        for(i=1;i<=n;i++) for(j=i+1;j<=n;j++) if (lens[j]>lens[i]) { t=lens[i]; lens[i]=lens[j]; lens[j]=t }
        half=total/2; acc=0
        for(i=1;i<=n;i++){ acc+=lens[i]; if (acc>=half){ n50=lens[i]; break } }
      }
      printf "%d\t%d\t%d\t%d\n", n,total,n50,maxlen
    }'
}

flagstat_to_row() {
  local bam="$1"
  samtools flagstat -@ 1 "$bam" 2>/dev/null |
  awk '
    /in total/ && !t {total=$1; t=1}
    /mapped \(/ && !m {mapped=$1; match($0,/\(([0-9.]+)%/,a); pct=a[1]; m=1}
    END{print total, mapped, (pct==""?0:pct)}'
}


shopt -s nullglob
for gdir in "$BASE"/*; do
  group="$(basename "$gdir")"

  asmfa="$gdir/asm/final.contigs.clean.fa"

  read c2 tot n50 maxlen < <(fa_stats_ge "$asmfa" 2000)

  bins_dir="$gdir/refine/${group}_DASTool_DASTool_bins"
  n_bins=$(find "$bins_dir" -maxdepth 1 -type f \( -name '*.fa' -o -name '*.fasta' -o -name '*.fna' \) | wc -l | tr -d ' ')

  mapdir="$gdir/map_concoct"

  for bam in "$mapdir"/*.sorted.bam; do
    s="$(basename "$bam")"
    s="${s%.sorted.bam}"

    echo -e "$s\t$c2\t$tot\t$n50\t$maxlen" >> "$ASM_TSV"

    read tot_reads mapped_reads pct <<<"$(flagstat_to_row "$bam")"
    echo -e "$s\t$tot_reads\t$mapped_reads\t$pct" >> "$MAP_TSV"

    echo -e "$s\t$n_bins\t$bins_dir" >> "$BIN_TSV"
  done
done
shopt -u nullglob
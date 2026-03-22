#!/usr/bin/env bash
set -euo pipefail

BASE="/lab/biohpc/Metagenomics_Huang/new_data_work"
WORK="/lab/biohpc/Metagenomics_Huang/new_data/work_new"
MASTER="$BASE/group_samples_master_newdata.tsv"
QUAST_ROOT="$BASE/quast"

OUTDIR="$BASE/summary"
OUTTSV="$OUTDIR/summary_like_pic.tsv"
mkdir -p "$OUTDIR"

MULTIQC_FASTQC="$WORK/qc/multiqc_data/multiqc_fastqc.txt"

TMP_DIR="$(mktemp -d)"
SAMPLE_READS_TSV="$TMP_DIR/sample_reads.tsv"
GROUP_MERGED_TSV="$TMP_DIR/group_merged_reads.tsv"
GROUP_METRICS_TSV="$TMP_DIR/group_metrics.tsv"
GROUP_ORDER_TSV="$TMP_DIR/group_order.tsv"
TMP_ROWS="$TMP_DIR/rows.tsv"

# -------- helpers --------

quast_metrics_for_group() {
  local group="$1"
  local rpt="$QUAST_ROOT/$group/report.tsv"

  awk '
    BEGIN{cont2k=""; cont0=""; n50="NA"; l50="NA"}
    {
      line=$0
      gsub(/\r/,"",line)
      gsub(/[[:space:]]{2,}/,"\t",line)
      split(line,a,"\t")
      key=a[1]; val=a[2]

      if(key ~ /^# contigs \(>= 2000 bp\)$/) cont2k=val
      if(key ~ /^# contigs \(>= 0 bp\)$/)    cont0=val
      if(key=="N50") n50=val
      if(key=="L50") l50=val
    }
    END{
      cont = (cont2k!="" ? cont2k : (cont0!="" ? cont0 : "NA"))
      print cont "\t" n50 "\t" l50
    }
  ' "$rpt"
}

mag_count_for_group() {
  local group="$1"
  local bins="$BASE/$group/refine/${group}_DASTool_DASTool_bins"

  shopt -s nullglob
  local arr=( "$bins"/*.fa "$bins"/*.fasta "$bins"/*.fna )
  shopt -u nullglob

  echo "${#arr[@]}"
}

# -------- Step 0 --------
awk -F'\t' 'NR>1 && $1!=""{print $1}' "$MASTER" | awk '!seen[$0]++' > "$GROUP_ORDER_TSV"

# -------- Step 1 --------
awk '
  BEGIN{FS=OFS="\t"}
  NR==1{
    for(i=1;i<=NF;i++){
      if($i=="Sample") s=i
      if($i=="Total Sequences") t=i
    }
    next
  }
  {
    samp=$s; tot=$t
    if(match(samp, /^MS[0-9]+/)){
      ms=substr(samp, RSTART, RLENGTH)
      if(samp ~ /_merged_R1_001$/){ r1[ms]=tot }
      if(samp ~ /_merged_1P$/){ p1[ms]=tot }
      seen[ms]=1
    }
  }
  END{
    for(ms in seen){
      val="NA"
      if(ms in r1) val=r1[ms]
      else if(ms in p1) val=p1[ms]
      print ms, val
    }
  }
' "$MULTIQC_FASTQC" | sort -t $'\t' -k1,1 > "$SAMPLE_READS_TSV"

# -------- Step 2 --------
awk -F'\t' '
  BEGIN{OFS="\t"}
  FNR==NR{reads[$1]=$2; next}
  NR==1{next}
  {
    g=$1; ms=$2
    v=reads[ms]
    if(v!="" && v!="NA") sum[g]+=v
    seen[g]=1
    if(!(g in first)) first[g]=++k
  }
  END{
    for(g in first){ ord[first[g]]=g }
    for(i=1;i<=k;i++){
      g=ord[i]
      if(g in sum) print g, sum[g]
      else print g, "NA"
    }
  }
' "$SAMPLE_READS_TSV" "$MASTER" > "$GROUP_MERGED_TSV"

# -------- Step 3 --------
echo -e "group\tcontigs_ge2kb\tn50\tl50\tmags" > "$GROUP_METRICS_TSV"
while read -r grp; do
  read cont n50 l50 < <(quast_metrics_for_group "$grp")
  mags="$(mag_count_for_group "$grp")"
  echo -e "${grp}\t${cont}\t${n50}\t${l50}\t${mags}" >> "$GROUP_METRICS_TSV"
done < "$GROUP_ORDER_TSV"

# -------- Step 4 --------
awk -F'\t' '
  BEGIN{OFS="\t"}
  FILENAME==ARGV[1]{ sreads[$1]=$2; next }
  FILENAME==ARGV[2]{ gmerge[$1]=$2; next }
  FILENAME==ARGV[3]{
    if(FNR==1) next
    g=$1
    cont[g]=$2; n50[g]=$3; l50[g]=$4; mags[g]=$5
    next
  }
  FILENAME==ARGV[4]{
    if(FNR==1) next
    g=$1; ms=$2

    pe = (ms in sreads ? sreads[ms] : "NA")
    mr = (g in gmerge ? gmerge[g] : "NA")
    c  = (g in cont  ? cont[g]  : "NA")
    n  = (g in n50   ? n50[g]   : "NA")
    l  = (g in l50   ? l50[g]   : "NA")
    m  = (g in mags  ? mags[g]  : "NA")

    print g, ms, pe, mr, c, n, l, m
  }
' "$SAMPLE_READS_TSV" "$GROUP_MERGED_TSV" "$GROUP_METRICS_TSV" "$MASTER" > "$TMP_ROWS"

echo -e "Co-assembly_group\tSample\tPaired-end_filtered_reads\tMerged_reads_len\tContigs >2000 bp\tN50\tL50\tMAGs" > "$OUTTSV"

LC_ALL=C sort -t $'\t' -k1,1 -k2,2V "$TMP_ROWS" >> "$OUTTSV"
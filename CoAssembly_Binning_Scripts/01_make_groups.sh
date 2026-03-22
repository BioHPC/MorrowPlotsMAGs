#!/usr/bin/env bash
set -euo pipefail

BASE="/lab/biohpc/Metagenomics_Huang"
TRIM="${BASE}/work/trim"
META="${BASE}/SingleAssembly/summary/sample_metadata.tsv"
OUT="${BASE}/CoAssemble"

MASTER="${OUT}/group_samples_master.tsv"
OUT0="${OUT}/group_samples_shard0.tsv"
OUT1="${OUT}/group_samples_shard1.tsv"

THREADS=32
N_FIRST=6

mkdir -p "${OUT}"

# Collect all R1 once
mapfile -t ALL_R1 < <(ls -1 "${TRIM}"/*_1P.fq.gz)

# Build sample-group pairs
mapfile -t SAMPLE_GROUPS < <(
  awk -F'\t' 'NR>1{printf "%s\t%s_%s\n", $1, $2, $3}' "${META}" | sort -u
)

# Function to collect reads
collect_reads() {
  local sample="$1"
  local pat="^${sample}(_|$)"
  local r1s=() r2s=()

  for r1 in "${ALL_R1[@]}"; do
    base=$(basename "$r1")
    base="${base%_1P.fq.gz}"
    if [[ "$base" =~ $pat ]]; then
      r2="${TRIM}/${base}_2P.fq.gz"
      r1s+=("${TRIM}/${base}_1P.fq.gz")
      r2s+=("${r2}")
    fi
  done

  printf "%s\t%s\n" "$(IFS=,; echo "${r1s[*]}")" "$(IFS=,; echo "${r2s[*]}")"
}

# Create master + per-group files
echo -e "group\tsample\tR1s\tR2s" > "$MASTER"

declare -A SEEN

while IFS=$'\t' read -r sample group; do
  read -r R1s R2s < <(collect_reads "$sample")

  echo -e "${group}\t${sample}\t${R1s}\t${R2s}" >> "$MASTER"

  if [[ -z "${SEEN[$group]:-}" ]]; then
    mkdir -p "${OUT}/${group}"
    echo -e "sample\tR1_count\tR1s\tR2s" > "${OUT}/${group}/group_samples.tsv"
    SEEN[$group]=1
  fi

  rcnt=$(tr ',' '\n' <<< "$R1s" | wc -l)
  echo -e "${sample}\t${rcnt}\t${R1s}\t${R2s}" >> "${OUT}/${group}/group_samples.tsv"

done < <(printf "%s\n" "${SAMPLE_GROUPS[@]}")

# Split into shards 

header=$(head -n1 "$MASTER")

mapfile -t GRPS < <(awk -F'\t' 'NR>1{print $1}' "$MASTER" | sort -u)

FIRST=("${GRPS[@]:0:N_FIRST}")
REST=("${GRPS[@]:N_FIRST}")

printf "%s\n" "$header" > "$OUT0"
printf "%s\n" "$header" > "$OUT1"

for g in "${FIRST[@]}"; do
  awk -F'\t' -v grp="$g" 'NR>1 && $1==grp' "$MASTER" >> "$OUT0"
done

for g in "${REST[@]}"; do
  awk -F'\t' -v grp="$g" 'NR>1 && $1==grp' "$MASTER" >> "$OUT1"
done
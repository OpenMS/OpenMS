#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCE_DIR="$(cd "$SCRIPT_DIR/../../.." && pwd)"
BUILD_DIR="${BUILD_DIR:-$SOURCE_DIR/../OpenMS-build}"
OPEN_DIA="${OPEN_DIA:-$BUILD_DIR/bin/OpenDIA}"
THREADS="${THREADS:-4}"
PERCOLATOR_THREADS="${PERCOLATOR_THREADS:-$THREADS}"
(( PERCOLATOR_THREADS > 3 )) && PERCOLATOR_THREADS=3
OPEN_DIA_SPEED_ARGS=(
  -TargetedDataExtraction:readOptions cacheWorkingInMemory
  -TargetedDataExtraction:batchSize 0
  -TargetedDataExtraction:innerBatchSize -1
  -TargetedDataExtraction:maxConcurrentSwaths -1
  -TargetedDataExtraction:outer_loop_threads -1
)
OPEN_DIA_ANALYSIS_ARGS=(
  -Rescoring:percolator:num_threads "$PERCOLATOR_THREADS"
  -Rescoring:percolator:seed 1
  -Inference:peptide:global true
  -Inference:protein:global true
  -Inference:peptide:run_specific false
  -Inference:protein:run_specific false
  -Inference:peptide:experiment_wide false
  -Inference:protein:experiment_wide false
  -Inference:gene:global false
  -Inference:gene:run_specific false
  -Inference:gene:experiment_wide false
  -Export:results:max_rs_peakgroup_qvalue 0.05
  # Native OSW annotation also filters with >=, so 1.01 retains every valid q-value.
  -Export:results:use_peptide_scores true
  -Export:results:use_protein_scores true
  -Export:results:max_global_peptide_qvalue 1.01
  -Export:results:max_global_protein_qvalue 1.01
  -Export:results:use_gene_scores false
  -Export:results:exclude_decoys false
  -Export:matrix:precursor true
  -Export:matrix:peptide true
  -Export:matrix:protein true
  -Export:matrix:gene false
  -Export:matrix:max_rs_peakgroup_qvalue 0.01
  -Export:matrix:use_peptide_scores true
  -Export:matrix:use_protein_scores true
  -Export:matrix:max_global_peptide_qvalue 0.01
  -Export:matrix:max_global_protein_qvalue 0.01
)

require_file()
{
  [[ -f "$1" ]] || { echo "ERROR: missing file: $1" >&2; exit 2; }
}

require_executable()
{
  [[ -x "$1" ]] || { echo "ERROR: missing executable: $1" >&2; exit 2; }
}

prepare_output()
{
  if [[ -e "$1" && -n "$(find "$1" -mindepth 1 -print -quit 2>/dev/null)" ]]; then
    echo "ERROR: output directory is not empty: $1" >&2
    exit 2
  fi
  mkdir -p "$1/logs" "$1/meta"
}

run_measured()
{
  local out="$1"
  shift
  local -a command=("$@")

  printf '%q ' "${command[@]}" > "$out/meta/command.sh"
  printf '\n' >> "$out/meta/command.sh"
  sha256sum "$OPEN_DIA" > "$out/meta/OpenDIA.sha256"
  git -C "$SOURCE_DIR" rev-parse HEAD > "$out/meta/source_revision.txt"
  git -C "$SOURCE_DIR" status --short > "$out/meta/source_status.txt"
  printf 'memory_max=%s\nmemory_swap_max=%s\n' "${BENCH_MEMORY_MAX:-unlimited}" "${BENCH_MEMORY_SWAP_MAX:-unlimited}" \
    > "$out/meta/resource_limits.txt"
  local -a limits=()
  if [[ -n "${BENCH_MEMORY_MAX:-}" ]]; then
    limits=(systemd-run --user --scope --quiet -p "MemoryMax=$BENCH_MEMORY_MAX")
    if [[ -n "${BENCH_MEMORY_SWAP_MAX:-}" ]]; then
      limits+=(-p "MemorySwapMax=$BENCH_MEMORY_SWAP_MAX")
    fi
    limits+=(--)
  fi
  printf 'Starting measured OpenDIA job: %s at %s\n' "$out" "$(date -Iseconds)"

  /usr/bin/time -v -o "$out/logs/resources.txt" \
    "${limits[@]}" "${command[@]}" > "$out/logs/stdout.log" 2> "$out/logs/stderr.log" &
  local pid=$!
  python3 "$SCRIPT_DIR/monitor_host.py" "$pid" "$out/logs" \
    --min-free-disk-gib "${MIN_FREE_DISK_GIB:-3}" &
  local host_monitor=$!

  (
    while kill -0 "$pid" 2>/dev/null; do
      du -sb "$out" 2>/dev/null | awk -v now="$(date -Iseconds)" '{print now "\t" $1}'
      sleep 1
    done
    du -sb "$out" 2>/dev/null | awk -v now="$(date -Iseconds)" '{print now "\t" $1}'
  ) > "$out/logs/disk_bytes.tsv" &
  local monitor=$!

  set +e
  wait "$pid"
  local status=$?
  wait "$monitor"
  wait "$host_monitor"
  local host_status=$?
  set -e

  printf '%s\n' "$status" > "$out/meta/exit_code.txt"
  awk 'BEGIN { max=0 } $2 > max { max=$2 } END { print max }' \
    "$out/logs/disk_bytes.tsv" > "$out/meta/peak_output_bytes.txt"
  find "$out" -type f -printf '%s\t%P\n' | sort -nr > "$out/meta/files.tsv"
  if (( host_status != 0 )); then
    echo "ERROR: host workload monitor failed for $out" >&2
    return "$host_status"
  fi
  return "$status"
}

#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCE_DIR="$(cd "$SCRIPT_DIR/../../.." && pwd)"
BUILD_DIR="${BUILD_DIR:-$SOURCE_DIR/../OpenMS-build}"
OPEN_DIA="${OPEN_DIA:-$BUILD_DIR/bin/OpenDIA}"
THREADS="${THREADS:-4}"
OPEN_DIA_SPEED_ARGS=(
  -TargetedDataExtraction:readOptions cacheWorkingInMemory
  -TargetedDataExtraction:batchSize 0
  -TargetedDataExtraction:innerBatchSize -1
  -TargetedDataExtraction:maxConcurrentSwaths -1
  -TargetedDataExtraction:outer_loop_threads -1
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

  /usr/bin/time -v -o "$out/logs/resources.txt" \
    "${command[@]}" > "$out/logs/stdout.log" 2> "$out/logs/stderr.log" &
  local pid=$!

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
  set -e

  printf '%s\n' "$status" > "$out/meta/exit_code.txt"
  awk 'BEGIN { max=0 } $2 > max { max=$2 } END { print max }' \
    "$out/logs/disk_bytes.tsv" > "$out/meta/peak_output_bytes.txt"
  find "$out" -type f -printf '%s\t%P\n' | sort -nr > "$out/meta/files.tsv"
  return "$status"
}

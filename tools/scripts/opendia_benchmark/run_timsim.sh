#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib.sh
source "$SCRIPT_DIR/lib.sh"

DATASET="${1:-/home/sing/Documents/datasets/OpenMS-TimSim-Fixture/packages/opendia_timsim_fixture_validation_3x3}"
OUT_ROOT="${2:-$PWD/opendia-benchmark/timsim}"
MODE="${3:-all}"

[[ "$MODE" == "all" || "$MODE" == "target" || "$MODE" == "entrapment" ]] || {
  echo "Usage: $0 [DATASET_DIR] [OUT_DIR] [all|target|entrapment]" >&2
  exit 2
}
require_executable "$OPEN_DIA"
require_file "$DATASET/metadata/fixture_manifest.json"
require_file "$DATASET/SHA256SUMS"

prepare_output "$OUT_ROOT"
if [[ "${VERIFY_INPUTS:-1}" == "1" ]]; then
  (
    cd "$DATASET"
    awk '$2 ~ /^\.\/(raw|libraries|truth|metadata)\//' SHA256SUMS | sha256sum --check -
  ) || { echo "ERROR: TimSim benchmark input checksum failed" >&2; exit 2; }
fi

RAW_ROOT="$OUT_ROOT/input"
mkdir -p "$RAW_ROOT"
for archive in "$DATASET"/raw/*.d.zip; do
  run="$(basename "$archive" .zip)"
  if [[ ! -f "$RAW_ROOT/$run/analysis.tdf" ]]; then
    unzip -q "$archive" -d "$RAW_ROOT"
  fi
done
mapfile -t RAW_FILES < <(find "$RAW_ROOT" -mindepth 1 -maxdepth 1 -type d -name '*.d' | sort)
(( ${#RAW_FILES[@]} == 6 )) || { echo "ERROR: expected 6 extracted .d inputs, found ${#RAW_FILES[@]}" >&2; exit 2; }
printf '%s\n' "${RAW_FILES[@]}" > "$OUT_ROOT/meta/raw_files.txt"

run_library()
{
  local name="$1"
  local library="$2"
  local out="$OUT_ROOT/$name"
  require_file "$library"
  prepare_output "$out"
  mkdir -p "$out/results" "$out/intermediates"
  sha256sum "$library" > "$out/meta/library.sha256"

  local -a command=(
    env "LD_LIBRARY_PATH=$BUILD_DIR/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
    "$OPEN_DIA"
    -in "${RAW_FILES[@]}"
    -tr "$library"
    -workflow:library_mode transition_list
    -workflow:working_format parquet
    -workflow:keep_intermediate_files true
    -workflow:intermediate_dir "$out/intermediates"
    -out_dir "$out/results"
    -TargetedDataExtraction:ion_mobility_mode auto
    "${OPEN_DIA_SPEED_ARGS[@]}"
    -TargetedDataExtraction:write_debug_files false
    -Export:results:max_rs_peakgroup_qvalue 0.01
    -Export:matrix:peptide true
    -Export:matrix:protein true
    -threads "$THREADS"
  )
  run_measured "$out" "${command[@]}"
}

if [[ "$MODE" == "all" || "$MODE" == "target" ]]; then
  run_library target_only "$DATASET/libraries/OpenSwathTimSim.target_only.transitions.tsv"
fi
if [[ "$MODE" == "all" || "$MODE" == "entrapment" ]]; then
  run_library entrapment "$DATASET/libraries/OpenSwathTimSim.transitions.tsv"
fi
if [[ -s "$OUT_ROOT/target_only/results/OpenDIA.results.tsv" && -s "$OUT_ROOT/entrapment/results/OpenDIA.results.tsv" ]]; then
  python3 "$SCRIPT_DIR/validate_timsim.py" "$DATASET" \
    "$OUT_ROOT/target_only/results/OpenDIA.results.tsv" \
    "$OUT_ROOT/entrapment/results/OpenDIA.results.tsv" \
    --output "$OUT_ROOT/validation.json"
fi

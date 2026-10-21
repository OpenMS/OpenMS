#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCE_DIR="$(cd "$SCRIPT_DIR/../../.." && pwd)"
BUILD_DIR="${1:-${BUILD_DIR:-$SOURCE_DIR/../OpenMS-build}}"
PRESET="${PRESET:-linux-x64-release}"
JOBS="${JOBS:-$(nproc)}"

if [[ -f "$BUILD_DIR/CMakeCache.txt" ]]; then
  cmake -S "$SOURCE_DIR" -B "$BUILD_DIR" \
    -DWITH_OPENTIMS=ON -DWITH_THERMO_RAW=ON
else
  cmake --preset "$PRESET" -S "$SOURCE_DIR" -B "$BUILD_DIR" \
    -DWITH_OPENTIMS=ON -DWITH_THERMO_RAW=ON
fi

cmake --build "$BUILD_DIR" --target OpenDIA OpenSwathExport -j "$JOBS"
for tool in OpenDIA OpenSwathExport; do
  set +e
  "$BUILD_DIR/bin/$tool" --help >/dev/null 2>&1
  help_status=$?
  set -e
  [[ "$help_status" == 0 || "$help_status" == 6 ]] || exit "$help_status"
done

{
  printf 'source_revision=%s\n' "$(git -C "$SOURCE_DIR" rev-parse HEAD)"
  if [[ -n "$(git -C "$SOURCE_DIR" status --porcelain)" ]]; then
    printf 'source_dirty=true\n'
  else
    printf 'source_dirty=false\n'
  fi
  printf 'build_dir=%s\n' "$BUILD_DIR"
  printf 'preset=%s\n' "$PRESET"
  printf 'binary_sha256=%s\n' "$(sha256sum "$BUILD_DIR/bin/OpenDIA" | cut -d' ' -f1)"
  printf 'export_binary_sha256=%s\n' "$(sha256sum "$BUILD_DIR/bin/OpenSwathExport" | cut -d' ' -f1)"
  grep -E '^(CMAKE_BUILD_TYPE|WITH_OPENTIMS|WITH_THERMO_RAW|CMAKE_GENERATOR):' \
    "$BUILD_DIR/CMakeCache.txt" || true
} | tee "$BUILD_DIR/OpenDIA.build-manifest.txt"

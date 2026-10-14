#!/usr/bin/env bash
# Installs the CI-only Windows tools (rclone, ccache) straight from their GitHub
# releases instead of the Chocolatey community feed, whose outages (e.g. 504
# Gateway Timeout) used to fail whole jobs. Each archive is pinned to a version
# and verified against its SHA256 before use.
#
# To update a tool, bump its version and replace the hash:
#   rclone: from the release's SHA256SUMS file
#   ccache: sha256sum of the downloaded zip (no checksum file is published)

set -eu
set -o pipefail

RCLONE_VERSION=1.75.2
RCLONE_SHA256=8a9db2bad1dc79a65b3348fd857069598e2aab7cf4d118b77b1df23055aff818
CCACHE_VERSION=4.14.1
CCACHE_SHA256=6219f3865ca59aec41ee4b678df171d5d35855ecb2b6dbbbd20690b3a68af7b4

install_dir="${RUNNER_TOOL_CACHE:-C:/tools}/openms-ci-tools"
mkdir -p "$install_dir"

# fetch <url> <sha256> <name of the directory inside the zip>
fetch() {
  local url=$1 sha=$2 dir=$3 zip
  zip="$install_dir/$(basename "$url")"
  curl --no-progress-meter --fail --retry 5 --retry-all-errors --retry-delay 10 -L -o "$zip" "$url"
  echo "$sha  $zip" | sha256sum --check --strict -
  7z x -y -o"$install_dir" "$zip" >/dev/null
  rm -f "$zip"
  local bin
  bin=$(cygpath -w "$install_dir/$dir" 2>/dev/null || echo "$install_dir/$dir")
  if [ -n "${GITHUB_PATH:-}" ]; then
    echo "$bin" >>"$GITHUB_PATH"
  fi
  export PATH="$install_dir/$dir:$PATH"
}

[ $# -gt 0 ] || set -- rclone ccache
for tool in "$@"; do
  case $tool in
    rclone)
      fetch "https://github.com/rclone/rclone/releases/download/v${RCLONE_VERSION}/rclone-v${RCLONE_VERSION}-windows-amd64.zip" \
        "$RCLONE_SHA256" "rclone-v${RCLONE_VERSION}-windows-amd64"
      rclone version | head -1 ;;
    ccache)
      fetch "https://github.com/ccache/ccache/releases/download/v${CCACHE_VERSION}/ccache-${CCACHE_VERSION}-windows-x86_64.zip" \
        "$CCACHE_SHA256" "ccache-${CCACHE_VERSION}-windows-x86_64"
      ccache --version | head -1 ;;
    *) echo "unknown tool: $tool" >&2; exit 1 ;;
  esac
done

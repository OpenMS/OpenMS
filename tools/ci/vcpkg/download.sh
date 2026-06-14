#!/bin/bash
set -e

REPO="castorNova2/openms-ci-artefacts"
TRIPLET="x64-linux"
BASELINE=$(jq -r '.["default-registry"].baseline' vcpkg-configuration.json)
RELEASE_TAG="vcpkg-cache-${TRIPLET}-${BASELINE:0:10}"
CACHE_DIR="${HOME}/.cache/vcpkg/archives"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
TMP_DIR="$SCRIPT_DIR/.vcpkg-cache-tmp"

mkdir -p "$TMP_DIR"

gh release download "$RELEASE_TAG" \
    --repo "$REPO" \
    --dir "$TMP_DIR" \
    --clobber

+shopt -s nullglob
for f in "$TMP_DIR"/*.zip; do
    filename=$(basename "$f")     
    prefix="${filename:0:2}"     
    mkdir -p "$CACHE_DIR/$prefix"
    mv "$f" "$CACHE_DIR/$prefix/$filename"
done
+shopt -s nullglob

rm -rf "$TMP_DIR"

echo "Done."
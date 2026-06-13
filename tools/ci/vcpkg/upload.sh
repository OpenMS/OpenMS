#!/bin/bash
set -e


VCPKG_CACHE_DIR="$HOME/.cache/vcpkg/archives"
TRIPLET="x64-linux"
BASELINE=$(jq -r '.["builtin-baseline"]' ../vcpkg.json)
RELEASE_TAG="vcpkg-cache-${TRIPLET}-${BASELINE:0:10}"
REPO="https://github.com/castorNova2/openms-ci-artefacts"

if ! gh release view "$RELEASE_TAG" --repo "$REPO" &>/dev/null; then
    gh release create "$RELEASE_TAG" \
        --repo "$REPO" \
        --title "vcpkg cache ($TRIPLET) baseline ${BASELINE:0:10}" \
        --notes " vcpkg cache for triplet $TRIPLET" \
        --prerelease
else
    echo "Release $RELEASE_TAG already exists"
fi

echo "Fetching already uploaded files"
UPLOADED=$(gh release view "$RELEASE_TAG" --repo "$REPO" --json assets --jq '.assets[].name')

COPIED=0
SKIPPED=0

while read -r archive; do
    FILENAME=$(basename "$archive")
    SUBDIR=$(basename "$(dirname "$archive")")
    CACHE_NAME="${SUBDIR}__${FILENAME}"

    if echo "$UPLOADED" | grep -q "^${CACHE_NAME}$"; then
        SKIPPED=$((SKIPPED + 1))
        continue
    fi

    echo "==> Uploading $CACHE_NAME..."
    gh release upload "$RELEASE_TAG" "$archive#$CACHE_NAME" \
        --repo "$REPO" \
        --clobber

    COPIED=$((COPIED + 1))
done < <(find "$VCPKG_CACHE_DIR" -name "*.zip")

echo "Upload complete"
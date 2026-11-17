#!/usr/bin/env bash
# Re-sync vendored Percolator from upstream.
# Run from src/openms/thirdparty/percolator/.
set -euo pipefail

UPSTREAM_URL="https://github.com/percolator/percolator.git"
UPSTREAM_REF="${1:-a9783ee8ded0af8148566bdff4d43d2ba444ce5e}"

here=$(pwd)
tmpdir=$(mktemp -d)
trap "rm -rf $tmpdir" EXIT

echo ">>> Cloning upstream Percolator at ref $UPSTREAM_REF..."
# --branch accepts branch/tag names but not bare SHAs. Clone master shallow,
# then fetch+checkout the requested ref (works for any commit / branch / tag).
git clone --depth 1 "$UPSTREAM_URL" "$tmpdir" 2>&1 | tail -2
cloned_sha=$(git -C "$tmpdir" rev-parse HEAD)
if [ "$cloned_sha" != "$UPSTREAM_REF" ]; then
  echo "  HEAD is $cloned_sha; fetching $UPSTREAM_REF..."
  git -C "$tmpdir" fetch --depth 1 origin "$UPSTREAM_REF" 2>&1 | tail -2
  git -C "$tmpdir" checkout FETCH_HEAD 2>&1 | tail -1
fi

# Copy and patch in a staging tree; the vendored files are only replaced once
# every patch has applied, so a failing patch leaves them untouched.
rel=$(git -C "$here" rev-parse --show-prefix)
stagedir=$(mktemp -d)
trap 'rm -rf "$tmpdir" "$stagedir"' EXIT
stage="$stagedir/$rel"
mkdir -p "$stage"
git init -q "$stagedir"

echo ">>> Copying whitelisted files..."
while IFS= read -r f; do
  [ -z "$f" ] && continue
  [[ "$f" =~ ^# ]] && continue
  if [ ! -f "$tmpdir/src/$f" ]; then
    echo "  WARN: $f not found in upstream; skipping"
    continue
  fi
  cp "$tmpdir/src/$f" "$stage/$f"
done < whitelist.txt

echo ">>> Applying local patches..."
# The patches are applied in name order to the copied upstream files, then the
# line endings are normalized to LF (upstream has none other since 3.09).
#
# All patches are `git diff` style with a/b/ prefixes (apply with -p1); a patch
# in `diff -urN` style with full paths (the former 01-namespace-wrap.patch) is
# applied with -p0. Auto-detect per-patch by inspecting the first line.
for p in "$here/patches/"*.patch; do
  [ -e "$p" ] || continue
  echo "  applying $p"
  if head -1 "$p" | grep -q "^diff --git a/"; then
    p_level=1
  else
    p_level=0
  fi
  if ! git -C "$stagedir" apply --whitespace=nowarn "-p$p_level" "$p" 2>&1; then
    echo "ERROR: patch $(basename "$p") did not apply; vendored files unchanged, upstream commit not recorded" >&2
    exit 1
  fi
done

echo ">>> Normalizing line endings (CRLF -> LF) post-patch..."
find "$stage" -maxdepth 1 -type f \( -name "*.h" -o -name "*.cpp" \) \
  -exec sed -i 's/\r$//' {} +

echo ">>> Installing patched files..."
cp "$stage"/* "$here/"

echo ">>> Recording upstream commit..."
commit_sha=$(git -C "$tmpdir" rev-parse HEAD)
echo "$commit_sha" > UPSTREAM_COMMIT
sed -i "s/^- \*\*Commit SHA\*\*:.*/- **Commit SHA**: $commit_sha/" UPSTREAM.md

echo ">>> Done. Review diff, then commit as:"
echo "   chore(percolator): sync from upstream $UPSTREAM_REF"

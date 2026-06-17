#!/usr/bin/env bash
# Re-sync vendored Percolator from upstream.
# Run from src/openms/thirdparty/percolator/.
set -euo pipefail

UPSTREAM_URL="https://github.com/percolator/percolator.git"
UPSTREAM_REF="${1:-eb157f74e963430e559e0d0bcd31291e4ad660ba}"

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

echo ">>> Copying whitelisted files..."
while IFS= read -r f; do
  [ -z "$f" ] && continue
  [[ "$f" =~ ^# ]] && continue
  if [ ! -f "$tmpdir/src/$f" ]; then
    echo "  WARN: $f not found in upstream; skipping"
    continue
  fi
  cp "$tmpdir/src/$f" "$here/$f"
done < whitelist.txt

echo ">>> Applying local patches..."
# The namespace-wrap patch was generated against upstream's CRLF-ending source,
# so apply FIRST (while copied files still have their original line endings),
# THEN normalize to LF for repo consistency.
#
# Patch format varies: 01-namespace-wrap.patch uses `diff -urN` with full
# paths (apply with -p0); 02-06 use `git format-patch` style with a/b/
# prefixes (apply with -p1). Auto-detect per-patch by inspecting first line.
repo_root=$(git -C "$here" rev-parse --show-toplevel)
for p in "$here/patches/"*.patch; do
  [ -e "$p" ] || continue
  echo "  applying $p"
  if head -1 "$p" | grep -q "^diff --git a/"; then
    p_level=1
  else
    p_level=0
  fi
  # --reject leaves partial failures visible in *.rej files; continue on fail.
  git -C "$repo_root" apply --reject --whitespace=nowarn "-p$p_level" "$p" 2>&1 || true
done

echo ">>> Normalizing line endings (CRLF -> LF) post-patch..."
find "$here" -maxdepth 1 -type f \( -name "*.h" -o -name "*.cpp" \) \
  -exec sed -i 's/\r$//' {} +

# Clean *.rej files; report any that remain.
rej_count=$(find "$here" -maxdepth 1 -name "*.rej" 2>/dev/null | wc -l)
if [ "$rej_count" -gt 0 ]; then
  echo "  WARN: $rej_count files had hunks that failed — see *.rej files"
  find "$here" -maxdepth 1 -name "*.rej" -printf "    %f\n"
fi

echo ">>> Recording upstream commit..."
commit_sha=$(git -C "$tmpdir" rev-parse HEAD)
echo "$commit_sha" > UPSTREAM_COMMIT
sed -i "s/^- \*\*Commit SHA\*\*:.*/- **Commit SHA**: $commit_sha/" UPSTREAM.md

echo ">>> Done. Review diff, then commit as:"
echo "   chore(percolator): sync from upstream $UPSTREAM_REF"

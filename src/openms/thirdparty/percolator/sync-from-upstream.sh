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
git clone --depth 1 --branch "$UPSTREAM_REF" "$UPSTREAM_URL" "$tmpdir" 2>&1 | tail -2

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
for p in patches/*.patch; do
  [ -e "$p" ] || continue
  echo "  applying $p"
  patch -p1 -d "$here" < "$p"
done

echo ">>> Recording upstream commit..."
commit_sha=$(git -C "$tmpdir" rev-parse HEAD)
echo "$commit_sha" > UPSTREAM_COMMIT
sed -i "s/^- \*\*Commit SHA\*\*:.*/- **Commit SHA**: $commit_sha/" UPSTREAM.md

echo ">>> Done. Review diff, then commit as:"
echo "   chore(percolator): sync from upstream $UPSTREAM_REF"

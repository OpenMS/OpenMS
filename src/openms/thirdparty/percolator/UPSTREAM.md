# Percolator vendored source tree

- **Upstream**: https://github.com/percolator/percolator
- **Pinned at**: commit eb157f7 (post-rel-3-08-01; PR #399 / I-spline PEP)
- **Commit SHA**: eb157f74e963430e559e0d0bcd31291e4ad660ba
- **License**: Apache 2.0 (see LICENSE-Apache-2.0.txt) plus the BSD 3-Clause licensed liblinear-derived TRON solver (see NOTICE-percolator.txt).

## What's here

A stripped-down copy of Percolator's `src/` tree, covering the PSM rescoring path only:
XML / tab / CLI / protein inference / peptide fragmentation code is excluded.

### Local replacement: TRON-based SVM solver

`ssl.{cpp,h}` and `tron.{cpp,h}` come from the
[`percolator-tron`](https://bitbucket.org/jthalloran/percolator_upgrade)
development branch (Halloran's TRON integration of liblinear v2.11). They
replace the legacy SVMlin-based L2-SVM-MFN solver from upstream
Percolator's `ssl.{cpp,h}`. The replacement:
- Sidesteps the SVMlin license question (its upstream headline is GPL v2+;
  the relicense to Apache-2.0 was scoped to Percolator itself, not to
  link-time inclusion in downstream non-Apache projects).
- Calls `extern "C"` BLAS (`dnrm2_`, `ddot_`, `daxpy_`, `dscal_`) which
  resolve at link time against OpenMS's existing libblas dependency.
- Is NOT in `whitelist.txt` — `sync-from-upstream.sh` won't overwrite
  these files. If a future upstream Percolator sync introduces a different
  ssl.cpp, this replacement will need to be re-applied manually.

## How to re-sync

Run `./sync-from-upstream.sh <tag-or-sha>`. The script copies files listed in
`whitelist.txt`, applies `patches/*.patch`, and records the new upstream commit
SHA in this file. Review the resulting diff and commit as
`chore(percolator): sync from upstream <ref>`.

A clean re-sync of the currently-pinned SHA reproduces the committed sources
bit-identically (verified 2026-10-09). If it does not, a downstream adaptation
was made by hand or a patch no longer applies — regenerate
`patches/01-namespace-wrap.patch` from the post-adaptation tree rather than
hand-editing sources and committing them, to keep sync self-healing.

## Patches

- `patches/01-namespace-wrap.patch` — covers all OpenMS adaptations except
  those in the later patches: namespace wrap in `OpenMS::Internal::Percolator`,
  `using namespace std;` placement, `std::` qualifications, additional
  includes (TabReader, Version.h), `::min/::max → std::min/std::max`,
  `parseOptions` body removal, enzyme / protein-inference / SQT feature
  drops, `#pragma once` on the new header-only regressors, etc.
  Regenerated 2026-04-24 (commit 9f227d7c03) from the current-tree-vs-upstream
  diff — replaces the former 01-06 chain which had accumulated drift.
  Regenerated again 2026-10-09: a hand edit had made it unparseable.
- `patches/02-atomic-include-negatives.patch` — makes the static
  `PosteriorEstimator::includeNegativesInResult` a `std::atomic<bool>`:
  `Scores::calcQvals` sets it from the OpenMP cross-validation loops, a data race
  on a plain `bool`.

### Regenerating the patch

Whenever you hand-fix a downstream issue in the vendored tree, capture it in
this patch rather than leaving it only in the source commit. The later patches
(`02-*.patch`, ...) are applied on top of patch 01 by `sync-from-upstream.sh`, so
they are reverse-applied on a staged copy of the tree first and stay out of
patch 01. To regenerate, run from this directory (it stops at the first error and
replaces patch 01 only when every step succeeded):

```bash
(
set -euo pipefail
new=patches/.01-namespace-wrap.patch.new
up=$(mktemp -d); stage=$(mktemp -d); trap 'rm -rf "$up" "$stage" "$new"' EXIT
whitelist=$(grep -v '^#' whitelist.txt | grep -v '^$' | LC_ALL=C sort)
for f in $whitelist; do
  gh api "repos/percolator/percolator/contents/src/$f?ref=$(cat UPSTREAM_COMMIT)" \
    --jq '.content' | base64 -d > "$up/$f"
  test -s "$up/$f"
done
mkdir -p "$stage/src/openms/thirdparty/percolator"
cp ./*.h ./*.cpp "$stage/src/openms/thirdparty/percolator/"
patches=$(ls -r patches/*.patch)  # an assignment, so set -e stops here if listing fails
for p in $patches; do
  [[ $p == patches/01-* ]] && continue
  pl=0; [[ $(head -1 "$p") == 'diff --git a/'* ]] && pl=1
  patch -s -R -p$pl -d "$stage" < "$p"
done
for f in $whitelist; do
  diff -urN --label "/tmp/percolator-preserved/$f" --label "src/openms/thirdparty/percolator/$f" \
    "$up/$f" "$stage/src/openms/thirdparty/percolator/$f" >> "$stage/01.patch" || [ $? -eq 1 ]
done
cp "$stage/01.patch" "$new"
mv "$new" patches/01-namespace-wrap.patch
)
```

Then check it: `./sync-from-upstream.sh "$(cat UPSTREAM_COMMIT)"` must exit 0 and
leave no change to the vendored `*.h`/`*.cpp` files. Pass the SHA: without an
argument the script syncs to its built-in default, which is not updated by a sync.

The `--label` flags are important — without them, `diff -urN` embeds
filesystem timestamps that make `git diff` noisy on every regeneration even
when the semantic patch content is unchanged.

## Licensing note

`ssl.{cpp,h}` and `tron.{cpp,h}` are BSD-3-Clause (liblinear v2.11 via
the percolator-tron fork). All other vendored files are Apache-2.0
(Percolator). See `NOTICE-percolator.txt` for the full attribution.

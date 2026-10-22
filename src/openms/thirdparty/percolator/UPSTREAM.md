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
- `patches/03-pep-padding-ties.patch` — keeps the PEP spline's points strictly
  increasing, which `BaseSpline::initiateQR` needs (it asserts `dx > 0` in Debug
  builds and divides by zero in Release builds, giving meaningless PEPs).
  `BaseSpline` keeps the distances in sparse vectors that treat |dx| <= 1e-15 as 0.
  The new `LogisticRegression::mergeTiesAndPad` (called from `setData`, after the
  log/logit transform is chosen from the real bin medians) merges points closer
  than 1e-12 (relative) in the transformed coordinates, and pads to the 4 points
  the spline needs, spaced like the data, in the same coordinates.
  `PosteriorEstimator::estimate` no longer pads the raw medians: upstream added at
  most 1e-20, lost to rounding unless the median is about 0, and padding raw
  scores in [0, 1] past 1 switched the spline from the logit to the log
  transform. Seen with ProSE `-rescore`. Still present in upstream master as of
  2026-10-12.
- `patches/04-cv-num-threads.patch` — gives the three OpenMP loops of
  `CrossValidation` a `num_threads` clause, so they run on at most `numThreads_`
  threads (OpenMS: `Percolator`'s `num_threads`) and never on more than OpenMP
  allows (a TOPP tool's `-threads`). Upstream sets the count with
  `omp_set_num_threads()` in its caller, which this tree does not include, so
  `numThreads_` was stored but unused.
  `numLoopThreads` is defined without OpenMP as well (it returns 1): Clang's
  `-fopenmp-simd` parses the `num_threads` clauses without defining `_OPENMP`.
- `patches/05-cv-folds-parallel.patch` — runs the remaining per-fold loops of
  cross-validation in parallel, on at most `numThreads_` threads like the loops of
  patch 04: the initial direction (`SanityCheck::getDefaultDirection`, which sorts
  each training set once per feature), the scoring of the test sets
  (`SanityCheck::getInitDirection`, `validateDirection`, `CrossValidation::train`)
  and of the training sets (`CrossValidation::doStep`, `mergeCpCnPairs`). Each fold
  works on its own `Scores`; the counts are added up in fold order, and the
  "Selected feature" reports (`Scores::getInitDirection` writes them to a stream
  it is given) are printed afterwards in fold order. None of these loops draws
  random numbers, so the results do not depend on the number of threads.
- `patches/06-packed-matrix-product-sparse.patch` — `PackedMatrix::packedMultiply`
  computes, per row, only the dot products with the columns that share an index
  with the row, instead of with every column: the others are 0 and are not
  stored. The products of the banded matrices of `BaseSpline` (the spline PEP,
  `pep_method` `nonparametric`) take linear instead of quadratic time, with the
  same result.

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
up=; stage=; new=; trap 'rm -rf "$up" "$stage" "$new"' EXIT
up=$(mktemp -d); stage=$(mktemp -d); new=$(mktemp -d patches/.regen.XXXXXX)
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
cp "$stage/01.patch" "$new/01.patch"
mv "$new/01.patch" patches/01-namespace-wrap.patch
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

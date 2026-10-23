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
  In `doStep`, the folds also fill their SVM inputs in parallel (not with nested
  cross-validation, which draws the nested sets at random), and
  `mergeCpCnPairs` counts the positives of all (Cpos, Cneg) pairs in parallel,
  each on its own copy of its test set: the count does not depend on the order
  of the set (PSMs that compare equal have the same label), and the test sets
  are discarded afterwards. Each of these loops catches the exceptions of its
  iterations and rethrows the first one, in iteration order, after the loop.
- `patches/06-spline-pep-sparse.patch` — speeds up the spline PEP (`BaseSpline`,
  `pep_method` `nonparametric`) without changing its result.
  `PackedMatrix::packedMultiply` computes, per row, only the dot products with
  the columns that share an index with the row, instead of with every column:
  the others are 0 and are not stored, so the products of the banded matrices
  take linear instead of quadratic time. `BaseSpline::solveInPlace` skips the
  rows whose index range does not include the column it looks for, and stops
  scanning a row at the first index past it (the indices are sorted and unique),
  instead of scanning every entry of every row.
- `patches/07-spline-pep-banded.patch` — speeds up the spline PEP further, again
  without changing its result. `BaseSpline::initiateQR` keeps `Q` and `R` as
  bands as well, and each iteration of `iterativeReweightedLeastSquares` runs
  `BaseSpline::bandedIteration`: the same floating-point operations, in the same
  order, as building `alpha / w`, `aWiQ`, `M = R + Qt * aWiQ`, `gamma = Qt * z`,
  `solveInPlace` and `gnew = z - aWiQ * gamma` with the sparse matrix code, on
  band rows in reused buffers and without allocations. The elimination makes the
  same pivot choices, row swaps and row operations as `solveInPlace` and visits
  only the rows that can hold the pivot column; the back substitution updates
  only the right-hand side (the row operations of that pass cannot change the
  result). Whenever dropping an entry (|value| <= 1e-15) would violate the band
  structure, or the elimination would leave it, `bandedIteration` returns before
  it commits the iteration's outputs and the matrix code runs instead. On FMA
  targets the compiler may contract `a * b + c` into a fused multiply-add; the
  band code keeps the expression shapes of the matrix code so that it is
  contracted the same way: the elimination stores its rounded product in a
  `volatile` (like `packedProd`, which stores it), and the dot products are
  loops `tot += a * b` like `packedDotProd` (under Clang with `clang loop
  unroll(disable) vectorize(disable)`: Clang's SLP vectoriser otherwise
  multiplied the terms of `Qt * z` separately). This is checked, not guaranteed:
  both paths ran in every IRLS iteration on 20 real and 300 synthetic inputs and
  gave the same bits with GCC 14.4 -O3 (no -march, `-march=haswell`,
  `-march=native` on Zen 4, `-march=haswell -ffp-contract=off`) and Clang 23 -O3
  (`-mfma` with the default `-ffp-contract=on`, `fast` and `off`, `-march=native
  -ffp-contract=fast`). Other compilers, versions and targets (AArch64, MSVC
  `/fp:contract`, x87, fast-math) are not covered.
- `patches/08-cv-sorts-tasks.patch` — speeds up cross-validation without changing
  its result. The sorts of `Scores` (`calcScores`, `generatePositiveTrainingSet`,
  `createXvalSetsBySpectrum`, `merge`, `postMergeStep`) call `sortLikeStdSort` from the
  OpenMS-local header `SortLikeStd.h` (not an upstream file, like `InProcessSampler.*`):
  inside an OpenMP parallel region it makes the comparisons and swaps of libstdc++'s
  `std::sort`, so it leaves the same permutation (also of PSMs that compare equal), but
  hands the right parts of the partitions to tasks (the algorithm of OpenMS'
  `FragmentIndex::sortPeptides_`); otherwise, and with a NaN score or mass (no strict weak
  ordering), it is `std::sort`.
  The loops over the folds (`CrossValidation`, `SanityCheck`) run their folds as tasks of a
  parallel region on up to 12 threads (`numTeamThreads`), whose other threads run the tasks
  of the sorts; `createXvalSetsBySpectrum` and `merge` are called in such a region as well.
  `Scores::getInitDirection` counts the positives of each feature in a task of its own, on
  (score, label) pairs sorted by score only (the q-values treat a group of equal scores as
  one, so the count does not depend on the order within a group), while the scores are
  sorted by each feature in turn as before: the order the last feature leaves is where the
  training starts from. The (Cpos, Cneg) pairs are trained as tasks, and each task then
  counts the positives of its pair on the pair's test set without sorting the set
  (`Scores::countPositives`, the count of `calcScoresAndQvals`, by the same argument), so
  the test sets are no longer copied. Nested cross-validation (which draws its nested sets at
  random) keeps its serial loop over the folds. Every loop, including those of
  `mergeCpCnPairs`, catches the exceptions of its iterations and rethrows the first one, in
  iteration order.
- `patches/10-cv-stages.patch` — speeds up the stages of cross-validation around the SVM
  training without changing its result (patch 09 is reserved for the SVM solver).
  `Scores::calcScores` scores the PSMs (each by `calcScore`, as before) in blocks of 2048, as
  tasks of the enclosing OpenMP parallel region if there is one (the fold tasks of
  `CrossValidation`), instead of one after the other.
  Each step of `CrossValidation::doStep` scores the training sets with the weights that the
  previous step's `mergeCpCnPairs` scored them with at its end; neither the weights nor the
  sets change in between, so the scores are the same: the step only sorts them again
  (`Scores::sortScoresAndCalcQvals`; the order is the same as after scoring them again) and
  sets the q-values. A flag of `CrossValidation` (`trainScoresOfWeights_`), set at the end of
  `mergeCpCnPairs` and cleared at the start of `mergeCpCnPairs`, `doStep`,
  `preIterationSetup` and `postIterationProcessing`, says when this holds.
  Without nested cross-validation, `doStep` trained and validated the SVMs of a fold on two
  copies of its training set; `Scores::generateTrainingSetsOfSortedSet` generates the SVM
  inputs from the set itself, without copying or sorting it, and the pairs are validated on
  the set itself (`countPositives` does not change it). The negatives are the decoys in their
  order. With `trainBestPositive` the copy was sorted by spectrum and label (`OrderScanLabel`),
  kept the first PSM of each group (`std::unique`) and sorted the kept PSMs by score again:
  the (file, scan, label) keys of the PSMs are sorted from the same order like `std::sort`
  (`sortLikeStdSort`, as before; the comparisons have the same outcomes, so the sort makes the
  same moves) to find the first PSM of each group, and if each positive (a kept target with
  q <= fdr) scores strictly higher (`operator>`) than the next one in the set's order, no two
  positives are equivalent and every sort by score leaves them in that order. Otherwise (and
  with a NaN score or mass) the copies as before.

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

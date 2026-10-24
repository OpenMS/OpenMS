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

OpenMS changes to the replacement (not in `patches/`: the sync neither copies
nor patches these files):
- Fast function object (OpenMS-local files `tron_kernels.h`,
  `tron_kernels_impl.h`, `tron_kernels_sse2.cpp`, `tron_kernels_avx2.cpp`,
  `cpu_avx2.h`). On x86-64 builds for the baseline instruction set (no AVX, no
  FMA: the build whose bits are the reference), `tron()` stores X in blocks of
  4 rows (`build_bt`) instead of a row-major Eigen matrix, and its function
  object (`FastTronFunction`) computes f, the gradient and Hv with SSE2 or AVX2
  kernels, chosen once at run time (`cpuHasAvx2`; `OPENMS_TRON_ISA=ref|sse2|avx2`
  overrides it, `tronSelectKernels` is the test hook). Every value is computed in
  the order of `TronFunction`: the margins in the order of Eigen's SSE2
  row-major gemv (even and odd columns in two lanes), the sums over rows in row
  order, every product rounded (the AVX2 file is compiled with `-mavx2 -mno-fma
  -ffp-contract=off`; compiled with FMA or without AVX2 it exports an empty table
  that is never selected); zero padding of rows and columns cannot change a sum
  that starts at +0, except in Hv with an infinite or NaN direction (0 * Inf),
  where `FastTronFunction::Hv` loops over the active rows like `TronFunction`.
  Further: the active rows are packed while `fun` has them in L1, f and g at
  w = 0 are computed once per solve (TRON and the gradient-ratio check ask three
  times), and the buffers are per thread and grow-only (freed at the end of a
  solve once they exceed 64 MB). The first use of each kernel table in a process
  solves two small problems both ways (about 1 ms) and falls back to
  `TronFunction` with a warning if any bit differs. Weights, outputs and the
  gradient ratio are bit-identical to `TronFunction` (class test
  `PercolatorTron_test`: whole solves for both kernels, n = 1..60, and fun /
  grad / Hv with infinite, NaN and huge values). Builds with AVX/FMA
  (`-march=...`), other architectures, problems with one row and problems with
  more than INT_MAX - 64 rows keep `TronFunction`. MSVC builds use the SSE2
  kernels only (no `/arch:AVX2`: its contraction is not verified). Checked with
  GCC 14.4 -O3 (conda-forge, Linux); MSVC and Apple Clang are covered only by
  the self-check. (With infinite or huge feature values, TRON's conjugate
  gradient does not terminate, in both function objects alike.)

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
  `CrossValidation`), instead of one after the other (the tasks share the weights).
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
  `Scores::getInitDirection` sorts the scores by each feature in turn (31 sorts per fold, from
  the order the previous sort left) on keys: the feature value, the rank of the PSM's (scan,
  mass, label) in the order `ScoreHolder::operator<` compares them after the score, and the
  position. The comparisons have the outcomes of `operator<`, so the sorts make the same
  moves, without reading `pPSM` and moving whole `ScoreHolder`s; the PSMs are put into the
  final order of the keys (and get the scores of the last feature) once, at the end. The
  positives of each feature are counted in a task on the (score, label) pairs in the order of
  its sort (sorted by score, as in the serial loop), instead of a task that sorted the pairs
  again; the chain hands the task a copy of its keys, the task builds the pairs. The ranks are
  sorted with `sortLikeStdSort` (any sort would do), the feature columns copied in parallel
  blocks.
  `PosteriorEstimator::estimatePi0` (called by `Scores::merge` and `postMergeStep` for each
  fold and for the merged set) runs its 100 bootstraps as tasks of the enclosing OpenMP
  parallel region (as in `CrossValidation::postIterationProcessing`). The bootstraps draw from
  the global `PseudoRandom` generator one after the other, min(n, 1000) numbers each, so each
  task draws its numbers from the generator's state at its first draw
  (`PseudoRandom::advance`: the generator multiplies by 279470273 modulo 4294967291, so it can
  jump ahead by a power of the multiplier; `PseudoRandom::lcg_rand(state)` draws like
  `lcg_rand` from a state of the caller's); the mean squared errors are added up in the order
  of the bootstraps, and the generator ends in the state the serial loop leaves. Outside a
  parallel region, the loop as before.
  `Scores::createXvalSetsBySpectrum` sorts the PSMs by spectrum hash on (hash, position) keys
  (the hashes `OrderScanHash` compares; same outcomes, same moves, see above) instead of
  hashing two PSMs for every comparison; inside an OpenMP parallel region it draws the folds
  (and RESET labels) in the order of the PSMs as before and then fills each fold's training
  and test set in a task of its own, in the same order. It then moved the feature rows of the
  test sets' PSMs to the start of the pool, one after the other by swaps (tracking the rows
  it displaced in a hash map); `FeatureMemoryPool::moveRowsInOrder` moves them there directly
  (copied out and back, in parallel blocks) when they are exactly the pool's rows (none free),
  each once, with the pool's row length: the swaps leave the same rows and pointers.
  Otherwise the swaps as before.
  `Scores::calcQvals` computes the q-values with `PosteriorEstimator::getQValuesInTasks` where
  it can: inside an OpenMP parallel region, for pi0 = 1 (no mix-max correction; the training
  sets of the cross-validation), a finite decoy factor >= 0 and at least 8192 PSMs. The
  PSMs are split into blocks of whole groups of equal scores at fixed positions; the targets
  and decoys of the blocks are counted first, so each block computes the FDRs of its groups
  from the same counts and by the same expression as `getQValues`; the q-values (minimum FDR
  from a position to the end, of numbers in [0, 1]) are the minima within each block and then
  with the minimum of the blocks after it. Otherwise `getQValues` as before.
  The sorts of `Scores` by score (`calcScores`, `merge`, `postMergeStep`) compare with
  `ScoreHolderGreater`, a copy of `operator>` (defined in another file, so each comparison was
  a function call) that the sort can inline: same outcomes, same moves.
- `patches/11-tdc-sort-tasks.patch` — `Scores::weedOutRedundantTDC` sorts the PSMs by
  spectrum (`OrderScanMassCharge`) with `sortLikeStdSort` (`SortLikeStd.h`) instead of
  `std::sort`, unless a score or mass is NaN (no strict weak ordering): the same permutation,
  with the tasks of an enclosing OpenMP parallel region. OpenMS' `Percolator::rescore` calls it
  in such a region (as the cross-validation calls `merge`), so that this sort, the sort of
  `postMergeStep` and the bootstraps of `estimatePi0` run as tasks; outside a region it is
  `std::sort` as before.
- `patches/12-cv-fold-chains.patch` — speeds up the training of the cross-validation without
  changing its result. Without quick validation (its stop test needs the sum over the folds after
  each step), nested cross-validation (which draws the nested sets at random) and the reports of
  each step's training sets and weights (`VERB > 2`, `reportPerformanceEachIteration`), no fold
  depends on another fold until the end of `CrossValidation::train`. `trainFoldChains` then runs
  all steps of a fold (scoring its training set, its SVM inputs, its (Cpos, Cneg) pairs as tasks,
  the choice of its best pair as in `mergeCpCnPairs`, the scoring of the training set with the
  chosen weights) and the scoring of its test set as one task, so that the folds no longer wait for
  each other at every stage of every step; the counts are added up in fold order afterwards and the
  reports printed as before. An error ends its fold's task; the error `train` would have thrown (the
  first by step, then by stage, then by fold) is rethrown after the reports of the steps before
  it. `CrossValidation::setupFoldInputs` is `doStep`'s setup of a fold's SVM inputs, now shared by
  both. Otherwise `train` steps through `doStep` as before.
  The (Cpos, Cneg) pairs of a fold count their positives (`Scores::countPositives`) with less
  work, to the same count. The fold's feature rows with their labels are collected once if
  the fold validates at least 3 pairs (`Scores::getFeatureRows`, sorted by address, so the scoring reads them in the order they are in
  memory; the training set keeps its PSMs and labels during the training), and the scores (each by
  `calcScore`, as before) are counted from the PSMs with the highest scores only where that gives
  the same count (`Scores::countPositivesOfTop`): with pi0 >= 1 and a finite decoy factor >= 0,
  `getQValues` makes no mix-max correction, the positives are the targets up to the end of the
  last group of equal scores whose min(FDR, 1) is below the threshold, and every group after the
  PSMs scoring at least as high as the k-th highest (k = n / 8, then 4 times more) has an FDR of
  at least (their decoys, all targets): each operation of the FDR's expression is monotone. If
  that bound is not below the threshold (checked before sorting), only those PSMs are sorted and
  walked; otherwise (also for a threshold above 1) all are sorted as before.
  The set-up of the folds (`preIterationSetup`: `createXvalSetsBySpectrum` and the initial
  directions) runs on up to 16 threads instead of 12 (`numSetupTeamThreads`).
- `patches/14-rowsscale-postprocessing.patch` — speeds up the steps before and after the training
  of the cross-validation on large sets, without changing their results.
  `CrossValidation::postIterationProcessing`: when the sanity check rejects the learned weights,
  the test sets are scored with the default weights as tasks of a parallel region (as
  `SanityCheck::validateDirection` does), so their scoring, sorts and q-values run as the region's
  tasks (`Scores::calcScores`, `SortLikeStd.h`, `PosteriorEstimator::getQValuesInTasks`): the same
  scores, order and q-values as the loop over the sets.

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

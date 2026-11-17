# Percolator vendored source tree

- **Upstream**: https://github.com/percolator/percolator
- **Pinned at**: commit a90106b (master after PR #405, "round three"; version 3.09.0 plus the unreleased changes of
  `src/ChangeLog`)
- **Commit SHA**: a90106b40928c01097f5bf0ddc09651877303241
- **License**: Apache 2.0 (see LICENSE-Apache-2.0.txt) plus the BSD 3-Clause licensed liblinear-derived TRON solver (see NOTICE-percolator.txt).

## What's here

A stripped-down copy of Percolator's `src/` tree, covering the PSM rescoring path only:
the CLI (`Caller`), multi-file input, protein inference, enzymes and the SVMlin solver are excluded.
`whitelist.txt` lists the upstream files; everything else here is OpenMS-local.

### Local replacement: TRON-based SVM solver

`ssl.{cpp,h}` and `tron.{cpp,h}` come from the
[`percolator-tron`](https://bitbucket.org/jthalloran/percolator_upgrade)
development branch (Halloran's TRON integration of liblinear v2.11). They
replace the SVMlin-based L2-SVM-MFN solver of upstream Percolator's `ssl.{cpp,h}`
(upstream still ships SVMlin there, with exact least-squares steps since 3.09). The replacement:
- Sidesteps the SVMlin license question (its upstream headline is GPL v2+;
  the relicense to Apache-2.0 was scoped to Percolator itself, not to
  link-time inclusion in downstream non-Apache projects). Nothing in this
  directory is copied from upstream's `ssl.{cpp,h}`.
- Is NOT in `whitelist.txt` — `sync-from-upstream.sh` won't overwrite
  these files.

`ssl.h` provides the interface that upstream's training code (CrossValidation,
`Scores::generate*TrainingSet`; Apache-2.0) calls, as an independent implementation:
`AlgIn` (the training set: owned copies of the rows, or rows in place; `row`, `rows`,
`resizeExamples`, `set`, `add`, `reserve`, `capacity`, `clear`), `options` (with upstream's
`directSolve`, ignored by TRON, and OpenMS' `weightsOnly`), an empty `SolverWorkspace` (TRON
pools its buffers itself), `prefetchRow` / `prefetchAddress`, and `L2_SVM_MFN` overloads that
forward to `tron()`. TRON starts from the weights it is given (a warm start, as upstream's
solver since 3.09; the outputs of the warm start are not read) and returns zero weights for an
empty training set, as upstream's solver does. The solvers read the n - 1 features of each row
and supply the bias column themselves, so rows in place need not end with 1.0 for them.

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
  The active rows are packed while `fun` has them in L1, f and g at w = 0 are computed once per
  solve, and the buffers are grow-only and come from a pool of workspaces that solves take and
  give back (a workspace above 64 MB is freed at the end of its solve, unless a
  `TronWorkspaceHold` exists: `CrossValidation::train` holds one, patch 05). The first use of each
  kernel table in a process solves two small problems both ways (about 1 ms) and falls back to
  `TronFunction` with a warning if any bit differs. Weights, outputs and the gradient ratio are
  bit-identical to `TronFunction` (class test `PercolatorTron_test`). Builds with AVX/FMA
  (`-march=...`), other architectures, problems with one row and problems with more than
  INT_MAX / 2 rows keep `TronFunction`. MSVC builds use the SSE2 kernels only. With
  `options::weightsOnly` tron() stops after TRON (no Outputs, no convergence check; returns 1);
  the cross-validation's solves ask for that (patch 05).
- `TronPairGroup`: solves on the same set that pass the same group share one BT4 copy of X, built
  once by the solves that start while it is built (bit-identical to solves alone). The
  cross-validation makes one group per training set and pass while several (Cpos, Cneg) pairs train
  on a set at once: in the first three iterations, until each fold's pair is fixed, and in every
  iteration of `classic_training` (patch 05). The memory of the pair solves is then one copy of X
  per training set, whatever the number of threads; the later single solves build their own.

### Local addition: Newton and LDA solvers (`newton_svm.{h,cpp}`)

OpenMS-local, not in `whitelist.txt`, used only when the cross-validation's solver is
not TRON (`CrossValidation::setSolver`, OpenMS `Percolator` parameter `solver`, ProSE's
`rescore_solver`; patch 05). `newtonSvm` minimizes the same objective as `tron()` with the modified finite
Newton method of Keerthi & DeCoste (2005, JMLR 6:341-361; L2-SVM-MFN), written from
the paper (not from SVMlin, see the licensing note above): per step one pass over the
active rows forms the explicit Hessian `I + 2 sum_I C_i x_i x_i^T` and the right-hand
side (rank-4 updates), a Cholesky factorization gives the Newton point, one pass over
all rows gives its margins, and the exact line search on the breakpoints of the
piecewise quadratic objective (selection by partitioning) gives the step. This is the
method of upstream's solver since 3.09 (exact least-squares steps, warm starts). It stops on TRON's criterion
(`||g(w)|| <= eps ||g(0)||`) and takes at most 128 columns. `ldaSvm` is Fisher's
discriminant (pooled within-class scatter, small ridge) from one pass, scaled and
shifted along its direction as in the paper's beta-seeding (eq. 8/9). Rows are cut
into fixed blocks of 8192 whose partial sums are added in block order, so the results
do not depend on the threads; within a parallel region the blocks of problems of at
least 8 blocks run as tasks (`OpenMPTasks.h`). Class test `PercolatorNewton_test`.

### Local additions: `OpenMPThreads.h`, `OpenMPTasks.h`, `InProcessSampler.*`, stubs

- `OpenMPThreads.h`: the thread cap of the vendored code (`maxThreads()`, `ThreadCap`, patch 04),
  `LoopErrors` (patch 08) and `TaskLoopScope` / `taskLoopThreads()` (the taskloops of the fold
  chains, patch 11).
- `OpenMPTasks.h`: `PERCOLATOR_OMP_TASKS` / `PERCOLATOR_TASK_PRAGMA` (OpenMP tasks only where the
  compiler has them; MSVC's OpenMP 2.0 has none). Used by `newton_svm.cpp` and the fold chains
  (patch 11); upstream uses no tasks.
- `InProcessSampler.{h,cpp}`: the reservoir and out-of-fold samples of OpenMS' capped training.
- `Globals.h`, `Version.h`, `config.h`: stubs for upstream's CMake-generated headers
  (`Globals.h` is upstream's `Globals.h.cmake` with the values of its `CommonCMake.txt`).

## How to re-sync

Run `./sync-from-upstream.sh <tag-or-sha>`. The script copies files listed in
`whitelist.txt`, applies `patches/*.patch` in name order, and records the new upstream commit
SHA in `UPSTREAM_COMMIT` and in this file. Review the resulting diff and commit as
`chore(percolator): sync from upstream <ref>`.

A clean re-sync of the pinned SHA reproduces the committed sources
bit-identically (verified 2026-11-17). If it does not, a downstream adaptation
was made by hand or a patch no longer applies — regenerate the patches (below)
rather than hand-editing sources and committing them, to keep sync self-healing.

## Patches

All patches are `git diff` style (a/b/ prefixes, applied with -p1) on the whitelisted files.

- `patches/01-namespace-wrap.patch` — the OpenMS adaptations of the upstream files: namespace
  `OpenMS::Internal::Percolator` around every file (opened after the includes, with the file's
  `using namespace std;` inside), `std::min`/`std::max` in `Vector.cpp`, includes of excluded headers
  removed (`Enzyme.h`, `Option.h`, `ProteinProbEstimator.h`, `TabFileValidator.h`, `ValidateTabFile.h`:
  `SetHandler.h` defines the `InputSegment` of the multi-file reader itself), `TabReader.h` and
  `<algorithm>` included where those had brought them, the protein-inference branch of
  `DataSet::readPsm` disabled, and `PosteriorEstimator::parseOptions` (qvality's command line)
  removed.
- `patches/02-pep-padding-ties.patch` (formerly 03) — keeps the PEP spline's points strictly
  increasing, which `BaseSpline::initiateQR` needs (it asserts `dx > 0` in Debug builds and divides by
  zero in Release builds, giving meaningless PEPs). The new `LogisticRegression::mergeTiesAndPad`
  (called from `setData`, after the log/logit transform is chosen from the real bin medians) merges
  points closer than 1e-12 (relative) in the transformed coordinates, and pads to the 4 points the
  spline needs, spaced like the data, in the same coordinates. `PosteriorEstimator::estimate` no longer
  pads the raw medians: upstream added at most 1e-20, lost to rounding unless the median is about 0,
  and padding raw scores in [0, 1] past 1 switched the spline from the logit to the log transform.
  Seen with ProSE `-rescore`. Still present in upstream a90106b.
- `patches/03-spline-fill-g.patch` (formerly 35) — `BaseSpline::iterativeReweightedLeastSquares`:
  after each step, `gnew` gets its dropped entries back as zeros (`fillPacked`, before `limitg`). The
  packed vector operations leave out entries of (almost) zero, but `evaluateSlope`, `splineEval` and
  `LogisticRegression::calcPZW` / `limitg` read g by position: a dropped entry shifted the later
  values, and an all-zero g (no entries) was read out of bounds (OpenMS issue #10880). A filled `gnew`
  is unchanged. Still present in upstream a90106b.
- `patches/04-omp-thread-cap.patch` (formerly 04) — the parallel regions of `CrossValidation.cpp`,
  `Scores.cpp` and `ParallelSort.h` size their teams by `maxThreads()` (`OpenMPThreads.h`):
  `omp_get_max_threads()`, at most the cap of the calling thread, which OpenMS' wrapper sets per call
  (`ThreadCap`: the `Percolator` parameter `num_threads`, at most `SysInfo::ompThreadsOfThisThread()`).
  Upstream sizes them by `omp_get_max_threads()`, which its executable sets once; OpenMS must not call
  `omp_set_num_threads()` (MSVC's vcomp keeps one number for the whole process). `inParallel()`
  replaces `omp_in_parallel()` in a clause, which did not compile with Clang's `-fopenmp-simd` (no
  `_OPENMP`).
- `patches/05-cv-solver-options.patch` (formerly 14, 15, 16 and the weights-only part of 18) —
  `CrossValidation::setSolver` (TRON, Newton or LDA; `trainSvm` dispatches every solve, Newton and LDA
  fall back to TRON from 0 where they cannot solve), `setTrainingOptions` (TRON's tolerance, a warm
  start of TRON with classic training, a single Cneg candidate), weights-only solves, a
  `TronWorkspaceHold` during `train()`, one `TronPairGroup` per training set for the solves of its
  pairs (formerly 18: one copy of X per set, not per solve and thread), and no outputs X w for the
  warm start (TRON does not read them). Newton stops at 1e-4 in the iterations before the last one
  of classic training without quick validation (which can stop after any iteration). With the
  defaults the training is upstream's.
- `patches/06-cv-fixed-folds.patch` (formerly 19) — out-of-fold capped training (OpenMS `subset_mode`
  `out_of_fold`, ProSE's default): `Scores::createXvalSetsBySpectrum` takes an optional function that
  gives each PSM's fold (no fold is drawn, the empty-fold rule does not apply); `setFixedFolds` passes
  it, nested sets are drawn per fold from a generator state of their own, and each fold's Cneg
  candidates come from its own training set. With fixed folds each fold stops by upstream's
  convergence rule on its own count: the sum over the folds that upstream uses would let a fold's
  held-out labels (in the other folds' training sets) decide when its model stops training.
  `setActiveFolds` (the learning curve's stopped folds), `getTrainingPsms` (before
  `postIterationProcessing`, which releases the training sets) and `Scores::setScoreHolders`.
- `patches/07-sanity-check-virtual.patch` — `SanityCheck::validateDirection` is virtual, so that the
  wrapper can check each fold's model on its own training set (out-of-fold training).
- `patches/08-omp-loop-errors.patch` — the parallel loops of the cross-validation and of
  `Scores::bestInitDirection` catch the exceptions of their iterations and rethrow the first after the
  loop (`LoopErrors`, `OpenMPThreads.h`): an exception that leaves an OpenMP region ends the process,
  and OpenMS runs the code as a library (e.g. a `std::bad_alloc` of a large run).
- `patches/09-spline-pep-sparse.patch` (formerly 06) — speeds up the spline PEP (`BaseSpline`,
  `pep_method` `nonparametric`) without changing its result. `PackedMatrix::packedMultiply` computes,
  per row, only the dot products with the columns that share an index with the row, instead of with
  every column: the others are 0 and are not stored, so the products of the banded matrices take
  linear instead of quadratic time. `BaseSpline::solveInPlace` skips the rows whose index range does
  not include the column it looks for, and stops scanning a row at the first index past it (the
  indices are sorted and unique).
- `patches/10-spline-pep-banded.patch` (formerly 07) — speeds up the spline PEP further, again without
  changing its result. `BaseSpline::initiateQR` keeps `Q` and `R` as bands as well, and each iteration
  of `iterativeReweightedLeastSquares` runs `BaseSpline::bandedIteration`: the same floating-point
  operations, in the same order, as building `alpha / w`, `aWiQ`, `M = R + Qt * aWiQ`, `gamma = Qt * z`,
  `solveInPlace` and `gnew = z - aWiQ * gamma` with the sparse matrix code, on band rows in reused
  buffers and without allocations. The elimination makes the same pivot choices, row swaps and row
  operations as `solveInPlace` and visits only the rows that can hold the pivot column. Whenever
  dropping an entry (|value| <= 1e-15) would violate the band structure, or the elimination would
  leave it, `bandedIteration` returns before it commits the iteration's outputs and the matrix code
  runs instead. On FMA targets the compiler may contract `a * b + c` into a fused multiply-add; the
  band code keeps the expression shapes of the matrix code so that it is contracted the same way: the
  elimination stores its rounded product in a `volatile` (like `packedProd`), and the dot products
  are loops `tot += a * b` like `packedDotProd` (under Clang with `clang loop unroll(disable)
  vectorize(disable)`: Clang's SLP vectoriser otherwise multiplied the terms of `Qt * z` separately).
  This is checked, not guaranteed: on the former tree (eb157f7) both paths ran in every IRLS
  iteration on 20 real and 300 synthetic inputs and gave the same bits with GCC 14.4 -O3 (no -march,
  `-march=haswell`, `-march=native` on Zen 4, `-march=haswell -ffp-contract=off`) and Clang 23 -O3
  (`-mfma` with the default `-ffp-contract=on`, `fast` and `off`, `-march=native -ffp-contract=fast`);
  on a90106b the `nonparametric` PEPs of `perc_bench` (100,000 PSMs, GCC 13.3 -O3) are bit-identical
  with and without patches 09 and 10. Other compilers, versions and targets (AArch64, MSVC
  `/fp:contract`, x87, fast-math) are not covered. (`isFilled`, which patch 10 uses, comes with
  patch 03.)
- `patches/11-cv-fold-chains.patch` (the fold chains of the former 12 and 33, on upstream's new
  training) — `CrossValidation::trainFoldChains`: the iterations of `train()` as one chain of OpenMP
  tasks per fold. A step of a fold is what `doStep` and `mergeCpCnPairs` do for it, in the same
  order: its training set scored with its weights, its SVM input, one task per (Cpos, Cneg) pair
  that solves the pair and counts its positives right after, the choice of the best pair and the
  count of the fold's positives. The last pair of a step runs the rest of it and starts the fold's
  next step, so a fold waits only for its own pairs, and the threads run whatever task is ready.
  No task waits for other tasks of the chains: a task waits only at the end of the taskloops it runs
  (the row loops of Scores, Newton's blocks), and the solves of a step's pairs while they build their
  shared copy of the SVM input (TronPairGroup); the per-thread buffers and `TaskLoopScope` rely on the
  tasks being tied. With fixed folds (each fold stops on its own count) and with classic
  training the chains never wait for each other; with upstream's stop rule on the sum of the folds
  the last fold to finish an iteration decides and starts the next steps. Inside the chains the
  scoring of a set (`Scores::calcScores`) and the copies of its training rows (`appendTrainingRows`)
  run as taskloops (`TaskLoopScope`, `OpenMPThreads.h`). The reports of the iterations are printed
  afterwards as the loop prints them; after an error no fold starts a step beyond the iteration of
  the first error, which is rethrown after the reports of the iterations before it. `trainsFoldChains()`: with OpenMP tasks and more than one
  thread, without nested cross-validation, `-R` or verbosity above 2, with all folds active and held
  at once; otherwise the loop runs. The results are those of the loop, bit for bit (scores, q-values
  and PEPs of `perc_bench` on 100,000 and 1,000,000 synthetic PSMs at 1, 2 and 4 threads: TRON,
  Newton and LDA, classic and new training, c_pos 0 and 0.1, `out_of_fold`, `out_of_fold_avg`,
  in-sample caps, learning curve, nested cross-validation, warm starts). Measured with `perc_bench`
  (1,000,000 synthetic PSMs, 30 features, 3 interleaved repetitions on 4 cores): 7-18% less time than
  the loop at 4 threads and 11-17% at 2 threads (100,000 PSMs: 22-29% at 4 threads); from 1 to 4
  threads 2.5-2.9x instead of 2.1-2.7x.

### The former patch set (01-35 on eb157f7)

Most of the former patches were speedups of the cross-validation, whose code upstream rewrote between
eb157f7 and a90106b (258 commits: total-order comparators and `parallelSort`, count-only grid evaluation,
no re-sorting of sorted sets, the SVM input built on all threads, PSMs in a block store, the I-spline fit
without a design matrix, and a training that fixes each fold's (Cpos, Cneg) pair after three iterations
and stops at a gain below 0.1%). Measured on synthetic PINs (`perc_bench`, 4 threads; time of
`Percolator::rescore`, peak RSS), the new tree with upstream's defaults is 3-4x faster than the former
tree with all 35 patches, with less memory, so they were not ported:

- superseded by upstream: 02 (the static flag is gone), 05, 08, 10 (folds and sorts), 09, 11, 12 (stages;
  its fold chains are re-implemented as 11), 13, 17, 28, 30 (initial direction, normalization), 20 and 32 (the
  best-scoring positive per spectrum: `OrderScanLabel` orders by score since 25cb5124; issue #10720),
  21 (PSM ownership: DataSet no longer owns its PSMs, the feature pool no longer zero-fills), 22-26
  (sort certificates, packed keys, set-up in tasks), 33 (OpenMP tasks: upstream uses none; the fold
  chains as 11), 34
  (comments of removed code);
- not ported: 27 (reuse of the validation scores), 29, 31 (code gone);
- ported: 01, 03, 04, 14, 15, 16, 18, 19, 35, 06, 07 (as 01-06, 09 and 10 above).

`SortLikeStd.h` (sorts with the moves of `std::sort`) is removed: with total-order comparators any
sort gives the same order.

### Regenerating the patches

The patches are the commits of a port branch: a commit with the pristine upstream files of the pinned
SHA (`whitelist.txt`) and the OpenMS-local files, then one commit per patch that touches only
whitelisted files (OpenMS-local files may change in separate commits). For each patch commit:

```bash
files=$(grep -v '^#' whitelist.txt | grep -v '^$' | sed 's#^#src/openms/thirdparty/percolator/#')
git diff --no-color --full-index <commit>~1 <commit> -- $files > patches/NN-name.patch
```

To rebuild such a branch from this directory: copy the whitelisted files of the pinned SHA into a fresh
repository, commit them with the OpenMS-local files, and apply and commit the patches one by one.
Then check: `./sync-from-upstream.sh "$(cat UPSTREAM_COMMIT)"` must exit 0 and leave no change to the
vendored `*.h`/`*.cpp` files. Pass the SHA: without an argument the script syncs to its built-in default.

## Licensing note

`ssl.{cpp,h}` and `tron.{cpp,h}` are BSD-3-Clause (liblinear v2.11 via
the percolator-tron fork); the `AlgIn` interface in `ssl.h` is OpenMS' own implementation.
All other vendored files are Apache-2.0 (Percolator). See `NOTICE-percolator.txt` for the full attribution.

# Percolator vendored source tree

- **Upstream**: https://github.com/percolator/percolator
- **Pinned at**: commit a9783ee8 (master of 2026-11-17, after PRs #407 and #408, the P-IRLS PEPs; version 3.09.0
  plus the unreleased changes of `src/ChangeLog`)
- **Commit SHA**: a9783ee8ded0af8148566bdff4d43d2ba444ce5e
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
  `TronWorkspaceHold` exists: `CrossValidation::train` holds one, patch 03). The first use of each
  kernel table in a process solves two small problems both ways (about 1 ms) and falls back to
  `TronFunction` with a warning if any bit differs. Weights, outputs and the gradient ratio are
  bit-identical to `TronFunction` (class test `PercolatorTron_test`). Builds with AVX/FMA
  (`-march=...`), other architectures, problems with one row and problems with more than
  INT_MAX / 2 rows keep `TronFunction`. MSVC builds use the SSE2 kernels only. With
  `options::weightsOnly` tron() stops after TRON (no Outputs, no convergence check; returns 1);
  the cross-validation's solves ask for that (patch 03).
- `TronPairGroup`: solves on the same set that pass the same group share one BT4 copy of X, built
  once by the solves that start while it is built (bit-identical to solves alone). The
  cross-validation makes one group per training set and pass while several (Cpos, Cneg) pairs train
  on a set at once: in the first three iterations, until each fold's pair is fixed, and in every
  iteration of `classic_training` (patch 03). The memory of the pair solves is then one copy of X
  per training set, whatever the number of threads; the later single solves build their own.

### Local addition: Newton and LDA solvers (`newton_svm.{h,cpp}`)

OpenMS-local, not in `whitelist.txt`, used when the cross-validation's solver is
not TRON (the default of OpenMS' `Percolator` and ProSE) (`CrossValidation::setSolver`, OpenMS `Percolator` parameter `solver`, ProSE's
`rescore_solver`; patch 03). `newtonSvm` minimizes the same objective as `tron()` with the modified finite
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

- `OpenMPThreads.h`: the thread cap of the vendored code (`maxThreads()`, `ThreadCap`, patch 02),
  `LoopErrors` (patch 06) and `TaskLoopScope` / `taskLoopThreads()` (the taskloops of the fold
  chains, patch 07).
- `OpenMPTasks.h`: `PERCOLATOR_OMP_TASKS` / `PERCOLATOR_TASK_PRAGMA` (OpenMP tasks only where the
  compiler has them; MSVC's OpenMP 2.0 has none). Used by `newton_svm.cpp` and the fold chains
  (patch 07); upstream uses no tasks.
- `InProcessSampler.{h,cpp}`: the reservoir and out-of-fold samples of OpenMS' capped training.
- `Globals.h`, `Version.h`, `config.h`: stubs for upstream's CMake-generated headers
  (`Globals.h` is upstream's `Globals.h.cmake` with the values of its `CommonCMake.txt`).

## How to re-sync

Run `./sync-from-upstream.sh <tag-or-sha>`. The script copies files listed in
`whitelist.txt`, applies `patches/*.patch` in name order, and records the new upstream commit
SHA in `UPSTREAM_COMMIT` and in this file. Review the resulting diff and commit as
`chore(percolator): sync from upstream <ref>`.

A clean re-sync of the pinned SHA reproduces the committed sources
bit-identically (verified 2026-11-18). If it does not, a downstream adaptation
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
- `patches/02-omp-thread-cap.patch` (former set: 04) — the parallel regions of `CrossValidation.cpp`,
  `Scores.cpp` and `ParallelSort.h` size their teams by `maxThreads()` (`OpenMPThreads.h`):
  `omp_get_max_threads()`, at most the cap of the calling thread, which OpenMS' wrapper sets per call
  (`ThreadCap`: the `Percolator` parameter `num_threads`, at most `SysInfo::ompThreadsOfThisThread()`).
  Upstream sizes them by `omp_get_max_threads()`, which its executable sets once; OpenMS must not call
  `omp_set_num_threads()` (MSVC's vcomp keeps one number for the whole process). `inParallel()`
  replaces `omp_in_parallel()` in a clause, which did not compile with Clang's `-fopenmp-simd` (no
  `_OPENMP`).
- `patches/03-cv-solver-options.patch` (former set: 14, 15, 16, 18) —
  `CrossValidation::setSolver` (TRON, Newton or LDA; `trainSvm` dispatches every solve, Newton and LDA
  fall back to TRON from 0 where they cannot solve), `setTrainingOptions` (TRON's tolerance, a warm
  start of TRON with classic training, a single Cneg candidate), weights-only solves, a
  `TronWorkspaceHold` during `train()`, one `TronPairGroup` per training set for the solves of its
  pairs (formerly 18: one copy of X per set, not per solve and thread), and no outputs X w for the
  warm start (TRON does not read them). Newton stops at 1e-4 in the iterations before the last one
  of classic training without quick validation (which can stop after any iteration). With the
  defaults of `CrossValidation` the training is upstream's (OpenMS' `Percolator` sets its own).
- `patches/04-cv-fixed-folds.patch` (former set: 19) — out-of-fold capped training (OpenMS `subset_mode`
  `out_of_fold`, ProSE's default): `Scores::createXvalSetsBySpectrum` takes an optional function that
  gives each PSM's fold (no fold is drawn, the empty-fold rule does not apply); `setFixedFolds` passes
  it, nested sets are drawn per fold from a generator state of their own, and each fold's Cneg
  candidates come from its own training set. With fixed folds, where the training stops early (new
  training or quick validation), each fold stops by upstream's convergence rule on its own count: the sum over the folds that upstream uses would let a fold's
  held-out labels (in the other folds' training sets) decide when its model stops training.
  `setActiveFolds` (the learning curve's stopped folds), `getTrainingPsms` (before
  `postIterationProcessing`, which releases the training sets) and `Scores::setScoreHolders`.
- `patches/05-sanity-check-virtual.patch` — `SanityCheck::validateDirection` is virtual, so that the
  wrapper can check each fold's model on its own training set (out-of-fold training).
- `patches/06-omp-loop-errors.patch` — the parallel loops of the cross-validation and of
  `Scores::bestInitDirection` catch the exceptions of their iterations and rethrow the first after the
  loop (`LoopErrors`, `OpenMPThreads.h`): an exception that leaves an OpenMP region ends the process,
  and OpenMS runs the code as a library (e.g. a `std::bad_alloc` of a large run).

- `patches/07-cv-fold-chains.patch` (the fold chains of the former set's 11, 12 and 33, on upstream's
  new training) — `CrossValidation::trainFoldChains`: the iterations of `train()` as one chain of OpenMP
  tasks per fold. A step of a fold is what `doStep` and `mergeCpCnPairs` do for it, in the same
  order: its training set scored with its weights, its SVM input, one task per (Cpos, Cneg) pair
  that solves the pair and counts its positives right after, the choice of the best pair and the
  count of the fold's positives. The last pair of a step runs the rest of it and starts the fold's
  next step, so a fold waits only for its own pairs, and the threads run whatever task is ready.
  No task waits for other tasks of the chains: a task waits only at the end of the taskloops it runs
  (the row loops of Scores, Newton's blocks), and the solves of a step's pairs while they build their
  shared copy of the SVM input (TronPairGroup); the per-thread buffers and `TaskLoopScope` rely on the
  tasks being tied. With fixed folds (each fold stops on its own count) and with classic
  training (no early stop without quick validation) the chains never wait for each other; with upstream's stop rule on the sum of the folds
  the last fold to finish an iteration decides and starts the next steps. Inside the chains the
  scoring of a set (`Scores::calcScores`) and the copies of its training rows (`appendTrainingRows`)
  run as taskloops (`TaskLoopScope`, `OpenMPThreads.h`). The reports of the iterations are printed
  afterwards as the loop prints them; after an error no fold starts a step beyond the iteration of
  the first error, which is rethrown after the reports of the iterations before it. The first error is the loop's:
  by iteration, then by stage in the order in which the loop runs the stages for all folds (e.g. the solves of all
  pairs before the feature rows and the counts), then fold and pair. The per-iteration bookkeeping grows with the
  iterations that run, not with the iteration cap. A step that a step starts is a task, which a runtime may run at once,
  nested in the starting task (e.g. libomp in a team of one thread): beyond four steps nested on a thread's stack the
  step waits in the thread's list until the thread's outermost step has finished its own, so the stack does not grow
  with the iterations. A fold's feature rows for the counts of its pairs are released when its grid is fixed.
  `trainsFoldChains()`: with OpenMP tasks and more than one thread, without nested cross-validation, `-R` or verbosity above 2, with all folds active and held
  at once; otherwise the loop runs. The results are those of the loop, bit for bit (scores, q-values
  and PEPs of `perc_bench` on 100,000 and 1,000,000 synthetic PSMs at 1, 2 and 4 threads: TRON,
  Newton and LDA, classic and new training, c_pos 0 and 0.1, `out_of_fold`, `out_of_fold_avg`,
  in-sample caps, learning curve, nested cross-validation, warm starts; on a90106b, whose training code
  is that of a9783ee8, and with OpenMS' defaults on a9783ee8). Measured with `perc_bench` (synthetic
  PSMs, 30 features, 3 interleaved repetitions on 4 cores): with OpenMS' defaults 16-17% less time
  than the loop at 4 threads on 1,000,000 PSMs and 37% on 100,000; in other configurations 7-29% at 4
  threads and 11-17% at 2; from 1 to 4 threads 2.5-2.9x instead of 2.1-2.7x.

Deliberate difference from the executable outside the patches (OpenMS wrapper, `Percolator.cpp`): with
`train_fdr` 0 the wrapper uses `test_fdr` also in the passes after a capped training (the `-N` second pass,
the out-of-fold merge, the learning curve), as the cross-validation resolves its selection FDR. The
executable of a9783ee8 normalizes the `-N` second pass at its unresolved `selectionFdr_` (`Caller.cpp`:
with `-F 0`, 0), so with `-trainFDR 0 -N` the two differ; with any `train_fdr` > 0 they agree.

### The spline PEP: fixed upstream (#408)

On a90106b the patch set had four more patches for the P-IRLS spline PEP (`pep_method` `nonparametric`,
ProSE's default; OpenMS' former 03, 35, 06 and 07): tied or padded bin medians made `BaseSpline::initiateQR`
assert `dx > 0` in Debug builds and divide by zero in Release builds (the padding of 1e-20 was lost to
rounding, and padding scores in [0, 1] past 1 switched the logit to the log transform); entries of g that
the packed vector code dropped shifted the later values, and an all-zero g was read out of bounds (OpenMS
#10880); and two speedups of the matrix code. Upstream's #408 (a9783ee8) replaced that code: the Reinsch
step is solved in double-double arithmetic (`BaseSpline::reinschStep`, which fills every entry of g and
refuses a step at equal knots), `LogisticRegression::setData` merges equal knots, `PosteriorEstimator`
pads distinctly, and `PackedMatrix` is gone. The class tests of the former patches
(`PercolatorPosteriorEstimator_test`: fewer than 4 bins, bin medians closer than the spline resolves,
scores between 0 and 1; `Percolator_test`: nonparametric PEP with fewer than 4 distinct scores) pass
without them. On synthetic PINs (`perc_bench`, ProSE's configuration, 100,000 and 1,000,000 PSMs) the
scores and q-values are bit-identical to those of a90106b with the four patches, the PEPs differ by at
most 1.3e-14, and `Percolator::rescore` is not slower (1M PSMs: 3.38-3.44 s against 3.65-3.74 s).

### The former patch set (01-35 on eb157f7)

Most of the former patches were speedups of the cross-validation, whose code upstream rewrote between
eb157f7 and a90106b (258 commits, 50 in `CrossValidation.{cpp,h}`: total-order comparators and `parallelSort`, count-only grid evaluation,
no re-sorting of sorted sets, the SVM input built on all threads, PSMs in a block store, the I-spline fit
without a design matrix, and a training that fixes each fold's (Cpos, Cneg) pair after three iterations
and stops at a gain below 0.1%). Measured on synthetic PINs (`perc_bench`, 4 threads; time of
`Percolator::rescore`, peak RSS), the new tree with upstream's defaults is 3-4x faster than the former
tree with all 35 patches, with less memory, so they were not ported:

- superseded by upstream: 02 (the static flag is gone), 05, 08, 10 (folds and sorts), 09, 11, 12 (stages;
  its fold chains are re-implemented as 07), 13, 17, 28, 30 (initial direction, normalization), 20 and 32 (the
  best-scoring positive per spectrum: `OrderScanLabel` orders by score since 25cb5124; issue #10720),
  21 (PSM ownership: DataSet no longer owns its PSMs, the feature pool no longer zero-fills), 22-26
  (sort certificates, packed keys, set-up in tasks), 33 (OpenMP tasks: upstream uses none; the fold
  chains as 07), 34
  (comments of removed code);
- not ported: 27 (reuse of the validation scores), 29, 31 (code gone);
- ported: 01, 04, 14, 15, 16, 18, 19 (as 01-04 above); 03, 35, 06 and 07 were ported to a90106b and are
  now covered by upstream's #408 (see above).

`SortLikeStd.h` (sorts with the moves of `std::sort`) is removed: with total-order comparators any
sort gives the same order.

### Regenerating the patches

The patches are the commits of a port branch: a commit with the pristine upstream files of the pinned
SHA (`whitelist.txt`) and the OpenMS-local files, then one commit per patch that touches only
whitelisted files (OpenMS-local files may change in separate commits). For each patch commit, run in
`src/openms/thirdparty/percolator/` of the port branch (the pathspecs are relative to it; `--no-relative`
keeps the paths in the patch relative to the repository root, as `sync-from-upstream.sh` applies them):

```bash
files=$(grep -v '^#' whitelist.txt | grep -v '^$')
git diff --no-color --no-relative --full-index <commit>~1 <commit> -- $files > patches/NN-name.patch
```

To rebuild such a branch from this directory: copy the whitelisted files of the pinned SHA into a fresh
repository, commit them with the OpenMS-local files, and apply and commit the patches one by one.
Then check: `./sync-from-upstream.sh "$(cat UPSTREAM_COMMIT)"` must exit 0 and leave no change to the
vendored `*.h`/`*.cpp` files. Pass the SHA: without an argument the script syncs to its built-in default.

## Licensing note

`ssl.{cpp,h}` and `tron.{cpp,h}` are BSD-3-Clause (liblinear v2.11 via
the percolator-tron fork); the `AlgIn` interface in `ssl.h` is OpenMS' own implementation.
All other vendored files are Apache-2.0 (Percolator). See `NOTICE-percolator.txt` for the full attribution.

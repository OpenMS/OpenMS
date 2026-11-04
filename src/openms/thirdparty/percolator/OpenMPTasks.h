/*
 * OpenMS addition (not part of upstream Percolator): the OpenMP task constructs of the
 * parallel cross-validation, where the compiler has them.
 *
 * The tasks (task, taskwait, taskgroup; taskloop is OpenMP 4.5) need _OPENMP >= 201511, as GCC
 * and Clang have it. MSVC's OpenMP (/openmp, /openmp:experimental) is OpenMP 2.0 and has no
 * tasks; /openmp:llvm has no taskloop. Without tasks, PERCOLATOR_TASK_PRAGMA(...) expands to
 * nothing: the code of a task runs in the thread that reaches it, and a parallel region that
 * only hosts tasks is not opened. The loops over the folds, and over the (cpos, cneg) pairs in
 * CrossValidation::doStep, are a parallel for instead (#if PERCOLATOR_OMP_TASKS ... #else), one
 * fold or pair per thread; in CrossValidation::trainFoldChains, the fold chains get one thread
 * each and the pairs of a chain run as a nested parallel for on the chain's share of the threads,
 * so that the training uses more threads than there are folds. The results are the same: each
 * iteration computes what its task computes, and those results do not depend on the threads.
 *
 * Defining PERCOLATOR_OMP_TASKS as 0 on the command line builds the variant without tasks with
 * any compiler (to test it with GCC or Clang).
 */
#ifndef PERCOLATOR_OPENMP_TASKS_H
#define PERCOLATOR_OPENMP_TASKS_H

#ifndef PERCOLATOR_OMP_TASKS
#if defined(_OPENMP) && _OPENMP >= 201511 && !(defined(_MSC_VER) && !defined(__clang__))
#define PERCOLATOR_OMP_TASKS 1
#else
#define PERCOLATOR_OMP_TASKS 0
#endif
#endif

#if PERCOLATOR_OMP_TASKS
#define PERCOLATOR_TASK_PRAGMA(...) _Pragma(#__VA_ARGS__)
#else
#define PERCOLATOR_TASK_PRAGMA(...)
#endif

#endif  // PERCOLATOR_OPENMP_TASKS_H

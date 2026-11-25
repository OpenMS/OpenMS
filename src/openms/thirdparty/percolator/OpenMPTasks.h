/*
 * OpenMS addition (not part of upstream Percolator): the OpenMP task constructs of the
 * parallel cross-validation, where the compiler has them.
 *
 * The tasks (task, taskwait, taskgroup; taskloop is OpenMP 4.5) need _OPENMP >= 201511, as GCC
 * and Clang have it. MSVC's OpenMP (/openmp, /openmp:experimental) is OpenMP 2.0 and has no
 * tasks; /openmp:llvm has no taskloop. Without tasks, PERCOLATOR_TASK_PRAGMA(...) expands to
 * nothing: CrossValidation::trainsFoldChains() is false and train() runs its loop over the
 * iterations (doStep, whose pairs run as a parallel for), the row loops of Scores run as parallel
 * loops, and the Newton/LDA solvers (newton_svm.cpp) run their row blocks one after the other on
 * the calling thread instead of as a taskloop. The results are the same (the partial sums of the
 * blocks are added in block order).
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

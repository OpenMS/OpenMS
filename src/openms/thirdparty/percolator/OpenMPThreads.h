// OpenMS-local (not upstream Percolator, not in whitelist.txt): the number of OpenMP threads of the vendored code.
//
// Upstream sizes its parallel regions by omp_get_max_threads(), which the percolator executable sets once for the
// process (Caller::run). OpenMS runs the vendored code in-process, as a library whose callers choose the threads per
// call: OpenMS' Percolator parameter num_threads, at most SysInfo::ompThreadsOfThisThread() (a TOPP tool's -threads,
// or a helper thread's share where the OpenMP runtime keeps one number of threads for the whole process, as MSVC's
// vcomp does, so that omp_set_num_threads() must not be called). The parallel regions therefore read maxThreads():
// omp_get_max_threads(), at most the cap that a ThreadCap of the calling thread sets. Without OpenMP (also Clang's
// -fopenmp-simd, which parses the clauses of every directive) maxThreads() is 1 and inParallel() is false.
#ifndef PERCOLATOR_OPENMP_THREADS_H
#define PERCOLATOR_OPENMP_THREADS_H

#include <algorithm>
#include <cstddef>
#include <exception>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace OpenMS { namespace Internal { namespace Percolator {

// the cap of the calling thread (0: none)
inline int& threadCapOfThisThread()
{
  static thread_local int cap = 0;
  return cap;
}

// The threads of a parallel region that the calling thread starts: omp_get_max_threads(), at most its cap.
inline int maxThreads()
{
#ifdef _OPENMP
  const int max = std::max(1, omp_get_max_threads());
  const int cap = threadCapOfThisThread();
  return cap > 0 ? std::min(cap, max) : max;
#else
  return 1;
#endif
}

// whether the calling thread runs in an active parallel region
inline bool inParallel()
{
#ifdef _OPENMP
  return omp_in_parallel() != 0;
#else
  return false;
#endif
}

// the threads for a loop of n iterations: min(n, maxThreads()), at least 1
inline int threadsFor(long long n)
{
  return static_cast<int>(std::max<long long>(1, std::min<long long>(n, maxThreads())));
}

// Caps the threads of the parallel regions that the calling thread starts while the object exists (threads <= 0:
// no cap). Nests: the destructor restores the previous cap.
class ThreadCap
{
 public:
  explicit ThreadCap(int threads) : previous_(threadCapOfThisThread()) { threadCapOfThisThread() = std::max(0, threads); }
  ~ThreadCap() { threadCapOfThisThread() = previous_; }
  ThreadCap(const ThreadCap&) = delete;
  ThreadCap& operator=(const ThreadCap&) = delete;

 private:
  int previous_;
};

// OpenMS: the fold chains of the cross-validation (CrossValidation::trainFoldChains) run their steps as tasks of a
// parallel region, where the parallel loops over rows (the scores of a set, the copies of its training rows) would run
// on one thread. A TaskLoopScope in a task lets the loops it calls run as taskloops on the threads of the region
// instead (taskLoopThreads). The flag belongs to the thread that runs the task: tasks are tied, so until the task ends
// only its descendant tasks can run on that thread (e.g. the tasks of a taskloop it waits for, or a pair task of the
// chains at a task scheduling point of the step that creates it); these see the flag set unless they open their own
// scope. The results do not depend on it.
inline bool& taskLoopsOfThisThread()
{
  static thread_local bool on = false;
  return on;
}

class TaskLoopScope
{
 public:
  TaskLoopScope() : previous_(taskLoopsOfThisThread()) { taskLoopsOfThisThread() = true; }
  ~TaskLoopScope() { taskLoopsOfThisThread() = previous_; }
  TaskLoopScope(const TaskLoopScope&) = delete;
  TaskLoopScope& operator=(const TaskLoopScope&) = delete;

 private:
  bool previous_;
};

// the threads of the region whose tasks may run the loops of the calling task (TaskLoopScope), or 1
inline int taskLoopThreads()
{
#ifdef _OPENMP
  return taskLoopsOfThisThread() && omp_in_parallel() ? omp_get_num_threads() : 1;
#else
  return 1;
#endif
}

// An exception must not leave an OpenMP region (the process would end; OpenMS calls the vendored code as a library). The
// iterations of a parallel loop catch theirs with capture(i), and rethrowFirst() throws the one of the first failed
// iteration after the loop.
class LoopErrors
{
 public:
  explicit LoopErrors(std::size_t numIterations) : errors_(numIterations) {}
  void capture(std::size_t i) noexcept { errors_[i] = std::current_exception(); }
  void rethrowFirst() const
  {
    for (const std::exception_ptr& e : errors_)
    {
      if (e) std::rethrow_exception(e);
    }
  }

 private:
  std::vector<std::exception_ptr> errors_;
};

}}}  // namespace OpenMS::Internal::Percolator

#endif

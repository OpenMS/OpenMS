/*
 * Percolator's TRON-based SVM solver.
 *
 * Wraps liblinear v2.11's TRON optimizer (src/tron.{cpp,h}) with a
 * function-object that encodes Percolator's L2-loss SVM with
 * per-class costs (C_i = c_i / (2*lambda)). References:
 *   - Halloran & Rocke (2018), "A Matter of Time: Faster Percolator
 *     Analysis...", JPR 17(5):1978-1982.
 *   - Lin, Weng, Keerthi (2008), "Trust Region Newton Method for
 *     Large-Scale Logistic Regression", JMLR 9:627-650.
 *
 * Apache-2.0 (Percolator codebase).
 */
#ifndef _svmlin_H
#define _svmlin_H

#include <algorithm>
#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

namespace OpenMS { namespace Internal { namespace Percolator {

/* Caller-facing iteration / convergence defaults. EPSILON is the outer
 * convergence tolerance for the TRON solver; tight enough to give
 * MFN-comparable weights at Percolator's default. CGITERMAX and
 * MFNITERMAX are passed through to options.cgitermax / mfnitermax —
 * cgitermax has no liblinear analogue and is silently ignored on the
 * TRON path; mfnitermax caps TRON's outer iterations. */
constexpr int    CGITERMAX  = 10000;
constexpr int    MFNITERMAX = 50;
constexpr double EPSILON    = 1e-7;

// AlgIn: the examples of one SVM solve. This is an OpenMS implementation, written for TRON, of the member functions
// that upstream Percolator's training code (CrossValidation, Scores::generate*TrainingSet, Apache-2.0) calls on its
// training sets; it shares no code with upstream's SVMlin-derived ssl.h.
//
// An example is a row of n doubles, its n - 1 feature values followed by the constant 1.0 of the bias, and a label
// (+1 or -1). With owned rows (the default), set() copies the features into storage_ and appends the 1.0 itself. With
// borrowed rows (rowsInPlace), set() keeps the caller's pointer; upstream passes rows that end with the 1.0
// (FeatureMemoryPool::canTrainInPlace). rowTable_ has the address of every example's row in either mode, and the
// solvers (tron(), newtonSvm()) walk it through rows(). They read only the n - 1 features of a row and supply the 1.0
// themselves, so borrowed rows without the bias value work for them too (OpenMS' RESET passes such rows).
class AlgIn {
 public:
  explicit AlgIn(const int numFeat, bool rowsInPlace = false) : n(numFeat), borrowed_(rowsInPlace) {}
  // (the AlgIn(size, numFeat) of the former interface; its second int would convert to the bool above)
  AlgIn(int, int) = delete;
  AlgIn(const AlgIn& src) : n(src.n), borrowed_(src.borrowed_) { copyExamplesOf(src); }
  AlgIn(AlgIn&& src) noexcept : n(src.n), borrowed_(src.borrowed_) { takeExamplesOf(src); }
  AlgIn& operator=(const AlgIn& src) {
    if (this != &src) {
      n = src.n;
      borrowed_ = src.borrowed_;
      storage_.reset();
      storageRows_ = 0;
      copyExamplesOf(src);
    }
    return *this;
  }
  AlgIn& operator=(AlgIn&& src) noexcept {
    if (this != &src) {
      n = src.n;
      borrowed_ = src.borrowed_;
      takeExamplesOf(src);
    }
    return *this;
  }

  int m = 0;              // examples
  int n;                  // values per row: the features and the bias
  int positives = 0;
  int negatives = 0;
  std::vector<double> Y;  // label of each example

  const double* row(int i) const { return rowTable_[static_cast<std::size_t>(i)]; }
  // OpenMS: the row of every example, for the solvers
  const double* const* rows() const { return rowTable_.data(); }
  bool rowsInPlace() const { return borrowed_; }

  // forgets all examples; the memory stays allocated for the next ones
  void clear() {
    rowTable_.clear();
    Y.clear();
    m = 0;
    positives = 0;
    negatives = 0;
  }
  // how many examples fit without another allocation
  std::size_t capacity() const { return borrowed_ ? rowTable_.capacity() : storageRows_; }
  // allocates room for numExamples examples (owned rows are not written until set() is called)
  void reserve(std::size_t numExamples) {
    if (!borrowed_) ensureStorage(numExamples);
    rowTable_.reserve(numExamples);
    Y.reserve(numExamples);
  }
  // appends one example (features of row, label)
  void add(const double* row, double label) {
    const std::size_t last = rowTable_.size();
    resizeExamples(last + 1);
    set(last, row, label);
  }
  // makes room for exactly numExamples examples; those below the previous count keep their rows and labels, and the
  // added ones must be written with set() before a solve
  void resizeExamples(std::size_t numExamples) {
    const std::size_t before = rowTable_.size();
    if (!borrowed_ && numExamples > storageRows_) {
      // grow by half at least, so that add() in a loop does not reallocate per example
      ensureStorage(std::max({numExamples, storageRows_ + storageRows_ / 2, kMinStorageRows}));
    }
    rowTable_.resize(numExamples, nullptr);
    if (!borrowed_) {
      for (std::size_t i = before; i < numExamples; ++i) rowTable_[i] = ownedRow(i);
    }
    Y.resize(numExamples, 0.0);
    m = static_cast<int>(numExamples);
  }
  // writes example i (< m): the features of row, then 1.0 for owned rows, and the label; different threads may
  // write different examples at once
  void set(std::size_t i, const double* row, double label) {
    if (borrowed_) {
      rowTable_[i] = row;
    } else {
      double* target = ownedRow(i);
      std::copy(row, row + (n - 1), target);
      target[n - 1] = 1.0;
    }
    Y[i] = label;
  }

 private:
  static constexpr std::size_t kMinStorageRows = 256;
  bool borrowed_;
  std::unique_ptr<double[]> storage_;  // owned rows, storageRows_ x n values (not value-initialised)
  std::size_t storageRows_ = 0;
  std::vector<const double*> rowTable_;  // the row of each example

  double* ownedRow(std::size_t i) const { return storage_.get() + i * static_cast<std::size_t>(n); }

  // owned rows: storage for at least rowsWanted rows, keeping the rows written so far
  void ensureStorage(std::size_t rowsWanted) {
    if (rowsWanted <= storageRows_) return;
    std::unique_ptr<double[]> larger = std::make_unique_for_overwrite<double[]>(rowsWanted * static_cast<std::size_t>(n));
    const std::size_t written = rowTable_.size() * static_cast<std::size_t>(n);
    if (written > 0) std::copy(storage_.get(), storage_.get() + written, larger.get());
    storage_ = std::move(larger);
    storageRows_ = rowsWanted;
    for (std::size_t i = 0; i < rowTable_.size(); ++i) rowTable_[i] = ownedRow(i);
  }
  void copyExamplesOf(const AlgIn& src) {
    m = src.m;
    positives = src.positives;
    negatives = src.negatives;
    Y = src.Y;
    if (borrowed_) {
      rowTable_ = src.rowTable_;
      return;
    }
    // no rows to keep (ensureStorage copies those of rowTable_, and storage_ may be gone: operator=)
    rowTable_.clear();
    ensureStorage(src.rowTable_.size());
    rowTable_.resize(src.rowTable_.size());
    for (std::size_t i = 0; i < src.rowTable_.size(); ++i) {
      rowTable_[i] = ownedRow(i);
      std::copy(src.rowTable_[i], src.rowTable_[i] + n, ownedRow(i));
    }
  }
  void takeExamplesOf(AlgIn& src) noexcept {
    m = src.m;
    positives = src.positives;
    negatives = src.negatives;
    Y = std::move(src.Y);
    storage_ = std::move(src.storage_);  // the row addresses in rowTable_ stay valid
    storageRows_ = src.storageRows_;
    rowTable_ = std::move(src.rowTable_);
    src.storageRows_ = 0;
    src.rowTable_.clear();
    src.Y.clear();
    src.m = src.positives = src.negatives = 0;
  }
};

// Cache hints for scans that follow pointers to rows or records (no-ops without a compiler prefetch builtin).
inline void prefetchAddress(const void* p) {
#if defined(__GNUC__) || defined(__clang__)
  __builtin_prefetch(p, 0, 3);
#else
  (void)p;
#endif
}

// one hint per 64-byte cache line of the n values of row
inline void prefetchRow(const double* row, int n) {
  const unsigned char* bytes = reinterpret_cast<const unsigned char*>(row);
  const std::size_t length = sizeof(double) * static_cast<std::size_t>(n);
  for (std::size_t offset = 0; offset < length; offset += 64) prefetchAddress(bytes + offset);
}

struct options {
  double lambda;
  double lambda_u;
  double epsilon;
  int cgitermax;
  int mfnitermax;
  // upstream: its SVMlin solver solves the least squares steps exactly; TRON has no such steps (ignored)
  bool directSolve = false;
  // OpenMS: the caller needs the weights only (tron(): the Outputs are not computed, the convergence check after
  // the solve is skipped and tron() returns 1); the weights are the same
  bool weightsOnly = false;
};

// The per-thread solver buffers in upstream's L2_SVM_MFN interface. tron() takes its buffers from a pool of
// workspaces (see TronWorkspaceHold), so this holds nothing.
struct SolverWorkspace {};

// TRON solver — the real implementation. Returns 1 if the post-solve
// relative gradient norm satisfies the TRON convergence criterion
// (||g(w)|| / ||g(w_0)|| <= opts.epsilon); 0 otherwise. Asserts
// opts.lambda > 0 — a precondition of the C_i scaling.
//
// opts.epsilon is passed directly to TRON as its outer convergence eps;
// at Percolator's default 1e-7 the TRON path converges to MFN-comparable
// weights. opts.cgitermax is silently ignored — vendored liblinear TRON
// has no cgitermax analogue (inner CG is an uncapped while-loop until
// cg-tolerance or trust-region exit).
//
// OpenMS: on x86-64 builds for the baseline instruction set (no AVX, no FMA), the
// function object runs on kernels for SSE2 or AVX2 (chosen at run time), with the
// same results bit for bit (tron_kernels.h). gradRatio (optional) receives the
// gradient-norm ratio (0 if ||g(w_0)|| = 0).
class TronPairGroup;
// Weights: n values, the start of the solve (0, or a warm start), the solution on return. Outputs: X w for the
// m examples on return (resized to at least m), unless opts.weightsOnly. With no examples (m = 0) the weights are 0
// and tron() returns 1.
int tron(const AlgIn& set, const options& opts,
         std::vector<double>& Weights, std::vector<double>& Outputs,
         double cpos, double cneg, double* gradRatio = nullptr,
         TronPairGroup* group = nullptr);

// OpenMS: solves of tron() on the same set (e.g. the (Cpos, Cneg) pairs of a cross-validation fold, which run at
// the same time) that pass the same group share one copy of the set's features in the solver's layout: it is
// built once, by the solves that start while it is built, instead of once per solve. Results are bit-identical to
// solves alone. A solve waits only while the copy is being built by running solves; solves with another set or
// size, or on the reference solver, build their own. The group must outlive its solves, and the set (its rows, their
// contents and the labels) must not change during the group's lifetime: a solve checks only the set's address, size,
// row array, first and last row and class counts, so one group per set and step, as CrossValidation
// makes them.
class TronPairGroup {
 public:
  explicit TronPairGroup(const AlgIn& set);
  ~TronPairGroup();
  TronPairGroup(const TronPairGroup&) = delete;
  TronPairGroup& operator=(const TronPairGroup&) = delete;
  struct State;
  State* state() { return state_.get(); }  // (ssl.cpp)

 private:
  std::unique_ptr<State> state_;
};

// OpenMS: asks Linux to back the 2 MiB pages that lie inside [p, p + bytes) with transparent huge pages (a hint, which
// the kernel follows in its "always" and "madvise" THP modes; nothing elsewhere, and nothing for blocks under 4 MiB).
// For large buffers that are written once and read many times (the training matrices of a large data set): 512
// times fewer page faults to fill them, far fewer TLB misses to read them. Call it before the first write. Values are
// not affected.
void adviseHugePages(const void* p, std::size_t bytes);

// OpenMS: whether a TronWorkspaceHold exists (the Newton solver keeps its buffers then as well)
bool tronWorkspacesHeld();

// OpenMS: while an object of this class exists, tron() keeps the buffers of its finished solves
// of any size for the next solves (otherwise those above 64 MB are freed after the solve); the
// last one to end frees them. CrossValidation holds one during its training.
class TronWorkspaceHold {
 public:
  TronWorkspaceHold();
  ~TronWorkspaceHold();
  TronWorkspaceHold(const TronWorkspaceHold&) = delete;
  TronWorkspaceHold& operator=(const TronWorkspaceHold&) = delete;
};

// OpenMS, for tests: selects the function object of tron() for all threads: "ref"
// (the reference TronFunction), "sse2" or "avx2" (the kernels); nullptr = the
// default (environment variable OPENMS_TRON_ISA = ref|sse2|avx2, else AVX2 if the
// CPU has it). Returns false, and changes nothing, if the path is not available
// in this build or on this CPU.
bool tronSelectKernels(const char* isa);

// OpenMS, for tests: the function object of tron() ("ref", "sse2" or "avx2") on the
// problem of tron(set, opts, ..., cpos, cneg) with opts.lambda = lambda, called as
// TRON calls it: f[0] = fun(0), g[0..n) = grad(0), Hs[0..n) = Hv(s), then f[1] =
// fun(w), g[n..2n) = grad(w), Hs[n..2n) = Hv(s) (n = set.n). Returns false if the
// kernels are not available.
bool tronEvaluateForTest(const char* isa, const AlgIn& set, double cpos, double cneg,
                         double lambda, const double* w, const double* s,
                         double* f, double* g, double* Hs);

// Drop-in replacement for the legacy SVMlin-based L2_SVM_MFN. The
// original implementation in this file was the modified-finite-Newton
// solver from Vikas Sindhwani's SVMlin (2006). Sindhwani granted
// Percolator a relicensing to redistribute the SVMlin code as part of
// Percolator under the project's Apache-2.0 terms, but the relicense
// was scoped to Percolator itself and does NOT cleanly cover link-time
// inclusion of perclibrary into downstream binaries that ship under
// non-Apache licenses — the upstream SVMlin headline is GPL v2+, so
// any third-party project that statically links perclibrary inherits
// an unresolved GPL-vs-license-X question. Replacing the body with a
// forwarder to TRON (BSD-3 vendored liblinear) sidesteps that question
// at the source. The public name is preserved so existing callers
// compile unchanged.
// The signatures are those of upstream's L2_SVM_MFN (3.09); the start is Weights (Outputs is not read).
inline int L2_SVM_MFN(const AlgIn& set, const options& opts,
                      std::vector<double>& Weights, std::vector<double>& Outputs,
                      double cpos, double cneg, SolverWorkspace& /*workspace*/) {
  return tron(set, opts, Weights, Outputs, cpos, cneg);
}

inline int L2_SVM_MFN(const AlgIn& set, const options& opts,
                      std::vector<double>& Weights, std::vector<double>& Outputs,
                      double cpos, double cneg, TronPairGroup* group = nullptr) {
  return tron(set, opts, Weights, Outputs, cpos, cneg, nullptr, group);
}

}}}  // namespace OpenMS::Internal::Percolator

#endif

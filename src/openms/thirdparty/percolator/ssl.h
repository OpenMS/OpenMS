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

#include <cstddef>
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

class AlgIn {
 public:
  AlgIn(const unsigned int size, const int numFeat);
  virtual ~AlgIn();
  // No copy / no move. AlgIn owns raw arrays and has no copy/move
  // semantics — any return-by-value path that doesn't get NRVO would
  // shallow-copy the pointers and double-free in the source dtor.
  // Deleting these turns any non-elided path into a compile error,
  // which is the safe failure mode.
  AlgIn(const AlgIn&) = delete;
  AlgIn& operator=(const AlgIn&) = delete;
  AlgIn(AlgIn&&) = delete;
  AlgIn& operator=(AlgIn&&) = delete;
  int m;          /* number of examples */
  int n;          /* number of features, including synthesized bias slot */
  int positives;
  int negatives;
  double** vals;  /* m pointers to external feature storage of length n-1 */
  double* Y;      /* labels, size m */
};

struct vector_double {
  int d;
  double* vec = nullptr;
  ~vector_double() { delete[] vec; }
};

struct vector_int {
  int d;
  int* vec = nullptr;
  ~vector_int() { delete[] vec; }
};

struct options {
  double lambda;
  double lambda_u;
  double epsilon;
  int cgitermax;
  int mfnitermax;
};

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
int tron(const AlgIn& set, options& opts,
         vector_double& Weights, vector_double& Outputs,
         double cpos, double cneg, double* gradRatio = nullptr);

// OpenMS: asks Linux to back the 2 MiB pages that lie inside [p, p + bytes) with transparent huge pages (a hint, which
// the kernel follows in its "always" and "madvise" THP modes; nothing elsewhere, and nothing for blocks under 4 MiB).
// For large buffers that are written once and read many times (the training matrices of a large data set): 512
// times fewer page faults to fill them, far fewer TLB misses to read them. Call it before the first write. Values are
// not affected.
void adviseHugePages(const void* p, std::size_t bytes);

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
inline int L2_SVM_MFN(const AlgIn& set, options& opts,
                      vector_double& Weights, vector_double& Outputs,
                      double cpos, double cneg) {
  return tron(set, opts, Weights, Outputs, cpos, cneg);
}

}}}  // namespace OpenMS::Internal::Percolator

#endif

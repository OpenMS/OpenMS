/* See ssl.h for the contract. */
#include "ssl.h"

#include <algorithm>
#include <atomic>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <memory>
#include <mutex>
#include <vector>

#include <Eigen/Dense>

#include "tron.h"

// OpenMS: the fast function object (FastTronFunction, tron_kernels.h) reproduces the bits of
// TronFunction as compiled for the x86-64 baseline (Eigen's SSE2 gemv, no FMA). Builds with AVX or
// FMA (-march=...) and other architectures keep TronFunction.
#if (defined(__x86_64__) || defined(_M_X64)) && !defined(__AVX__) && !defined(__FMA__)
#define PERCOLATOR_TRON_KERNELS 1
#include "cpu_avx2.h"
#include "tron_kernels.h"
#endif

namespace OpenMS { namespace Internal { namespace Percolator {

// AlgIn lifecycle — was previously defined in the legacy SVMlin
// ssl.cpp; carried forward unchanged so existing callers (CrossValidation,
// tests) keep their construction sites.
AlgIn::AlgIn(const unsigned int size, const int numFeat) {
  vals = new double*[size];
  Y = new double[size];
  n = numFeat;
  positives = 0;
  negatives = 0;
}

AlgIn::~AlgIn() {
  delete[] vals;
  delete[] Y;
}

namespace {

// Row-major dense matrix — TRON's hot paths (X*w in fun, row scatter
// in grad/Hv) are row-oriented; column-major Eigen::MatrixXd makes
// each X.row(i) a stride-m gather, which is catastrophic on
// m≈10000 problems. Row-major fixes that.
using DenseRowMatrix =
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

// Silence liblinear's per-iteration logging.
void silent_print(const char* /*buf*/) {}

// Inner-CG tolerance — hardcoded per Halloran & Rocke (2018). Vendored
// liblinear TRON has no cgitermax knob, so inner CG is uncapped (while-loop
// until cg-tolerance or trust-region exit).
//
// TRON's outer convergence eps is taken from opts.epsilon — at Percolator's
// default 1e-7 this gives tight (MFN-comparable) weights; users wanting
// Halloran's published speedup can pass a looser --epsilon (e.g. 0.01).
constexpr double kTronEpsCg = 0.1;

// TronFunction implements liblinear's `function` interface for
// Percolator's L2-loss SVM with per-class cost weighting.
//
// Math (proof in the design spec): Percolator MFN minimizes
//   F_MFN(w) = (lambda/2)||w||^2 + (1/2) sum c_i (1 - t_i w.x_i)^2
// liblinear's L2-loss SVC minimizes
//   F_TRON(w) = (1/2)||w||^2 + sum C_i (1 - t_i w.x_i)^2
// To make TRON's optimum coincide with MFN's, set
//   C_i = c_i / (2 * lambda)
// (this fixed scaling makes F_TRON = F_MFN / lambda, a positive
//  multiplicative reparameterization that preserves the minimizer).
class TronFunction : public function {
 public:
  TronFunction(const DenseRowMatrix& X,
               const Eigen::VectorXd& t,
               const Eigen::VectorXd& C)
      : X_(X), t_(t), C_(C),
        m_(static_cast<int>(X.rows())), n_(static_cast<int>(X.cols())),
        margin_(m_), I_active_(m_), z_active_(m_), sizeI_(0) {
    assert(t_.size() == m_);
    assert(C_.size() == m_);
  }

  // f(w) = (1/2)||w||^2 + sum_{i:margin_i<1} C_i (1 - margin_i)^2
  double fun(double* w_raw) override {
    Eigen::Map<const Eigen::VectorXd> w(w_raw, n_);
    margin_ = (X_ * w).cwiseProduct(t_);  // margin_i = t_i * (X w)_i

    double f = 0.5 * w.squaredNorm();
    for (int i = 0; i < m_; ++i) {
      if (margin_(i) < 1.0) {
        const double d = 1.0 - margin_(i);
        f += C_(i) * d * d;
      }
    }
    return f;
  }

  // grad: g = w + 2 * sum_{i in I} C_i * t_i * (margin_i - 1) * x_i
  // where I = {i : margin_i < 1}.
  // Single-pass over examples: scan + accumulate fused; raw pointer
  // arithmetic on row-major X_ to give the compiler simple inner loops.
  // Caches the active-set indices for the subsequent Hv call.
  void grad(double* w_raw, double* g_raw) override {
    sizeI_ = 0;
    for (int j = 0; j < n_; ++j) g_raw[j] = 0.0;

    const double* X_data = X_.data();
    for (int i = 0; i < m_; ++i) {
      if (margin_(i) < 1.0) {
        I_active_[sizeI_++] = i;
        const double z = 2.0 * C_(i) * t_(i) * (margin_(i) - 1.0);
        const double* row = X_data + (size_t)i * n_;
        for (int j = 0; j < n_; ++j) g_raw[j] += z * row[j];
      }
    }
    for (int j = 0; j < n_; ++j) g_raw[j] += w_raw[j];
  }

  // Hv(s) = s + 2 * sum_{i in I} C_i (x_i . s) x_i  (Gauss-Newton Hessian)
  // Same raw-pointer / row-major access pattern as grad.
  void Hv(double* s_raw, double* Hs_raw) override {
    for (int j = 0; j < n_; ++j) Hs_raw[j] = 0.0;

    const double* X_data = X_.data();
    for (int k = 0; k < sizeI_; ++k) {
      const int i = I_active_[k];
      const double* row = X_data + (size_t)i * n_;
      double xTs = 0.0;
      for (int j = 0; j < n_; ++j) xTs += row[j] * s_raw[j];
      const double scale = 2.0 * C_(i) * xTs;
      for (int j = 0; j < n_; ++j) Hs_raw[j] += scale * row[j];
    }
    for (int j = 0; j < n_; ++j) Hs_raw[j] += s_raw[j];
  }

  int get_nr_variable(void) override { return n_; }

  // Public for the post-solve gradient-ratio check.
  void compute_grad_at(const Eigen::VectorXd& w_eig,
                       Eigen::VectorXd& g_out) {
    g_out.resize(n_);
    fun(const_cast<double*>(w_eig.data()));   // refresh margin_
    grad(const_cast<double*>(w_eig.data()), g_out.data());
  }

 private:
  const DenseRowMatrix& X_;
  const Eigen::VectorXd& t_;
  const Eigen::VectorXd& C_;
  const int m_;
  const int n_;
  Eigen::VectorXd margin_;       // margin_i = t_i * (X w)_i, refreshed by fun()
  std::vector<int> I_active_;    // active-set indices, refreshed by grad()
  Eigen::VectorXd z_active_;     // z_active[k] = C_i t_i (margin_i - 1)
  int sizeI_;                    // number of active examples
};

// ---- OpenMS: fast function object (bit-identical to TronFunction) ----

enum TronKernelChoice { kTronRef = 0, kTronSse2 = 1, kTronAvx2 = 2 };

#ifdef PERCOLATOR_TRON_KERNELS
// the AVX2 table is empty if tron_kernels_avx2.cpp was not compiled for AVX2 without FMA
bool avx2Available() { return TronKernels::avx2::ops.hv != nullptr && cpuHasAvx2(); }
#endif

// OPENMS_TRON_ISA=ref|sse2|avx2 overrides the choice by CPU (for tests); avx2 needs a CPU with AVX2.
int defaultTronKernels() {
#ifdef PERCOLATOR_TRON_KERNELS
  const char* e = std::getenv("OPENMS_TRON_ISA");
  if (e != nullptr && std::strcmp(e, "ref") == 0) return kTronRef;
  if (e != nullptr && std::strcmp(e, "sse2") == 0) return kTronSse2;
  if (e != nullptr && *e != '\0' && std::strcmp(e, "avx2") != 0) {
    static std::once_flag warned;
    std::call_once(warned, [e] {
      std::fprintf(stderr, "Percolator TRON: ignoring unknown OPENMS_TRON_ISA=%s (ref, sse2 or avx2)\n", e);
    });
  }
  return avx2Available() ? kTronAvx2 : kTronSse2;
#else
  return kTronRef;
#endif
}

std::atomic<int>& tronKernels() {
  static std::atomic<int> k{defaultTronKernels()};
  return k;
}

#ifdef PERCOLATOR_TRON_KERNELS
// Grow-only buffer of doubles, 64-byte aligned (BT4 blocks are loaded with aligned loads).
struct AlignedDoubles {
  std::vector<double> v;
  double* p = nullptr;
  size_t cap = 0;
  // the old contents are not kept; if the allocation throws, the buffer is empty
  double* get(size_t n) {
    if (n > cap) {
      const size_t c = std::max(n, cap + cap / 2);
      std::vector<double>().swap(v);
      p = nullptr;
      cap = 0;
      v.resize(c + 8);
      p = reinterpret_cast<double*>((reinterpret_cast<std::uintptr_t>(v.data()) + 63) & ~std::uintptr_t(63));
      cap = c;
    }
    return p;
  }
  size_t bytes() const { return v.capacity() * sizeof(double); }
};

// Workspace of a solve of tron(). The cross-validation runs many solves one after the other (on
// the threads of an OpenMP team), so the buffers are taken from a pool and given back, and are
// allocated (and page-faulted) once, not per solve.
struct TronWorkspace {
  AlignedDoubles xb;    // X in BT4 layout
  AlignedDoubles a[2];  // active rows in BT4 layout (current, pending)
  std::vector<double> t, C, ones, out, margin, c2, z, sc, h;
  std::vector<int> idx;

  size_t bytes() const {
    size_t b = xb.bytes() + a[0].bytes() + a[1].bytes() + idx.capacity() * sizeof(int);
    for (const std::vector<double>* v : {&t, &C, &ones, &out, &margin, &c2, &z, &sc, &h}) b += v->capacity() * sizeof(double);
    return b;
  }
};

// Workspaces above this size are freed at the end of tron(), so that a huge problem does not keep
// its buffers in the pool, unless a TronWorkspaceHold exists (then at the end of the last one).
constexpr size_t kTronWorkspaceKeepBytes = size_t(64) << 20;

struct TronWorkspacePool {
  std::mutex mutex;
  std::vector<std::unique_ptr<TronWorkspace> > idle;
  int holds = 0;
};

TronWorkspacePool& tronWorkspacePool() {
  static TronWorkspacePool pool;
  return pool;
}

// A workspace from the pool for one solve, given back at the end (freed if it is large and no
// TronWorkspaceHold exists)
class TronWorkspaceLease {
 public:
  TronWorkspaceLease() {
    TronWorkspacePool& pool = tronWorkspacePool();
    {
      std::lock_guard<std::mutex> lock(pool.mutex);
      if (!pool.idle.empty()) {
        ws_ = std::move(pool.idle.back());
        pool.idle.pop_back();
      }
    }
    if (!ws_) ws_ = std::make_unique<TronWorkspace>();
  }
  ~TronWorkspaceLease() {
    TronWorkspacePool& pool = tronWorkspacePool();
    std::lock_guard<std::mutex> lock(pool.mutex);
    if (pool.holds == 0 && ws_->bytes() > kTronWorkspaceKeepBytes) return;  // freed
    try {
      pool.idle.push_back(std::move(ws_));
    } catch (...) {  // (freed)
    }
  }
  TronWorkspaceLease(const TronWorkspaceLease&) = delete;
  TronWorkspaceLease& operator=(const TronWorkspaceLease&) = delete;
  TronWorkspace& get() { return *ws_; }

 private:
  std::unique_ptr<TronWorkspace> ws_;
};

template <class T> void growTo(std::vector<T>& v, size_t n) {
  if (v.size() < n) v.resize(n);
}

// TronFunction on X in BT4 layout (tron_kernels.h), same bits. Differences in what is computed:
//  - fun keeps the active rows (in order) and packs them into the pending buffer while the chunk
//    of X is in L1 (capacity: 1.25 x the last active count; if too small, grad packs them);
//  - f and g at w = 0 (all bits zero) are remembered: TRON evaluates them twice, the gradient-ratio
//    check of tron() a third time;
//  - with all rows active, grad and Hv use X itself.
class FastTronFunction : public function {
 public:
  FastTronFunction(const TronKernels::Ops& ops, TronWorkspace& ws, const double* XB,
                   int m, int n, const double* t, const double* C)
      : ops_(ops), ws_(ws), XB_(XB), t_(t), C_(C), m_(m), mp_((m + 31) & ~31),
        n_(n), np_((n + 3) & ~3), g0_(n) {
    growTo(ws_.margin, mp_);
    growTo(ws_.c2, mp_);
    growTo(ws_.z, mp_);
    growTo(ws_.sc, mp_);
    growTo(ws_.h, np_);
    growTo(ws_.idx, m_ + 1);
  }

  double fun(double* w_raw) override {
    bool zero = true;
    for (int j = 0; j < n_; ++j) {
      std::uint64_t b;
      std::memcpy(&b, w_raw + j, sizeof(b));
      zero &= b == 0;
    }
    if (zero && have0_) {
      memo_ = true;
      return f0_;
    }
    memo_ = false;
    zero_ = zero;
    Eigen::Map<const Eigen::VectorXd> w(w_raw, n_);
    double f = 0.5 * w.squaredNorm();
    const std::int64_t est = std::max<std::int64_t>(est_, 1024) * 5 / 4;
    const int capRows = (int)std::min<std::int64_t>(mp_, (est + 31) & ~std::int64_t(31));
    double* Ap = zero ? nullptr : ws_.a[pend_].get((size_t)capRows * np_);
    int packed = 0;
    nfa_ = ops_.fun_bt(XB_, n_, np_, w_raw, t_, C_, m_, ws_.margin.data(), ws_.idx.data(),
                       Ap, capRows, &f, &packed);
    packed_ = packed != 0;
    if (zero) f0_ = f;
    else est_ = nfa_;
    return f;
  }

  void grad(double* w_raw, double* g_raw) override {
    double* c2 = ws_.c2.data();
    if (memo_) {  // g at w = 0 known; all rows are active
      for (int j = 0; j < n_; ++j) g_raw[j] = g0_[j];
      na_ = m_;
      nq_ = mp_;
      A_ = XB_;
      for (int q = 0; q < mp_; ++q) c2[q] = 2.0 * C_[q];
      return;
    }
    na_ = nfa_;
    const bool all = na_ == m_;
    nq_ = all ? mp_ : (na_ + 31) & ~31;
    const double* mg = ws_.margin.data();
    const int* idx = ws_.idx.data();
    double* z = ws_.z.data();
    for (int q = 0; q < na_; ++q) {
      const int i = idx[q];
      c2[q] = 2.0 * C_[i];
      z[q] = 2.0 * C_[i] * t_[i] * (mg[i] - 1.0);
    }
    for (int q = na_; q < nq_; ++q) {
      c2[q] = 0.0;
      z[q] = 0.0;
    }
    if (all) {
      A_ = XB_;
    } else if (packed_) {
      std::swap(cur_, pend_);
      A_ = ws_.a[cur_].p;
    } else {
      double* a = ws_.a[cur_].get((size_t)nq_ * np_);
      ops_.pack_bt(XB_, np_, idx, na_, 0, nq_, a);
      A_ = a;
    }
    packed_ = false;
    double* gt = ws_.h.data();
    ops_.scatter(A_, np_, z, nq_, 0, np_ / 4, gt);
    for (int j = 0; j < n_; ++j) g_raw[j] = gt[j] + w_raw[j];
    if (zero_ && all) {  // the memo assumes all rows active (a NaN margin is not)
      for (int j = 0; j < n_; ++j) g0_[j] = g_raw[j];
      have0_ = true;
    }
  }

  void Hv(double* s_raw, double* Hs_raw) override {
    bool finite = true;
    for (int j = 0; j < n_; ++j) finite &= std::isfinite(s_raw[j]);
    if (!finite) return HvActiveRows(s_raw, Hs_raw);
    double* h = ws_.h.data();
    ops_.hv(A_, n_, np_, s_raw, ws_.c2.data(), nq_, h, ws_.sc.data());
    for (int j = 0; j < n_; ++j) Hs_raw[j] = h[j] + s_raw[j];
  }

  int get_nr_variable(void) override { return n_; }

 private:
  // Hv as TronFunction::Hv, over the active rows only: with an infinite or NaN s (a gradient that
  // overflowed), the zero rows that pad the kernels' blocks would add 0 * Inf = NaN.
  void HvActiveRows(const double* s, double* Hs) const {
    for (int j = 0; j < n_; ++j) Hs[j] = 0.0;
    const double* c2 = ws_.c2.data();
    for (int q = 0; q < na_; ++q) {
      const double* x = A_ + (std::ptrdiff_t)(q >> 2) * np_ * 4 + (q & 3);
      double xTs = 0.0;
      for (int j = 0; j < n_; ++j) xTs += x[4 * j] * s[j];
      const double scale = c2[q] * xTs;
      for (int j = 0; j < n_; ++j) Hs[j] += scale * x[4 * j];
    }
    for (int j = 0; j < n_; ++j) Hs[j] += s[j];
  }

  const TronKernels::Ops& ops_;
  TronWorkspace& ws_;
  const double* XB_;  // mp x np, BT4
  const double* t_;   // labels, zero-padded to mp
  const double* C_;   // costs, zero-padded to mp
  const int m_, mp_, n_, np_;
  const double* A_ = nullptr;  // rows used by grad's scatter and Hv: XB_ or the packed active rows
  int nfa_ = 0, na_ = 0, nq_ = 0, est_ = 0, cur_ = 0, pend_ = 1;
  bool packed_ = false, zero_ = false, memo_ = false, have0_ = false;
  double f0_ = 0.0;
  std::vector<double> g0_;
};
#endif

// Post-solve convergence indicator: ||grad(w_final)|| / ||grad(w_0)||.
int gradRatioConverged(const Eigen::VectorXd& g0, const Eigen::VectorXd& g_final,
                       double epsilon, double* gradRatio) {
  const double g0_norm = g0.norm();
  if (g0_norm == 0.0) {  // pathological: no signal — declare convergence
    if (gradRatio != nullptr) *gradRatio = 0.0;
    return 1;
  }
  const double ratio = g_final.norm() / g0_norm;
  if (gradRatio != nullptr) *gradRatio = ratio;
  return (ratio <= epsilon) ? 1 : 0;
}

// the reference: TronFunction on a row-major Eigen matrix
int tronRef(const AlgIn& set, options& opts, vector_double& Weights, vector_double& Outputs,
            double cpos, double cneg, double* gradRatio) {
  const int m = set.m;
  const int n = Weights.d;
  const int num_real = n - 1;

  // Build dense matrix X of shape (m, n) with synthesized 1.0 bias column.
  // Row-major so X.row(i) is contiguous in memory (the hot pattern in
  // grad and Hv).
  DenseRowMatrix X(m, n);
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < num_real; ++j) {
      X(i, j) = set.vals[i][j];
    }
    X(i, num_real) = 1.0;
  }

  // Labels
  Eigen::VectorXd t_eig(m);
  for (int i = 0; i < m; ++i) t_eig(i) = set.Y[i];

  // Per-example costs: C_i = c_i / (2 * lambda).
  const double inv_2lambda = 1.0 / (2.0 * opts.lambda);
  Eigen::VectorXd C_eig(m);
  for (int i = 0; i < m; ++i) {
    C_eig(i) = (set.Y[i] > 0.0 ? cpos : cneg) * inv_2lambda;
  }

  // Solve.
  TronFunction fn(X, t_eig, C_eig);
  TRON tron_solver(&fn, opts.epsilon, kTronEpsCg, opts.mfnitermax);
  tron_solver.set_print_string(silent_print);
  tron_solver.tron(Weights.vec);

  // Reconstruct Outputs[i] = X.row(i) . w  (bias is already inside w[n-1]).
  Eigen::Map<const Eigen::VectorXd> w_eig(Weights.vec, n);
  Eigen::Map<Eigen::VectorXd> y_out(Outputs.vec, m);
  y_out = X * w_eig;

  // Post-solve convergence indicator: ||grad(w_final)|| / ||grad(w_0)||.
  // w_0 = 0 here (caller passes zero-initialized Weights).
  Eigen::VectorXd g0, g_final;
  Eigen::VectorXd w_zero = Eigen::VectorXd::Zero(n);
  fn.compute_grad_at(w_zero, g0);
  fn.compute_grad_at(w_eig, g_final);
  return gradRatioConverged(g0, g_final, opts.epsilon, gradRatio);
}

#ifdef PERCOLATOR_TRON_KERNELS
const TronKernels::Ops& tronOps(int kernels) {
  return kernels == kTronAvx2 ? TronKernels::avx2::ops : TronKernels::sse2::ops;
}

// the same with FastTronFunction on X in BT4 layout (m >= 2)
int tronFast(const TronKernels::Ops& ops, const AlgIn& set, options& opts, vector_double& Weights,
             vector_double& Outputs, double cpos, double cneg, double* gradRatio) {
  const int m = set.m;
  const int n = Weights.d;
  const int num_real = n - 1;
  const double inv_2lambda = 1.0 / (2.0 * opts.lambda);
  TronWorkspaceLease lease;
  TronWorkspace& ws = lease.get();
  const int np = (n + 3) & ~3;
  const int mp = (m + 31) & ~31;
  int ret;
  {
    // X (m, n) with the 1.0 bias column, in BT4 layout, zero-padded to (mp, np)
    double* XB = ws.xb.get((size_t)mp * (size_t)np);
    ops.build_bt(set.vals, m, num_real, np, XB);
    ws.t.assign(mp, 0.0);
    ws.C.assign(mp, 0.0);
    for (int i = 0; i < m; ++i) {
      ws.t[i] = set.Y[i];
      ws.C[i] = (set.Y[i] > 0.0 ? cpos : cneg) * inv_2lambda;
    }

    FastTronFunction fn(ops, ws, XB, m, n, ws.t.data(), ws.C.data());
    TRON tron_solver(&fn, opts.epsilon, kTronEpsCg, opts.mfnitermax);
    tron_solver.set_print_string(silent_print);
    tron_solver.tron(Weights.vec);

    // Outputs = X * w: the margins with t = 1 ((0 + cc) * 1, as Eigen's 0 + 1 * cc)
    ws.ones.assign(mp, 1.0);
    growTo(ws.out, mp);
    ops.margins_bt(XB, n, np, Weights.vec, ws.ones.data(), 0, mp, ws.out.data());
    std::copy(ws.out.begin(), ws.out.begin() + m, Outputs.vec);

    Eigen::VectorXd g0(n), g_final(n);
    Eigen::VectorXd w_zero = Eigen::VectorXd::Zero(n);
    fn.fun(w_zero.data());
    fn.grad(w_zero.data(), g0.data());
    fn.fun(Weights.vec);
    fn.grad(Weights.vec, g_final.data());
    ret = gradRatioConverged(g0, g_final, opts.epsilon, gradRatio);
  }
  return ret;
}

// One-time check per kernel table, at its first use: two small problems (odd and even n, mixed
// scales, rows entering and leaving the active set) solved with TronFunction and with the kernels
// must give the same bits; otherwise (another Eigen or compiler than the ones checked) tron() uses
// TronFunction.
bool selfCheckOnce(const TronKernels::Ops& ops) {
  std::uint64_t r = 0x9E3779B97F4A7C15ULL;
  auto uniform = [&r]() {
    r = r * 6364136223846793005ULL + 1442695040888963407ULL;
    return double(r >> 11) * (1.0 / 9007199254740992.0);
  };
  const int sizes[2][2] = {{160, 8}, {64, 45}};  // (m, features): n = 9 (one kernel), 46 (generic)
  for (const auto& sz : sizes) {
    const int m = sz[0], nf = sz[1];
    std::vector<double> vals((size_t)m * nf);
    AlgIn in(m, nf + 1);
    in.m = m;
    for (int i = 0; i < m; ++i) {
      const double label = uniform() < 0.4 ? 1.0 : -1.0;
      in.Y[i] = label;
      for (int j = 0; j < nf; ++j) {
        const double v = uniform() + uniform() - 1.0 + (j % 3 == 0 ? 0.3 * label : 0.0);
        vals[(size_t)i * nf + j] = v * (j % 4 == 1 ? 1e3 : j % 4 == 2 ? 1e-3 : 1.0);
      }
      in.vals[i] = vals.data() + (size_t)i * nf;
    }
    options opts{1.0, 1.0, EPSILON, CGITERMAX, MFNITERMAX};
    vector_double w[2], out[2];
    double ratio[2];
    int conv[2];
    for (int k = 0; k < 2; ++k) {
      w[k].d = nf + 1;
      w[k].vec = new double[nf + 1]();
      out[k].d = m;
      out[k].vec = new double[m]();
    }
    conv[0] = tronRef(in, opts, w[0], out[0], 1.0, 3.0, &ratio[0]);
    conv[1] = tronFast(ops, in, opts, w[1], out[1], 1.0, 3.0, &ratio[1]);
    if (conv[0] != conv[1] || std::memcmp(&ratio[0], &ratio[1], sizeof(double)) != 0 ||
        std::memcmp(w[0].vec, w[1].vec, sizeof(double) * (nf + 1)) != 0 ||
        std::memcmp(out[0].vec, out[1].vec, sizeof(double) * m) != 0) {
      return false;
    }
  }
  return true;
}

bool kernelsVerified(int kernels) {
  static std::once_flag once[3];
  static bool ok[3] = {true, false, false};
  std::call_once(once[kernels], [kernels] {
    ok[kernels] = selfCheckOnce(tronOps(kernels));
    if (!ok[kernels]) {
      std::fprintf(stderr, "Percolator TRON: the %s kernels do not reproduce the reference solver in this "
                           "build; using the reference solver.\n", tronOps(kernels).name);
    }
  });
  return ok[kernels];
}
#endif

}  // namespace

TronWorkspaceHold::TronWorkspaceHold() {
#ifdef PERCOLATOR_TRON_KERNELS
  TronWorkspacePool& pool = tronWorkspacePool();
  std::lock_guard<std::mutex> lock(pool.mutex);
  ++pool.holds;
#endif
}

TronWorkspaceHold::~TronWorkspaceHold() {
#ifdef PERCOLATOR_TRON_KERNELS
  TronWorkspacePool& pool = tronWorkspacePool();
  std::lock_guard<std::mutex> lock(pool.mutex);
  if (--pool.holds == 0) {  // the large workspaces are freed
    pool.idle.erase(std::remove_if(pool.idle.begin(), pool.idle.end(),
                                   [](const std::unique_ptr<TronWorkspace>& ws) {
                                     return ws->bytes() > kTronWorkspaceKeepBytes;
                                   }),
                    pool.idle.end());
  }
#endif
}

bool tronSelectKernels(const char* isa) {
  int k;
  if (isa == nullptr) k = defaultTronKernels();
  else if (std::strcmp(isa, "ref") == 0) k = kTronRef;
#ifdef PERCOLATOR_TRON_KERNELS
  else if (std::strcmp(isa, "sse2") == 0) k = kTronSse2;
  else if (std::strcmp(isa, "avx2") == 0 && avx2Available()) k = kTronAvx2;
#endif
  else return false;
#ifdef PERCOLATOR_TRON_KERNELS
  if (isa != nullptr && k != kTronRef && !kernelsVerified(k)) return false;
#endif
  tronKernels().store(k, std::memory_order_relaxed);
  return true;
}

#ifdef PERCOLATOR_TRON_KERNELS
// Whether the kernels take a problem of m rows and n columns: with one row, Eigen computes X * w as a dot
// product (another order); the kernels index rows and columns with int, padded to multiples of 32 and 4, and
// form row offsets of a few superblocks ahead (prefetch), so both stay far below INT_MAX.
static bool fastDimensions(int m, int n) {
  return m >= 2 && n >= 1 && m <= std::numeric_limits<int>::max() / 2 && n <= (1 << 16);
}
#endif

bool tronEvaluateForTest(const char* isa, const AlgIn& set, double cpos, double cneg,
                         double lambda, const double* w, const double* s,
                         double* f, double* g, double* Hs) {
  const int m = set.m, n = set.n;
  std::vector<double> zero(n, 0.0), wv(w, w + n), sv(s, s + n);
  auto evaluate = [&](function& fn) {
    f[0] = fn.fun(zero.data());
    fn.grad(zero.data(), g);
    fn.Hv(sv.data(), Hs);
    f[1] = fn.fun(wv.data());
    fn.grad(wv.data(), g + n);
    fn.Hv(sv.data(), Hs + n);
  };
  const double inv_2lambda = 1.0 / (2.0 * lambda);
  if (std::strcmp(isa, "ref") == 0) {
    DenseRowMatrix X(m, n);
    Eigen::VectorXd t(m), C(m);
    for (int i = 0; i < m; ++i) {
      for (int j = 0; j < n - 1; ++j) X(i, j) = set.vals[i][j];
      X(i, n - 1) = 1.0;
      t(i) = set.Y[i];
      C(i) = (set.Y[i] > 0.0 ? cpos : cneg) * inv_2lambda;
    }
    TronFunction fn(X, t, C);
    evaluate(fn);
    return true;
  }
#ifdef PERCOLATOR_TRON_KERNELS
  int k = -1;
  if (std::strcmp(isa, "sse2") == 0) k = kTronSse2;
  else if (std::strcmp(isa, "avx2") == 0 && avx2Available()) k = kTronAvx2;
  if (k < 0 || !fastDimensions(m, n)) return false;
  const TronKernels::Ops& ops = tronOps(k);
  TronWorkspace ws;
  const int np = (n + 3) & ~3, mp = (m + 31) & ~31;
  double* XB = ws.xb.get((size_t)mp * (size_t)np);
  ops.build_bt(set.vals, m, n - 1, np, XB);
  ws.t.assign(mp, 0.0);
  ws.C.assign(mp, 0.0);
  for (int i = 0; i < m; ++i) {
    ws.t[i] = set.Y[i];
    ws.C[i] = (set.Y[i] > 0.0 ? cpos : cneg) * inv_2lambda;
  }
  FastTronFunction fn(ops, ws, XB, m, n, ws.t.data(), ws.C.data());
  evaluate(fn);
  return true;
#else
  return false;
#endif
}

int tron(const AlgIn& set, options& opts,
         vector_double& Weights, vector_double& Outputs,
         double cpos, double cneg, double* gradRatio) {
  if (opts.lambda <= 0.0) {
    std::fprintf(stderr,
                 "TRON adapter: opts.lambda must be > 0 (got %g); aborting.\n",
                 opts.lambda);
    return 0;
  }

#ifdef PERCOLATOR_TRON_KERNELS
  const int kernels = tronKernels().load(std::memory_order_relaxed);
  const int m = set.m, n = Weights.d;
  if (kernels != kTronRef && fastDimensions(m, n) && kernelsVerified(kernels)) {
    return tronFast(tronOps(kernels), set, opts, Weights, Outputs, cpos, cneg, gradRatio);
  }
#endif
  return tronRef(set, opts, Weights, Outputs, cpos, cneg, gradRatio);
}

}}}  // namespace OpenMS::Internal::Percolator

/* See ssl.h for the contract. */
#include "ssl.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <vector>

#include <Eigen/Dense>

#include "tron.h"

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

}  // namespace

int tron(const AlgIn& set, options& opts,
         vector_double& Weights, vector_double& Outputs,
         double cpos, double cneg) {
  if (opts.lambda <= 0.0) {
    std::fprintf(stderr,
                 "TRON adapter: opts.lambda must be > 0 (got %g); aborting.\n",
                 opts.lambda);
    return 0;
  }

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
  const double g0_norm = g0.norm();
  if (g0_norm == 0.0) return 1;  // pathological: no signal — declare convergence
  const double ratio = g_final.norm() / g0_norm;
  return (ratio <= opts.epsilon) ? 1 : 0;
}

}}}  // namespace OpenMS::Internal::Percolator

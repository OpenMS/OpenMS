/*
 * OpenMS-local (not upstream Percolator, not in whitelist.txt): solvers for the
 * problem of tron() (ssl.h) that the cross-validation can choose instead of TRON
 * (CrossValidation::setSolver).
 *
 *  - newtonSvm: the modified finite Newton method for the L2-loss linear SVM of
 *    Keerthi & DeCoste (2005), "A Modified Finite Newton Method for Fast Solution
 *    of Large Scale Linear SVMs", JMLR 6:341-361 (L2-SVM-MFN), written from the
 *    paper for few features (n <= kNewtonMaxColumns): the regularized least squares
 *    problem of each step is solved with the explicit Hessian (one pass over the
 *    active rows forms X_I^T D X_I) and a Cholesky factorization instead of CG, and
 *    the step length comes from the exact line search on the sorted breakpoints of
 *    the piecewise quadratic objective (section 4.4 of the paper).
 *  - ldaSvm: Fisher's linear discriminant of the two classes (pooled within-class
 *    scatter), scaled and shifted to minimize the objective of tron() along it with
 *    all rows active (the beta-seeding of section 4.1 of the paper), from one pass.
 *
 * Both minimize, as tron():
 *   F(w) = (1/2)||w||^2 + sum_i C_i max(0, 1 - t_i w.x_i)^2,  C_i = c_i / (2 lambda),
 * with x_i = (features of row i, 1) (the bias is the last weight and regularized).
 * The rows are cut into fixed blocks whose partial sums are added in block order:
 * the result does not depend on the number of threads. Within a parallel region
 * (OpenMP), the blocks of large problems run as tasks of the team.
 *
 * Apache-2.0 (Percolator codebase).
 */
#ifndef PERCOLATOR_NEWTON_SVM_H
#define PERCOLATOR_NEWTON_SVM_H

#include "ssl.h"

namespace OpenMS { namespace Internal { namespace Percolator {

// the solver of the cross-validation's SVMs (CrossValidation::setSolver)
enum SvmSolver {
  kSolverTron = 0,    // tron() (default)
  kSolverNewton = 1,  // newtonSvm(), warm-started from the pair's weights of the step before
  kSolverLda = 2      // ldaSvm(), one pass per fold and step
};

// the largest number of columns (features + bias) newtonSvm() and ldaSvm() take
constexpr int kNewtonMaxColumns = 128;

struct NewtonStats {
  int iterations = 0;     // Newton steps
  int passes = 0;         // passes over all rows (margins; the Hessian passes read the active rows only)
  long long hessianRows = 0;  // active rows summed over the Hessian passes
  double gradRatio = 0.0;     // ||g(w)|| / ||g(0)|| at the end
};

// Minimizes F (above) for the costs cpos, cneg and opts.lambda by Newton's method,
// from w (set.n weights; all zero: the first step solves the least squares problem
// on all rows, i.e. the LS-SVM / Fisher discriminant). Ends when ||g(w)|| <=
// eps * ||g(0)|| (TRON's criterion) or after opts.mfnitermax steps. On return w holds
// the solution and outputs (if not null, set.m values) X w. Returns 1 if the
// criterion was met, 0 if not, and -1 (w unchanged) if the problem is not supported
// (n > kNewtonMaxColumns, m < 1, opts.lambda <= 0) or the factorization failed.
int newtonSvm(const AlgIn& set, const options& opts, double cpos, double cneg, double eps,
              double* w, double* outputs, NewtonStats* stats = nullptr);

// Fisher's discriminant of the rows with Y > 0 against those with Y < 0, scaled and
// shifted for F (above) with the costs cpos, cneg and opts.lambda. Returns false
// (w unchanged) if a class has fewer than two rows, n > kNewtonMaxColumns or
// opts.lambda <= 0; w then needs another solver.
bool ldaSvm(const AlgIn& set, const options& opts, double cpos, double cneg, double* w);

// F(w) (above), for tests and diagnostics (one pass, in row order)
double svmObjective(const AlgIn& set, double lambda, double cpos, double cneg, const double* w);

}}}  // namespace OpenMS::Internal::Percolator

#endif

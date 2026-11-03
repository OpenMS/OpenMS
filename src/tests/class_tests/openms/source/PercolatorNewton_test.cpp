// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Oliver Kohlbacher $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

// Vendored Percolator (OpenMS::Internal::Percolator): this test links the
// static OpenMS_Percolator library directly, its symbols are hidden in libOpenMS.
#include "newton_svm.h"

#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace OpenMS;
using namespace std;
using namespace OpenMS::Internal::Percolator;

namespace
{
  // deterministic pseudo-random numbers (the same on every platform)
  struct Rng
  {
    uint64_t s;
    double uniform() // [0, 1)
    {
      s = s * 6364136223846793005ULL + 1442695040888963407ULL;
      return double(s >> 11) * (1.0 / 9007199254740992.0);
    }
    double normal() { return uniform() + uniform() + uniform() + uniform() - 2.0; }
  };

  // kind 0: overlapping classes; 1: separable, large margins; 2: a constant and a zero feature,
  // duplicated rows, mixed scales; 3: negatives only
  struct Problem
  {
    int m, nf;
    vector<double> vals;
    vector<double> y;
    vector<double*> rows;
    AlgIn in;
    Problem(int m_, int nf_, int kind, uint64_t seed)
      : m(m_), nf(nf_), vals(size_t(m_) * nf_), y(m_), rows(m_), in(m_, nf_ + 1)
    {
      Rng r{seed};
      for (int i = 0; i < m; ++i)
      {
        const double label = (kind != 3 && r.uniform() < 0.4) ? 1.0 : -1.0;
        y[i] = label;
        for (int j = 0; j < nf; ++j)
        {
          double v = r.normal();
          if (kind == 0 || kind == 3) v += (j % 3 == 0 ? 0.4 : 0.0) * label;
          else if (kind == 1) v = (j == 0 ? 50.0 * label : 0.01 * v);
          else
          {
            if (j == 1) v = 1.0;
            else if (j == 2 % nf) v = 0.0;
            else v *= (j % 2 == 0 ? 1e3 : 1e-3);
            v += 0.3 * label;
          }
          vals[size_t(i) * nf + j] = v;
        }
        if (kind == 2 && i % 7 == 6)
        {
          y[i] = y[i - 1];
          for (int j = 0; j < nf; ++j) vals[size_t(i) * nf + j] = vals[size_t(i - 1) * nf + j];
        }
      }
      in.m = m;
      for (int i = 0; i < m; ++i)
      {
        in.vals[i] = vals.data() + size_t(i) * nf;
        in.Y[i] = y[i];
        if (y[i] > 0) ++in.positives;
        else ++in.negatives;
      }
    }
  };

  options makeOptions(double lambda = 1.0)
  {
    options opts;
    opts.lambda = lambda;
    opts.lambda_u = 1.0;
    opts.epsilon = EPSILON;
    opts.cgitermax = CGITERMAX;
    opts.mfnitermax = MFNITERMAX;
    return opts;
  }

  vector<double> tronWeights(const Problem& p, double cpos, double cneg, double lambda = 1.0)
  {
    options opts = makeOptions(lambda);
    vector_double w, out;
    w.d = p.nf + 1;
    w.vec = new double[w.d]();
    out.d = p.m;
    out.vec = new double[out.d]();
    tron(p.in, opts, w, out, cpos, cneg);
    return vector<double>(w.vec, w.vec + w.d);
  }

  // Newton from start (zero: cold) on threads threads (blocks as tasks of a parallel region if > 1)
  int newton(const Problem& p, double cpos, double cneg, vector<double>& w, int threads = 1, double eps = EPSILON,
             vector<double>* out = nullptr, NewtonStats* st = nullptr, double lambda = 1.0)
  {
    options opts = makeOptions(lambda);
    w.resize(p.nf + 1, 0.0);
    if (out) out->assign(p.m, 0.0);
    int ret = -2;
#ifdef _OPENMP
    if (threads > 1)
    {
#pragma omp parallel num_threads(threads)
#pragma omp single
      ret = newtonSvm(p.in, opts, cpos, cneg, eps, w.data(), out ? out->data() : nullptr, st);
      return ret;
    }
#endif
    (void)threads;
    return newtonSvm(p.in, opts, cpos, cneg, eps, w.data(), out ? out->data() : nullptr, st);
  }

  // a problem with the given rows of one feature
  struct Tiny
  {
    vector<double> x, y;
    AlgIn in;
    Tiny(vector<double> x_, vector<double> y_) : x(std::move(x_)), y(std::move(y_)), in(int(x.size()), 2)
    {
      in.m = int(x.size());
      for (int i = 0; i < in.m; ++i)
      {
        in.vals[i] = &x[i];
        in.Y[i] = y[i];
      }
    }
  };

  bool sameBits(const vector<double>& a, const vector<double>& b)
  {
    return a.size() == b.size() && std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0;
  }
}

START_TEST(PercolatorNewton, "$Id$")

/////////////////////////////////////////////////////////////

START_SECTION((int newtonSvm(const AlgIn& set, const options& opts, double cpos, double cneg, double eps, double* w, double* outputs, NewtonStats* stats)))
{
  // the optimum of TRON (to EPSILON) within 1e-6 of the objective, on small to medium problems of every kind
  const int sizes[][2] = {{1, 3}, {2, 1}, {40, 4}, {500, 9}, {3000, 30}, {20000, 31}};
  const double costs[][2] = {{0.1, 0.1}, {0.1, 0.4}, {1.0, 10.0}, {10.0, 0.3}};
  int checked = 0, better = 0;
  double maxRel = -1.0;
  for (const auto& sz : sizes)
  {
    for (int kind = 0; kind < 4; ++kind)
    {
      Problem p(sz[0], sz[1], kind, 17 + 31 * kind + sz[0]);
      for (const auto& c : costs)
      {
        for (double lambda : {1.0, 0.01})
        {
          const vector<double> wt = tronWeights(p, c[0], c[1], lambda);
          vector<double> wn, out;
          NewtonStats st;
          const int ret = newton(p, c[0], c[1], wn, 1, EPSILON, &out, &st, lambda);
          TEST_EQUAL(ret >= 0, true)
          const double ft = svmObjective(p.in, lambda, c[0], c[1], wt.data());
          const double fn = svmObjective(p.in, lambda, c[0], c[1], wn.data());
          const double rel = (fn - ft) / ft;
          maxRel = std::max(maxRel, rel);
          if (rel < -1e-6) ++better;
          // never worse than TRON (beyond 1e-6); on the well-conditioned problems the same objective
          TEST_EQUAL(rel <= 1e-6, true)
          if (kind == 0) TEST_EQUAL(std::fabs(rel) <= 1e-6, true)
          TEST_EQUAL(ret == 1 && st.gradRatio <= EPSILON, true)
          // and within 1e-6 of the optimum (Newton to a tolerance of 1e-14 from there); with features
          // scaled by 1e3 and 1e-3 or one class only, the gradient criterion leaves up to 1e-4 (as TRON's)
          vector<double> wo = wn;
          options tight = makeOptions(lambda);
          newtonSvm(p.in, tight, c[0], c[1], 1e-14, wo.data(), nullptr);
          const double fo = svmObjective(p.in, lambda, c[0], c[1], wo.data());
          TEST_EQUAL((fn - fo) / fo <= (kind >= 2 ? 1e-4 : 1e-6), true)
          // outputs = X w
          double maxErr = 0.0;
          for (int i = 0; i < p.m; ++i)
          {
            double s = wn[p.nf];
            for (int j = 0; j < p.nf; ++j) s += p.vals[size_t(i) * p.nf + j] * wn[j];
            maxErr = std::max(maxErr, std::fabs(s - out[i]) / (1.0 + std::fabs(s)));
          }
          TEST_EQUAL(maxErr < 1e-9, true)
          ++checked;
        }
      }
    }
  }
  STATUS("solves " << checked << ", max rel. objective difference to TRON " << maxRel << ", Newton better by more than 1e-6: " << better)

  // warm start from the solution: no step, the same weights
  {
    Problem p(3000, 30, 0, 5);
    vector<double> w;
    newton(p, 0.1, 0.3, w);
    vector<double> w2 = w;
    NewtonStats st;
    options opts = makeOptions();
    TEST_EQUAL(newtonSvm(p.in, opts, 0.1, 0.3, EPSILON, w2.data(), nullptr, &st), 1)
    TEST_EQUAL(st.iterations, 0)
    TEST_EQUAL(sameBits(w, w2), true)
    // warm start from the solution of other costs: few steps, the optimum of TRON
    vector<double> w3 = w;
    TEST_EQUAL(newtonSvm(p.in, opts, 0.1, 1.0, EPSILON, w3.data(), nullptr, &st), 1)
    const vector<double> wt = tronWeights(p, 0.1, 1.0);
    const double ft = svmObjective(p.in, 1.0, 0.1, 1.0, wt.data());
    TEST_EQUAL(std::fabs(svmObjective(p.in, 1.0, 0.1, 1.0, w3.data()) - ft) <= 1e-6 * ft, true)
  }

  // the result does not depend on the threads: blocks as tasks (m >= 8 blocks of 8192 rows)
  {
    Problem p(70001, 31, 0, 99);
    vector<double> w1, w4, w7, o1, o4;
    TEST_EQUAL(newton(p, 0.1, 0.3, w1, 1, EPSILON, &o1) >= 0, true)
    TEST_EQUAL(newton(p, 0.1, 0.3, w4, 4, EPSILON, &o4) >= 0, true)
    TEST_EQUAL(newton(p, 0.1, 0.3, w7, 7) >= 0, true)
    TEST_EQUAL(sameBits(w1, w4), true)
    TEST_EQUAL(sameBits(w1, w7), true)
    TEST_EQUAL(sameBits(o1, o4), true)
    const vector<double> wt = tronWeights(p, 0.1, 0.3);
    const double ft = svmObjective(p.in, 1.0, 0.1, 0.3, wt.data());
    TEST_EQUAL(std::fabs(svmObjective(p.in, 1.0, 0.1, 0.3, w4.data()) - ft) <= 1e-6 * ft, true)
  }

  // the iteration limit: 0 (not converged), the last iterate
  {
    Tiny t({1.0, -1.0}, {1.0, -1.0});
    options opts = makeOptions();
    opts.mfnitermax = 0;
    vector<double> w(2, 0.0);
    TEST_EQUAL(newtonSvm(t.in, opts, 1.0, 1.0, EPSILON, w.data(), nullptr), 0)
    opts.mfnitermax = MFNITERMAX;
    TEST_EQUAL(newtonSvm(t.in, opts, 1.0, 1.0, EPSILON, w.data(), nullptr), 1)
    TEST_REAL_SIMILAR(w[0], 2.0 / 3.0)
    // huge costs (squared gradient components overflow): not "converged" at 0, the optimum (1, 0)
    vector<double> v(2, 0.0);
    TEST_EQUAL(newtonSvm(t.in, opts, 1e160, 1e160, EPSILON, v.data(), nullptr), 1)
    TEST_REAL_SIMILAR(v[0], 1.0)
    // huge costs and a warm start that I + H would hide in the gradient: the optimum all the same
    Tiny q({1.0, 1.0, 1.0, 1.0}, {1.0, -1.0, 1.0, -1.0});
    vector<double> z(2, 0.0), h = {std::ldexp(1.0, 27), -std::ldexp(1.0, 27)};
    const double huge = std::ldexp(1.0, 52);
    TEST_EQUAL(newtonSvm(q.in, opts, huge, huge, EPSILON, z.data(), nullptr) >= 0, true)
    const int rh = newtonSvm(q.in, opts, huge, huge, EPSILON, h.data(), nullptr);
    if (rh >= 0)
    {
      const double fz = svmObjective(q.in, 1.0, huge, huge, z.data());
      TEST_EQUAL(svmObjective(q.in, 1.0, huge, huge, h.data()) <= fz * (1.0 + 1e-6), true)
    }
    // infinite costs: -1, w unchanged
    vector<double> u(2, 0.0);
    TEST_EQUAL(newtonSvm(t.in, opts, std::numeric_limits<double>::infinity(), 1.0, EPSILON, u.data(), nullptr), -1)
    TEST_EQUAL(u[0], 0.0)
  }

  // not supported: w unchanged, -1
  {
    Problem p(50, 130, 0, 3);
    vector<double> w(131, 0.5);
    options opts = makeOptions();
    TEST_EQUAL(newtonSvm(p.in, opts, 1.0, 1.0, EPSILON, w.data(), nullptr), -1)
    TEST_EQUAL(w[0], 0.5)
    Problem q(50, 4, 0, 3);
    vector<double> v(5, 0.0);
    opts.lambda = 0.0;
    TEST_EQUAL(newtonSvm(q.in, opts, 1.0, 1.0, EPSILON, v.data(), nullptr), -1)
  }
}
END_SECTION

START_SECTION((bool ldaSvm(const AlgIn& set, const options& opts, double cpos, double cneg, double* w)))
{
  options opts = makeOptions();
  // overlapping classes: Fisher's direction ranks like the SVM (most predictions agree in sign)
  Problem p(20000, 30, 0, 11);
  vector<double> w(31, 0.0);
  TEST_EQUAL(ldaSvm(p.in, opts, 0.1, 0.3, w.data()), true)
  const vector<double> wt = tronWeights(p, 0.1, 0.3);
  double dotLT = 0.0, nl = 0.0, nt = 0.0;
  for (int j = 0; j < 30; ++j)
  {
    dotLT += w[j] * wt[j];
    nl += w[j] * w[j];
    nt += wt[j] * wt[j];
  }
  TEST_EQUAL(dotLT / std::sqrt(nl * nt) > 0.9, true)
  // the same with the blocks as tasks
#ifdef _OPENMP
  vector<double> w4(31, 0.0);
  bool ok4 = false;
#pragma omp parallel num_threads(4)
#pragma omp single
  ok4 = ldaSvm(p.in, opts, 0.1, 0.3, w4.data());
  TEST_EQUAL(ok4, true)
  TEST_EQUAL(sameBits(w, w4), true)
#endif
  // a constant and a zero feature (ridge)
  Problem q(3000, 9, 2, 4);
  vector<double> v(10, 0.0);
  TEST_EQUAL(ldaSvm(q.in, opts, 1.0, 1.0, v.data()), true)
  for (double x : v) TEST_EQUAL(std::isfinite(x), true)
  // the scale is positive (Fisher's ranking) for any costs, or LDA declines
  {
    Tiny t({-2.1, -1.9, -3.1, -2.9}, {1.0, 1.0, -1.0, -1.0});
    for (double cneg : {0.05, 0.1, 1.0})
    {
      vector<double> v(2, 0.0);
      if (ldaSvm(t.in, opts, 0.1, cneg, v.data())) TEST_EQUAL(v[0] > 0.0, true)
    }
  }
  // fewer than two rows of a class: false, w unchanged
  Problem r(500, 9, 3, 4);
  vector<double> u(10, 0.25);
  TEST_EQUAL(ldaSvm(r.in, opts, 1.0, 1.0, u.data()), false)
  TEST_EQUAL(u[3], 0.25)
}
END_SECTION

START_SECTION((double svmObjective(const AlgIn& set, double lambda, double cpos, double cneg, const double* w)))
{
  Problem p(3, 1, 0, 1);
  // rows x = (v, 1); w = 0: every margin 0, F = sum C_i with C_i = c_i / (2 lambda)
  vector<double> w(2, 0.0);
  double expected = 0.0;
  for (int i = 0; i < 3; ++i) expected += (p.y[i] > 0 ? 0.5 : 2.0) / 2.0;
  TEST_REAL_SIMILAR(svmObjective(p.in, 1.0, 0.5, 2.0, w.data()), expected)
}
END_SECTION

/////////////////////////////////////////////////////////////
END_TEST

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
#include "ssl.h"

#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

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

  struct Problem
  {
    int m, nf;
    vector<double> vals; // m x nf
    vector<double> y;
  };

  // kind 0: overlapping classes; 1: separable, large margins (rows leave the active set);
  // 2: a constant and a zero feature, duplicated rows, mixed scales
  Problem makeProblem(int m, int nf, int kind, uint64_t seed)
  {
    Rng r{seed};
    Problem p{m, nf, vector<double>(size_t(m) * nf), vector<double>(m)};
    for (int i = 0; i < m; ++i)
    {
      const double label = (r.uniform() < 0.4) ? 1.0 : -1.0;
      p.y[i] = label;
      for (int j = 0; j < nf; ++j)
      {
        double v = r.normal();
        if (kind == 0) v += (j % 3 == 0 ? 0.4 : 0.0) * label;
        else if (kind == 1) v = (j == 0 ? 50.0 * label : 0.01 * v);
        else
        {
          if (j == 1) v = 1.0;
          else if (j == 2 % nf) v = 0.0;
          else v *= (j % 2 == 0 ? 1e3 : 1e-3);
          v += 0.3 * label;
        }
        p.vals[size_t(i) * nf + j] = v;
      }
      if (kind == 2 && i % 7 == 6) // duplicate of the previous row
      {
        p.y[i] = p.y[i - 1];
        for (int j = 0; j < nf; ++j) p.vals[size_t(i) * nf + j] = p.vals[size_t(i - 1) * nf + j];
      }
    }
    return p;
  }

  struct Result
  {
    vector<double> w, out;
    int converged;
    double ratio;
  };

  Result solve(const Problem& p, double cpos, double cneg, double lambda = 1.0)
  {
    AlgIn in(p.m, p.nf + 1);
    in.m = p.m;
    for (int i = 0; i < p.m; ++i)
    {
      in.vals[i] = const_cast<double*>(p.vals.data() + size_t(i) * p.nf);
      in.Y[i] = p.y[i];
    }
    options opts;
    opts.lambda = lambda;
    opts.lambda_u = 1.0;
    opts.epsilon = EPSILON;
    opts.cgitermax = CGITERMAX;
    opts.mfnitermax = MFNITERMAX;
    vector_double w, out;
    w.d = p.nf + 1;
    w.vec = new double[w.d]();
    out.d = p.m;
    out.vec = new double[out.d]();
    Result res;
    res.ratio = -1.0;
    res.converged = tron(in, opts, w, out, cpos, cneg, &res.ratio);
    res.w.assign(w.vec, w.vec + w.d);
    res.out.assign(out.vec, out.vec + out.d);
    return res;
  }

  bool sameBits(const vector<double>& a, const vector<double>& b)
  {
    return a.size() == b.size() && (a.empty() || std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0);
  }

  bool sameResult(const Result& a, const Result& b)
  {
    return sameBits(a.w, b.w) && sameBits(a.out, b.out) && a.converged == b.converged
           && std::memcmp(&a.ratio, &b.ratio, sizeof(double)) == 0;
  }

  // restores the default kernel selection, also when a check fails
  struct KernelGuard
  {
    ~KernelGuard() { tronSelectKernels(nullptr); }
  };

  // the kernels available here; on x86-64 baseline builds (where ssl.cpp uses them) SSE2 must be
  vector<string> availableKernels()
  {
    KernelGuard guard;
    vector<string> isas;
    for (const char* isa : {"sse2", "avx2"})
    {
      if (tronSelectKernels(isa)) isas.push_back(isa);
    }
    return isas;
  }
}

START_TEST(PercolatorTron, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION([EXTRA] tronSelectKernels)
{
  TEST_TRUE(tronSelectKernels("ref"))
  TEST_FALSE(tronSelectKernels("no-such-isa"))
  TEST_TRUE(tronSelectKernels(nullptr))
}
END_SECTION

START_SECTION([EXTRA] kernels available)
{
  const vector<string> isas = availableKernels();
  for (const string& isa : isas)
  {
    STATUS("kernels " << isa << " available (self-check passed)")
  }
#if (defined(__x86_64__) || defined(_M_X64)) && !defined(__AVX__) && !defined(__FMA__)
  // the build in which tron() uses the kernels: SSE2 must be there and pass the self-check
  TEST_TRUE(!isas.empty() && isas[0] == "sse2")
#endif
}
END_SECTION

START_SECTION([EXTRA] tron: the kernels give the bits of the reference function object)
{
  // On x86-64 baseline builds, tron() runs on SSE2 or AVX2 kernels (tron_kernels.h); weights,
  // outputs and the gradient-norm ratio must equal those of the reference TronFunction bit for bit.
  KernelGuard guard;
  const vector<string> isas = availableKernels();
  const vector<int> ns = {1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 13, 16, 17, 23, 32, 33, 41, 44, 60}; // incl. the bias column
  const vector<int> ms = {1, 2, 3, 5, 7, 31, 32, 33, 63, 100, 257, 1001, 2049, 4999};
  const double costs[][2] = {{0.1, 0.3}, {1.0, 3.0}, {30.0, 10.0}, {0.0, 0.0}};
  int problems = 0, identical = 0;
  for (int kind = 0; kind < 3; ++kind)
  {
    for (size_t a = 0; a < ns.size(); ++a)
    {
      for (size_t b = 0; b < ms.size(); ++b)
      {
        if ((a + b + kind) % 2 == 1 && ms[b] > 300) continue; // keep the run time short
        const Problem p = makeProblem(ms[b], ns[a] - 1, kind, 1000003ULL * a + 7919ULL * b + kind);
        const double* c = costs[(a + b + kind) % 4];
        tronSelectKernels("ref");
        const Result ref = solve(p, c[0], c[1]);
        for (const string& isa : isas)
        {
          tronSelectKernels(isa.c_str());
          const Result r = solve(p, c[0], c[1]);
          ++problems;
          const bool same = sameResult(r, ref);
          identical += same;
          if (!same)
          {
            STATUS("differs: " << isa << " kind " << kind << " n " << ns[a] << " m " << ms[b])
          }
        }
      }
    }
  }
  STATUS(identical << " of " << problems << " problems identical")
  TEST_EQUAL(identical, problems)
}
END_SECTION

START_SECTION([EXTRA] tron: extreme costs and NaN or subnormal values)
{
  // Whole solves with extreme or NaN values must give the bits of the reference as well. (With an
  // infinite value, huge values or an overflowing gradient, e.g. 3 rows x = 1e10, lambda = 1e-300,
  // the reference's conjugate gradient never terminates; those are checked function by function
  // in the next section.)
  KernelGuard guard;
  const vector<string> isas = availableKernels();
  struct Case { const char* name; Problem p; double cpos, cneg, lambda; };
  vector<Case> cases;
  for (int nf : {1, 4, 9})
  {
    for (int m : {3, 40, 333})
    {
      Problem p = makeProblem(m, nf, 0, 77 + m + nf);
      cases.push_back({"huge costs", p, 1.0, 2.0, 1e-300});
      cases.push_back({"tiny costs", p, 1e-300, 1e-300, 1.0});
      Problem nan = p;
      nan.vals[size_t(m - 1) * nf] = std::numeric_limits<double>::quiet_NaN();
      cases.push_back({"NaN value", nan, 1.0, 3.0, 1.0});
      Problem subn = p;
      for (double& v : subn.vals) v *= 1e-310;
      cases.push_back({"subnormal values", subn, 1.0, 3.0, 1.0});
    }
  }
  int problems = 0, identical = 0;
  for (const Case& c : cases)
  {
    tronSelectKernels("ref");
    const Result ref = solve(c.p, c.cpos, c.cneg, c.lambda);
    for (const string& isa : isas)
    {
      tronSelectKernels(isa.c_str());
      const Result r = solve(c.p, c.cpos, c.cneg, c.lambda);
      ++problems;
      const bool same = sameResult(r, ref);
      identical += same;
      if (!same)
      {
        STATUS("differs: " << isa << " " << c.name << " m " << c.p.m << " features " << c.p.nf)
      }
    }
  }
  STATUS(identical << " of " << problems << " problems identical")
  TEST_EQUAL(identical, problems)
}
END_SECTION

START_SECTION([EXTRA] function object: f and gradient and Hv with infinite or NaN or huge values)
{
  // fun, grad and Hv as TRON calls them, kernels against the reference, with directions s that hold
  // Inf or NaN (an overflowing gradient) and weights that leave no row, some rows or all rows active
  const vector<string> isas = availableKernels();
  const double inf = std::numeric_limits<double>::infinity(), nan = std::numeric_limits<double>::quiet_NaN();
  int evaluations = 0, identical = 0;
  for (int n = 1; n <= 13; ++n)
  {
    for (int m : {2, 3, 5, 8, 33, 70})
    {
      for (int kind = 0; kind < 3; ++kind)
      {
        Problem p = makeProblem(m, n - 1, kind, 31 * n + m + kind);
        if (kind == 2 && n > 1) p.vals[0] = 1e300; // a huge value
        AlgIn in(p.m, n);
        in.m = p.m;
        for (int i = 0; i < p.m; ++i)
        {
          in.vals[i] = p.vals.data() + size_t(i) * (n - 1);
          in.Y[i] = p.y[i];
        }
        Rng r{uint64_t(7 * n + 13 * m + kind)};
        for (int variant = 0; variant < 6; ++variant)
        {
          vector<double> w(n), s(n);
          for (int j = 0; j < n; ++j)
          {
            w[j] = r.normal() * (variant < 2 ? 0.1 : variant < 4 ? 1.0 : 1e3); // 1e3: few or no active rows
            s[j] = r.normal();
          }
          if (variant % 2 == 1) s[(variant / 2) % n] = (variant == 1 ? inf : variant == 3 ? -inf : nan);
          if (variant == 4) s[n - 1] = 1e308;
          vector<double> f_ref(2), g_ref(2 * n), h_ref(2 * n);
          tronEvaluateForTest("ref", in, 1.0, 3.0, 1.0, w.data(), s.data(), f_ref.data(), g_ref.data(), h_ref.data());
          for (const string& isa : isas)
          {
            vector<double> f(2), g(2 * n), h(2 * n);
            TEST_TRUE(tronEvaluateForTest(isa.c_str(), in, 1.0, 3.0, 1.0, w.data(), s.data(), f.data(), g.data(), h.data()))
            ++evaluations;
            const bool same = sameBits(f, f_ref) && sameBits(g, g_ref) && sameBits(h, h_ref);
            identical += same;
            if (!same)
            {
              STATUS("differs: " << isa << " n " << n << " m " << m << " kind " << kind << " variant " << variant)
            }
          }
        }
      }
    }
  }
  // the overflowing gradient of 3 rows x = 1e10 with C = 5e299 (TRON's direction is then -g = +Inf)
  {
    vector<double> vals = {1e10, 1e10, 1e10};
    AlgIn in(3, 2);
    in.m = 3;
    const double y[3] = {1.0, -1.0, 1.0};
    for (int i = 0; i < 3; ++i) { in.vals[i] = vals.data() + i; in.Y[i] = y[i]; }
    vector<double> w = {0.0, 0.0}, s = {inf, inf}, f_ref(2), g_ref(4), h_ref(4);
    tronEvaluateForTest("ref", in, 1.0, 1.0, 1e-300, w.data(), s.data(), f_ref.data(), g_ref.data(), h_ref.data());
    TEST_FALSE(std::isfinite(g_ref[0]))
    for (const string& isa : isas)
    {
      vector<double> f(2), g(4), h(4);
      tronEvaluateForTest(isa.c_str(), in, 1.0, 1.0, 1e-300, w.data(), s.data(), f.data(), g.data(), h.data());
      ++evaluations;
      identical += sameBits(f, f_ref) && sameBits(g, g_ref) && sameBits(h, h_ref);
    }
  }
  STATUS(identical << " of " << evaluations << " evaluations identical")
  TEST_EQUAL(identical, evaluations)
}
END_SECTION

START_SECTION([EXTRA] tron: concurrent solves)
{
  // per-thread workspaces: solves on several threads at once give the bits of the serial reference
  KernelGuard guard;
  const vector<string> isas = availableKernels();
  vector<Problem> problems;
  for (int k = 0; k < 24; ++k) problems.push_back(makeProblem(200 + 397 * (k % 7), 4 + (k % 11) * 3, k % 3, 9000 + k));
  tronSelectKernels("ref");
  vector<Result> ref(problems.size());
  for (size_t k = 0; k < problems.size(); ++k) ref[k] = solve(problems[k], 1.0, 3.0);
  for (const string& isa : isas)
  {
    tronSelectKernels(isa.c_str());
    vector<Result> res(problems.size());
#pragma omp parallel for schedule(dynamic, 1)
    for (int k = 0; k < (int)problems.size(); ++k) res[k] = solve(problems[k], 1.0, 3.0);
    int identical = 0;
    for (size_t k = 0; k < problems.size(); ++k) identical += sameResult(res[k], ref[k]);
    TEST_EQUAL(identical, (int)problems.size())
  }
}
END_SECTION

START_SECTION([EXTRA] tron: kernels after a larger and before a smaller problem on the same thread)
{
  KernelGuard guard;
  // the per-thread buffers grow and are reused; results must not depend on what ran before
  const Problem big = makeProblem(3000, 30, 0, 42);
  const Problem small = makeProblem(77, 11, 0, 43);
  tronSelectKernels("ref");
  const Result ref_small = solve(small, 1.0, 3.0);
  const Result ref_big = solve(big, 1.0, 3.0);
  for (const char* isa : {"sse2", "avx2"})
  {
    if (!tronSelectKernels(isa)) continue;
    const Result s1 = solve(small, 1.0, 3.0);
    const Result b1 = solve(big, 1.0, 3.0);
    const Result s2 = solve(small, 1.0, 3.0);
    TEST_TRUE(sameBits(s1.w, ref_small.w) && sameBits(s1.out, ref_small.out))
    TEST_TRUE(sameBits(b1.w, ref_big.w) && sameBits(b1.out, ref_big.out))
    TEST_TRUE(sameBits(s2.w, ref_small.w) && sameBits(s2.out, ref_small.out))
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

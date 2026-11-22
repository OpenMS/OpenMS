// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenDIAGBM.h>

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <random>

using namespace OpenMS::OpenDIA;

namespace
{
  int failures = 0;
  void check(bool ok, const std::string& what)
  {
    if (!ok)
    {
      ++failures;
      std::cerr << "FAIL: " << what << '\n';
    }
  }

  double auc(const std::vector<double>& s, const std::vector<std::uint8_t>& y)
  {
    std::vector<std::size_t> o(s.size());
    std::iota(o.begin(), o.end(), 0);
    std::sort(o.begin(), o.end(), [&](std::size_t a, std::size_t b) { return s[a] < s[b]; });
    double rank_sum = 0.0, np = 0.0;
    for (std::size_t i = 0; i < o.size(); ++i)
      if (y[o[i]])
      {
        rank_sum += static_cast<double>(i + 1);
        np += 1.0;
      }
    const double nn = static_cast<double>(s.size()) - np;
    return (rank_sum - np * (np + 1.0) / 2.0) / (np * nn);
  }

  struct Problem
  {
    std::size_t d = 0;
    std::vector<float> x;
    std::vector<std::uint8_t> y;
    std::vector<double> w;
  };

  // Nonlinear separable problem: y = 1 iff x0 * x1 > 0 or x2 > 1.5, plus noise features.
  Problem makeProblem(std::size_t n, std::size_t d, unsigned seed)
  {
    Problem p;
    p.d = d;
    p.x.resize(n * d);
    p.y.resize(n);
    p.w.assign(n, 1.0);
    std::mt19937_64 rng(seed);
    std::normal_distribution<float> nd(0.0f, 1.0f);
    for (std::size_t i = 0; i < n; ++i)
    {
      for (std::size_t j = 0; j < d; ++j) p.x[i * d + j] = nd(rng);
      p.y[i] = (p.x[i * d] * p.x[i * d + 1] > 0.0f) || p.x[i * d + 2] > 1.5f;
    }
    return p;
  }
}

int main()
{
  // (a) Separable synthetic problem.
  {
    const std::size_t n = 20000;
    auto p = makeProblem(n, 8, 1);
    std::vector<std::uint32_t> tr, va, te;
    for (std::uint32_t i = 0; i < n; ++i) (i % 10 < 7 ? tr : i % 10 < 8 ? va : te).push_back(i);
    HistGBM m;
    HistGBM::Options o;
    o.threads = 2;
    m.fit(p.x.data(), p.d, tr, va, p.y, p.w, o);
    std::vector<double> s;
    std::vector<std::uint8_t> y;
    for (auto i : te)
    {
      s.push_back(m.predict(p.x.data() + i * p.d));
      y.push_back(p.y[i]);
    }
    const double a = auc(s, y);
    std::cout << "separable AUC " << a << " trees " << m.trees.size() << '\n';
    check(a > 0.99, "separable AUC > 0.99");
  }

  // (b) Thread-count independence.
  {
    const std::size_t n = 30000;
    auto p = makeProblem(n, 12, 2);
    for (std::size_t i = 0; i < n; ++i) p.w[i] = p.y[i] ? 2.5 : 1.0;
    std::vector<std::uint32_t> tr, va;
    for (std::uint32_t i = 0; i < n; ++i) (i % 10 ? tr : va).push_back(i);
    HistGBM a, b;
    HistGBM::Options o;
    o.max_iter = 60;
    o.threads = 1;
    a.fit(p.x.data(), p.d, tr, va, p.y, p.w, o);
    o.threads = 4;
    b.fit(p.x.data(), p.d, tr, va, p.y, p.w, o);
    bool same = a.trees.size() == b.trees.size() && a.init == b.init;
    for (std::size_t i = 0; same && i < n; ++i) same = a.predict(p.x.data() + i * p.d) == b.predict(p.x.data() + i * p.d);
    check(same, "identical model for threads 1 and 4");
  }

  // (c) End-to-end gbmRescore on synthetic windows.
  {
    CrossFitScorer::Input in;
    const std::size_t pairs = 3000, wpp = 4, d = 6;
    in.dimension = d;
    std::mt19937_64 rng(7);
    std::normal_distribution<float> nd(0.0f, 1.0f);
    std::vector<float> weak;
    for (std::uint32_t pr = 0; pr < pairs; ++pr)
      for (int dec = 0; dec < 2; ++dec)
        for (std::size_t w = 0; w < wpp; ++w)
        {
          const bool signal = !dec && w == 0 && pr % 2 == 0;  // half the targets are present
          for (std::size_t j = 0; j < d; ++j) in.features.push_back(nd(rng));
          if (signal)
          {
            in.features[in.features.size() - d] += 3.0f;
            in.features[in.features.size() - d + 1] += 2.0f;
          }
          in.precursor.push_back(2 * pr + dec);
          in.pair.push_back(pr);
          in.decoy.push_back(static_cast<std::uint8_t>(dec));
          in.fold.push_back(static_cast<std::uint8_t>(pr % 3));
          in.pair_hash.push_back(CrossFitScorer::fnv1a("p" + std::to_string(pr)));
          // Weakly informative incoming score: a noisy view of feature 0.
          weak.push_back(in.features[in.features.size() - d] + 0.7f * nd(rng));
        }
    const auto count = [&](const std::vector<double>& q)
    {
      std::size_t c = 0;
      for (std::size_t r = 0; r < q.size(); ++r) c += !in.decoy[r] && q[r] <= 0.01;
      return c;
    };
    std::vector<std::uint32_t> all(in.rows());
    std::iota(all.begin(), all.end(), 0u);
    const auto comp = CrossFitScorer::pooledPairedTDC(in, weak, all);
    std::vector<double> q0(in.rows(), 1.0);
    for (std::size_t i = 0; i < comp.ranked.size(); ++i) q0[comp.ranked[i]] = comp.q[i];
    std::vector<float> flat(in.rows());
    for (auto& v : flat) v = nd(rng);
    const auto compf = CrossFitScorer::pooledPairedTDC(in, flat, all);
    std::vector<double> qf(in.rows(), 1.0);
    for (std::size_t i = 0; i < compf.ranked.size(); ++i) qf[compf.ranked[i]] = compf.q[i];
    for (double frac : {0.02, 0.0})
    {
      const auto res = gbmRescore(in, weak, frac, 1, 2, [](const std::string& l) { std::cout << l << '\n'; });
      const std::size_t c = count(res.q_value);
      std::cout << "gbmRescore frac " << frac << ": targets at q<=0.01 " << c << " (weak input " << count(q0) << ", uninformative "
                << count(qf) << ")\n";
      check(c > count(q0) && c > count(qf), "gbmRescore gains targets");
      std::size_t kept = 0;
      for (float v : res.score) kept += v > -1e29f;
      check(kept == 2 * pairs, "one kept window per precursor");
    }
  }

  // Optional timing run: OPENDIA_GBM_BENCH_ROWS=<rows> (49 features, OPENDIA_GBM_BENCH_THREADS threads, default 8).
  if (const char* env = std::getenv("OPENDIA_GBM_BENCH_ROWS"))
  {
    const std::size_t n = std::strtoull(env, nullptr, 10);
    const char* te = std::getenv("OPENDIA_GBM_BENCH_THREADS");
    auto p = makeProblem(n, 49, 3);
    std::mt19937_64 rng(9);
    std::normal_distribution<float> nd(0.0f, 1.0f);
    for (auto& v : p.x) v += 0.8f * nd(rng);  // label noise via feature noise, so boosting runs long
    std::vector<std::uint32_t> tr, va;
    for (std::uint32_t i = 0; i < n; ++i) (i % 10 ? tr : va).push_back(i);
    HistGBM m;
    HistGBM::Options o;
    o.threads = te ? std::atoi(te) : 8;
    const auto t0 = std::chrono::steady_clock::now();
    m.fit(p.x.data(), p.d, tr, va, p.y, p.w, o);
    const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::cout << "bench " << tr.size() << " x 49 threads " << o.threads << ": " << m.trees.size() << " trees in " << sec << " s\n";
  }

  if (failures) std::cerr << failures << " failure(s)\n";
  else std::cout << "all OpenDIAGBM tests passed\n";
  return failures ? EXIT_FAILURE : EXIT_SUCCESS;
}

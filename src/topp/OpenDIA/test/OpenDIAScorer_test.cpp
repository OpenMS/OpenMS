// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenDIAScorer.h>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <random>
#include <stdexcept>

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

  /// Synthetic population: `pairs` target-decoy pairs, `windows` windows per precursor, 3 folds. Targets have one true
  /// window (feature mean shift `effect` in the first two columns) if `present`; all other windows are noise.
  CrossFitScorer::Input synthetic(std::size_t pairs, std::size_t windows, double effect, std::uint32_t seed, std::vector<bool>& truth)
  {
    std::mt19937_64 rng(seed);
    std::normal_distribution<double> noise(0.0, 1.0);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    CrossFitScorer::Input in;
    in.dimension = 4;
    for (std::uint32_t p = 0; p < pairs; ++p)
    {
      const bool present = u(rng) < 0.5;
      for (int decoy = 0; decoy < 2; ++decoy)
      {
        for (std::size_t w = 0; w < windows; ++w)
        {
          const bool real = !decoy && present && w == 0;
          for (std::size_t j = 0; j < 4; ++j) in.features.push_back(static_cast<float>(noise(rng) + (real && j < 2 ? effect : 0.0)));
          in.precursor.push_back(2 * p + decoy);
          in.pair.push_back(p);
          in.decoy.push_back(static_cast<std::uint8_t>(decoy));
          in.fold.push_back(static_cast<std::uint8_t>(p % 3));
          in.pair_hash.push_back(CrossFitScorer::fnv1a("v" + std::to_string(p)));
          truth.push_back(real);
          const std::uint32_t row = static_cast<std::uint32_t>(in.precursor.size() - 1);
          // Seeds: a fifth of the real windows as positives, decoy windows as negatives.
          if (real && p % 5 == 0) in.seeds.push_back({row, true, CrossFitScorer::fnv1a("v" + std::to_string(p))});
          if (decoy && w == 0 && p % 2 == 0) in.seeds.push_back({row, false, CrossFitScorer::fnv1a("v" + std::to_string(p))});
        }
      }
    }
    return in;
  }
}

int main()
{
  // 1. Gradient of the readout objective against central finite differences.
  {
    WhitenedNetwork net(5, 7, 3, 11);
    std::mt19937_64 rng(3);
    std::normal_distribution<float> nd(0.0f, 1.0f);
    std::vector<float> data(40 * 5);
    for (auto& v : data) v = nd(rng);
    for (std::size_t i = 0; i < 10 * 5; ++i) data[i] += 0.8f;
    std::vector<const float*> pos, neg;
    for (std::size_t r = 0; r < 10; ++r) pos.push_back(data.data() + r * 5);
    for (std::size_t r = 10; r < 40; ++r) neg.push_back(data.data() + r * 5);
    std::vector<double> g;
    net.objective(pos, neg, 1e-2, &g, 1);
    double worst = 0.0;
    for (std::size_t p = 0; p < net.parameters().size(); ++p)
    {
      const double keep = net.parameters()[p], h = 1e-6;
      net.parameters()[p] = keep + h;
      const double up = net.objective(pos, neg, 1e-2, nullptr, 1);
      net.parameters()[p] = keep - h;
      const double down = net.objective(pos, neg, 1e-2, nullptr, 1);
      net.parameters()[p] = keep;
      const double fd = (up - down) / (2 * h);
      worst = std::max(worst, std::abs(fd - g[p]) / std::max(1.0, std::abs(fd)));
    }
    check(worst < 1e-5, "gradient matches finite differences (worst rel err " + std::to_string(worst) + ")");
    std::vector<double> g4;
    net.objective(pos, neg, 1e-2, &g4, 4);
    check(g4 == g, "gradient independent of thread count");
  }

  // 2. Readout in one dimension: w = (mu+ - mu-) / (var- + lambda).
  {
    const std::vector<double> p = {3.0, 5.0}, n = {0.0, 1.0, 2.0};
    double obj = 0.0;
    const auto w = WhitenedNetwork::readout(p, n, 1, 0.5, obj);
    check(std::abs(w[0] - 3.0 / 1.5) < 1e-12 && std::abs(obj - 9.0 / 1.5) < 1e-12, "1-d readout closed form");
  }

  // 3. Pooled paired TDC on a hand example: one fold, decoy winners at 0, 4, 1, 2.
  {
    CrossFitScorer::Input in;
    in.dimension = 1;
    // pair 0: target 5 vs decoy 0; pair 1: target 3 (two windows, 3 and -1) vs decoy 4; pair 2: target 1.5 vs decoy 1;
    // pair 3: decoy 2 only.
    const std::vector<std::tuple<std::uint32_t, std::uint32_t, int, float>> rows = {
      {0, 0, 0, 5.0f}, {1, 0, 1, 0.0f}, {2, 1, 0, 3.0f}, {2, 1, 0, -1.0f}, {3, 1, 1, 4.0f}, {4, 2, 0, 1.5f}, {5, 2, 1, 1.0f}, {6, 3, 1, 2.0f}};
    std::vector<float> score;
    for (const auto& [pre, pair, dec, s] : rows)
    {
      in.precursor.push_back(pre);
      in.pair.push_back(pair);
      in.decoy.push_back(static_cast<std::uint8_t>(dec));
      in.fold.push_back(0);
      in.pair_hash.push_back(0);
      in.features.push_back(0.0f);
      score.push_back(s);
    }
    const auto c = CrossFitScorer::pooledPairedTDC(in, score, {0, 1, 2, 3, 4, 5, 6, 7});
    check(c.winners.size() == 7, "one winner per precursor");
    // Ranked: target 5 (T1), decoy 4 (D1), decoy 2 (D2), target 1.5 (T2).
    check(c.ranked == std::vector<std::uint32_t>({0, 4, 7, 5}), "paired competition order");
    // raw q: 1/1, 2/1 -> 1, 3/1 -> 1, 3/2 -> 1 ; reverse min -> 1, 1, 1, 1; first is min(1, 1) = 1.
    check(std::abs(c.q[0] - 1.0) < 1e-12 && std::abs(c.q[3] - 1.0) < 1e-12, "q values with (D + 1) / T");
  }

  // 4. Tail calibration is monotone within a fold and keeps NOT_KEPT rows.
  {
    std::vector<bool> truth;
    auto in = synthetic(900, 1, 0.0, 5, truth);
    std::vector<float> s(in.rows());
    std::mt19937_64 rng(9);
    std::normal_distribution<float> nd(0.0f, 1.0f);
    for (auto& v : s) v = nd(rng);
    s[0] = CrossFitScorer::NOT_KEPT;
    const auto c = CrossFitScorer::calibrateTail(in, s);
    check(c[0] == CrossFitScorer::NOT_KEPT, "NOT_KEPT stays");
    bool monotone = true;
    for (std::size_t a = 1; a < s.size(); ++a)
      for (std::size_t b = 1; b < s.size(); b += 37)
        if (in.fold[a] == in.fold[b] && s[a] < s[b] && c[a] > c[b] + 1e-6f) monotone = false;
    check(monotone, "tail calibration is monotone within a fold");
  }

  // 5. End to end on a separable synthetic population: finds the true windows, keeps decoy q honest, thread independent.
  {
    std::vector<bool> truth;
    const auto in = synthetic(3000, 3, 3.0, 21, truth);
    CrossFitScorer::Options o;
    o.max_epochs = 6;
    o.steps_per_epoch = 20;
    o.batch_positive = 64;
    o.batch_negative = 256;
    o.hidden = 16;
    o.embed = 4;
    o.seed = 3;
    o.threads = 1;
    const auto a = CrossFitScorer::run(in, o);
    o.threads = 4;
    const auto b = CrossFitScorer::run(in, o);
    check(a.score == b.score && a.q_value == b.q_value, "run independent of thread count");
    std::size_t real = 0, passing = 0, false_hits = 0;
    for (std::size_t r = 0; r < in.rows(); ++r)
    {
      real += truth[r];
      if (!in.decoy[r] && a.q_value[r] <= 0.01)
      {
        ++passing;
        false_hits += !truth[r];
      }
    }
    std::cerr << "synthetic: real " << real << " passing " << passing << " false " << false_hits << '\n';
    check(passing > real * 8 / 10, "finds most true windows at q <= 0.01");
    check(false_hits <= passing / 50 + 2, "false share at q <= 0.01 near the nominal level");
  }

  // 6. Input validation.
  {
    std::vector<bool> truth;
    auto in = synthetic(30, 1, 1.0, 1, truth);
    in.fold[1] = (in.fold[0] + 1) % 3; // decoy of pair 0 in another fold
    bool threw = false;
    try { CrossFitScorer::run(in, {}); } catch (const std::invalid_argument&) { threw = true; }
    check(threw, "rejects a pair spanning two folds");
  }

  // 7. RT width selection: true windows sit within +-10 s of the library RT, every other window (decoys, targets'
  // chance windows) is uniform over +-100 s. Narrowing the window removes chance competition, so the chosen width must
  // be below the full tolerance, cover the true windows, and keep at least the full-window target count.
  {
    std::vector<bool> truth;
    auto in = synthetic(3000, 3, 2.0, 23, truth);
    std::mt19937_64 rng(4);
    std::uniform_real_distribution<double> wide(-100.0, 100.0), near(-10.0, 10.0);
    std::vector<double> rt(in.rows());
    for (std::size_t r = 0; r < in.rows(); ++r) rt[r] = truth[r] ? near(rng) : wide(rng);
    CrossFitScorer::Options o;
    o.max_epochs = 6;
    o.steps_per_epoch = 20;
    o.batch_positive = 64;
    o.batch_negative = 256;
    o.hidden = 16;
    o.embed = 4;
    o.seed = 5;
    auto res = CrossFitScorer::run(in, o);
    std::size_t full = 0;
    for (std::size_t r = 0; r < in.rows(); ++r) full += !in.decoy[r] && res.q_value[r] <= 0.01;
    std::vector<double> grid;
    for (int k = 8; k <= 20; ++k) grid.push_back(0.05 * k);
    auto selected = res;
    const auto width = CrossFitScorer::selectWidth(in, selected, rt, 100.0, grid);
    std::size_t narrow = 0, outside = 0;
    for (std::size_t r = 0; r < in.rows(); ++r)
    {
      narrow += !in.decoy[r] && selected.q_value[r] <= 0.01;
      outside += selected.score[r] > CrossFitScorer::NOT_KEPT && std::abs(rt[r]) > width[in.fold[r]];
    }
    std::cerr << "width selection: widths " << width[0] << " " << width[1] << " " << width[2] << ", targets " << full << " -> " << narrow << '\n';
    check(width.size() == 3, "one width per fold");
    bool below = true;
    for (double w : width) below = below && w < 100.0 && w >= 10.0;
    check(below, "chosen widths below the tolerance and covering the true windows");
    check(outside == 0, "no kept window outside its fold's width");
    check(narrow >= full, "width selection keeps at least the full-window target count");
    auto again = res;
    check(CrossFitScorer::selectWidth(in, again, rt, 100.0, grid) == width && again.q_value == selected.q_value, "width selection is deterministic");
  }

  if (failures) return EXIT_FAILURE;
  std::cerr << "OpenDIAScorer_test passed\n";
  return EXIT_SUCCESS;
}

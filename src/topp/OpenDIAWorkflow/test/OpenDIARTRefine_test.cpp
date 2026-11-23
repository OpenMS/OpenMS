// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include "OpenDIARTRefine.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>

using namespace OpenMS::OpenDIA;

namespace
{
  void require(bool condition, const char* message)
  {
    if (!condition) throw std::runtime_error(message);
  }
  bool near(double a, double b, double tol) { return std::abs(a - b) <= tol; }
}

int main()
{
  try
  {
    // Tokeniser.
    const auto t = RTRefineModel::tokenize("AAC(UniMod:4)M(UniMod:35)K");
    require(t == std::vector<std::string>{"A", "A", "C(UniMod:4)", "M(UniMod:35)", "K"}, "tokens of a modified sequence");
    const auto nt = RTRefineModel::tokenize(".(UniMod:1)PEPTIDEK");
    require(nt.size() == 9 && nt[0] == "(UniMod:1)" && nt[1] == "P" && nt[8] == "K", "N-terminal modification is its own token");

    // Features: vocabulary {(UniMod:1), A, C(UniMod:4), K, M(UniMod:35)}, alphabet {A, C, K, M}.
    const RTRefineModel model({"AAC(UniMod:4)M(UniMod:35)K", "(UniMod:1)AK"});
    require(model.vocabulary().size() == 5 && model.alphabet().size() == 4, "vocabulary and alphabet");
    require(model.vocabulary()[0] == "(UniMod:1)" && model.vocabulary()[1] == "A", "vocabulary sorted");
    require(model.dimension() == 5 + 8 + 2 + 4 + 1, "dimension");
    const auto f = model.features("AAC(UniMod:4)M(UniMod:35)K", 3, 500.5);
    require(f.size() == model.dimension(), "feature length");
    require(f[0] == 0 && f[1] == 2 && f[2] == 1 && f[3] == 1 && f[4] == 1, "token counts");
    require(f[5 + 0] == 1 && f[5 + 4 + 2] == 1, "terminal one-hots (A first, K last)");
    require(f[13] == 5 && near(f[14], std::log(5.0), 1e-12), "residue count and log");
    require(f[15] == 0 && f[16] == 0 && f[17] == 1 && f[18] == 0 && f[19] == 500.5, "charge one-hot and m/z");
    const auto g = model.features("(UniMod:1)AK", 7, 1.0);
    require(g[0] == 1 && g[13] == 2 && g[18] == 1, "N-terminal mod counted, charge clipped to 4");

    // Cholesky: [[4,2],[2,3]] x = [2,1] -> x = [0.5, 0]; [[4,2],[2,3]] x = [6,5] -> x = [1, 1].
    std::vector<double> x;
    require(choleskySolve({4, 2, 2, 3}, {2, 1}, 2, x) && near(x[0], 0.5, 1e-12) && near(x[1], 0.0, 1e-12), "cholesky 1");
    require(choleskySolve({4, 2, 2, 3}, {6, 5}, 2, x) && near(x[0], 1.0, 1e-12) && near(x[1], 1.0, 1e-12), "cholesky 2");
    require(!choleskySolve({1, 2, 2, 1}, {1, 1}, 2, x), "not SPD rejected");

    // Ridge recovers y = 3 count(K) - 2 length on random sequences.
    std::mt19937 rng(7);
    const std::string letters = "ACDEKLPRSTV";
    std::vector<std::string> seqs;
    for (int i = 0; i < 400; ++i)
    {
      std::string s;
      const int len = 6 + static_cast<int>(rng() % 15);
      for (int k = 0; k < len; ++k) s += letters[rng() % letters.size()];
      seqs.push_back(s);
    }
    const RTRefineModel m2(seqs);
    const std::size_t d = m2.dimension();
    std::vector<double> xs, ys;
    std::vector<std::size_t> rows;
    std::vector<std::uint32_t> folds;
    for (std::size_t i = 0; i < seqs.size(); ++i)
    {
      const auto fi = m2.features(seqs[i], 2 + static_cast<int>(i % 2), 400.0 + static_cast<double>(i));
      xs.insert(xs.end(), fi.begin(), fi.end());
      ys.push_back(3.0 * static_cast<double>(std::count(seqs[i].begin(), seqs[i].end(), 'K')) - 2.0 * static_cast<double>(seqs[i].size()));
      rows.push_back(i);
      folds.push_back(static_cast<std::uint32_t>(i % 5));
    }
    const auto fit = fitRidge(xs, d, rows, ys, 1e-8);
    double worst = 0.0;
    for (std::size_t i = 0; i < seqs.size(); ++i) worst = std::max(worst, std::abs(fit.predict(xs.data() + i * d) - ys[i]));
    require(worst < 1e-6, "ridge recovers the linear relation");

    // CV: deterministic, picks the smallest alpha on noise-free data, out-of-fold predictions close.
    const auto a = fitRTRefine(xs, d, rows, ys, folds);
    const auto b = fitRTRefine(xs, d, rows, ys, folds);
    require(a.alpha == b.alpha && a.out_of_fold == b.out_of_fold && a.correction == b.correction, "CV deterministic");
    require(near(a.alpha, 0.01, 1e-12), "smallest alpha chosen on noise-free data");
    require(a.training == seqs.size() && a.correction.size() == seqs.size() && a.out_of_fold.size() == seqs.size(), "result sizes");
    for (std::size_t i = 0; i < seqs.size(); ++i) require(near(a.out_of_fold[i], ys[i], 0.1), "out-of-fold close");
  }
  catch (const std::exception& e)
  {
    std::cerr << "FAILED: " << e.what() << "\n";
    return 1;
  }
  std::cout << "OpenDIARTRefine_test passed\n";
  return 0;
}

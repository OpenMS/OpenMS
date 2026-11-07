// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenDIACNNRTWindow.h>

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace OpenMS::OpenDIA;

namespace
{
  void check(bool condition, const char* what)
  {
    if (!condition)
    {
      throw std::runtime_error(std::string("CNN RT window test failed: ") + what);
    }
  }

  template <class F>
  void throws(F f, const char* what)
  {
    try
    {
      f();
    }
    catch (const std::invalid_argument&)
    {
      return;
    }
    check(false, what);
  }
}

int main()
{
  // A noiseless line y = 100 + 2x is reproduced exactly between knots and held constant outside.
  std::vector<std::pair<double, double>> line;
  for (int i = 0; i < 100; ++i)
  {
    line.emplace_back(i, 100.0 + 2.0 * i);
  }
  const auto exact = fitRTMap(line);
  // Endpoint anchors beyond the outer knots see constant extrapolation, so a few are rejected at 3 s.
  check(exact.library_rt.size() == 20 && exact.sigma == 0.0 && exact.anchors >= 96, "knot count and zero sigma");
  check(std::fabs(exact(50.25) - 200.5) < 1e-12, "interpolation");
  check(exact(-10.0) == exact.run_rt.front() && exact(1000.0) == exact.run_rt.back(), "constant extrapolation");

  // Gross outliers are rejected: the map ignores them and they do not inflate sigma.
  auto noisy = line;
  noisy.emplace_back(40.0, 5000.0);
  noisy.emplace_back(60.0, -5000.0);
  const auto robust = fitRTMap(noisy);
  check(robust.anchors <= exact.anchors && std::fabs(robust(50.25) - 200.5) < 1e-12 && robust.sigma == 0.0, "outlier rejection");

  // Input order does not matter.
  std::vector<std::pair<double, double>> reversed(noisy.rbegin(), noisy.rend());
  const auto again = fitRTMap(reversed);
  check(again.library_rt == robust.library_rt && again.run_rt == robust.run_rt, "order invariance");

  throws([] { fitRTMap(std::vector<std::pair<double, double>>(19, {1.0, 1.0})); }, "too few anchors");
  throws([] { fitRTMap(std::vector<std::pair<double, double>>(30, {NAN, 1.0})); }, "non-finite anchors");

  // Two folds of 30 targets on y = 100 + 2x, plus one far-off window per fold and one decoy.
  std::vector<CNNQValueInput> rows;
  std::vector<double> apex, library;
  for (int fold = 0; fold < 2; ++fold)
  {
    for (int i = 0; i < 30; ++i)
    {
      const std::string id = std::to_string(fold) + "_" + std::to_string(i);
      rows.push_back({"s" + id, "p" + id, "g" + id, false, fold, 5.0 + i});
      library.push_back(i * 3.0);
      apex.push_back(100.0 + 6.0 * i);
    }
    const std::string far = std::to_string(fold) + "_far";
    rows.push_back({"s" + far, "p" + std::to_string(fold) + "_0", "g" + std::to_string(fold) + "_0", false, fold, 1.0});
    library.push_back(0.0);
    apex.push_back(400.0);
    rows.push_back({"d" + far, "d" + far, "gd" + far, true, fold, 0.5});
    library.push_back(30.0);
    apex.push_back(160.0);
  }
  std::vector<RTWindowFold> report;
  const auto mask = crossfitRTWindowMask(rows, apex, library, report, 30);
  check(report.size() == 2 && report[0].map.anchors >= 20 && report[0].half_window == 30.0, "fold reports");
  check(report[0].precursors == 31 && report[0].kept_precursors == 31, "kept precursors");
  for (std::size_t i = 0; i < rows.size(); ++i)
  {
    const bool far_window = rows[i].sample_id.find("_far") != std::string::npos && !rows[i].is_decoy;
    check(mask[i] == !far_window, "window mask");
  }
  // Fold 0's map comes only from fold 1 anchors: corrupting fold 0's own top targets leaves it unchanged.
  auto shifted = apex;
  for (std::size_t i = 0; i < 30; ++i)
  {
    shifted[i] += 1000.0;
  }
  std::vector<RTWindowFold> shifted_report;
  crossfitRTWindowMask(rows, shifted, library, shifted_report, 30);
  check(shifted_report[0].map.run_rt == report[0].map.run_rt, "out-of-fold anchors");
  check(shifted_report[1].map.run_rt != report[1].map.run_rt, "other fold uses these anchors");

  throws([&] { std::vector<RTWindowFold> r; crossfitRTWindowMask(rows, {1.0}, library, r); }, "size mismatch");
  std::vector<CNNQValueInput> one_fold(rows.begin(), rows.begin() + 31);
  throws([&] { std::vector<RTWindowFold> r; crossfitRTWindowMask(one_fold, std::vector<double>(apex.begin(), apex.begin() + 31),
                                                                 std::vector<double>(library.begin(), library.begin() + 31), r); },
         "single fold");
  std::cout << "OpenDIACNNRTWindow tests passed\n";
  return 0;
}

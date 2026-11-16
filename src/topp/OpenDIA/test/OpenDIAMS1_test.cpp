// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenDIAMS1.h>

#include <cmath>
#include <cstdlib>
#include <iostream>

using namespace OpenMS::OpenDIA;

namespace
{
#include "OpenDIAMS1_cases.inc"

  int failures = 0;
  void check(bool ok, const char* what)
  {
    if (!ok)
    {
      std::cerr << "FAILED: " << what << '\n';
      ++failures;
    }
  }
}

int main()
{
  // 1. Isotope envelope against the reference (ms1_precheck.py) for plain, carbamidomethyl and oxidised sequences.
  {
    const auto e = MS1Evidence::envelope("PEPTIDEK");
    const double ref[4] = {0.6034096066741153, 0.28908116292464126, 0.08763753628556781, 0.019871694115675718};
    bool ok = true;
    for (int i = 0; i < 4; ++i) ok = ok && std::abs(e[i] - ref[i]) < 1e-9;
    check(ok, "PEPTIDEK envelope matches the reference");
    const auto m = MS1Evidence::envelope("AAAAEAALKAQLNC(UniMod:4)M(UniMod:35)LR");
    const double refm[4] = {0.36006641091481706, 0.33809235332889226, 0.20704397862731033, 0.09479725712898034};
    ok = true;
    for (int i = 0; i < 4; ++i) ok = ok && std::abs(m[i] - refm[i]) < 1e-9;
    check(ok, "modified sequence envelope matches the reference");
    const auto u = MS1Evidence::envelope("PEPT(UniMod:21)IDEK");
    check(u[0] == 0.0 && u[3] == 0.0, "unknown modification gives a zero envelope");
    const auto x = MS1Evidence::envelope("PEPTXDEK");
    check(x[0] == 0.0, "unknown residue gives a zero envelope");
  }

  // 2. Features of real arm-E raw01 windows against the Python precheck (float32 arrays, so 1e-4 relative).
  for (const auto& c : kCases)
  {
    IsotopeTraces iso{};
    for (int i = 0; i < 4; ++i)
      for (int b = 0; b < 32; ++b) iso[i][b] = static_cast<float>(c.ms1[i][b]);
    std::array<double, 32> xs{};
    for (int b = 0; b < 32; ++b) xs[b] = c.xs[b];
    const auto env = MS1Evidence::envelope(c.seq);
    bool env_ok = true;
    for (int i = 0; i < 4; ++i) env_ok = env_ok && std::abs(env[i] - c.env[i]) < 1e-9;
    check(env_ok, "case envelope matches");
    const auto f = MS1Evidence::features(iso, xs, env);
    for (int k = 0; k < 5; ++k)
    {
      const bool ok = std::abs(f[k] - c.feat[k]) <= 1e-4 * std::max(1.0, std::abs(c.feat[k]));
      if (!ok) std::cerr << c.seq << " feature " << k << ": " << f[k] << " vs " << c.feat[k] << '\n';
      check(ok, "feature matches the precheck");
    }
  }

  // 3. A clean synthetic peak: M+0..M+3 in the envelope ratio, co-eluting with the fragments.
  {
    const auto env = MS1Evidence::envelope("PEPTIDEK");
    IsotopeTraces iso{};
    std::array<double, 32> xs{};
    for (int b = 0; b < 32; ++b)
    {
      const double g = std::exp(-0.5 * (b - 16.0) * (b - 16.0) / 4.0);
      for (int i = 0; i < 4; ++i) iso[i][b] = static_cast<float>(1000.0 * env[i] * g);
      xs[b] = 500.0 * g;
    }
    const auto f = MS1Evidence::features(iso, xs, env);
    check(f[0] > 0.999 && f[1] > 0.999 && f[2] > 0.999 && f[3] == 0.0 && f[4] > 0.0, "clean co-eluting isotope peak scores high");
  }

  if (failures) return EXIT_FAILURE;
  std::cerr << "OpenDIAMS1_test passed\n";
  return EXIT_SUCCESS;
}

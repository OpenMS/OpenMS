// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include "NuXLSearchIndex.h"
///////////////////////////

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>

#include <map>
#include <random>

using namespace OpenMS;
using namespace std;

namespace
{
  // spectrum with the given (sorted) peak m/z and the charge array NuXLSpectrumIndex reads
  MSSpectrum spectrumWith(const vector<double>& mzs)
  {
    MSSpectrum s;
    s.getIntegerDataArrays().resize(1);
    for (Size i = 0; i < mzs.size(); ++i)
    {
      s.push_back(Peak1D(mzs[i], static_cast<float>(i + 1)));
      s.getIntegerDataArrays()[NuXLConstants::IA_CHARGE_INDEX].push_back(1 + static_cast<int>(i % 3));
    }
    return s;
  }

  // queries around every peak, on every bucket edge and outside the m/z range
  vector<double> queries(const vector<double>& mzs, std::mt19937& rng)
  {
    vector<double> q;
    for (double mz : mzs)
    {
      q.push_back(mz);
      q.push_back(std::nextafter(mz, 0.0));
      q.push_back(std::nextafter(mz, 1e9));
    }
    for (Size i = 1; i + 1 < mzs.size(); ++i) { q.push_back(0.5 * (mzs[i - 1] + mzs[i])); } // equidistant: ties
    for (double edge = 0.0; edge < (mzs.empty() ? 10.0 : mzs.back() + 10.0); edge += NuXLSpectrumIndex::BIN)
    {
      q.push_back(edge);
      q.push_back(std::nextafter(edge, 0.0));
    }
    std::uniform_real_distribution<double> any(-5.0, (mzs.empty() ? 10.0 : mzs.back()) + 50.0);
    for (int i = 0; i < 2000; ++i) { q.push_back(any(rng)); }
    return q;
  }
}

START_TEST(NuXLSearchIndex, "$Id$")

START_SECTION([NuXLSpectrumIndex] lowerBound and nearest match std::lower_bound and MSSpectrum::findNearest)
{
  std::mt19937 rng(42);
  vector<vector<double>> peak_lists;
  // random spectra with duplicate m/z values, closely spaced peaks and peaks exactly on bucket edges
  for (int s = 0; s < 20; ++s)
  {
    std::uniform_real_distribution<double> mz(100.0, 2000.0);
    vector<double> mzs;
    const int n = 1 + static_cast<int>(rng() % 400);
    for (int i = 0; i < n; ++i) { mzs.push_back(mz(rng)); }
    mzs.push_back(mzs.front()); // duplicate
    mzs.push_back(400.0);       // on a bucket edge
    mzs.push_back(400.0 + 1e-9);
    std::sort(mzs.begin(), mzs.end());
    peak_lists.push_back(mzs);
  }
  peak_lists.push_back({500.0}); // single peak
  // more peaks than the 16 bit bucket index holds: binary search fallback
  vector<double> many;
  for (int i = 0; i < 70000; ++i) { many.push_back(100.0 + 0.025 * i); }
  peak_lists.push_back(many);

  Size checked = 0, mismatches = 0;
  for (const auto& mzs : peak_lists)
  {
    const MSSpectrum s = spectrumWith(mzs);
    NuXLSpectrumIndex si;
    si.build(s);
    TEST_EQUAL(si.mz.size(), s.size())
    for (double x : queries(mzs, rng))
    {
      const Size expected_lb = std::lower_bound(si.mz.begin(), si.mz.end(), x) - si.mz.begin();
      if (si.lowerBound(x) != expected_lb) ++mismatches;
      if (si.nearest(x) != s.findNearest(x)) ++mismatches;
      ++checked;
    }
  }
  TEST_EQUAL(checked > 50000, true)
  TEST_EQUAL(mismatches, 0)

  // the copied peak data
  const MSSpectrum s = spectrumWith({100.0, 200.0, 300.0});
  NuXLSpectrumIndex si;
  si.build(s);
  TEST_REAL_SIMILAR(si.mz[1], 200.0)
  TEST_REAL_SIMILAR(si.intensity[2], 3.0)
  TEST_EQUAL(si.charge[1], 2)
}
END_SECTION

START_SECTION([NuXLPrecursorIndex] lowerBound and upperBound match the multimap)
{
  std::mt19937 rng(7);
  std::uniform_real_distribution<double> mass(300.0, 6000.0);
  multimap<double, pair<Size, int>> mm;
  for (Size i = 0; i < 20000; ++i)
  {
    const double m = mass(rng);
    mm.insert({m, {i, 0}});
    if (i % 50 == 0) mm.insert({m, {i, 1}}); // same mass twice
  }
  mm.insert({1000.0, {99999, 0}}); // on a bucket edge
  NuXLPrecursorIndex pi;
  pi.build(mm);
  TEST_EQUAL(pi.mass.size(), mm.size())

  // the index order is the multimap order
  Size pos = 0, order_mismatches = 0;
  for (const auto& e : mm)
  {
    if (pi.mass[pos] != e.first || pi.scan_and_isotope[pos] != e.second) ++order_mismatches;
    ++pos;
  }
  TEST_EQUAL(order_mismatches, 0)

  vector<double> q{0.0, 299.0, 1000.0, std::nextafter(1000.0, 0.0), std::nextafter(1000.0, 1e9), 7000.0};
  for (const auto& e : mm) { q.push_back(e.first); q.push_back(std::nextafter(e.first, 0.0)); }
  for (int i = 0; i < 20000; ++i) { q.push_back(mass(rng)); }
  // reference: the multimap's keys in its order (std::distance on the multimap itself would be linear per query)
  vector<double> keys;
  for (const auto& e : mm) { keys.push_back(e.first); }
  Size mismatches = 0;
  for (double x : q)
  {
    if (pi.lowerBound(x) != static_cast<Size>(std::lower_bound(keys.begin(), keys.end(), x) - keys.begin())) ++mismatches;
    if (pi.upperBound(x) != static_cast<Size>(std::upper_bound(keys.begin(), keys.end(), x) - keys.begin())) ++mismatches;
  }
  // spot checks against the multimap itself
  for (double x : {299.0, 1000.0, std::nextafter(1000.0, 0.0), 7000.0})
  {
    if (pi.lowerBound(x) != static_cast<Size>(std::distance(mm.begin(), mm.lower_bound(x)))) ++mismatches;
    if (pi.upperBound(x) != static_cast<Size>(std::distance(mm.begin(), mm.upper_bound(x)))) ++mismatches;
  }
  TEST_EQUAL(mismatches, 0)

  NuXLPrecursorIndex empty;
  empty.build(multimap<double, pair<Size, int>>());
  TEST_EQUAL(empty.lowerBound(500.0), 0)
  TEST_EQUAL(empty.upperBound(500.0), 0)
}
END_SECTION

START_SECTION([NuXLLogFactorialTable] same values as the summation it replaces)
{
  const NuXLLogFactorialTable table;
  Size mismatches = 0;
  for (UInt x = 0; x < 5000; ++x) // beyond the table size, too
  {
    double z = 0;
    for (double y = 2; y <= static_cast<double>(x); ++y) { z += std::log(y); }
    if (table(x) != z) ++mismatches; // bitwise identical
  }
  TEST_EQUAL(mismatches, 0)
}
END_SECTION

START_SECTION([NuXLScoringScratch] matchOdds returns the memoized value of the function)
{
  NuXLScoringScratch scratch;
  Size calls = 0;
  auto value = [](size_t n, size_t m, double p) { return static_cast<double>(n) + static_cast<double>(m) / 100.0 + p; };
  auto f = [&calls, &value](size_t n, size_t m, double p) { ++calls; return value(n, m, p); };
  // (TEST_REAL_SIMILAR evaluates its arguments more than once, so the calls are made outside of it)
  double v = scratch.matchOdds(2000, 3, 0.5, f);
  TEST_REAL_SIMILAR(v, value(2000, 3, 0.5))
  TEST_EQUAL(calls, 1)
  v = scratch.matchOdds(2000, 3, 0.5, f); // cached
  TEST_REAL_SIMILAR(v, value(2000, 3, 0.5))
  TEST_EQUAL(calls, 1)
  v = scratch.matchOdds(2000, 7, 0.5, f);
  TEST_REAL_SIMILAR(v, value(2000, 7, 0.5))
  TEST_EQUAL(calls, 2)
  v = scratch.matchOdds(2000, 1, 0.5, f); // below the largest count seen
  TEST_REAL_SIMILAR(v, value(2000, 1, 0.5))
  TEST_EQUAL(calls, 3)
  v = scratch.matchOdds(2000, 7, 0.5, f); // still cached after the row grew
  TEST_REAL_SIMILAR(v, value(2000, 7, 0.5))
  TEST_EQUAL(calls, 3)
  v = scratch.matchOdds(5, 6, 0.5, f); // more matches than peaks: not cached
  TEST_REAL_SIMILAR(v, value(5, 6, 0.5))
  v = scratch.matchOdds(5, 6, 0.5, f);
  TEST_EQUAL(calls, 5)
  // a row only grows to the largest match count seen, not to N + 1
  TEST_EQUAL(scratch.match_odds[2000].size(), 8)
  TEST_EQUAL(scratch.match_odds[5].size(), 0)
}
END_SECTION

START_SECTION([NuXLFragmentAdductRegistry] ranks follow NuXLFragmentAdductDefinition::operator<)
{
  vector<NuXLFragmentAdductDefinition> defs{
    {EmpiricalFormula("C9H10N2O5"), "U'", 226.0594},
    {EmpiricalFormula("C4H4N2O2"), "U'", 112.0273},
    {EmpiricalFormula("C9H12N2O6"), "U", 244.0695},
    {EmpiricalFormula("C9H10N2O5"), "U-H2O", 226.0594}, // same mass and formula, other name
    {EmpiricalFormula("C4H4N2O2"), "U'", 112.0273},     // equal to the second one
  };
  NuXLFragmentAdductRegistry reg;
  for (const auto& d : defs) { reg.add(&d); }
  map<string, vector<vector<double>>> block{{"U", {{1.0}}}};
  reg.finalize(block);

  Size inconsistent = 0;
  for (Size a = 0; a < defs.size(); ++a)
  {
    for (Size b = 0; b < defs.size(); ++b)
    {
      if ((reg.rank[a] < reg.rank[b]) != (defs[a] < defs[b])) ++inconsistent;
    }
  }
  TEST_EQUAL(inconsistent, 0)
  TEST_EQUAL(reg.rank[1], reg.rank[4])
  // names share ids
  TEST_EQUAL(reg.name_id.size(), 3)
  TEST_EQUAL(reg.id_to_name_id[0], reg.id_to_name_id[1])
  TEST_EQUAL(reg.block_lists_by_name_id[reg.name_id.at("U")] == &block.at("U"), true)
  TEST_EQUAL(reg.block_lists_by_name_id[reg.name_id.at("U'")] == nullptr, true)

  // nucleotide tags of a spectrum by name id
  MSSpectrum s = spectrumWith({100.0});
  s.getStringDataArrays().resize(1);
  s.getStringDataArrays()[0].push_back("U-H2O");
  NuXLSpectrumIndex si;
  si.build(s);
  reg.annotateSpectrum(si);
  TEST_EQUAL(si.has_tag_for_fa_id[reg.name_id.at("U-H2O")], 1)
  TEST_EQUAL(si.has_tag_for_fa_id[reg.name_id.at("U")], 0)
}
END_SECTION

END_TEST

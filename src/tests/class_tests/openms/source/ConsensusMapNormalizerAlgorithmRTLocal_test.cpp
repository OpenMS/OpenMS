// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmRTLocal.h>
///////////////////////////

#include <cmath>
#include <functional>
#include <limits>

using namespace OpenMS;
using namespace std;

namespace
{
  using Settings = ConsensusMapNormalizerAlgorithmRTLocal::Settings;

  /// A consensus map with @p n_maps maps and @p n_features features at RT 10, 20, 30, ...
  /// intensity(feature f, map m) = base(f) * factor(f, m); a factor <= 0 leaves the handle out.
  /// The m-th map gets the map index @p map_ids[m] (default: m).
  ConsensusMap makeMap(Size n_maps, Size n_features, const function<double(Size, Size)>& factor, const vector<UInt64>& map_ids = {})
  {
    auto id = [&map_ids](Size m) { return map_ids.empty() ? UInt64(m) : map_ids[m]; };
    ConsensusMap map;
    for (Size m = 0; m < n_maps; ++m)
    {
      map.getColumnHeaders()[id(m)].filename = "map" + std::to_string(id(m)) + ".mzML";
      map.getColumnHeaders()[id(m)].size = n_features;
    }
    UInt64 uid = 0;
    for (Size f = 0; f < n_features; ++f)
    {
      const double rt = 10.0 * double(f + 1);
      const double base = 1000.0 * double(1 + (f * 37) % 101); // spread over two orders of magnitude
      ConsensusFeature cf;
      cf.setRT(rt);
      cf.setMZ(400.0 + double(f));
      cf.setIntensity(123456.0);
      for (Size m = 0; m < n_maps; ++m)
      {
        const double fac = factor(f, m);
        if (fac <= 0.0) continue;
        Peak2D p;
        p.setRT(rt);
        p.setMZ(400.0 + double(f));
        p.setIntensity(base * fac);
        cf.insert(id(m), p, ++uid);
      }
      map.push_back(cf);
    }
    return map;
  }

  /// Intensity of feature @p f in map @p m (0 if absent).
  double intensity(const ConsensusMap& map, Size f, Size m)
  {
    for (const auto& fh : map[f].getFeatures())
    {
      if (fh.getMapIndex() == m) return fh.getIntensity();
    }
    return 0.0;
  }
}

START_TEST(ConsensusMapNormalizerAlgorithmRTLocal, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((static std::vector<std::vector<double>> computeOffsets(const ConsensusMap& map, const Settings& settings, const std::string& acc_filter, const std::string& desc_filter)))
{
  // constant loading factors 1, 2, 4, 0.5: the offsets are log2(factor) minus their mean (0.5)
  {
    const vector<double> factors = {1.0, 2.0, 4.0, 0.5};
    ConsensusMap map = makeMap(4, 60, [&factors](Size, Size m) { return factors[m]; });
    vector<vector<double>> offsets = ConsensusMapNormalizerAlgorithmRTLocal::computeOffsets(map, Settings(), "", "");
    TEST_EQUAL(offsets.size(), 60)
    TEST_EQUAL(offsets[0].size(), 4)
    TEST_REAL_SIMILAR(offsets[0][0], -0.5)
    TEST_REAL_SIMILAR(offsets[0][1], 0.5)
    TEST_REAL_SIMILAR(offsets[30][2], 1.5)
    TEST_REAL_SIMILAR(offsets[59][3], -1.5)

    // a window larger than the map covers the whole map
    Settings s;
    s.neighbors = std::numeric_limits<Size>::max();
    offsets = ConsensusMapNormalizerAlgorithmRTLocal::computeOffsets(map, s, "", "");
    TEST_REAL_SIMILAR(offsets[0][1], 0.5)
    TEST_REAL_SIMILAR(offsets[30][2], 1.5)
    TEST_REAL_SIMILAR(offsets[59][3], -1.5)
  }

  // a feature in fewer than min_maps maps is not a reference, but takes the offset of its map's neighbours
  {
    ConsensusMap map = makeMap(3, 60, [](Size f, Size m) { return (f == 30 && m != 1) ? 0.0 : (m == 1 ? 4.0 : 1.0); });
    TEST_EQUAL(map[30].getFeatures().size(), 1)
    vector<vector<double>> offsets = ConsensusMapNormalizerAlgorithmRTLocal::computeOffsets(map, Settings(), "", "");
    TEST_REAL_SIMILAR(offsets[29][1], 4.0 / 3.0)
    TEST_REAL_SIMILAR(offsets[30][0], 4.0 / 3.0)  // its only handle, in map 1

    // the same near the start of the map, with a window larger than the map
    map = makeMap(3, 60, [](Size f, Size m) { return (f == 1 && m != 1) ? 0.0 : (m == 1 ? 4.0 : 1.0); });
    Settings s;
    s.neighbors = std::numeric_limits<Size>::max();
    offsets = ConsensusMapNormalizerAlgorithmRTLocal::computeOffsets(map, s, "", "");
    TEST_REAL_SIMILAR(offsets[1][0], 4.0 / 3.0)
  }

  // a map with fewer than min_reference reference features is not normalized; the others are
  {
    ConsensusMap map = makeMap(4, 60, [](Size f, Size m) { return (m == 3 && f >= 10) ? 0.0 : (m == 1 ? 2.0 : 1.0); });
    vector<vector<double>> offsets = ConsensusMapNormalizerAlgorithmRTLocal::computeOffsets(map, Settings(), "", "");
    TEST_REAL_SIMILAR(offsets[0][3], 0.0) // map 3 holds 10 reference features, fewer than 20
    TEST_EQUAL(offsets[0][1] > 0.4, true)
    TEST_EQUAL(offsets[20][1] > 0.4, true)
  }

  // fewer than two maps qualify: nothing is normalized
  {
    ConsensusMap map = makeMap(3, 60, [](Size, Size m) { return m == 1 ? 2.0 : 1.0; });
    Settings s;
    s.min_reference = 100;
    vector<vector<double>> offsets = ConsensusMapNormalizerAlgorithmRTLocal::computeOffsets(map, s, "", "");
    for (const auto& per_feature : offsets)
    {
      for (double o : per_feature) TEST_REAL_SIMILAR(o, 0.0)
    }
  }

  // map indices with gaps, as FileFilter leaves them when it extracts maps 1, 3 and 4 of five: same offsets as maps 0, 1, 2
  {
    const vector<double> factors = {1.0, 4.0, 0.5};
    ConsensusMap dense = makeMap(3, 60, [&factors](Size, Size m) { return factors[m]; });
    ConsensusMap filtered = makeMap(3, 60, [&factors](Size, Size m) { return factors[m]; }, {1, 3, 4});
    TEST_EQUAL(filtered.isMapConsistent(), true)
    vector<vector<double>> offsets_dense = ConsensusMapNormalizerAlgorithmRTLocal::computeOffsets(dense, Settings(), "", "");
    vector<vector<double>> offsets = ConsensusMapNormalizerAlgorithmRTLocal::computeOffsets(filtered, Settings(), "", "");
    TEST_REAL_SIMILAR(offsets[0][0], -1.0 / 3.0)
    TEST_REAL_SIMILAR(offsets[0][1], 5.0 / 3.0)
    TEST_REAL_SIMILAR(offsets[0][2], -4.0 / 3.0)
    TEST_EQUAL(offsets == offsets_dense, true)

    // a handle whose map has no column header
    Peak2D p;
    p.setRT(10.0);
    p.setMZ(400.0);
    p.setIntensity(5000.0);
    filtered[0].insert(2, p, 999999);
    TEST_EXCEPTION(Exception::ElementNotFound, ConsensusMapNormalizerAlgorithmRTLocal::computeOffsets(filtered, Settings(), "", ""))
  }

  // invalid settings
  {
    ConsensusMap map = makeMap(3, 10, [](Size, Size) { return 1.0; });
    Settings s;
    s.neighbors = 0;
    TEST_EXCEPTION(Exception::InvalidParameter, ConsensusMapNormalizerAlgorithmRTLocal::computeOffsets(map, s, "", ""))
    s = Settings();
    s.rounds = 0;
    TEST_EXCEPTION(Exception::InvalidParameter, ConsensusMapNormalizerAlgorithmRTLocal::computeOffsets(map, s, "", ""))
    s = Settings();
    s.min_maps = 0;
    TEST_EXCEPTION(Exception::InvalidParameter, ConsensusMapNormalizerAlgorithmRTLocal::computeOffsets(map, s, "", ""))
    s = Settings();
    s.min_reference = 0;
    TEST_EXCEPTION(Exception::InvalidParameter, ConsensusMapNormalizerAlgorithmRTLocal::computeOffsets(map, s, "", ""))
  }
}
END_SECTION

START_SECTION((static void normalizeMaps(ConsensusMap& map, const Settings& settings, const std::string& acc_filter, const std::string& desc_filter)))
{
  // constant loading factors are removed completely
  {
    const vector<double> factors = {1.0, 2.0, 4.0, 0.5};
    ConsensusMap map = makeMap(4, 60, [&factors](Size, Size m) { return factors[m]; });
    ConsensusMapNormalizerAlgorithmRTLocal::normalizeMaps(map, Settings(), "", "");
    for (Size f = 0; f < map.size(); f += 7)
    {
      TEST_REAL_SIMILAR(intensity(map, f, 1), intensity(map, f, 0))
      TEST_REAL_SIMILAR(intensity(map, f, 2), intensity(map, f, 0))
      TEST_REAL_SIMILAR(intensity(map, f, 3), intensity(map, f, 0))
    }
    // the consensus intensity is left alone
    TEST_REAL_SIMILAR(map[0].getIntensity(), 123456.0)
  }

  // a gain drifting along RT in one map is followed: map 1 rises linearly from 0 to 2 log2 units
  {
    const Size n = 200;
    ConsensusMap map = makeMap(3, n, [n](Size f, Size m) { return m == 1 ? std::exp2(2.0 * double(f) / double(n - 1)) : 1.0; });
    Settings s;
    s.neighbors = 5;
    ConsensusMapNormalizerAlgorithmRTLocal::normalizeMaps(map, s, "", "");
    double max_error = 0.0;
    for (Size f = 2 * s.neighbors; f + 2 * s.neighbors < n; ++f) // away from the ends, where the window is one-sided
    {
      max_error = std::max(max_error, std::fabs(std::log2(intensity(map, f, 1) / intensity(map, f, 0))));
      max_error = std::max(max_error, std::fabs(std::log2(intensity(map, f, 2) / intensity(map, f, 0))));
    }
    TEST_EQUAL(max_error < 1e-6, true)
  }

  // a minority of changed features does not move the maps: every 5th feature is 4x higher in map 2
  {
    ConsensusMap map = makeMap(4, 100, [](Size f, Size m) { return (m == 2 && f % 5 == 0) ? 4.0 : 1.0; });
    const ConsensusMap before = map;
    ConsensusMapNormalizerAlgorithmRTLocal::normalizeMaps(map, Settings(), "", "");
    for (Size f = 0; f < map.size(); ++f)
    {
      for (Size m = 0; m < 4; ++m)
      {
        TEST_REAL_SIMILAR(intensity(map, f, m), intensity(before, f, m))
      }
    }
  }

  // zero intensities stay zero and are not used
  {
    ConsensusMap map = makeMap(3, 60, [](Size, Size m) { return m == 1 ? 2.0 : 1.0; });
    map.getColumnHeaders()[3].filename = "map3.mzML";
    Peak2D zero;
    zero.setRT(60.0);
    zero.setMZ(405.0);
    zero.setIntensity(0.0);
    map[5].insert(FeatureHandle(3, zero, 999999));
    vector<vector<double>> offsets = ConsensusMapNormalizerAlgorithmRTLocal::computeOffsets(map, Settings(), "", "");
    TEST_REAL_SIMILAR(offsets[5][3], 0.0)
    ConsensusMapNormalizerAlgorithmRTLocal::normalizeMaps(map, Settings(), "", "");
    TEST_REAL_SIMILAR(intensity(map, 5, 3), 0.0)
    TEST_REAL_SIMILAR(intensity(map, 5, 1), intensity(map, 5, 0))
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

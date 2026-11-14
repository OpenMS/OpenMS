// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse$
// $Authors: Hendrik Brauer, Oliver Kohlbacher, Johannes Junker$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmThreshold.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

namespace
{
  /// Two maps with the map indices @p ids and three features; map ids[1] is twice as intense as map ids[0].
  /// The column header sizes are @p sizes (the map with the larger size is the reference).
  ConsensusMap makeMap(const vector<UInt64>& ids, const vector<Size>& sizes)
  {
    ConsensusMap map;
    for (Size i = 0; i < ids.size(); ++i)
    {
      map.getColumnHeaders()[ids[i]].filename = "map" + std::to_string(ids[i]) + ".mzML";
      map.getColumnHeaders()[ids[i]].size = sizes[i];
    }
    UInt64 uid = 0;
    for (double intensity : {100.0, 200.0, 300.0})
    {
      ConsensusFeature cf;
      cf.setRT(intensity);
      cf.setMZ(500.0);
      Peak2D p;
      p.setRT(intensity);
      p.setMZ(500.0);
      p.setIntensity(intensity);
      cf.insert(ids[0], p, ++uid);
      p.setIntensity(2.0 * intensity);
      cf.insert(ids[1], p, ++uid);
      map.push_back(cf);
    }
    return map;
  }

  /// All handle intensities of one map, in consensus feature order.
  vector<double> intensitiesOf(const ConsensusMap& map, UInt64 map_index)
  {
    vector<double> out;
    for (const auto& cf : map)
    {
      for (const auto& fh : cf.getFeatures())
      {
        if (fh.getMapIndex() == map_index) { out.push_back(fh.getIntensity()); }
      }
    }
    return out;
  }
}

START_TEST(ConsensusMapNormalizerAlgorithmThreshold, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ConsensusMapNormalizerAlgorithmThreshold* ptr = nullptr;
ConsensusMapNormalizerAlgorithmThreshold* nullPointer = nullptr;
START_SECTION(ConsensusMapNormalizerAlgorithmThreshold())
{
  ptr = new ConsensusMapNormalizerAlgorithmThreshold();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~ConsensusMapNormalizerAlgorithmThreshold())
{
	delete ptr;
}
END_SECTION

START_SECTION((static std::vector<double> computeCorrelation(const ConsensusMap & map, const double & ratio_threshold, const std::string & acc_filter, const std::string & desc_filter)))
{
  // one factor per map, in the order of the column headers; the reference is the map with the larger size.
  // Map indices need not be consecutive (FileFilter -consensus:map keeps them).
  for (const vector<UInt64>& ids : {vector<UInt64>{0, 1}, vector<UInt64>{1, 3}})
  {
    vector<double> ratios = ConsensusMapNormalizerAlgorithmThreshold::computeCorrelation(makeMap(ids, {3, 3}), 0.4, "", "");
    TEST_EQUAL(ratios.size(), 2)
    TEST_REAL_SIMILAR(ratios[0], 1.0) // ties: the first map is the reference
    TEST_REAL_SIMILAR(ratios[1], 0.5)
    ratios = ConsensusMapNormalizerAlgorithmThreshold::computeCorrelation(makeMap(ids, {3, 9}), 0.4, "", "");
    TEST_REAL_SIMILAR(ratios[0], 2.0)
    TEST_REAL_SIMILAR(ratios[1], 1.0)
  }

  // a feature handle whose map has no column header
  ConsensusMap map = makeMap({1, 3}, {3, 3});
  map.getColumnHeaders().erase(3);
  TEST_EXCEPTION(Exception::ElementNotFound, ConsensusMapNormalizerAlgorithmThreshold::computeCorrelation(map, 0.4, "", ""))
}
END_SECTION

START_SECTION((static void normalizeMaps(ConsensusMap & map, const std::vector<double> & ratios)))
{
  for (const vector<UInt64>& ids : {vector<UInt64>{0, 1}, vector<UInt64>{1, 3}})
  {
    ConsensusMap map = makeMap(ids, {3, 3});
    ConsensusMapNormalizerAlgorithmThreshold::normalizeMaps(map, {1.0, 0.5});
    TEST_REAL_SIMILAR(intensitiesOf(map, ids[0])[2], 300.0)
    TEST_REAL_SIMILAR(intensitiesOf(map, ids[1])[2], 300.0)
  }
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


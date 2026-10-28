// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmQuantile.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

namespace
{
  /// Two maps with the map indices @p ids, three features: map ids[0] holds 1, 2, 3 and map ids[1] holds 6, 2, 4
  /// (in feature order). The quantile reference distribution is 1.5, 3, 4.5.
  ConsensusMap makeMap(const vector<UInt64>& ids)
  {
    ConsensusMap map;
    for (UInt64 id : ids)
    {
      map.getColumnHeaders()[id].filename = "map" + std::to_string(id) + ".mzML";
      map.getColumnHeaders()[id].size = 3;
    }
    const vector<pair<double, double>> values = {{1.0, 6.0}, {2.0, 2.0}, {3.0, 4.0}};
    UInt64 uid = 0;
    for (Size f = 0; f < values.size(); ++f)
    {
      ConsensusFeature cf;
      cf.setRT(100.0 * double(f + 1));
      cf.setMZ(500.0);
      Peak2D p;
      p.setRT(cf.getRT());
      p.setMZ(500.0);
      p.setIntensity(values[f].first);
      cf.insert(ids[0], p, ++uid);
      p.setIntensity(values[f].second);
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

START_TEST(ConsensusMapNormalizerAlgorithmQuantile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ConsensusMapNormalizerAlgorithmQuantile* ptr = nullptr;
ConsensusMapNormalizerAlgorithmQuantile* null_ptr = nullptr;
START_SECTION(ConsensusMapNormalizerAlgorithmQuantile())
{
	ptr = new ConsensusMapNormalizerAlgorithmQuantile();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~ConsensusMapNormalizerAlgorithmQuantile())
{
	delete ptr;
}
END_SECTION

START_SECTION((virtual ~ConsensusMapNormalizerAlgorithmQuantile()))
{
  // TODO
}
END_SECTION

START_SECTION((static void normalizeMaps(ConsensusMap &map)))
{
  // both maps get the reference distribution 1.5, 3, 4.5, in the rank order of their own values;
  // map indices need not be consecutive (FileFilter -consensus:map keeps them)
  for (const vector<UInt64>& ids : {vector<UInt64>{0, 1}, vector<UInt64>{1, 3}})
  {
    ConsensusMap map = makeMap(ids);
    ConsensusMapNormalizerAlgorithmQuantile::normalizeMaps(map);
    const vector<double> a = intensitiesOf(map, ids[0]);
    const vector<double> b = intensitiesOf(map, ids[1]);
    TEST_REAL_SIMILAR(a[0], 1.5)
    TEST_REAL_SIMILAR(a[1], 3.0)
    TEST_REAL_SIMILAR(a[2], 4.5)
    TEST_REAL_SIMILAR(b[0], 4.5)
    TEST_REAL_SIMILAR(b[1], 1.5)
    TEST_REAL_SIMILAR(b[2], 3.0)
  }
}
END_SECTION

START_SECTION((static void resample(const std::vector< double > &data_in, std::vector< double > &data_out, UInt n_resampling_points)))
{
  // TODO
}
END_SECTION

START_SECTION((static void extractIntensityVectors(const ConsensusMap &map, std::vector< std::vector< double > > &out_intensities)))
{
  // one row per column header, in map index order (maps 1 and 3 -> rows 0 and 1)
  vector<vector<double>> ints;
  ConsensusMapNormalizerAlgorithmQuantile::extractIntensityVectors(makeMap({1, 3}), ints);
  TEST_EQUAL(ints.size(), 2)
  TEST_EQUAL(ints[0].size(), 3)
  TEST_REAL_SIMILAR(ints[0][2], 3.0)
  TEST_REAL_SIMILAR(ints[1][0], 6.0)

  // a feature handle whose map has no column header
  ConsensusMap map = makeMap({1, 3});
  map.getColumnHeaders().erase(3);
  TEST_EXCEPTION(Exception::ElementNotFound, ConsensusMapNormalizerAlgorithmQuantile::extractIntensityVectors(map, ints))
}
END_SECTION

START_SECTION((static void setNormalizedIntensityValues(const std::vector< std::vector< double > > &feature_ints, ConsensusMap &map)))
{
  ConsensusMap map = makeMap({1, 3});
  ConsensusMapNormalizerAlgorithmQuantile::setNormalizedIntensityValues({{10.0, 20.0, 30.0}, {40.0, 50.0, 60.0}}, map);
  TEST_REAL_SIMILAR(intensitiesOf(map, 1)[1], 20.0)
  TEST_REAL_SIMILAR(intensitiesOf(map, 3)[2], 60.0)
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST




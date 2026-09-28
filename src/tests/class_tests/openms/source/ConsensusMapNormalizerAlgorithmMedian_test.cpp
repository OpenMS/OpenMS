// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmMedian.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

namespace
{
  /// One consensus feature holding one handle per (map index, intensity) pair.
  ConsensusFeature makeCF(const vector<pair<UInt64, double>>& handles, double rt, double mz, UInt64& uid)
  {
    ConsensusFeature cf;
    cf.setRT(rt);
    cf.setMZ(mz);
    // Deliberately not the sum of the handles: the sections below assert that normalization leaves
    // this value alone, so it has to be recognisable.
    cf.setIntensity(123456.0);
    for (const auto& h : handles)
    {
      Peak2D p;
      p.setRT(rt);
      p.setMZ(mz);
      p.setIntensity(h.second);
      cf.insert(h.first, p, ++uid);
    }
    return cf;
  }

  /// Three maps, three features. Per-map medians are 200 / 400 / 800.
  ConsensusMap makeMap(const vector<Size>& header_sizes)
  {
    ConsensusMap map;
    for (Size i = 0; i < header_sizes.size(); ++i)
    {
      map.getColumnHeaders()[i].filename = "map" + std::to_string(i) + ".mzML";
      map.getColumnHeaders()[i].size = header_sizes[i];
    }
    UInt64 uid = 0;
    map.push_back(makeCF({{0, 100.0}, {1, 200.0}, {2, 400.0}}, 100.0, 500.0, uid));
    map.push_back(makeCF({{0, 300.0}, {1, 600.0}, {2, 1200.0}}, 200.0, 600.0, uid));
    map.push_back(makeCF({{0, 200.0}, {1, 400.0}, {2, 800.0}}, 300.0, 700.0, uid));
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

START_TEST(ConsensusMapNormalizerAlgorithmMedian, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ConsensusMapNormalizerAlgorithmMedian* ptr = nullptr;
ConsensusMapNormalizerAlgorithmMedian* null_ptr = nullptr;
START_SECTION(ConsensusMapNormalizerAlgorithmMedian())
{
	ptr = new ConsensusMapNormalizerAlgorithmMedian();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~ConsensusMapNormalizerAlgorithmMedian())
{
	delete ptr;
	NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual ~ConsensusMapNormalizerAlgorithmMedian()))
{
  ConsensusMapNormalizerAlgorithmMedian* p = new ConsensusMapNormalizerAlgorithmMedian();
  delete p;
  NOT_TESTABLE
}
END_SECTION

START_SECTION((static Size computeMedians(const ConsensusMap & map, std::vector<double> & medians, const std::string& acc_filter, const std::string& desc_filter)))
{
  // The medians themselves: one per map, over that map's own detected intensities.
  {
    ConsensusMap map = makeMap({3, 9, 3});
    vector<double> medians;
    ConsensusMapNormalizerAlgorithmMedian::computeMedians(map, medians, "", "");
    TEST_EQUAL(medians.size(), 3)
    TEST_REAL_SIMILAR(medians[0], 200.0)
    TEST_REAL_SIMILAR(medians[1], 400.0)
    TEST_REAL_SIMILAR(medians[2], 800.0)
  }

  // The returned reference index is the argmax of ColumnHeader::size -- NOT of the number of
  // handles actually present. All three maps hold three handles here; only the stored sizes differ.
  {
    ConsensusMap map = makeMap({3, 9, 3});
    vector<double> medians;
    TEST_EQUAL(ConsensusMapNormalizerAlgorithmMedian::computeMedians(map, medians, "", ""), 1)
  }

  // Ties keep the lowest index: the comparison is '>', so an equal size does not displace.
  {
    ConsensusMap map = makeMap({9, 9, 9});
    vector<double> medians;
    TEST_EQUAL(ConsensusMapNormalizerAlgorithmMedian::computeMedians(map, medians, "", ""), 0)
  }

  // Degenerate reference (OpenMS issue #9864, finding A3). ColumnHeader::size is left at its 0
  // initializer by some linking paths -- FeatureGroupingAlgorithmQT never touches column headers --
  // and then '0 > 0' never fires and the reference silently collapses to map index 0, whatever the
  // data says. Here map 2 holds the most handles by a clear margin and is still not selected.
  // ProteomicsLFQ used to reach this through every fractionated run; it no longer normalizes, so
  // the only remaining caller is the ConsensusMapNormalizer TOPP tool, whose consensusXML input
  // restores the sizes on read. This section pins the defect so it stays visible while it is latent.
  {
    ConsensusMap map = makeMap({0, 0, 0});
    UInt64 uid = 1000;
    map.push_back(makeCF({{2, 1600.0}}, 400.0, 800.0, uid));
    map.push_back(makeCF({{2, 3200.0}}, 500.0, 900.0, uid));
    TEST_EQUAL(intensitiesOf(map, 0).size(), 3)
    TEST_EQUAL(intensitiesOf(map, 2).size(), 5)
    vector<double> medians;
    TEST_EQUAL(ConsensusMapNormalizerAlgorithmMedian::computeMedians(map, medians, "", ""), 0)
  }

  // Silent bail-out: if ANY map has no feature passing the filters, every median is set to 1.0 and
  // the function returns 0, which makes the whole normalization a multiply-by-one. Only a warning
  // is logged, so a caller cannot tell this apart from "normalized successfully" by return value.
  {
    ConsensusMap map = makeMap({3, 3, 3});
    map.getColumnHeaders()[3].filename = "map3.mzML"; // declared, but no handle ever references it
    map.getColumnHeaders()[3].size = 0;
    vector<double> medians;
    TEST_EQUAL(ConsensusMapNormalizerAlgorithmMedian::computeMedians(map, medians, "", ""), 0)
    TEST_EQUAL(medians.size(), 4)
    for (Size i = 0; i < medians.size(); ++i) { TEST_REAL_SIMILAR(medians[i], 1.0) }
  }

  // A column header missing from the map is an error, not a silently skipped map.
  {
    ConsensusMap map = makeMap({3, 3, 3});
    map.getColumnHeaders().erase(1);
    map.getColumnHeaders()[7].filename = "map7.mzML"; // keeps the count at 3, but index 1 is gone
    vector<double> medians;
    TEST_EXCEPTION(Exception::ElementNotFound,
                   ConsensusMapNormalizerAlgorithmMedian::computeMedians(map, medians, "", ""))
  }
}
END_SECTION

START_SECTION((static void normalizeMaps(ConsensusMap & map, NormalizationMethod method, const std::string& acc_filter, const std::string& desc_filter)))
{
  // NM_SCALE multiplies every handle of map j by medians[reference] / medians[j].
  // Reference is map 1 (largest stored size), median 400: map 0 doubles, map 1 is unchanged,
  // map 2 halves.
  {
    ConsensusMap map = makeMap({3, 9, 3});
    ConsensusMapNormalizerAlgorithmMedian::normalizeMaps(map, ConsensusMapNormalizerAlgorithmMedian::NM_SCALE, "", "");
    const vector<double> m0 = intensitiesOf(map, 0);
    const vector<double> m1 = intensitiesOf(map, 1);
    const vector<double> m2 = intensitiesOf(map, 2);
    TEST_REAL_SIMILAR(m0[0], 200.0)
    TEST_REAL_SIMILAR(m0[1], 600.0)
    TEST_REAL_SIMILAR(m0[2], 400.0)
    TEST_REAL_SIMILAR(m1[0], 200.0)
    TEST_REAL_SIMILAR(m1[1], 600.0)
    TEST_REAL_SIMILAR(m1[2], 400.0)
    TEST_REAL_SIMILAR(m2[0], 200.0)
    TEST_REAL_SIMILAR(m2[1], 600.0)
    TEST_REAL_SIMILAR(m2[2], 400.0)
  }

  // Only the handles are rewritten. The consensus feature's own intensity is left untouched, so
  // after normalization it no longer agrees with the handles it aggregates.
  {
    ConsensusMap map = makeMap({3, 9, 3});
    ConsensusMapNormalizerAlgorithmMedian::normalizeMaps(map, ConsensusMapNormalizerAlgorithmMedian::NM_SCALE, "", "");
    for (const auto& cf : map) { TEST_REAL_SIMILAR(cf.getIntensity(), 123456.0) }
  }

  // The scale factor shares a common numerator, so within-map ratios are invariant to which map is
  // chosen as the reference: only the absolute scale of each map moves. Same fixture, reference
  // forced to map 0 by the stored sizes.
  {
    ConsensusMap map = makeMap({9, 3, 3});
    ConsensusMapNormalizerAlgorithmMedian::normalizeMaps(map, ConsensusMapNormalizerAlgorithmMedian::NM_SCALE, "", "");
    const vector<double> m2 = intensitiesOf(map, 2);
    TEST_REAL_SIMILAR(m2[0], 100.0)
    TEST_REAL_SIMILAR(m2[1], 300.0)
    TEST_REAL_SIMILAR(m2[2], 200.0)
    TEST_REAL_SIMILAR(m2[1] / m2[0], 3.0) // unchanged by the reference choice
  }

  // NM_SHIFT ignores the size-derived reference entirely: it scans the medians for the largest and
  // shifts towards that, which is why the degenerate reference above cannot affect it. Map 2 has
  // the largest median (800), so map 0 shifts by +600 and map 2 not at all -- even though the
  // stored sizes name map 1 as the reference.
  {
    ConsensusMap map = makeMap({3, 9, 3});
    ConsensusMapNormalizerAlgorithmMedian::normalizeMaps(map, ConsensusMapNormalizerAlgorithmMedian::NM_SHIFT, "", "");
    const vector<double> m0 = intensitiesOf(map, 0);
    const vector<double> m2 = intensitiesOf(map, 2);
    TEST_REAL_SIMILAR(m0[0], 700.0)
    TEST_REAL_SIMILAR(m0[1], 900.0)
    TEST_REAL_SIMILAR(m2[0], 400.0)
    TEST_REAL_SIMILAR(m2[1], 1200.0)
  }

  // The bail-out is a no-op, not an error: unchanged intensities, no exception.
  {
    ConsensusMap map = makeMap({3, 3, 3});
    map.getColumnHeaders()[3].filename = "map3.mzML";
    ConsensusMapNormalizerAlgorithmMedian::normalizeMaps(map, ConsensusMapNormalizerAlgorithmMedian::NM_SCALE, "", "");
    const vector<double> m2 = intensitiesOf(map, 2);
    TEST_REAL_SIMILAR(m2[0], 400.0)
    TEST_REAL_SIMILAR(m2[1], 1200.0)
    TEST_REAL_SIMILAR(m2[2], 800.0)
  }
}
END_SECTION

START_SECTION((static bool passesFilters_(ConsensusMap::ConstIterator cf_it, const ConsensusMap& map, const std::string& acc_filter, const std::string& desc_filter)))
{
  ConsensusMap map = makeMap({3, 3, 3});

  // Empty filters pass everything, including features carrying no identification at all.
  TEST_TRUE(ConsensusMapNormalizerAlgorithmMedian::passesFilters_(map.begin(), map, "", ""))

  // A non-empty accession filter needs an identification to match against; an unidentified feature
  // is dropped rather than passed through.
  TEST_FALSE(ConsensusMapNormalizerAlgorithmMedian::passesFilters_(map.begin(), map, "ALBU", ""))

  PeptideHit hit;
  hit.setSequence(AASequence::fromString("PEPTIDER"));
  PeptideEvidence ev;
  ev.setProteinAccession("P02769|ALBU_BOVIN");
  hit.addPeptideEvidence(ev);
  PeptideIdentification pid;
  pid.setHits({hit});
  map[0].setPeptideIdentifications({pid});

  TEST_TRUE(ConsensusMapNormalizerAlgorithmMedian::passesFilters_(map.begin(), map, "ALBU", ""))
  TEST_FALSE(ConsensusMapNormalizerAlgorithmMedian::passesFilters_(map.begin(), map, "TRYP", ""))
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Michal Startek $
// $Authors: Michal Startek $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/WNetMatcher.h>

using namespace OpenMS;
using namespace std;

START_TEST(WNetMatcher, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((static DistanceMetric metricFromString(const std::string& s)))
  TEST_EQUAL(WNetMatcher::metricFromString("L1") == WNetMatcher::DistanceMetric::L1, true)
  TEST_EQUAL(WNetMatcher::metricFromString("L2") == WNetMatcher::DistanceMetric::L2, true)
  TEST_EQUAL(WNetMatcher::metricFromString("LINF") == WNetMatcher::DistanceMetric::LINF, true)
  TEST_EQUAL(WNetMatcher::metricFromString("unknown") == WNetMatcher::DistanceMetric::LINF, true)
END_SECTION

START_SECTION((static MatchResult match(const std::vector<std::array<double, 2>>& positions_a, const std::vector<double>& intensities_a, const std::vector<std::array<double, 2>>& positions_b, const std::vector<double>& intensities_b, DistanceMetric metric, double max_distance, double trash_cost)))
{
  // Empty inputs should return empty result
  {
    vector<array<double, 2>> empty_pos;
    vector<double> empty_int;
    vector<array<double, 2>> pos_b = {{1.0, 2.0}};
    vector<double> int_b = {1.0};

    auto result = WNetMatcher::match(empty_pos, empty_int, pos_b, int_b);
    TEST_EQUAL(result.matched_pairs.size(), 0)
    TEST_REAL_SIMILAR(result.cost, 0.0)

    result = WNetMatcher::match(pos_b, int_b, empty_pos, empty_int);
    TEST_EQUAL(result.matched_pairs.size(), 0)
    TEST_REAL_SIMILAR(result.cost, 0.0)
  }

  // Two identical point sets should match perfectly
  {
    vector<array<double, 2>> positions = {{1.0, 10.0}, {2.0, 20.0}, {3.0, 30.0}};
    vector<double> intensities = {100.0, 200.0, 150.0};

    auto result = WNetMatcher::match(positions, intensities, positions, intensities);
    TEST_EQUAL(result.matched_pairs.size(), 3)
  }

  // Close points should match, distant points should not
  {
    vector<array<double, 2>> pos_a = {{1.0, 10.0}, {2.0, 20.0}, {3.0, 30.0}};
    vector<double> int_a = {100.0, 200.0, 150.0};
    vector<array<double, 2>> pos_b = {{1.1, 10.1}, {2.1, 20.1}, {50.0, 500.0}};
    vector<double> int_b = {110.0, 190.0, 50.0};

    auto result = WNetMatcher::match(pos_a, int_a, pos_b, int_b,
                                     WNetMatcher::DistanceMetric::LINF,
                                     5.0, 5.0);
    // First two pairs are close (distance ~0.1), third is far (distance ~470)
    TEST_EQUAL(result.matched_pairs.size(), 2)
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

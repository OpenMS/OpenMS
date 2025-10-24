// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing$
// $Authors: Justin Sing$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/MATH/STATISTICS/RankData.h>
#include <vector>
#include <cmath>
#include <limits>
#include <iostream>

using namespace OpenMS;
using namespace OpenMS::Math;
using namespace std;

///////////////////////////

START_TEST(RankData, "$Id$")

/////////////////////////////////////////////////////////////

// Helper: check all values in a vector are NaN
auto all_nan = [](const std::vector<double>& v) -> bool
{
  for (double x : v) { if (!std::isnan(x)) return false; }
  return true;
};

// Helper: compare two vectors element-wise with TEST_REAL_SIMILAR
auto expect_equal_real_vectors = [](const std::vector<double>& a, const std::vector<double>& b)
{
  TEST_EQUAL(a.size(), b.size())
  for (Size i = 0; i < a.size(); ++i)
  {
    TEST_REAL_SIMILAR(a[i], b[i])
  }
};

START_SECTION((RankData::rankdata - basic tie methods on doubles))
{
  using M = RankData::Method;

  const std::vector<double> a{0.0, 2.0, 3.0, 2.0};

  auto r_avg   = RankData::rankdata(a, M::Average, RankData::NaNPolicy::Propagate);
  auto r_min   = RankData::rankdata(a, M::Min,     RankData::NaNPolicy::Propagate);
  auto r_max   = RankData::rankdata(a, M::Max,     RankData::NaNPolicy::Propagate);
  auto r_dense = RankData::rankdata(a, M::Dense,   RankData::NaNPolicy::Propagate);
  auto r_ord   = RankData::rankdata(a, M::Ordinal, RankData::NaNPolicy::Propagate);

  expect_equal_real_vectors(r_avg,   {1.0, 2.5, 4.0, 2.5});
  expect_equal_real_vectors(r_min,   {1.0, 2.0, 4.0, 2.0});
  expect_equal_real_vectors(r_max,   {1.0, 3.0, 4.0, 3.0});
  expect_equal_real_vectors(r_dense, {1.0, 2.0, 3.0, 2.0});
  expect_equal_real_vectors(r_ord,   {1.0, 2.0, 4.0, 3.0});
}
END_SECTION

START_SECTION((RankData::rankdata - all-equal values and ordinal stability))
{
  using M = RankData::Method;
  const std::vector<double> a{5.0, 5.0, 5.0};

  auto r_avg   = RankData::rankdata(a, M::Average);
  auto r_min   = RankData::rankdata(a, M::Min);
  auto r_max   = RankData::rankdata(a, M::Max);
  auto r_dense = RankData::rankdata(a, M::Dense);
  auto r_ord   = RankData::rankdata(a, M::Ordinal);

  expect_equal_real_vectors(r_avg,   {2.0, 2.0, 2.0});
  expect_equal_real_vectors(r_min,   {1.0, 1.0, 1.0});
  expect_equal_real_vectors(r_max,   {3.0, 3.0, 3.0});
  expect_equal_real_vectors(r_dense, {1.0, 1.0, 1.0});
  // ordinal should break ties by original order (stable)
  expect_equal_real_vectors(r_ord,   {1.0, 2.0, 3.0});
}
END_SECTION

START_SECTION((RankData::rankdata - NaN policy: propagate))
{
  using M = RankData::Method;
  using P = RankData::NaNPolicy;

  const double NaN = std::numeric_limits<double>::quiet_NaN();
  const std::vector<double> a{0.0, 2.0, 3.0, NaN, -2.0, NaN};

  // default method=Average; explicit nan_policy=Propagate
  auto r = RankData::rankdata(a, M::Average, P::Propagate);

  TEST_EQUAL(r.size(), a.size())
  TEST_TRUE(all_nan(r))
}
END_SECTION

START_SECTION((RankData::rankdata - NaN policy: omit))
{
  using M = RankData::Method;
  using P = RankData::NaNPolicy;

  const double NaN = std::numeric_limits<double>::quiet_NaN();
  const std::vector<double> a{0.0, 2.0, 3.0, NaN, -2.0, NaN};

  auto r = RankData::rankdata(a, M::Average, P::Omit);

  TEST_EQUAL(r.size(), a.size())
  // Expected (ignoring NaNs, rank among values [-2,0,2,3] -> [1,2,3,4]):
  // indices: 0->2, 1->3, 2->4, 3->NaN, 4->1, 5->NaN
  TEST_REAL_SIMILAR(r[0], 2.0)
  TEST_REAL_SIMILAR(r[1], 3.0)
  TEST_REAL_SIMILAR(r[2], 4.0)
  TEST_TRUE(std::isnan(r[3]))
  TEST_REAL_SIMILAR(r[4], 1.0)
  TEST_TRUE(std::isnan(r[5]))
}
END_SECTION

START_SECTION((RankData::rankdata - NaN policy: raise))
{
  using M = RankData::Method;
  using P = RankData::NaNPolicy;

  const double NaN = std::numeric_limits<double>::quiet_NaN();
  const std::vector<double> a{1.0, NaN, 2.0};

  bool threw = false;
  try
  {
    (void)RankData::rankdata(a, M::Average, P::Raise);
  }
  catch (const std::invalid_argument&)
  {
    threw = true;
  }
  TEST_EQUAL(threw, true)
}
END_SECTION

START_SECTION((RankData::rankdata - integer input support))
{
  using M = RankData::Method;

  const std::vector<int> ai{3, 3, 1};

  auto r_avg = RankData::rankdata(ai, M::Average);
  auto r_min = RankData::rankdata(ai, M::Min);
  auto r_max = RankData::rankdata(ai, M::Max);
  auto r_den = RankData::rankdata(ai, M::Dense);
  auto r_ord = RankData::rankdata(ai, M::Ordinal);

  // Values sorted: [1, 3, 3]; ranks: [1, 2, 3]
  // average for the two '3's => (2+3)/2 = 2.5
  expect_equal_real_vectors(r_avg, {2.5, 2.5, 1.0});
  expect_equal_real_vectors(r_min, {2.0, 2.0, 1.0});
  expect_equal_real_vectors(r_max, {3.0, 3.0, 1.0});
  expect_equal_real_vectors(r_den, {2.0, 2.0, 1.0});
  // ordinal should keep original order for ties (first '3' before second '3')
  expect_equal_real_vectors(r_ord, {2.0, 3.0, 1.0});
}
END_SECTION

START_SECTION((RankData free-function overloads by value))
{
  using M = RankData::Method;
  using P = RankData::NaNPolicy;

  // double (by value)
  std::vector<double> ad {0.0, 2.0, 2.0};
  auto rd = rankdata_double(ad, M::Min, P::Propagate);
  std::vector<double> expd {1.0, 2.0, 2.0};
  TEST_EQUAL(rd.size(), expd.size());
  for (Size i = 0; i < rd.size(); ++i) TEST_REAL_SIMILAR(rd[i], expd[i]);

  // float (by value)
  std::vector<float> af {0.0f, 2.0f, 3.0f, 2.0f};
  auto rf = rankdata_float(af, M::Dense, P::Propagate);
  std::vector<double> expf {1.0, 2.0, 3.0, 2.0};
  TEST_EQUAL(rf.size(), expf.size());
  for (Size i = 0; i < rf.size(); ++i) TEST_REAL_SIMILAR(rf[i], expf[i]);

  // int (by value)
  std::vector<int> ai {5, 5, 1};
  auto ri = rankdata_int(ai, M::Average, P::Propagate);
  std::vector<double> expi {2.5, 2.5, 1.0};
  TEST_EQUAL(ri.size(), expi.size());
  for (Size i = 0; i < ri.size(); ++i) TEST_REAL_SIMILAR(ri[i], expi[i]);
}
END_SECTION

START_SECTION((RankData::rankdata - NaN policy: omit with all-NaN input))
{
  using M = RankData::Method;
  using P = RankData::NaNPolicy;

  const double NaN = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> a {NaN, NaN, NaN};

  auto r = RankData::rankdata(a, M::Average, P::Omit);
  TEST_EQUAL(r.size(), a.size());
  for (double x : r) TEST_EQUAL(std::isnan(x), true);
}
END_SECTION

START_SECTION((RankData::rankdata - ordinal stability across two tie blocks))
{
  using M = RankData::Method;

  // Two tie groups: {1,1} then {3,3}; ordinal assigns by stable sorted order.
  // Input indices:   0   1    2   3
  std::vector<double> a {3.0, 3.0, 1.0, 1.0};

  auto r = RankData::rankdata(a, M::Ordinal);
  // Sorted stable indices by value: [2,3,0,1] -> ranks map back to [3,4,1,2]
  std::vector<double> exp {3.0, 4.0, 1.0, 2.0};

  TEST_EQUAL(r.size(), exp.size());
  for (Size i = 0; i < r.size(); ++i) TEST_REAL_SIMILAR(r[i], exp[i]);
}
END_SECTION

START_SECTION((RankData::rankdata - empty input))
{
  std::vector<double> empty;
  auto r = RankData::rankdata(empty);
  TEST_EQUAL(r.size(), 0)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

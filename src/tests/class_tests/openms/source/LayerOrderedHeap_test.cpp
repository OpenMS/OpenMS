// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Christie Mathews $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/DATASTRUCTURES/LayerOrderedHeap.h>
///////////////////////////

#include <algorithm>
#include <functional>
#include <numeric>
#include <random>
#include <vector>

using namespace OpenMS;
using namespace std;

START_TEST(LayerOrderedHeap, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
START_SECTION((LayerOrderedHeap()))
{
  LayerOrderedHeap<double> loh;
  TEST_EQUAL(loh.empty(), true)
  TEST_EQUAL(loh.size(), 0)
  TEST_EQUAL(loh.numLayers(), 0)
}
END_SECTION

START_SECTION((explicit LayerOrderedHeap(std::vector<T>& data, Compare comp)))
{
  vector<double> data = {5.0, 3.0, 8.0, 1.0, 9.0, 2.0, 7.0};
  LayerOrderedHeap<double> loh(data);
  TEST_EQUAL(loh.size(), 7)
  TEST_EQUAL(loh.empty(), false)
  TEST_NOT_EQUAL(loh.numLayers(), 0)
}
END_SECTION

START_SECTION(([EXTRA] layer ordering invariant))
{
  mt19937 rng(42);
  uniform_real_distribution<double> dist(0.0, 1000.0);

  vector<double> data(1000);
  for (auto& x : data) x = dist(rng);

  LayerOrderedHeap<double> loh(data);

  for (size_t i = 0; i + 1 < loh.numLayers(); ++i)
  {
    const auto& cur = loh.layer(i);
    const auto& nxt = loh.layer(i + 1);
    double cur_max = *max_element(data.begin() + cur.begin, data.begin() + cur.end);
    double nxt_min = *min_element(data.begin() + nxt.begin, data.begin() + nxt.end);
    TEST_EQUAL(cur_max <= nxt_min, true)
  }
}
END_SECTION

START_SECTION(([EXTRA] descending layer ordering))
{
  mt19937 rng(123);
  uniform_real_distribution<double> dist(0.0, 1000.0);

  vector<double> data(500);
  for (auto& x : data) x = dist(rng);

  LayerOrderedHeap<double, greater<double>> loh(data, greater<double>());

  for (size_t i = 0; i + 1 < loh.numLayers(); ++i)
  {
    const auto& cur = loh.layer(i);
    const auto& nxt = loh.layer(i + 1);
    // with greater<>, layer[i] has LARGER values than layer[i+1]
    double cur_min = *min_element(data.begin() + cur.begin, data.begin() + cur.end);
    double nxt_max = *max_element(data.begin() + nxt.begin, data.begin() + nxt.end);
    TEST_EQUAL(cur_min >= nxt_max, true)
  }
}
END_SECTION

START_SECTION(([EXTRA] geometric layer sizes))
{
  vector<double> data(100);
  iota(data.begin(), data.end(), 0.0);

  LayerOrderedHeap<double> loh(data);

  size_t expected_size = 1;
  size_t remaining = 100;
  for (size_t i = 0; i < loh.numLayers(); ++i)
  {
    size_t actual = std::min(expected_size, remaining);
    TEST_EQUAL(loh.layer(i).size(), actual)
    remaining -= actual;
    expected_size *= 2;
  }
}
END_SECTION

START_SECTION(([EXTRA] single element))
{
  vector<double> data = {42.0};
  LayerOrderedHeap<double> loh(data);
  TEST_EQUAL(loh.size(), 1)
  TEST_EQUAL(loh.numLayers(), 1)
  TEST_EQUAL(loh.layer(0).min_val, 42.0)
  TEST_EQUAL(loh.layer(0).max_val, 42.0)
}
END_SECTION

START_SECTION(([EXTRA] empty input))
{
  vector<double> data;
  LayerOrderedHeap<double> loh(data);
  TEST_EQUAL(loh.empty(), true)
  TEST_EQUAL(loh.numLayers(), 0)
}
END_SECTION

START_SECTION((size_t layerForRank(size_t k) const))
{
  vector<double> data = {10.0, 20.0, 30.0};
  LayerOrderedHeap<double> loh(data);

  // layer 0 has 1 element (rank 0), layer 1 has 2 elements (ranks 1-2)
  TEST_EQUAL(loh.layerForRank(0), 0)
  TEST_EQUAL(loh.layerForRank(1), 1)
  TEST_EQUAL(loh.layerForRank(2), 1)
}
END_SECTION

START_SECTION(([EXTRA] layer min/max metadata))
{
  mt19937 rng(99);
  uniform_real_distribution<double> dist(0.0, 1000.0);

  vector<double> data(200);
  for (auto& x : data) x = dist(rng);

  LayerOrderedHeap<double> loh(data);

  for (size_t i = 0; i < loh.numLayers(); ++i)
  {
    const auto& lay = loh.layer(i);
    double actual_min = *min_element(data.begin() + lay.begin, data.begin() + lay.end);
    double actual_max = *max_element(data.begin() + lay.begin, data.begin() + lay.end);
    TEST_REAL_SIMILAR(lay.min_val, actual_min)
    TEST_REAL_SIMILAR(lay.max_val, actual_max)
  }
}
END_SECTION

END_TEST

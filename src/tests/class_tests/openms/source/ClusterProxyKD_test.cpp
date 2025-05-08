// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Veit $
// $Authors: Johannes Veit $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmKD.h>

using namespace OpenMS;
using namespace std;

START_TEST(ClusterProxyKD, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ClusterProxyKD* ptr = nullptr;
ClusterProxyKD* nullPointer = nullptr;

START_SECTION((ClusterProxyKD()))
  ptr = new ClusterProxyKD();
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~ClusterProxyKD()))
  delete ptr;
END_SECTION

START_SECTION((ClusterProxyKD(Size size, double avg_distance, Size center_index)))
  ptr = new ClusterProxyKD(1, 0.2, 3);
  TEST_NOT_EQUAL(ptr, nullPointer);
  delete ptr;
END_SECTION

ClusterProxyKD proxy_0;
ClusterProxyKD proxy_1(10, 0.01, 4);
ClusterProxyKD proxy_2(9, 0.001, 3);
ClusterProxyKD proxy_3(9, 0.01, 2);
ClusterProxyKD proxy_4(9, 0.01, 1);
ClusterProxyKD proxy_5 = proxy_1;

START_SECTION((ClusterProxyKD(const ClusterProxyKD& rhs)))
  ptr = new ClusterProxyKD(proxy_1);
  TEST_NOT_EQUAL(ptr, nullPointer);
  TEST_EQUAL(ptr->getSize(), proxy_1.getSize());
  TEST_REAL_SIMILAR(ptr->getAvgDistance(), proxy_1.getAvgDistance());
  TEST_EQUAL(ptr->getCenterIndex(), proxy_1.getCenterIndex());
  delete ptr;
END_SECTION

START_SECTION((ClusterProxyKD& operator=(const ClusterProxyKD& rhs)))
  ClusterProxyKD proxy_5 = proxy_1;
  TEST_EQUAL(proxy_5.getSize(), proxy_1.getSize());
  TEST_REAL_SIMILAR(proxy_5.getAvgDistance(), proxy_1.getAvgDistance());
  TEST_EQUAL(proxy_5.getCenterIndex(), proxy_1.getCenterIndex());
END_SECTION

START_SECTION((bool operator<(const ClusterProxyKD& rhs) const))
  TEST_EQUAL(proxy_1 < proxy_2, true)
  TEST_EQUAL(proxy_1 < proxy_3, true)
  TEST_EQUAL(proxy_1 < proxy_4, true)
  TEST_EQUAL(proxy_2 < proxy_3, true)
  TEST_EQUAL(proxy_2 < proxy_4, true)
  TEST_EQUAL(proxy_3 < proxy_4, true)
  TEST_EQUAL(proxy_2 < proxy_1, false)
  TEST_EQUAL(proxy_3 < proxy_1, false)
  TEST_EQUAL(proxy_4 < proxy_1, false)
  TEST_EQUAL(proxy_3 < proxy_2, false)
  TEST_EQUAL(proxy_4 < proxy_2, false)
  TEST_EQUAL(proxy_4 < proxy_3, false)
  TEST_EQUAL(proxy_1 < proxy_1, false)
END_SECTION

START_SECTION((bool operator!=(const ClusterProxyKD& rhs) const))
  TEST_EQUAL(proxy_0 != proxy_0, false)
  TEST_EQUAL(proxy_1 != proxy_1, false)
  TEST_EQUAL(proxy_1 != proxy_5, false)
  TEST_FALSE(proxy_0 == proxy_1)
  TEST_FALSE(proxy_1 == proxy_2)
END_SECTION

START_SECTION((bool operator==(const ClusterProxyKD& rhs) const))
  TEST_TRUE(proxy_0 == proxy_0)
  TEST_TRUE(proxy_1 == proxy_1)
  TEST_TRUE(proxy_1 == proxy_5)
  TEST_EQUAL(proxy_0 == proxy_1, false)
  TEST_EQUAL(proxy_1 == proxy_2, false)
END_SECTION

START_SECTION((Size getSize() const))
  TEST_EQUAL(proxy_0.getSize(), 0)
  TEST_EQUAL(proxy_1.getSize(), 10)
END_SECTION

START_SECTION((bool isValid() const))
  TEST_EQUAL(proxy_0.isValid(), false)
  TEST_EQUAL(proxy_1.isValid(), true)
END_SECTION

START_SECTION((double getAvgDistance() const))
  TEST_REAL_SIMILAR(proxy_0.getAvgDistance(), 0.0)
  TEST_REAL_SIMILAR(proxy_1.getAvgDistance(), 0.01)
END_SECTION

START_SECTION((Size getCenterIndex() const))
  TEST_EQUAL(proxy_0.getCenterIndex(), 0)
  TEST_EQUAL(proxy_1.getCenterIndex(), 4)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST

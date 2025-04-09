// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Veit $
// $Authors: Johannes Veit $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/ANALYSIS/QUANTITATION/KDTreeFeatureNode.h>
#include <OpenMS/ANALYSIS/QUANTITATION/KDTreeFeatureMaps.h>
#include <OpenMS/KERNEL/FeatureMap.h>

using namespace OpenMS;
using namespace std;

START_TEST(KDTreeFeatureNode, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

Feature f1;
f1.setCharge(2);
f1.setIntensity(100);
f1.setMZ(400);
f1.setRT(1000);

Feature f2;
f2.setCharge(3);
f2.setIntensity(1000);
f2.setMZ(500);
f2.setRT(2000);

FeatureMap fmap;
fmap.push_back(f1);
fmap.push_back(f2);

vector<FeatureMap> fmaps;
fmaps.push_back(fmap);

Param p;
p.setValue("rt_tol", 100);
p.setValue("mz_tol", 10);
p.setValue("mz_unit", "ppm");

KDTreeFeatureMaps* kd_data_ptr = new KDTreeFeatureMaps(fmaps, p);

KDTreeFeatureNode* ptr = nullptr;
KDTreeFeatureNode* nullPointer = nullptr;

START_SECTION((KDTreeFeatureNode(KDTreeFeatureMaps* data, Size idx)))
  ptr = new KDTreeFeatureNode(kd_data_ptr, 0);
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~KDTreeFeatureNode()))
  delete ptr;
END_SECTION

KDTreeFeatureNode node_1(kd_data_ptr, 1);

START_SECTION((KDTreeFeatureNode(const KDTreeFeatureNode& rhs)))
  ptr = new KDTreeFeatureNode(node_1);
  TEST_NOT_EQUAL(ptr, nullPointer)
  TEST_EQUAL(ptr->getIndex(), node_1.getIndex())
  TEST_REAL_SIMILAR((*ptr)[0], node_1[0])
  TEST_REAL_SIMILAR((*ptr)[1], node_1[1])
  delete ptr;
END_SECTION

START_SECTION((KDTreeFeatureNode& operator=(KDTreeFeatureNode const& rhs)))
  KDTreeFeatureNode node_2 = node_1;
  TEST_EQUAL(node_2.getIndex(), node_1.getIndex())
  TEST_REAL_SIMILAR(node_2[0], node_1[0])
  TEST_REAL_SIMILAR(node_2[1], node_1[1])
END_SECTION

START_SECTION((Size getIndex() const))
  TEST_EQUAL(node_1.getIndex(), 1)
END_SECTION

START_SECTION((value_type operator[](Size i) const))
  TEST_REAL_SIMILAR(node_1[0], 2000)
  TEST_REAL_SIMILAR(node_1[1], 500)
END_SECTION

delete kd_data_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST

// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Axel Walter $
// $Authors: Axel Walter $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/QC/FeatureSummary.h>
#include <OpenMS/KERNEL/FeatureMap.h>

///////////////////////////

START_TEST(FeatureSummary, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

FeatureSummary* ptr = nullptr;
FeatureSummary* null_ptr = nullptr;

START_SECTION((FeatureSummary()))
{
  ptr = new FeatureSummary();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION((virtual ~FeatureSummary()))
{
  delete ptr;
}
END_SECTION

START_SECTION((Result compute(const FeatureMap& feature_map)))
{
  FeatureSummary fs;

  // three features, two carry an "rt_deviation" meta value (10 and 20 sec)
  FeatureMap fm;
  Feature f1; f1.setMetaValue("rt_deviation", 10.0); fm.push_back(f1);
  Feature f2; f2.setMetaValue("rt_deviation", 20.0); fm.push_back(f2);
  Feature f3; fm.push_back(f3); // no rt_deviation -> ignored in the mean
  FeatureSummary::Result r = fs.compute(fm);
  TEST_EQUAL(r.feature_count, 3)
  TEST_REAL_SIMILAR(r.rt_shift_mean, 15.0) // (10 + 20) / 2

  // empty map -> no features, mean defaults to 0
  FeatureMap empty;
  FeatureSummary::Result r0 = fs.compute(empty);
  TEST_EQUAL(r0.feature_count, 0)
  TEST_REAL_SIMILAR(r0.rt_shift_mean, 0.0)

  // features without any rt_deviation -> count set, mean stays 0
  FeatureMap no_dev;
  no_dev.push_back(Feature());
  no_dev.push_back(Feature());
  FeatureSummary::Result rn = fs.compute(no_dev);
  TEST_EQUAL(rn.feature_count, 2)
  TEST_REAL_SIMILAR(rn.rt_shift_mean, 0.0)
}
END_SECTION

START_SECTION(([FeatureSummary::Result] bool operator==(const Result& rhs) const))
{
  FeatureSummary::Result a, b;
  TEST_EQUAL(a == b, true) // default-constructed: count 0, mean 0

  a.feature_count = 5;
  TEST_EQUAL(a == b, false)
  b.feature_count = 5;
  TEST_EQUAL(a == b, true)

  a.rt_shift_mean = 3.5;
  TEST_EQUAL(a == b, false)
  b.rt_shift_mean = 3.5;
  TEST_EQUAL(a == b, true)
}
END_SECTION

START_SECTION((const std::string& getName() const override))
{
  FeatureSummary fs;
  TEST_STRING_EQUAL(fs.getName(), "Summary of features from featureXML file")
}
END_SECTION

START_SECTION((QCBase::Status requirements() const override))
{
  FeatureSummary fs;
  TEST_EQUAL(fs.requirements() == QCBase::Status(QCBase::Requires::PREFDRFEAT), true)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

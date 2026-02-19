// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FORMAT/FeatureMapArrowIO.h>
///////////////////////////

#include <OpenMS/config.h>

#ifdef WITH_PARQUET

#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/DATASTRUCTURES/ConvexHull2D.h>

#include <arrow/api.h>

using namespace OpenMS;
using namespace std;

START_TEST(FeatureMapArrowIO, "$Id$")

/////////////////////////////////////////////////////////////
// Export tests
/////////////////////////////////////////////////////////////

START_SECTION(exportFeaturesToArrow - empty FeatureMap)
{
  FeatureMap fm;
  auto table = FeatureMapArrowIO::exportFeaturesToArrow(fm);
  TEST_NOT_EQUAL(table, nullptr)
  TEST_EQUAL(table->num_rows(), 0)
  TEST_EQUAL(table->num_columns(), 17)
}
END_SECTION

START_SECTION(exportFeaturesToArrow - single feature with convex hulls and metavalues)
{
  FeatureMap fm;
  Feature f;
  f.setRT(100.5);
  f.setMZ(500.25);
  f.setIntensity(1000.0f);
  f.setCharge(2);
  f.setOverallQuality(0.95f);
  f.setQuality(0, 0.9f);   // RT quality
  f.setQuality(1, 0.85f);  // MZ quality
  f.setWidth(1.5f);
  f.setUniqueId(12345);

  // Add convex hull with 4 points
  ConvexHull2D hull;
  ConvexHull2D::PointArrayType points;
  points.push_back(DPosition<2>(99.0, 499.5));
  points.push_back(DPosition<2>(99.0, 501.0));
  points.push_back(DPosition<2>(102.0, 501.0));
  points.push_back(DPosition<2>(102.0, 499.5));
  hull.setHullPoints(points);
  f.getConvexHulls().push_back(hull);

  // Add metavalues
  f.setMetaValue("my_int", 42);
  f.setMetaValue("my_float", 3.14);
  f.setMetaValue("my_string", String("hello"));

  fm.push_back(f);

  auto table = FeatureMapArrowIO::exportFeaturesToArrow(fm);
  TEST_NOT_EQUAL(table, nullptr)
  TEST_EQUAL(table->num_rows(), 1)
  TEST_EQUAL(table->num_columns(), 17)

  // Verify scalar columns
  auto col_unique_id = std::static_pointer_cast<arrow::Int64Array>(table->GetColumnByName("unique_id")->chunk(0));
  TEST_EQUAL(col_unique_id->Value(0), 12345)

  auto col_parent = table->GetColumnByName("parent_feature_id")->chunk(0);
  TEST_EQUAL(col_parent->IsNull(0), true)

  auto col_depth = std::static_pointer_cast<arrow::Int32Array>(table->GetColumnByName("depth")->chunk(0));
  TEST_EQUAL(col_depth->Value(0), 0)

  auto col_rt = std::static_pointer_cast<arrow::DoubleArray>(table->GetColumnByName("rt")->chunk(0));
  TEST_REAL_SIMILAR(col_rt->Value(0), 100.5)

  auto col_mz = std::static_pointer_cast<arrow::DoubleArray>(table->GetColumnByName("mz")->chunk(0));
  TEST_REAL_SIMILAR(col_mz->Value(0), 500.25)

  auto col_intensity = std::static_pointer_cast<arrow::FloatArray>(table->GetColumnByName("intensity")->chunk(0));
  TEST_REAL_SIMILAR(col_intensity->Value(0), 1000.0f)

  auto col_charge = std::static_pointer_cast<arrow::Int32Array>(table->GetColumnByName("charge")->chunk(0));
  TEST_EQUAL(col_charge->Value(0), 2)

  auto col_oq = std::static_pointer_cast<arrow::FloatArray>(table->GetColumnByName("overall_quality")->chunk(0));
  TEST_REAL_SIMILAR(col_oq->Value(0), 0.95f)

  auto col_qrt = std::static_pointer_cast<arrow::FloatArray>(table->GetColumnByName("quality_rt")->chunk(0));
  TEST_REAL_SIMILAR(col_qrt->Value(0), 0.9f)

  auto col_qmz = std::static_pointer_cast<arrow::FloatArray>(table->GetColumnByName("quality_mz")->chunk(0));
  TEST_REAL_SIMILAR(col_qmz->Value(0), 0.85f)

  auto col_width = std::static_pointer_cast<arrow::FloatArray>(table->GetColumnByName("width")->chunk(0));
  TEST_REAL_SIMILAR(col_width->Value(0), 1.5f)

  // Verify bounding box
  auto col_rt_bb_min = std::static_pointer_cast<arrow::DoubleArray>(table->GetColumnByName("rt_bb_min")->chunk(0));
  auto col_rt_bb_max = std::static_pointer_cast<arrow::DoubleArray>(table->GetColumnByName("rt_bb_max")->chunk(0));
  auto col_mz_bb_min = std::static_pointer_cast<arrow::DoubleArray>(table->GetColumnByName("mz_bb_min")->chunk(0));
  auto col_mz_bb_max = std::static_pointer_cast<arrow::DoubleArray>(table->GetColumnByName("mz_bb_max")->chunk(0));
  TEST_REAL_SIMILAR(col_rt_bb_min->Value(0), 99.0)
  TEST_REAL_SIMILAR(col_rt_bb_max->Value(0), 102.0)
  TEST_REAL_SIMILAR(col_mz_bb_min->Value(0), 499.5)
  TEST_REAL_SIMILAR(col_mz_bb_max->Value(0), 501.0)

  // Verify convex hulls (list<struct{hull_index, points}>)
  auto col_hulls = std::static_pointer_cast<arrow::ListArray>(table->GetColumnByName("convex_hulls")->chunk(0));
  TEST_EQUAL(col_hulls->value_length(0), 1)  // 1 hull
  auto hull_struct = std::static_pointer_cast<arrow::StructArray>(col_hulls->value_slice(0));
  auto hull_idx_arr = std::static_pointer_cast<arrow::Int32Array>(hull_struct->field(0));
  TEST_EQUAL(hull_idx_arr->Value(0), 0)
  auto points_list = std::static_pointer_cast<arrow::ListArray>(hull_struct->field(1));
  TEST_EQUAL(points_list->value_length(0), 4)  // 4 points

  // Verify metavalues (list<struct{name, value, value_type}>)
  auto col_mv = std::static_pointer_cast<arrow::ListArray>(table->GetColumnByName("metavalues")->chunk(0));
  TEST_EQUAL(col_mv->value_length(0), 4)  // 4 metavalues (my_int, my_float, my_string, + FWHM from setWidth)
}
END_SECTION

START_SECTION(exportFeaturesToArrow - subordinate features)
{
  FeatureMap fm;

  // Create top-level feature (unique_id=1) with 2 subordinates
  Feature f1;
  f1.setRT(100.0);
  f1.setMZ(500.0);
  f1.setIntensity(1000.0f);
  f1.setCharge(2);
  f1.setUniqueId(1);

  Feature sub1;
  sub1.setRT(100.1);
  sub1.setMZ(500.1);
  sub1.setIntensity(600.0f);
  sub1.setCharge(2);
  sub1.setUniqueId(2);

  Feature sub2;
  sub2.setRT(100.2);
  sub2.setMZ(500.2);
  sub2.setIntensity(400.0f);
  sub2.setCharge(2);
  sub2.setUniqueId(3);

  // Give subordinate 2 (unique_id=2) its own subordinate (unique_id=4)
  Feature subsub;
  subsub.setRT(100.15);
  subsub.setMZ(500.15);
  subsub.setIntensity(200.0f);
  subsub.setCharge(2);
  subsub.setUniqueId(4);

  sub1.getSubordinates().push_back(subsub);
  f1.getSubordinates().push_back(sub1);
  f1.getSubordinates().push_back(sub2);
  fm.push_back(f1);

  auto table = FeatureMapArrowIO::exportFeaturesToArrow(fm);
  TEST_NOT_EQUAL(table, nullptr)
  TEST_EQUAL(table->num_rows(), 4)

  auto col_uid = std::static_pointer_cast<arrow::Int64Array>(table->GetColumnByName("unique_id")->chunk(0));
  auto col_parent = std::static_pointer_cast<arrow::Int64Array>(table->GetColumnByName("parent_feature_id")->chunk(0));
  auto col_depth = std::static_pointer_cast<arrow::Int32Array>(table->GetColumnByName("depth")->chunk(0));

  // Row 0: top-level feature (uid=1), parent=null, depth=0
  TEST_EQUAL(col_uid->Value(0), 1)
  TEST_EQUAL(col_parent->IsNull(0), true)
  TEST_EQUAL(col_depth->Value(0), 0)

  // Row 1: subordinate (uid=2), parent=1, depth=1
  TEST_EQUAL(col_uid->Value(1), 2)
  TEST_EQUAL(col_parent->Value(1), 1)
  TEST_EQUAL(col_depth->Value(1), 1)

  // Row 2: sub-subordinate (uid=4), parent=2, depth=2
  TEST_EQUAL(col_uid->Value(2), 4)
  TEST_EQUAL(col_parent->Value(2), 2)
  TEST_EQUAL(col_depth->Value(2), 2)

  // Row 3: subordinate (uid=3), parent=1, depth=1
  TEST_EQUAL(col_uid->Value(3), 3)
  TEST_EQUAL(col_parent->Value(3), 1)
  TEST_EQUAL(col_depth->Value(3), 1)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST

#else // WITH_PARQUET

START_TEST(FeatureMapArrowIO, "$Id$")
END_TEST

#endif // WITH_PARQUET

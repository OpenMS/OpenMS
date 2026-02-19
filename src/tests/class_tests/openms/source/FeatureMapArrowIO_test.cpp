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
// Round-trip tests (export -> import)
/////////////////////////////////////////////////////////////

START_SECTION(importFeaturesFromArrow - round-trip with subordinates and hulls and metavalues)
{
  // === Build source FeatureMap ===
  FeatureMap fm_in;

  // Feature 1 (uid=100): with subordinates, hull, and metavalues
  Feature f1;
  f1.setRT(50.0);
  f1.setMZ(400.0);
  f1.setIntensity(500.0f);
  f1.setCharge(3);
  f1.setOverallQuality(0.8f);
  f1.setQuality(0, 0.7f);
  f1.setQuality(1, 0.6f);
  f1.setWidth(2.0f);
  f1.setUniqueId(100);

  // Convex hull with 4 points
  ConvexHull2D hull1;
  ConvexHull2D::PointArrayType pts1;
  pts1.push_back(DPosition<2>(49.0, 399.0));
  pts1.push_back(DPosition<2>(49.0, 401.0));
  pts1.push_back(DPosition<2>(51.0, 401.0));
  pts1.push_back(DPosition<2>(51.0, 399.0));
  hull1.setHullPoints(pts1);
  f1.getConvexHulls().push_back(hull1);

  // Metavalues (note: setWidth also adds FWHM metavalue)
  f1.setMetaValue("my_int", 42);
  f1.setMetaValue("my_float", 3.14);
  f1.setMetaValue("my_string", String("hello"));

  // Subordinate A (uid=101)
  Feature subA;
  subA.setRT(50.5);
  subA.setMZ(400.1);
  subA.setIntensity(200.0f);
  subA.setCharge(3);
  subA.setUniqueId(101);

  // Subordinate B (uid=102) with its own sub-subordinate
  Feature subB;
  subB.setRT(49.5);
  subB.setMZ(399.9);
  subB.setIntensity(300.0f);
  subB.setCharge(3);
  subB.setUniqueId(102);

  // Sub-subordinate C (uid=103) under B
  Feature subC;
  subC.setRT(49.0);
  subC.setMZ(399.8);
  subC.setIntensity(100.0f);
  subC.setCharge(3);
  subC.setUniqueId(103);

  subB.getSubordinates().push_back(subC);
  f1.getSubordinates().push_back(subA);
  f1.getSubordinates().push_back(subB);

  // Feature 2 (uid=200): no subordinates, no hulls
  Feature f2;
  f2.setRT(150.0);
  f2.setMZ(600.0);
  f2.setIntensity(2000.0f);
  f2.setCharge(2);
  f2.setUniqueId(200);

  fm_in.push_back(f1);
  fm_in.push_back(f2);

  // === Export to Arrow ===
  auto table = FeatureMapArrowIO::exportFeaturesToArrow(fm_in);
  TEST_NOT_EQUAL(table, nullptr)

  // Total rows: f1 + subA + subB + subC + f2 = 5
  TEST_EQUAL(table->num_rows(), 5)

  // === Import back ===
  FeatureMap fm_out;
  bool ok = FeatureMapArrowIO::importFeaturesFromArrow(table, fm_out);
  TEST_EQUAL(ok, true)

  // === Verify structure ===
  TEST_EQUAL(fm_out.size(), 2)  // 2 top-level features

  // Feature 1 checks
  const Feature& out_f1 = fm_out[0];
  TEST_REAL_SIMILAR(out_f1.getRT(), 50.0)
  TEST_REAL_SIMILAR(out_f1.getMZ(), 400.0)
  TEST_REAL_SIMILAR(out_f1.getIntensity(), 500.0f)
  TEST_EQUAL(out_f1.getCharge(), 3)
  TEST_REAL_SIMILAR(out_f1.getOverallQuality(), 0.8f)
  TEST_REAL_SIMILAR(out_f1.getQuality(0), 0.7f)
  TEST_REAL_SIMILAR(out_f1.getQuality(1), 0.6f)
  TEST_REAL_SIMILAR(out_f1.getWidth(), 2.0f)
  TEST_EQUAL(out_f1.getUniqueId(), 100)

  // Check convex hull
  TEST_EQUAL(out_f1.getConvexHulls().size(), 1)
  const auto& out_hull = out_f1.getConvexHulls()[0];
  const auto& out_pts = out_hull.getHullPoints();
  TEST_EQUAL(out_pts.size(), 4)
  TEST_REAL_SIMILAR(out_pts[0][0], 49.0)
  TEST_REAL_SIMILAR(out_pts[0][1], 399.0)
  TEST_REAL_SIMILAR(out_pts[2][0], 51.0)
  TEST_REAL_SIMILAR(out_pts[2][1], 401.0)

  // Check metavalues - int stays int, float stays float, string stays string
  TEST_EQUAL(out_f1.metaValueExists("my_int"), true)
  TEST_EQUAL(out_f1.metaValueExists("my_float"), true)
  TEST_EQUAL(out_f1.metaValueExists("my_string"), true)
  TEST_EQUAL(static_cast<int>(out_f1.getMetaValue("my_int")), 42)
  TEST_REAL_SIMILAR(static_cast<double>(out_f1.getMetaValue("my_float")), 3.14)
  TEST_EQUAL(out_f1.getMetaValue("my_string").toString(), "hello")

  // Check metavalue types are preserved
  TEST_EQUAL(out_f1.getMetaValue("my_int").valueType(), DataValue::INT_VALUE)
  TEST_EQUAL(out_f1.getMetaValue("my_float").valueType(), DataValue::DOUBLE_VALUE)
  TEST_EQUAL(out_f1.getMetaValue("my_string").valueType(), DataValue::STRING_VALUE)

  // Check subordinates of feature 1
  TEST_EQUAL(out_f1.getSubordinates().size(), 2)

  const Feature& out_subA = out_f1.getSubordinates()[0];
  TEST_REAL_SIMILAR(out_subA.getRT(), 50.5)
  TEST_REAL_SIMILAR(out_subA.getMZ(), 400.1)
  TEST_REAL_SIMILAR(out_subA.getIntensity(), 200.0f)
  TEST_EQUAL(out_subA.getCharge(), 3)
  TEST_EQUAL(out_subA.getUniqueId(), 101)
  TEST_EQUAL(out_subA.getSubordinates().size(), 0)

  const Feature& out_subB = out_f1.getSubordinates()[1];
  TEST_REAL_SIMILAR(out_subB.getRT(), 49.5)
  TEST_REAL_SIMILAR(out_subB.getMZ(), 399.9)
  TEST_REAL_SIMILAR(out_subB.getIntensity(), 300.0f)
  TEST_EQUAL(out_subB.getCharge(), 3)
  TEST_EQUAL(out_subB.getUniqueId(), 102)
  TEST_EQUAL(out_subB.getSubordinates().size(), 1)

  // Sub-subordinate C under B
  const Feature& out_subC = out_subB.getSubordinates()[0];
  TEST_REAL_SIMILAR(out_subC.getRT(), 49.0)
  TEST_REAL_SIMILAR(out_subC.getMZ(), 399.8)
  TEST_REAL_SIMILAR(out_subC.getIntensity(), 100.0f)
  TEST_EQUAL(out_subC.getCharge(), 3)
  TEST_EQUAL(out_subC.getUniqueId(), 103)
  TEST_EQUAL(out_subC.getSubordinates().size(), 0)

  // Feature 2 checks
  const Feature& out_f2 = fm_out[1];
  TEST_REAL_SIMILAR(out_f2.getRT(), 150.0)
  TEST_REAL_SIMILAR(out_f2.getMZ(), 600.0)
  TEST_REAL_SIMILAR(out_f2.getIntensity(), 2000.0f)
  TEST_EQUAL(out_f2.getCharge(), 2)
  TEST_EQUAL(out_f2.getUniqueId(), 200)
  TEST_EQUAL(out_f2.getSubordinates().size(), 0)
  TEST_EQUAL(out_f2.getConvexHulls().size(), 0)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST

#else // WITH_PARQUET

START_TEST(FeatureMapArrowIO, "$Id$")
END_TEST

#endif // WITH_PARQUET

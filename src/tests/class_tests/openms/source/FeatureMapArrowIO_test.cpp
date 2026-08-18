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
#include <OpenMS/SYSTEM/File.h>
///////////////////////////

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/DATASTRUCTURES/ConvexHull2D.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/DataProcessing.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinHit.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/config.h>
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
  f.setQuality(0, 0.9f);  // RT quality
  f.setQuality(1, 0.85f); // MZ quality
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
  f.setMetaValue("my_string", std::string("hello"));

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

  auto col_oq = std::static_pointer_cast<arrow::FloatArray>(table->GetColumnByName("quality")->chunk(0));
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
  TEST_EQUAL(col_hulls->value_length(0), 1) // 1 hull
  auto hull_struct = std::static_pointer_cast<arrow::StructArray>(col_hulls->value_slice(0));
  auto hull_idx_arr = std::static_pointer_cast<arrow::Int32Array>(hull_struct->field(0));
  TEST_EQUAL(hull_idx_arr->Value(0), 0)
  auto points_list = std::static_pointer_cast<arrow::ListArray>(hull_struct->field(1));
  TEST_EQUAL(points_list->value_length(0), 4) // 4 points

  // Verify metavalues (list<struct{name, value, value_type}>)
  auto col_mv = std::static_pointer_cast<arrow::ListArray>(table->GetColumnByName("metavalues")->chunk(0));
  TEST_EQUAL(col_mv->value_length(0), 3) // 3 metavalues (my_int, my_float, my_string; FWHM is excluded)
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

START_SECTION(importFeaturesFromArrow - round - trip with subordinates and hulls and metavalues)
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
  f1.setMetaValue("my_string", std::string("hello"));
  f1.setMetaValue("test_int_list", DataValue(IntList {1, 2, 3}));
  f1.setMetaValue("test_double_list", DataValue(DoubleList {1.5, 2.5}));
  f1.setMetaValue("test_string_list", DataValue(StringList {"a", "b", "c"}));

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
  TEST_EQUAL(fm_out.size(), 2) // 2 top-level features

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

  // Check list metavalue types are preserved
  TEST_EQUAL(out_f1.getMetaValue("test_int_list").valueType(), DataValue::INT_LIST)
  TEST_EQUAL(out_f1.getMetaValue("test_int_list") == DataValue(IntList {1, 2, 3}), true)
  TEST_EQUAL(out_f1.getMetaValue("test_double_list").valueType(), DataValue::DOUBLE_LIST)
  TEST_EQUAL(out_f1.getMetaValue("test_double_list") == DataValue(DoubleList {1.5, 2.5}), true)
  TEST_EQUAL(out_f1.getMetaValue("test_string_list").valueType(), DataValue::STRING_LIST)
  TEST_EQUAL(out_f1.getMetaValue("test_string_list") == DataValue(StringList {"a", "b", "c"}), true)

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

START_SECTION([EXTRA] importFeaturesFromArrow - large integer metavalue (> INT32_MAX) survives round-trip)
{
  // Regression for ArrowIOHelpers::readMetaValues integer parsing (finding #19):
  // an "int"-typed meta-value whose value exceeds INT32_MAX must be parsed with
  // std::stoll and stored as SignedSize. Before the fix std::stol -> static_cast<int>
  // truncated 5000000000 down to 705032704.
  FeatureMap fm_in;
  Feature f;
  f.setRT(10.0);
  f.setMZ(200.0);
  f.setIntensity(100.0f);
  f.setCharge(1);
  f.setUniqueId(555);
  f.setMetaValue("big", (OpenMS::SignedSize)5000000000LL); // > INT32_MAX (2147483647)
  fm_in.push_back(f);

  // Round-trip through the Arrow table (export then import).
  auto table = FeatureMapArrowIO::exportFeaturesToArrow(fm_in);
  TEST_NOT_EQUAL(table, nullptr)
  TEST_EQUAL(table->num_rows(), 1)

  FeatureMap fm_out;
  bool ok = FeatureMapArrowIO::importFeaturesFromArrow(table, fm_out);
  TEST_EQUAL(ok, true)
  TEST_EQUAL(fm_out.size(), 1)

  const Feature& out_f = fm_out[0];
  TEST_EQUAL(out_f.metaValueExists("big"), true)
  // Stored as an integer meta-value, not truncated to a smaller type on import.
  TEST_EQUAL(out_f.getMetaValue("big").valueType(), DataValue::INT_VALUE)
  TEST_EQUAL((OpenMS::SignedSize)out_f.getMetaValue("big"), 5000000000LL)
}
END_SECTION

/////////////////////////////////////////////////////////////
// PSM export and round-trip tests
/////////////////////////////////////////////////////////////

START_SECTION(exportPSMsToArrow - empty FeatureMap)
{
  FeatureMap fm;
  auto table = FeatureMapArrowIO::exportPSMsToArrow(fm);
  TEST_NOT_EQUAL(table, nullptr)
  TEST_EQUAL(table->num_rows(), 0)
  // Should have feature_unique_id column + all QPX PSM columns
  auto feature_id_col = table->GetColumnByName("feature_unique_id");
  TEST_NOT_EQUAL(feature_id_col, nullptr)
}
END_SECTION

START_SECTION(exportPSMsToArrow - feature and unassigned PSMs)
{
  FeatureMap fm;

  // Set up ProteinIdentification
  ProteinIdentification prot_id;
  prot_id.setIdentifier("run_1");
  prot_id.setSearchEngine("TestEngine");
  prot_id.setScoreType("TestScore");
  fm.getProteinIdentifications().push_back(prot_id);

  // Feature (uid=1000) with 1 PeptideIdentification containing 1 PeptideHit
  Feature f1;
  f1.setRT(100.0);
  f1.setMZ(500.0);
  f1.setIntensity(1000.0f);
  f1.setCharge(2);
  f1.setUniqueId(1000);

  PeptideIdentification pep_id1;
  pep_id1.setIdentifier("run_1");
  pep_id1.setScoreType("TestScore");
  pep_id1.setHigherScoreBetter(true);
  pep_id1.setRT(100.0);
  pep_id1.setMZ(500.0);
  PeptideHit hit1;
  hit1.setSequence(AASequence::fromString("PEPTIDER"));
  hit1.setScore(0.95);
  hit1.setCharge(2);
  pep_id1.getHits().push_back(hit1);
  f1.getPeptideIdentifications().push_back(pep_id1);

  fm.push_back(f1);

  // 1 unassigned PeptideIdentification with 1 PeptideHit
  PeptideIdentification pep_id2;
  pep_id2.setIdentifier("run_1");
  pep_id2.setScoreType("TestScore");
  pep_id2.setHigherScoreBetter(true);
  pep_id2.setRT(200.0);
  pep_id2.setMZ(600.0);
  PeptideHit hit2;
  hit2.setSequence(AASequence::fromString("ACDEFGHIK"));
  hit2.setScore(0.80);
  hit2.setCharge(3);
  pep_id2.getHits().push_back(hit2);
  fm.getUnassignedPeptideIdentifications().push_back(pep_id2);

  // Export PSMs to Arrow
  auto table = FeatureMapArrowIO::exportPSMsToArrow(fm);
  TEST_NOT_EQUAL(table, nullptr)
  TEST_EQUAL(table->num_rows(), 2)

  // Verify feature_unique_id column
  auto feature_id_chunked = table->GetColumnByName("feature_unique_id");
  TEST_NOT_EQUAL(feature_id_chunked, nullptr)
  auto feature_id_arr = std::static_pointer_cast<arrow::Int64Array>(feature_id_chunked->chunk(0));

  // Row 0: feature PSM -> feature_unique_id = 1000
  TEST_EQUAL(feature_id_arr->IsNull(0), false)
  TEST_EQUAL(feature_id_arr->Value(0), 1000)

  // Row 1: unassigned PSM -> feature_unique_id = null
  TEST_EQUAL(feature_id_arr->IsNull(1), true)

  // Verify sequence column
  auto seq_col = std::static_pointer_cast<arrow::StringArray>(table->GetColumnByName("sequence")->chunk(0));
  TEST_EQUAL(seq_col->GetString(0), "PEPTIDER")
  TEST_EQUAL(seq_col->GetString(1), "ACDEFGHIK")

  // Verify score column
  auto score_col = std::static_pointer_cast<arrow::DoubleArray>(table->GetColumnByName("score")->chunk(0));
  TEST_REAL_SIMILAR(score_col->Value(0), 0.95)
  TEST_REAL_SIMILAR(score_col->Value(1), 0.80)

  // Verify charge column
  auto charge_col = std::static_pointer_cast<arrow::Int32Array>(table->GetColumnByName("precursor_charge")->chunk(0));
  TEST_EQUAL(charge_col->Value(0), 2)
  TEST_EQUAL(charge_col->Value(1), 3)
}
END_SECTION

START_SECTION(importPSMsFromArrow - PSM round - trip)
{
  // === Build source FeatureMap ===
  FeatureMap fm_in;

  // Set up ProteinIdentification
  ProteinIdentification prot_id;
  prot_id.setIdentifier("run_1");
  prot_id.setSearchEngine("TestEngine");
  prot_id.setScoreType("TestScore");
  fm_in.getProteinIdentifications().push_back(prot_id);

  // Feature (uid=1000) with 1 PeptideIdentification
  Feature f1;
  f1.setRT(100.0);
  f1.setMZ(500.0);
  f1.setIntensity(1000.0f);
  f1.setCharge(2);
  f1.setUniqueId(1000);

  PeptideIdentification pep_id1;
  pep_id1.setIdentifier("run_1");
  pep_id1.setScoreType("TestScore");
  pep_id1.setHigherScoreBetter(true);
  pep_id1.setRT(100.0);
  pep_id1.setMZ(500.0);
  PeptideHit hit1;
  hit1.setSequence(AASequence::fromString("PEPTIDER"));
  hit1.setScore(0.95);
  hit1.setCharge(2);
  pep_id1.getHits().push_back(hit1);
  f1.getPeptideIdentifications().push_back(pep_id1);

  fm_in.push_back(f1);

  // Unassigned PeptideIdentification
  PeptideIdentification pep_id2;
  pep_id2.setIdentifier("run_1");
  pep_id2.setScoreType("TestScore");
  pep_id2.setHigherScoreBetter(true);
  pep_id2.setRT(200.0);
  pep_id2.setMZ(600.0);
  PeptideHit hit2;
  hit2.setSequence(AASequence::fromString("ACDEFGHIK"));
  hit2.setScore(0.80);
  hit2.setCharge(3);
  pep_id2.getHits().push_back(hit2);
  fm_in.getUnassignedPeptideIdentifications().push_back(pep_id2);

  // === Export features and PSMs ===
  auto features_table = FeatureMapArrowIO::exportFeaturesToArrow(fm_in);
  TEST_NOT_EQUAL(features_table, nullptr)

  auto psm_table = FeatureMapArrowIO::exportPSMsToArrow(fm_in);
  TEST_NOT_EQUAL(psm_table, nullptr)
  TEST_EQUAL(psm_table->num_rows(), 2)

  // === Import: first features, then PSMs ===
  FeatureMap fm_out;
  fm_out.setProteinIdentifications(fm_in.getProteinIdentifications());

  bool ok = FeatureMapArrowIO::importFeaturesFromArrow(features_table, fm_out);
  TEST_EQUAL(ok, true)
  TEST_EQUAL(fm_out.size(), 1)

  ok = FeatureMapArrowIO::importPSMsFromArrow(psm_table, fm_out);
  TEST_EQUAL(ok, true)

  // === Verify feature PSMs ===
  const Feature& out_f1 = fm_out[0];
  TEST_EQUAL(out_f1.getUniqueId(), 1000)
  TEST_EQUAL(out_f1.getPeptideIdentifications().size(), 1)

  const PeptideIdentification& out_pid1 = out_f1.getPeptideIdentifications()[0];
  TEST_EQUAL(out_pid1.getHits().size(), 1)
  TEST_EQUAL(out_pid1.getScoreType(), "TestScore")
  TEST_REAL_SIMILAR(out_pid1.getRT(), 100.0)

  const PeptideHit& out_hit1 = out_pid1.getHits()[0];
  TEST_EQUAL(out_hit1.getSequence().toString(), "PEPTIDER")
  TEST_REAL_SIMILAR(out_hit1.getScore(), 0.95)
  TEST_EQUAL(out_hit1.getCharge(), 2)

  // === Verify unassigned PSMs ===
  TEST_EQUAL(fm_out.getUnassignedPeptideIdentifications().size(), 1)

  const PeptideIdentification& out_pid2 = fm_out.getUnassignedPeptideIdentifications()[0];
  TEST_EQUAL(out_pid2.getHits().size(), 1)
  TEST_EQUAL(out_pid2.getScoreType(), "TestScore")

  const PeptideHit& out_hit2 = out_pid2.getHits()[0];
  TEST_EQUAL(out_hit2.getSequence().toString(), "ACDEFGHIK")
  TEST_REAL_SIMILAR(out_hit2.getScore(), 0.80)
  TEST_EQUAL(out_hit2.getCharge(), 3)
}
END_SECTION

START_SECTION([EXTRA] importPSMsFromArrow - interleaved P_ID rows group by value not adjacency)
{
  // Regression for FeatureMapArrowIO::importPSMsFromArrow grouping (finding #20):
  // PSM rows must be grouped into PeptideIdentifications by their P_ID VALUE, not by
  // row adjacency. An externally sorted/filtered Parquet file can interleave rows of
  // different P_IDs; the buggy adjacency-based grouping split each interleaved
  // PeptideIdentification into several single-hit PeptideIdentifications.

  // === Build a FeatureMap with 2 features, each with one 2-hit PeptideIdentification ===
  FeatureMap fm_in;

  ProteinIdentification prot_id;
  prot_id.setIdentifier("run_ilv");
  prot_id.setScoreType("TestScore");
  prot_id.setHigherScoreBetter(true);
  fm_in.getProteinIdentifications().push_back(prot_id);

  auto make_feature = [](UInt64 uid, double rt, double mz, const std::string& seqA,
                         double scoreA, const std::string& seqB, double scoreB) -> Feature
  {
    Feature f;
    f.setRT(rt);
    f.setMZ(mz);
    f.setIntensity(1000.0f);
    f.setCharge(2);
    f.setUniqueId(uid);

    PeptideIdentification pep_id;
    pep_id.setIdentifier("run_ilv");
    pep_id.setScoreType("TestScore");
    pep_id.setHigherScoreBetter(true);
    pep_id.setRT(rt);
    pep_id.setMZ(mz);

    PeptideHit h1;
    h1.setSequence(AASequence::fromString(seqA));
    h1.setScore(scoreA);
    h1.setCharge(2);
    pep_id.getHits().push_back(h1);

    PeptideHit h2;
    h2.setSequence(AASequence::fromString(seqB));
    h2.setScore(scoreB);
    h2.setCharge(2);
    pep_id.getHits().push_back(h2);

    f.getPeptideIdentifications().push_back(pep_id);
    return f;
  };

  fm_in.push_back(make_feature(9001, 100.0, 500.0, "PEPTIDER", 0.95, "PEPTIDEK", 0.80));
  fm_in.push_back(make_feature(9002, 200.0, 600.0, "ACDEFGHIK", 0.90, "ACDEFGHIR", 0.70));

  // === Export PSMs: 2 pep-ids x 2 hits = 4 rows ===
  // Row layout: [0]=f9001 hit0 (P_ID 0), [1]=f9001 hit1 (P_ID 0),
  //             [2]=f9002 hit0 (P_ID 1), [3]=f9002 hit1 (P_ID 1)
  auto psm_table = FeatureMapArrowIO::exportPSMsToArrow(fm_in);
  TEST_NOT_EQUAL(psm_table, nullptr)
  TEST_EQUAL(psm_table->num_rows(), 4)

  // === Permute rows into interleaved P_ID order {0, 2, 1, 3} ===
  // Result P_ID sequence becomes 0, 1, 0, 1 -- the adjacency-based grouping would
  // create four single-hit PeptideIdentifications; grouping by value keeps two 2-hit
  // PeptideIdentifications. Built by concatenating single-row table slices (equivalent
  // to arrow::compute::Take with indices {0,2,1,3}, without the compute dependency).
  auto interleaved_result = arrow::ConcatenateTables({
    psm_table->Slice(0, 1),
    psm_table->Slice(2, 1),
    psm_table->Slice(1, 1),
    psm_table->Slice(3, 1)});
  TEST_TRUE(interleaved_result.ok())
  auto interleaved = *interleaved_result;
  TEST_EQUAL(interleaved->num_rows(), 4)

  // === Import into a fresh FeatureMap whose 2 features carry the matching unique ids ===
  FeatureMap fm_out;
  fm_out.setProteinIdentifications(fm_in.getProteinIdentifications());
  Feature of1;
  of1.setUniqueId(9001);
  Feature of2;
  of2.setUniqueId(9002);
  fm_out.push_back(of1);
  fm_out.push_back(of2);

  bool ok = FeatureMapArrowIO::importPSMsFromArrow(interleaved, fm_out);
  TEST_EQUAL(ok, true)

  // === Each feature must receive exactly ONE PeptideIdentification with TWO hits ===
  // (buggy adjacency grouping yields two 1-hit PeptideIdentifications per feature).
  TEST_EQUAL(fm_out.size(), 2)
  TEST_EQUAL(fm_out[0].getUniqueId(), 9001)
  TEST_EQUAL(fm_out[0].getPeptideIdentifications().size(), 1)
  TEST_EQUAL(fm_out[0].getPeptideIdentifications()[0].getHits().size(), 2)

  TEST_EQUAL(fm_out[1].getUniqueId(), 9002)
  TEST_EQUAL(fm_out[1].getPeptideIdentifications().size(), 1)
  TEST_EQUAL(fm_out[1].getPeptideIdentifications()[0].getHits().size(), 2)
}
END_SECTION

/////////////////////////////////////////////////////////////
// Full Parquet directory round-trip test
/////////////////////////////////////////////////////////////

START_SECTION(exportToParquet / importFromParquet - full round - trip)
{
  FeatureMap fm;

  // --- ProteinIdentification ---
  ProteinIdentification prot_id;
  prot_id.setIdentifier("run_full_1");
  prot_id.setSearchEngine("Comet");
  prot_id.setSearchEngineVersion("2023.01");
  prot_id.setScoreType("expect");
  prot_id.setHigherScoreBetter(false);

  ProteinIdentification::SearchParameters sp;
  sp.db = "uniprot_human";
  sp.mass_type = ProteinIdentification::PeakMassType::MONOISOTOPIC;
  sp.precursor_mass_tolerance = 10.0;
  sp.precursor_mass_tolerance_ppm = true;
  sp.fragment_mass_tolerance = 0.02;
  sp.fragment_mass_tolerance_ppm = false;
  sp.digestion_enzyme = *ProteaseDB::getInstance()->getEnzyme("Trypsin");
  sp.missed_cleavages = 2;
  prot_id.setSearchParameters(sp);

  ProteinHit ph;
  ph.setAccession("P12345");
  ph.setScore(0.001);
  prot_id.insertHit(ph);

  ProteinIdentification::ProteinGroup pg;
  pg.probability = 0.99;
  pg.accessions = {"P12345"};
  prot_id.insertProteinGroup(pg);

  fm.setProteinIdentifications({prot_id});

  // --- Feature 1: with subordinate and PSM ---
  Feature f1;
  f1.setRT(100.0);
  f1.setMZ(500.0);
  f1.setIntensity(1000.0f);
  f1.setCharge(2);
  f1.setOverallQuality(0.9f);
  f1.setUniqueId(1001);
  f1.setMetaValue("custom_val", 42);

  // Convex hull
  ConvexHull2D hull;
  std::vector<DPosition<2>> pts = {{99.0, 499.0}, {99.0, 501.0}, {101.0, 501.0}, {101.0, 499.0}};
  hull.setHullPoints(pts);
  f1.getConvexHulls().push_back(hull);

  // Subordinate
  Feature sub;
  sub.setRT(100.5);
  sub.setMZ(500.1);
  sub.setIntensity(500.0f);
  sub.setCharge(2);
  sub.setUniqueId(1002);
  f1.getSubordinates().push_back(sub);

  // PeptideIdentification on feature 1
  PeptideIdentification pep_id1;
  pep_id1.setScoreType("expect");
  pep_id1.setHigherScoreBetter(false);
  pep_id1.setIdentifier("run_full_1");
  PeptideHit pep_hit1;
  pep_hit1.setSequence(AASequence::fromString("PEPTIDER"));
  pep_hit1.setScore(0.001);
  pep_hit1.setCharge(2);
  pep_id1.insertHit(pep_hit1);
  f1.setPeptideIdentifications({pep_id1});

  fm.push_back(f1);

  // --- Feature 2: simple, no PSMs ---
  Feature f2;
  f2.setRT(200.0);
  f2.setMZ(600.0);
  f2.setIntensity(2000.0f);
  f2.setCharge(3);
  f2.setUniqueId(2001);
  fm.push_back(f2);

  // --- Unassigned PeptideIdentification ---
  PeptideIdentification unassigned;
  unassigned.setScoreType("expect");
  unassigned.setHigherScoreBetter(false);
  unassigned.setIdentifier("run_full_1");
  PeptideHit unassigned_hit;
  unassigned_hit.setSequence(AASequence::fromString("ACDEFGHIK"));
  unassigned_hit.setScore(0.5);
  unassigned_hit.setCharge(3);
  unassigned.insertHit(unassigned_hit);
  fm.setUnassignedPeptideIdentifications({unassigned});

  // --- Export to temp directory ---
  std::string tmp_dir;
  NEW_TMP_FILE(tmp_dir)
  tmp_dir += ".fmd";

  TEST_EQUAL(FeatureMapArrowIO::exportToParquet(fm, tmp_dir), true)

  // --- Import back ---
  FeatureMap imported;
  TEST_EQUAL(FeatureMapArrowIO::importFromParquet(tmp_dir, imported), true)

  // --- Verify features ---
  TEST_EQUAL(imported.size(), 2)

  // Feature 1
  TEST_REAL_SIMILAR(imported[0].getRT(), 100.0)
  TEST_REAL_SIMILAR(imported[0].getMZ(), 500.0)
  TEST_REAL_SIMILAR(imported[0].getIntensity(), 1000.0)
  TEST_EQUAL(imported[0].getCharge(), 2)
  TEST_EQUAL(imported[0].getUniqueId(), 1001)
  TEST_EQUAL(imported[0].getSubordinates().size(), 1)
  TEST_EQUAL(imported[0].getConvexHulls().size(), 1)
  TEST_EQUAL(imported[0].getConvexHulls()[0].getHullPoints().size(), 4)
  TEST_EQUAL(int(imported[0].getMetaValue("custom_val")), 42)

  // Feature 1 subordinate
  TEST_REAL_SIMILAR(imported[0].getSubordinates()[0].getRT(), 100.5)
  TEST_EQUAL(imported[0].getSubordinates()[0].getUniqueId(), 1002)

  // Feature 2
  TEST_REAL_SIMILAR(imported[1].getRT(), 200.0)
  TEST_EQUAL(imported[1].getCharge(), 3)
  TEST_EQUAL(imported[1].getUniqueId(), 2001)
  TEST_EQUAL(imported[1].getSubordinates().size(), 0)

  // --- Verify PSMs ---
  TEST_EQUAL(imported[0].getPeptideIdentifications().size(), 1)
  TEST_EQUAL(imported[0].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "PEPTIDER")
  TEST_REAL_SIMILAR(imported[0].getPeptideIdentifications()[0].getHits()[0].getScore(), 0.001)

  TEST_EQUAL(imported.getUnassignedPeptideIdentifications().size(), 1)
  TEST_EQUAL(imported.getUnassignedPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "ACDEFGHIK")

  // --- Verify protein identifications ---
  TEST_EQUAL(imported.getProteinIdentifications().size(), 1)
  // Identifier synthesized on load per IdXMLFile.cpp:530 parity — stored "run_full_1"
  // becomes `<search_engine>_<date>_<UniqueIdGenerator>`. All pep_id collections
  // (per-feature + unassigned) are re-stamped in lock-step.
  const std::string& fm_synth_id = imported.getProteinIdentifications()[0].getIdentifier();
  TEST_NOT_EQUAL(fm_synth_id, "")
  TEST_NOT_EQUAL(fm_synth_id, "run_full_1")
  TEST_STRING_EQUAL(imported[0].getPeptideIdentifications()[0].getIdentifier(), fm_synth_id);
  TEST_STRING_EQUAL(imported.getUnassignedPeptideIdentifications()[0].getIdentifier(), fm_synth_id);
  TEST_EQUAL(imported.getProteinIdentifications()[0].getSearchEngine(), "Comet")
  TEST_EQUAL(imported.getProteinIdentifications()[0].getHits().size(), 1)
  TEST_EQUAL(imported.getProteinIdentifications()[0].getHits()[0].getAccession(), "P12345")
  TEST_EQUAL(imported.getProteinIdentifications()[0].getProteinGroups().size(), 1)
}
END_SECTION

START_SECTION([EXTRA] exportToParquet / importFromParquet - parquet->parquet content idempotency)
{
  // The sections above check memory -> parquet -> memory fidelity. This pins the complementary contract
  // (issue #9460 §3): a *second* serialization is idempotent -- re-exporting an imported bundle and
  // re-importing it yields field-identical content. It catches asymmetric export/import bugs where one
  // side silently drops or rewrites a field. The run identifier is re-synthesized on every load, so we
  // compare content (features, hulls, subordinates, metavalues, PSMs, protein hits/groups) and assert
  // identifier consistency *within* each map, not the volatile run-id string across maps.
  FeatureMap fm;

  ProteinIdentification prot_id;
  prot_id.setIdentifier("run_idem");
  prot_id.setSearchEngine("Comet");
  prot_id.setScoreType("expect");
  prot_id.setHigherScoreBetter(false);
  ProteinHit ph;
  ph.setAccession("Q99999");
  ph.setScore(0.002);
  prot_id.insertHit(ph);
  ProteinIdentification::ProteinGroup pg;
  pg.probability = 0.95;
  pg.accessions = {"Q99999"};
  prot_id.insertProteinGroup(pg);
  fm.setProteinIdentifications({prot_id});

  Feature f1;
  f1.setRT(111.0);
  f1.setMZ(555.5);
  f1.setIntensity(1234.0f);
  f1.setCharge(2);
  f1.setUniqueId(7001);
  f1.setMetaValue("label", "alpha");
  f1.setMetaValue("count", 7);
  ConvexHull2D hull;
  hull.setHullPoints({{110.0, 554.0}, {110.0, 557.0}, {112.0, 557.0}, {112.0, 554.0}});
  f1.getConvexHulls().push_back(hull);
  Feature sub;
  sub.setRT(111.2);
  sub.setMZ(555.6);
  sub.setIntensity(600.0f);
  sub.setCharge(2);
  sub.setUniqueId(7002);
  f1.getSubordinates().push_back(sub);
  PeptideIdentification pep;
  pep.setScoreType("expect");
  pep.setHigherScoreBetter(false);
  pep.setIdentifier("run_idem");
  PeptideHit pep_hit;
  pep_hit.setSequence(AASequence::fromString("PEPTIDEK"));
  pep_hit.setScore(0.002);
  pep_hit.setCharge(2);
  pep.insertHit(pep_hit);
  f1.setPeptideIdentifications({pep});
  fm.push_back(f1);

  Feature f2;
  f2.setRT(222.0);
  f2.setMZ(666.6);
  f2.setIntensity(4321.0f);
  f2.setCharge(3);
  f2.setUniqueId(8001);
  fm.push_back(f2);

  // first serialization round: memory -> parquet -> fm1
  std::string dir1;
  NEW_TMP_FILE(dir1) dir1 += ".fmd";
  TEST_EQUAL(FeatureMapArrowIO::exportToParquet(fm, dir1), true)
  FeatureMap fm1;
  TEST_EQUAL(FeatureMapArrowIO::importFromParquet(dir1, fm1), true)

  // second serialization round: re-export fm1 -> parquet -> fm2
  std::string dir2;
  NEW_TMP_FILE(dir2) dir2 += ".fmd";
  TEST_EQUAL(FeatureMapArrowIO::exportToParquet(fm1, dir2), true)
  FeatureMap fm2;
  TEST_EQUAL(FeatureMapArrowIO::importFromParquet(dir2, fm2), true)

  // the two imported maps must be field-identical
  TEST_EQUAL(fm1.size(), fm2.size())
  TEST_EQUAL(fm1.size(), 2)
  for (Size i = 0; i < fm1.size(); ++i)
  {
    TEST_REAL_SIMILAR(fm1[i].getRT(), fm2[i].getRT())
    TEST_REAL_SIMILAR(fm1[i].getMZ(), fm2[i].getMZ())
    TEST_REAL_SIMILAR(fm1[i].getIntensity(), fm2[i].getIntensity())
    TEST_EQUAL(fm1[i].getCharge(), fm2[i].getCharge())
    TEST_EQUAL(fm1[i].getUniqueId(), fm2[i].getUniqueId())
    TEST_EQUAL(fm1[i].getSubordinates().size(), fm2[i].getSubordinates().size())
    TEST_EQUAL(fm1[i].getConvexHulls().size(), fm2[i].getConvexHulls().size())
  }
  TEST_EQUAL(fm2[0].getMetaValue("label").toString(), "alpha")
  TEST_EQUAL(int(fm2[0].getMetaValue("count")), 7)
  TEST_EQUAL(fm2[0].getConvexHulls()[0].getHullPoints().size(), 4)
  TEST_EQUAL(fm2[0].getSubordinates()[0].getUniqueId(), 7002)
  TEST_EQUAL(fm1[0].getMetaValue("label").toString(), fm2[0].getMetaValue("label").toString())
  TEST_EQUAL(int(fm1[0].getMetaValue("count")), int(fm2[0].getMetaValue("count")))
  const auto& hull_points1 = fm1[0].getConvexHulls()[0].getHullPoints();
  const auto& hull_points2 = fm2[0].getConvexHulls()[0].getHullPoints();
  TEST_EQUAL(hull_points1.size(), hull_points2.size())
  for (Size i = 0; i < hull_points1.size(); ++i)
  {
    TEST_REAL_SIMILAR(hull_points1[i][0], hull_points2[i][0])
    TEST_REAL_SIMILAR(hull_points1[i][1], hull_points2[i][1])
  }
  TEST_REAL_SIMILAR(fm1[0].getSubordinates()[0].getRT(), fm2[0].getSubordinates()[0].getRT())
  TEST_REAL_SIMILAR(fm1[0].getSubordinates()[0].getMZ(), fm2[0].getSubordinates()[0].getMZ())
  TEST_REAL_SIMILAR(fm1[0].getSubordinates()[0].getIntensity(), fm2[0].getSubordinates()[0].getIntensity())
  TEST_EQUAL(fm1[0].getSubordinates()[0].getCharge(), fm2[0].getSubordinates()[0].getCharge())
  TEST_EQUAL(fm1[0].getSubordinates()[0].getUniqueId(), fm2[0].getSubordinates()[0].getUniqueId())
  TEST_EQUAL(fm1[0].getPeptideIdentifications().size(), fm2[0].getPeptideIdentifications().size())
  TEST_EQUAL(fm1[0].getPeptideIdentifications()[0].getHits().size(), fm2[0].getPeptideIdentifications()[0].getHits().size())
  TEST_STRING_EQUAL(fm1[0].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(),
                    fm2[0].getPeptideIdentifications()[0].getHits()[0].getSequence().toString())
  TEST_REAL_SIMILAR(fm1[0].getPeptideIdentifications()[0].getHits()[0].getScore(), fm2[0].getPeptideIdentifications()[0].getHits()[0].getScore())
  TEST_EQUAL(fm1[0].getPeptideIdentifications()[0].getHits()[0].getCharge(), fm2[0].getPeptideIdentifications()[0].getHits()[0].getCharge())
  TEST_STRING_EQUAL(fm2[0].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "PEPTIDEK")
  TEST_REAL_SIMILAR(fm2[0].getPeptideIdentifications()[0].getHits()[0].getScore(), 0.002)
  TEST_EQUAL(fm1.getProteinIdentifications().size(), fm2.getProteinIdentifications().size())
  TEST_EQUAL(fm1.getProteinIdentifications()[0].getHits().size(), fm2.getProteinIdentifications()[0].getHits().size())
  TEST_EQUAL(fm1.getProteinIdentifications()[0].getHits()[0].getAccession(), fm2.getProteinIdentifications()[0].getHits()[0].getAccession())
  TEST_REAL_SIMILAR(fm1.getProteinIdentifications()[0].getHits()[0].getScore(), fm2.getProteinIdentifications()[0].getHits()[0].getScore())
  TEST_EQUAL(fm2.getProteinIdentifications()[0].getHits()[0].getAccession(), "Q99999")
  TEST_EQUAL(fm1.getProteinIdentifications()[0].getProteinGroups().size(), fm2.getProteinIdentifications()[0].getProteinGroups().size())
  const auto& group1 = fm1.getProteinIdentifications()[0].getProteinGroups()[0];
  const auto& group2 = fm2.getProteinIdentifications()[0].getProteinGroups()[0];
  TEST_REAL_SIMILAR(group1.probability, group2.probability)
  TEST_EQUAL(group1.accessions.size(), group2.accessions.size())
  TEST_EQUAL(group1.accessions[0], group2.accessions[0])
  TEST_EQUAL(group2.accessions[0], "Q99999")
  // identifier consistency within each imported map: the stamped pep-id identifier equals the protein run identifier
  TEST_STRING_EQUAL(fm1[0].getPeptideIdentifications()[0].getIdentifier(), fm1.getProteinIdentifications()[0].getIdentifier())
  TEST_STRING_EQUAL(fm2[0].getPeptideIdentifications()[0].getIdentifier(), fm2.getProteinIdentifications()[0].getIdentifier())
}
END_SECTION

/////////////////////////////////////////////////////////////
// FeatureMap-level metadata round-trip tests
/////////////////////////////////////////////////////////////

START_SECTION(exportToParquet / importFromParquet - FeatureMap metadata round - trip(DocumentIdentifier + DataProcessing))
{
  FeatureMap fm;

  // --- Set DocumentIdentifier fields ---
  fm.setIdentifier("test_document_lsid_123");
  fm.setLoadedFilePath("/data/experiments/test_run.featureXML");
  // Note: setLoadedFileType requires an actual file path and reads file content,
  // so we don't set it here. The loaded_file_type is stored informatively but
  // cannot be restored without the original file.

  // --- Set DataProcessing ---
  DataProcessing dp1;
  dp1.getSoftware().setName("FeatureFinderCentroided");
  dp1.getSoftware().setVersion("3.2.0");
  dp1.setCompletionTime(DateTime::fromString("2025-06-15T14:30:00", "yyyy-MM-ddThh:mm:ss"));
  dp1.getProcessingActions().insert(DataProcessing::PEAK_PICKING);
  dp1.getProcessingActions().insert(DataProcessing::FILTERING);
  dp1.setMetaValue("parameter_file", std::string("params.ini"));
  dp1.setMetaValue("num_threads", 8);

  DataProcessing dp2;
  dp2.getSoftware().setName("MapAlignerPoseClustering");
  dp2.getSoftware().setVersion("3.2.0");
  dp2.setCompletionTime(DateTime::fromString("2025-06-15T15:00:00", "yyyy-MM-ddThh:mm:ss"));
  dp2.getProcessingActions().insert(DataProcessing::ALIGNMENT);
  dp2.setMetaValue("max_rt_shift", 300.5);

  fm.setDataProcessing({dp1, dp2});

  // --- Add a simple feature so the FeatureMap is non-empty ---
  Feature f;
  f.setRT(100.0);
  f.setMZ(500.0);
  f.setIntensity(1000.0f);
  f.setCharge(2);
  f.setUniqueId(42);
  fm.push_back(f);

  // --- Set up minimal ProteinIdentification (required for parquet export) ---
  ProteinIdentification prot_id;
  prot_id.setIdentifier("run_1");
  fm.setProteinIdentifications({prot_id});

  // --- Export to temp directory ---
  std::string tmp_dir;
  NEW_TMP_FILE(tmp_dir)
  tmp_dir += ".fmd";

  TEST_EQUAL(FeatureMapArrowIO::exportToParquet(fm, tmp_dir), true)

  // --- Import back ---
  FeatureMap imported;
  TEST_EQUAL(FeatureMapArrowIO::importFromParquet(tmp_dir, imported), true)

  // --- Verify DocumentIdentifier ---
  TEST_EQUAL(imported.getIdentifier(), "test_document_lsid_123")
  TEST_EQUAL(imported.getLoadedFilePath(), "/data/experiments/test_run.featureXML")

  // --- Verify DataProcessing ---
  TEST_EQUAL(imported.getDataProcessing().size(), 2)

  const DataProcessing& out_dp1 = imported.getDataProcessing()[0];
  TEST_EQUAL(out_dp1.getSoftware().getName(), "FeatureFinderCentroided")
  TEST_EQUAL(out_dp1.getSoftware().getVersion(), "3.2.0")
  TEST_EQUAL(out_dp1.getCompletionTime().isValid(), true)
  TEST_EQUAL(out_dp1.getCompletionTime().toString("yyyy-MM-ddThh:mm:ss"), "2025-06-15T14:30:00")
  TEST_EQUAL(out_dp1.getProcessingActions().count(DataProcessing::PEAK_PICKING), 1)
  TEST_EQUAL(out_dp1.getProcessingActions().count(DataProcessing::FILTERING), 1)
  TEST_EQUAL(out_dp1.getProcessingActions().size(), 2)
  TEST_EQUAL(out_dp1.metaValueExists("parameter_file"), true)
  TEST_EQUAL(out_dp1.getMetaValue("parameter_file").toString(), "params.ini")
  TEST_EQUAL(static_cast<int>(out_dp1.getMetaValue("num_threads")), 8)

  const DataProcessing& out_dp2 = imported.getDataProcessing()[1];
  TEST_EQUAL(out_dp2.getSoftware().getName(), "MapAlignerPoseClustering")
  TEST_EQUAL(out_dp2.getSoftware().getVersion(), "3.2.0")
  TEST_EQUAL(out_dp2.getCompletionTime().isValid(), true)
  TEST_EQUAL(out_dp2.getCompletionTime().toString("yyyy-MM-ddThh:mm:ss"), "2025-06-15T15:00:00")
  TEST_EQUAL(out_dp2.getProcessingActions().count(DataProcessing::ALIGNMENT), 1)
  TEST_EQUAL(out_dp2.getProcessingActions().size(), 1)
  TEST_REAL_SIMILAR(static_cast<double>(out_dp2.getMetaValue("max_rt_shift")), 300.5)

  // --- Verify feature data is still correct ---
  TEST_EQUAL(imported.size(), 1)
  TEST_REAL_SIMILAR(imported[0].getRT(), 100.0)
  TEST_EQUAL(imported[0].getUniqueId(), 42)
}
END_SECTION

/////////////////////////////////////////////////////////////
// PSM completeness round-trip test
/////////////////////////////////////////////////////////////

START_SECTION(exportToParquet / importFromParquet - PSM completeness round
              - trip(metavalues, additional_scores, is_decoy, protein_accessions, higher_score_better))
{
  FeatureMap fm;

  // --- ProteinIdentification with higher_score_better=false ---
  ProteinIdentification prot_id;
  prot_id.setIdentifier("run_psm_test");
  prot_id.setSearchEngine("Comet");
  prot_id.setScoreType("expect");
  prot_id.setHigherScoreBetter(false);

  ProteinHit ph1;
  ph1.setAccession("P12345");
  ph1.setScore(0.001);
  prot_id.insertHit(ph1);

  ProteinHit ph2;
  ph2.setAccession("Q67890");
  ph2.setScore(0.01);
  prot_id.insertHit(ph2);

  fm.setProteinIdentifications({prot_id});

  // --- Feature with PeptideIdentification ---
  Feature f1;
  f1.setRT(100.0);
  f1.setMZ(500.0);
  f1.setIntensity(1000.0f);
  f1.setCharge(2);
  f1.setUniqueId(5001);

  PeptideIdentification pep_id1;
  pep_id1.setIdentifier("run_psm_test");
  pep_id1.setScoreType("expect");
  pep_id1.setHigherScoreBetter(false);
  pep_id1.setRT(100.0);
  pep_id1.setMZ(500.25);
  pep_id1.setSpectrumReference("spectrum=42");

  // Hit 1: target with protein accessions, additional scores, metavalues
  PeptideHit hit1;
  hit1.setSequence(AASequence::fromString("PEPTIDER"));
  hit1.setScore(0.001);
  hit1.setCharge(2);
  hit1.setRank(1);
  hit1.setMetaValue("target_decoy", "target");
  hit1.setMetaValue("MS:1002252", 0.95); // additional score: xcorr
  hit1.setMetaValue("MS:1002253", 12.5); // additional score: deltacn
  hit1.setMetaValue("predicted_RT", 99.5);
  hit1.setMetaValue("ion_mobility", 0.85);
  hit1.setMetaValue("custom_psm_int", 42);
  hit1.setMetaValue("custom_psm_str", "test_value");

  // Add PeptideEvidence (protein accessions)
  PeptideEvidence ev1;
  ev1.setProteinAccession("P12345");
  hit1.addPeptideEvidence(ev1);
  PeptideEvidence ev2;
  ev2.setProteinAccession("Q67890");
  hit1.addPeptideEvidence(ev2);

  // Hit 2: decoy, single protein accession
  PeptideHit hit2;
  hit2.setSequence(AASequence::fromString("ACDEFGHIK"));
  hit2.setScore(0.5);
  hit2.setCharge(2);
  hit2.setRank(2);
  hit2.setMetaValue("target_decoy", "decoy");

  PeptideEvidence ev3;
  ev3.setProteinAccession("P12345");
  hit2.addPeptideEvidence(ev3);

  pep_id1.insertHit(hit1);
  pep_id1.insertHit(hit2);
  f1.setPeptideIdentifications({pep_id1});

  fm.push_back(f1);

  // --- Unassigned PeptideIdentification with spectrum-level metavalue ---
  PeptideIdentification pep_id2;
  pep_id2.setIdentifier("run_psm_test");
  pep_id2.setScoreType("expect");
  pep_id2.setHigherScoreBetter(false);
  pep_id2.setRT(200.0);
  pep_id2.setMZ(600.0);
  pep_id2.setMetaValue("spectrum_custom", "spec_meta_val");

  PeptideHit hit3;
  hit3.setSequence(AASequence::fromString("ACDEFGHIK"));
  hit3.setScore(0.1);
  hit3.setCharge(3);
  hit3.setRank(1);
  hit3.setMetaValue("target_decoy", "target");
  pep_id2.insertHit(hit3);
  fm.setUnassignedPeptideIdentifications({pep_id2});

  // --- Export and import ---
  std::string tmp_dir;
  NEW_TMP_FILE(tmp_dir)
  tmp_dir += ".fmd";

  TEST_EQUAL(FeatureMapArrowIO::exportToParquet(fm, tmp_dir), true)

  FeatureMap imported;
  TEST_EQUAL(FeatureMapArrowIO::importFromParquet(tmp_dir, imported), true)

  // --- Verify higher_score_better from ProteinIdentification ---
  TEST_EQUAL(imported[0].getPeptideIdentifications().size(), 1)
  const PeptideIdentification& out_pid1 = imported[0].getPeptideIdentifications()[0];
  TEST_EQUAL(out_pid1.isHigherScoreBetter(), false)
  TEST_EQUAL(out_pid1.getScoreType(), "expect")
  TEST_REAL_SIMILAR(out_pid1.getRT(), 100.0)
  TEST_REAL_SIMILAR(out_pid1.getMZ(), 500.25)

  // --- Verify hit 1 (target with accessions, scores, metavalues) ---
  TEST_EQUAL(out_pid1.getHits().size(), 2)

  // Find the hit with score ~0.001 (hit1)
  const PeptideHit* target_hit = nullptr;
  const PeptideHit* decoy_hit = nullptr;
  for (const auto& h : out_pid1.getHits())
  {
    if (h.getScore() < 0.01) target_hit = &h;
    else
      decoy_hit = &h;
  }
  TEST_NOT_EQUAL(target_hit, nullptr)
  TEST_NOT_EQUAL(decoy_hit, nullptr)

  // target hit: sequence, score, charge
  TEST_EQUAL(target_hit->getSequence().toString(), "PEPTIDER")
  TEST_REAL_SIMILAR(target_hit->getScore(), 0.001)
  TEST_EQUAL(target_hit->getCharge(), 2)

  // target hit: is_decoy -> target_decoy metavalue
  TEST_EQUAL(static_cast<std::string>(target_hit->getMetaValue("target_decoy")), "target")

  // target hit: protein accessions (PeptideEvidence)
  TEST_EQUAL(target_hit->getPeptideEvidences().size(), 2)
  std::set<std::string> accs;
  for (const auto& ev : target_hit->getPeptideEvidences())
  {
    accs.insert(ev.getProteinAccession());
  }
  TEST_EQUAL(accs.count("P12345"), 1)
  TEST_EQUAL(accs.count("Q67890"), 1)

  // decoy hit
  TEST_EQUAL(decoy_hit->getSequence().toString(), "ACDEFGHIK")
  TEST_REAL_SIMILAR(decoy_hit->getScore(), 0.5)
  TEST_EQUAL(static_cast<std::string>(decoy_hit->getMetaValue("target_decoy")), "decoy")
  TEST_EQUAL(decoy_hit->getPeptideEvidences().size(), 1)
  TEST_EQUAL(decoy_hit->getPeptideEvidences()[0].getProteinAccession(), "P12345")

  // --- Verify unassigned PSM ---
  TEST_EQUAL(imported.getUnassignedPeptideIdentifications().size(), 1)
  const PeptideIdentification& out_pid2 = imported.getUnassignedPeptideIdentifications()[0];
  TEST_EQUAL(out_pid2.isHigherScoreBetter(), false)
  TEST_REAL_SIMILAR(out_pid2.getRT(), 200.0)
  TEST_EQUAL(out_pid2.getHits().size(), 1)
  TEST_EQUAL(out_pid2.getHits()[0].getSequence().toString(), "ACDEFGHIK")
  TEST_REAL_SIMILAR(out_pid2.getHits()[0].getScore(), 0.1)
  TEST_EQUAL(out_pid2.getHits()[0].getCharge(), 3)
}
END_SECTION

/////////////////////////////////////////////////////////////
// PSM per-PSM higher_score_better + scan/reference_file_name round-trip
/////////////////////////////////////////////////////////////

START_SECTION(exportToParquet / importFromParquet - per - PSM higher_score_better independent from ProteinIdentification)
{
  FeatureMap fm;

  // ProteinIdentification with higher_score_better=false (run-level)
  ProteinIdentification prot_id;
  prot_id.setIdentifier("run_hsb_test");
  prot_id.setSearchEngine("Comet");
  prot_id.setScoreType("expect");
  prot_id.setHigherScoreBetter(false);
  fm.setProteinIdentifications({prot_id});

  // Feature with PeptideIdentification where higher_score_better=true (differs from run-level)
  Feature f1;
  f1.setRT(100.0);
  f1.setMZ(500.0);
  f1.setIntensity(1000.0f);
  f1.setCharge(2);
  f1.setUniqueId(6001);

  PeptideIdentification pep_id1;
  pep_id1.setIdentifier("run_hsb_test");
  pep_id1.setScoreType("xcorr");
  pep_id1.setHigherScoreBetter(true); // different from run-level false
  pep_id1.setRT(100.0);
  pep_id1.setMZ(500.25);
  pep_id1.setSpectrumReference("spectrum=99");

  PeptideHit hit1;
  hit1.setSequence(AASequence::fromString("PEPTIDER"));
  hit1.setScore(2.5);
  hit1.setCharge(2);
  hit1.setRank(1);
  hit1.setMetaValue("target_decoy", "target");
  pep_id1.insertHit(hit1);
  f1.setPeptideIdentifications({pep_id1});
  fm.push_back(f1);

  // Unassigned PeptideIdentification with higher_score_better=false (same as run-level)
  PeptideIdentification pep_id2;
  pep_id2.setIdentifier("run_hsb_test");
  pep_id2.setScoreType("expect");
  pep_id2.setHigherScoreBetter(false);
  pep_id2.setRT(200.0);
  pep_id2.setMZ(600.0);

  PeptideHit hit2;
  hit2.setSequence(AASequence::fromString("ACDEFGHIK"));
  hit2.setScore(0.05);
  hit2.setCharge(3);
  hit2.setRank(1);
  hit2.setMetaValue("target_decoy", "target");
  pep_id2.insertHit(hit2);
  fm.setUnassignedPeptideIdentifications({pep_id2});

  // --- Export and import ---
  std::string tmp_dir;
  NEW_TMP_FILE(tmp_dir)
  tmp_dir += ".fmd";

  TEST_EQUAL(FeatureMapArrowIO::exportToParquet(fm, tmp_dir), true)

  FeatureMap imported;
  TEST_EQUAL(FeatureMapArrowIO::importFromParquet(tmp_dir, imported), true)

  // Verify per-PSM higher_score_better: pep_id1 should be true (not run-level false)
  TEST_EQUAL(imported[0].getPeptideIdentifications().size(), 1)
  const PeptideIdentification& out_pid1 = imported[0].getPeptideIdentifications()[0];
  TEST_EQUAL(out_pid1.isHigherScoreBetter(), true) // per-PSM value, NOT run-level false
  TEST_EQUAL(out_pid1.getScoreType(), "xcorr")

  // Verify unassigned PSM: higher_score_better should be false
  TEST_EQUAL(imported.getUnassignedPeptideIdentifications().size(), 1)
  const PeptideIdentification& out_pid2 = imported.getUnassignedPeptideIdentifications()[0];
  TEST_EQUAL(out_pid2.isHigherScoreBetter(), false)
  TEST_EQUAL(out_pid2.getScoreType(), "expect")

  // Verify scan and reference_file_name survive round-trip via metavalue
  // (scan is derived from spectrum_reference on export, reference_file_name from ProteinIdentification)
  const PeptideHit& out_hit1 = out_pid1.getHits()[0];
  TEST_EQUAL(out_hit1.metaValueExists("scan"), true)
  TEST_EQUAL(static_cast<int>(out_hit1.getMetaValue("scan")), 99)
}
END_SECTION

/////////////////////////////////////////////////////////////
// FeatureMap-level MetaValue round-trip (spectra_data + scalars)
/////////////////////////////////////////////////////////////

START_SECTION(exportToParquet / importFromParquet - FeatureMap - level MetaValue round - trip)
{
  FeatureMap fm;

  // setPrimaryMSRunPath stores its argument into the `spectra_data` meta-value
  // on the FeatureMap itself (FeatureMap.cpp:415). This is the path that
  // pre-fix FeatureMapArrowIO dropped on store.
  fm.setPrimaryMSRunPath(StringList {"run_A.mzML", "run_B.mzML"});

  // Plus a few scalar / typed meta-values to exercise the typed deserializer.
  fm.setMetaValue("custom_int", 1234);
  fm.setMetaValue("custom_double", 2.71828);
  fm.setMetaValue("custom_string", std::string("free-form text"));
  fm.setMetaValue("custom_int_list", DataValue(IntList {10, 20, 30}));
  fm.setMetaValue("custom_double_list", DataValue(DoubleList {1.5, 2.5}));

  // Minimal ProteinIdentification so exportToParquet succeeds.
  ProteinIdentification prot_id;
  prot_id.setIdentifier("run_mv_test");
  fm.setProteinIdentifications({prot_id});

  // One feature so the map is non-empty.
  Feature f;
  f.setRT(50.0);
  f.setMZ(400.0);
  f.setIntensity(500.0f);
  f.setCharge(1);
  f.setUniqueId(101);
  fm.push_back(f);

  std::string tmp_dir;
  NEW_TMP_FILE(tmp_dir)
  tmp_dir += ".fmd";

  TEST_EQUAL(FeatureMapArrowIO::exportToParquet(fm, tmp_dir), true)

  FeatureMap imported;
  TEST_EQUAL(FeatureMapArrowIO::importFromParquet(tmp_dir, imported), true)

  // Verify spectra_data round-trips and getPrimaryMSRunPath reads back the paths.
  StringList ms_runs;
  imported.getPrimaryMSRunPath(ms_runs);
  TEST_EQUAL(ms_runs.size(), 2)
  TEST_EQUAL(ms_runs[0], "run_A.mzML")
  TEST_EQUAL(ms_runs[1], "run_B.mzML")

  // Verify scalar meta-values restore with the correct type.
  TEST_EQUAL(imported.metaValueExists("custom_int"), true)
  TEST_EQUAL(static_cast<int>(imported.getMetaValue("custom_int")), 1234)
  TEST_REAL_SIMILAR(static_cast<double>(imported.getMetaValue("custom_double")), 2.71828)
  TEST_EQUAL(StringUtils::toStr(imported.getMetaValue("custom_string")), "free-form text")

  // List-typed meta-values restore as their original types.
  IntList out_il = imported.getMetaValue("custom_int_list");
  TEST_EQUAL(out_il.size(), 3)
  TEST_EQUAL(out_il[0], 10)
  TEST_EQUAL(out_il[2], 30)
  DoubleList out_dl = imported.getMetaValue("custom_double_list");
  TEST_EQUAL(out_dl.size(), 2)
  TEST_REAL_SIMILAR(out_dl[0], 1.5)
  TEST_REAL_SIMILAR(out_dl[1], 2.5)
}
END_SECTION

/////////////////////////////////////////////////////////////
// Fix #2b: exportToParquet rejects duplicate ProtID identifiers (XML-lane parity)
/////////////////////////////////////////////////////////////

START_SECTION(exportToParquet - duplicate ProteinIdentification identifiers throw Exception::InvalidValue)
{
  FeatureMap fm;

  ProteinIdentification p1;
  p1.setIdentifier("dup");
  ProteinIdentification p2;
  p2.setIdentifier("dup");
  fm.setProteinIdentifications({p1, p2});

  Feature f;
  f.setRT(50.0);
  f.setMZ(400.0);
  f.setIntensity(500.0f);
  f.setCharge(1);
  f.setUniqueId(101);
  fm.push_back(f);

  std::string tmp_dir;
  NEW_TMP_FILE(tmp_dir)
  tmp_dir += ".fmd";

  // The store-side check fires before any Arrow builder is allocated, so no
  // partial .featureparquet exists on disk after the throw.
  TEST_EXCEPTION(Exception::InvalidValue, FeatureMapArrowIO::exportToParquet(fm, tmp_dir))
}
END_SECTION


START_SECTION(([EXTRA] a failed write leaves no partial .parquet behind))
{
  FeatureMap fmap;
  Feature f;
  f.setRT(100.0);
  f.setMZ(500.0);
  f.setIntensity(1000.0f);
  fmap.push_back(f);

  const std::string dir = File::getTempDirectory() + "/" + File::getUniqueName() + "_fmio";
  TEST_TRUE(File::makeDir(dir))

  TEST_TRUE(FeatureMapArrowIO::exportToParquet(fmap, dir))
  TEST_TRUE(File::exists(dir + "/features.parquet"))
  File::remove(dir + "/features.parquet");

  // arrow::io::FileOutputStream::Open creates and truncates the file before the table is written,
  // so a failure afterwards leaves a fragment with no Parquet footer, which a reader reports as
  // corrupt. A row group size of 0 is refused by Parquet for a non-empty table, which reaches
  // that failure deterministically on every platform. This is the FIRST file of the collection,
  // so nothing else has been written yet - collection-level atomicity is a separate concern.
  ParquetWriteConfig no_row_group;
  no_row_group.row_group_size = 0;
  TEST_FALSE(FeatureMapArrowIO::exportToParquet(fmap, dir, no_row_group))
  TEST_FALSE(File::exists(dir + "/features.parquet"))

  File::removeDirRecursively(dir);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST

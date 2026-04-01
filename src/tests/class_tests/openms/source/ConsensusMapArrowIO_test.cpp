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
#include <OpenMS/FORMAT/ConsensusMapArrowIO.h>
///////////////////////////

#include <OpenMS/config.h>

#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/ConsensusFeature.h>
#include <OpenMS/KERNEL/FeatureHandle.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/METADATA/DataProcessing.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/ProteinHit.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>

#include <arrow/api.h>

using namespace OpenMS;
using namespace std;

START_TEST(ConsensusMapArrowIO, "$Id$")

/////////////////////////////////////////////////////////////
// Export tests
/////////////////////////////////////////////////////////////

START_SECTION(exportFeaturesToArrow - empty ConsensusMap)
{
  ConsensusMap cmap;
  auto table = ConsensusMapArrowIO::exportFeaturesToArrow(cmap);
  TEST_NOT_EQUAL(table, nullptr)
  TEST_EQUAL(table->num_rows(), 0)
  TEST_EQUAL(table->num_columns(), 9)
}
END_SECTION

START_SECTION(exportFeaturesToArrow - single feature with handles and metavalues)
{
  ConsensusMap cmap;

  // Set up column headers
  ConsensusMap::ColumnHeader ch0;
  ch0.filename = "file0.mzML";
  ch0.label = "light";
  ch0.size = 100;
  ch0.unique_id = 111;
  cmap.getColumnHeaders()[0] = ch0;

  ConsensusMap::ColumnHeader ch1;
  ch1.filename = "file1.mzML";
  ch1.label = "heavy";
  ch1.size = 200;
  ch1.unique_id = 222;
  cmap.getColumnHeaders()[1] = ch1;

  ConsensusFeature cf;
  cf.setRT(100.5);
  cf.setMZ(500.25);
  cf.setIntensity(1000.0f);
  cf.setCharge(2);
  cf.setQuality(0.95f);
  cf.setWidth(1.5f);
  cf.setUniqueId(12345);

  // Add 2 FeatureHandles
  FeatureHandle h0;
  h0.setMapIndex(0);
  h0.setUniqueId(100);
  h0.setRT(100.0);
  h0.setMZ(500.2);
  h0.setIntensity(800.0f);
  h0.setCharge(2);
  h0.setWidth(1.2f);
  cf.insert(h0);

  FeatureHandle h1;
  h1.setMapIndex(1);
  h1.setUniqueId(200);
  h1.setRT(101.0);
  h1.setMZ(500.3);
  h1.setIntensity(1200.0f);
  h1.setCharge(2);
  h1.setWidth(1.8f);
  cf.insert(h1);

  // Add metavalues
  cf.setMetaValue("my_int", 42);
  cf.setMetaValue("my_float", 3.14);
  cf.setMetaValue("my_string", String("hello"));

  cmap.push_back(cf);

  auto table = ConsensusMapArrowIO::exportFeaturesToArrow(cmap);
  TEST_NOT_EQUAL(table, nullptr)
  TEST_EQUAL(table->num_rows(), 1)
  TEST_EQUAL(table->num_columns(), 9)

  // Verify scalar columns
  auto col_unique_id = std::static_pointer_cast<arrow::Int64Array>(table->GetColumnByName("unique_id")->chunk(0));
  TEST_EQUAL(col_unique_id->Value(0), 12345)

  auto col_rt = std::static_pointer_cast<arrow::DoubleArray>(table->GetColumnByName("rt")->chunk(0));
  TEST_REAL_SIMILAR(col_rt->Value(0), 100.5)

  auto col_mz = std::static_pointer_cast<arrow::DoubleArray>(table->GetColumnByName("mz")->chunk(0));
  TEST_REAL_SIMILAR(col_mz->Value(0), 500.25)

  auto col_intensity = std::static_pointer_cast<arrow::FloatArray>(table->GetColumnByName("intensity")->chunk(0));
  TEST_REAL_SIMILAR(col_intensity->Value(0), 1000.0)

  auto col_charge = std::static_pointer_cast<arrow::Int32Array>(table->GetColumnByName("charge")->chunk(0));
  TEST_EQUAL(col_charge->Value(0), 2)

  auto col_quality = std::static_pointer_cast<arrow::FloatArray>(table->GetColumnByName("quality")->chunk(0));
  TEST_REAL_SIMILAR(col_quality->Value(0), 0.95)

  auto col_width = std::static_pointer_cast<arrow::FloatArray>(table->GetColumnByName("width")->chunk(0));
  TEST_EQUAL(col_width->IsNull(0), false)
  TEST_REAL_SIMILAR(col_width->Value(0), 1.5)

  // Verify handles column (list of structs)
  auto handles_chunked = table->GetColumnByName("handles");
  TEST_NOT_EQUAL(handles_chunked, nullptr)
  auto handles_list = std::static_pointer_cast<arrow::ListArray>(handles_chunked->chunk(0));
  TEST_EQUAL(handles_list->value_length(0), 2)

  auto handle_structs = std::static_pointer_cast<arrow::StructArray>(handles_list->value_slice(0));
  auto h_map_idx = std::static_pointer_cast<arrow::Int64Array>(handle_structs->field(0));
  auto h_uid = std::static_pointer_cast<arrow::Int64Array>(handle_structs->field(1));
  auto h_rt = std::static_pointer_cast<arrow::DoubleArray>(handle_structs->field(2));
  auto h_intensity = std::static_pointer_cast<arrow::FloatArray>(handle_structs->field(4));

  // Handles are stored sorted by (map_index, unique_id)
  TEST_EQUAL(h_map_idx->Value(0), 0)
  TEST_EQUAL(h_uid->Value(0), 100)
  TEST_REAL_SIMILAR(h_rt->Value(0), 100.0)
  TEST_REAL_SIMILAR(h_intensity->Value(0), 800.0)

  TEST_EQUAL(h_map_idx->Value(1), 1)
  TEST_EQUAL(h_uid->Value(1), 200)
  TEST_REAL_SIMILAR(h_rt->Value(1), 101.0)
  TEST_REAL_SIMILAR(h_intensity->Value(1), 1200.0)

  // Verify metavalues column
  auto mv_chunked = table->GetColumnByName("metavalues");
  TEST_NOT_EQUAL(mv_chunked, nullptr)
  auto mv_list = std::static_pointer_cast<arrow::ListArray>(mv_chunked->chunk(0));
  TEST_EQUAL(mv_list->value_length(0), 3)
}
END_SECTION

START_SECTION(exportPSMsToArrow - empty ConsensusMap)
{
  ConsensusMap cmap;
  // Need at least an empty ProteinIdentification for QPXFile
  ProteinIdentification prot_id;
  prot_id.setIdentifier("run_1");
  cmap.setProteinIdentifications({prot_id});

  auto table = ConsensusMapArrowIO::exportPSMsToArrow(cmap);
  TEST_NOT_EQUAL(table, nullptr)
  TEST_EQUAL(table->num_rows(), 0)
}
END_SECTION

START_SECTION(exportPSMsToArrow - feature and unassigned PSMs)
{
  ConsensusMap cmap;

  ProteinIdentification prot_id;
  prot_id.setIdentifier("run_1");
  prot_id.setSearchEngine("TestEngine");
  prot_id.setScoreType("TestScore");
  cmap.setProteinIdentifications({prot_id});

  // ConsensusFeature (uid=1000) with 1 PeptideIdentification
  ConsensusFeature cf;
  cf.setRT(100.0);
  cf.setMZ(500.0);
  cf.setIntensity(1000.0f);
  cf.setCharge(2);
  cf.setUniqueId(1000);

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
  cf.getPeptideIdentifications().push_back(pep_id1);

  cmap.push_back(cf);

  // 1 unassigned PeptideIdentification
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
  cmap.getUnassignedPeptideIdentifications().push_back(pep_id2);

  auto table = ConsensusMapArrowIO::exportPSMsToArrow(cmap);
  TEST_NOT_EQUAL(table, nullptr)
  TEST_EQUAL(table->num_rows(), 2)

  // Verify consensus_feature_unique_id column
  auto cf_id_chunked = table->GetColumnByName("consensus_feature_unique_id");
  TEST_NOT_EQUAL(cf_id_chunked, nullptr)
  auto cf_id_arr = std::static_pointer_cast<arrow::Int64Array>(cf_id_chunked->chunk(0));

  // Row 0: feature PSM -> consensus_feature_unique_id = 1000
  TEST_EQUAL(cf_id_arr->IsNull(0), false)
  TEST_EQUAL(cf_id_arr->Value(0), 1000)

  // Row 1: unassigned PSM -> consensus_feature_unique_id = null
  TEST_EQUAL(cf_id_arr->IsNull(1), true)

  // Verify sequence column
  auto seq_col = std::static_pointer_cast<arrow::StringArray>(table->GetColumnByName("sequence")->chunk(0));
  TEST_EQUAL(seq_col->GetString(0), "PEPTIDER")
  TEST_EQUAL(seq_col->GetString(1), "ACDEFGHIK")

  // Verify score column
  auto score_col = std::static_pointer_cast<arrow::DoubleArray>(table->GetColumnByName("score")->chunk(0));
  TEST_REAL_SIMILAR(score_col->Value(0), 0.95)
  TEST_REAL_SIMILAR(score_col->Value(1), 0.80)
}
END_SECTION

/////////////////////////////////////////////////////////////
// Feature round-trip tests
/////////////////////////////////////////////////////////////

START_SECTION(importFeaturesFromArrow - feature round-trip with handles and metavalues)
{
  ConsensusMap cmap_in;

  // Feature 1: with handles and metavalues
  ConsensusFeature cf1;
  cf1.setRT(100.5);
  cf1.setMZ(500.25);
  cf1.setIntensity(1000.0f);
  cf1.setCharge(2);
  cf1.setQuality(0.95f);
  cf1.setWidth(1.5f);
  cf1.setUniqueId(1001);
  cf1.setMetaValue("my_int", 42);
  cf1.setMetaValue("my_float", 3.14);
  cf1.setMetaValue("my_string", String("hello"));
  cf1.setMetaValue("test_int_list", DataValue(IntList{1, 2, 3}));
  cf1.setMetaValue("test_double_list", DataValue(DoubleList{1.5, 2.5}));
  cf1.setMetaValue("test_string_list", DataValue(StringList{"a", "b", "c"}));

  FeatureHandle h0;
  h0.setMapIndex(0);
  h0.setUniqueId(100);
  h0.setRT(100.0);
  h0.setMZ(500.2);
  h0.setIntensity(800.0f);
  h0.setCharge(2);
  h0.setWidth(1.2f);
  cf1.insert(h0);

  FeatureHandle h1;
  h1.setMapIndex(1);
  h1.setUniqueId(200);
  h1.setRT(101.0);
  h1.setMZ(500.3);
  h1.setIntensity(1200.0f);
  h1.setCharge(3);
  h1.setWidth(1.8f);
  cf1.insert(h1);

  cmap_in.push_back(cf1);

  // Feature 2: minimal, no handles, no metavalues
  ConsensusFeature cf2;
  cf2.setRT(200.0);
  cf2.setMZ(600.0);
  cf2.setIntensity(2000.0f);
  cf2.setCharge(3);
  cf2.setQuality(0.5f);
  cf2.setUniqueId(2001);
  cmap_in.push_back(cf2);

  // Export
  auto table = ConsensusMapArrowIO::exportFeaturesToArrow(cmap_in);
  TEST_NOT_EQUAL(table, nullptr)
  TEST_EQUAL(table->num_rows(), 2)

  // Import
  ConsensusMap cmap_out;
  bool ok = ConsensusMapArrowIO::importFeaturesFromArrow(table, cmap_out);
  TEST_EQUAL(ok, true)
  TEST_EQUAL(cmap_out.size(), 2)

  // Verify feature 1
  const auto& out1 = cmap_out[0];
  TEST_REAL_SIMILAR(out1.getRT(), 100.5)
  TEST_REAL_SIMILAR(out1.getMZ(), 500.25)
  TEST_REAL_SIMILAR(out1.getIntensity(), 1000.0)
  TEST_EQUAL(out1.getCharge(), 2)
  TEST_REAL_SIMILAR(out1.getQuality(), 0.95)
  TEST_REAL_SIMILAR(out1.getWidth(), 1.5)
  TEST_EQUAL(out1.getUniqueId(), 1001)

  // Verify handles
  const auto& handles = out1.getFeatures();
  TEST_EQUAL(handles.size(), 2)

  auto it = handles.begin();
  TEST_EQUAL(it->getMapIndex(), 0)
  TEST_EQUAL(it->getUniqueId(), 100)
  TEST_REAL_SIMILAR(it->getRT(), 100.0)
  TEST_REAL_SIMILAR(it->getMZ(), 500.2)
  TEST_REAL_SIMILAR(it->getIntensity(), 800.0)
  TEST_EQUAL(it->getCharge(), 2)
  TEST_REAL_SIMILAR(it->getWidth(), 1.2)

  ++it;
  TEST_EQUAL(it->getMapIndex(), 1)
  TEST_EQUAL(it->getUniqueId(), 200)
  TEST_REAL_SIMILAR(it->getRT(), 101.0)
  TEST_REAL_SIMILAR(it->getMZ(), 500.3)
  TEST_REAL_SIMILAR(it->getIntensity(), 1200.0)
  TEST_EQUAL(it->getCharge(), 3)
  TEST_REAL_SIMILAR(it->getWidth(), 1.8)

  // Verify metavalues
  TEST_EQUAL(int(out1.getMetaValue("my_int")), 42)
  TEST_REAL_SIMILAR(double(out1.getMetaValue("my_float")), 3.14)
  TEST_EQUAL(String(out1.getMetaValue("my_string")), "hello")

  // Check list metavalue types are preserved
  TEST_EQUAL(out1.getMetaValue("test_int_list").valueType(), DataValue::INT_LIST)
  TEST_EQUAL(out1.getMetaValue("test_int_list") == DataValue(IntList{1, 2, 3}), true)
  TEST_EQUAL(out1.getMetaValue("test_double_list").valueType(), DataValue::DOUBLE_LIST)
  TEST_EQUAL(out1.getMetaValue("test_double_list") == DataValue(DoubleList{1.5, 2.5}), true)
  TEST_EQUAL(out1.getMetaValue("test_string_list").valueType(), DataValue::STRING_LIST)
  TEST_EQUAL(out1.getMetaValue("test_string_list") == DataValue(StringList{"a", "b", "c"}), true)

  // Verify feature 2
  const auto& out2 = cmap_out[1];
  TEST_REAL_SIMILAR(out2.getRT(), 200.0)
  TEST_REAL_SIMILAR(out2.getMZ(), 600.0)
  TEST_REAL_SIMILAR(out2.getIntensity(), 2000.0)
  TEST_EQUAL(out2.getCharge(), 3)
  TEST_REAL_SIMILAR(out2.getQuality(), 0.5)
  TEST_EQUAL(out2.getUniqueId(), 2001)
  TEST_EQUAL(out2.getFeatures().size(), 0)
}
END_SECTION

/////////////////////////////////////////////////////////////
// PSM round-trip tests
/////////////////////////////////////////////////////////////

START_SECTION(importPSMsFromArrow - PSM round-trip)
{
  ConsensusMap cmap_in;

  ProteinIdentification prot_id;
  prot_id.setIdentifier("run_1");
  prot_id.setSearchEngine("TestEngine");
  prot_id.setScoreType("TestScore");
  cmap_in.setProteinIdentifications({prot_id});

  // ConsensusFeature (uid=1000) with 1 PeptideIdentification
  ConsensusFeature cf;
  cf.setRT(100.0);
  cf.setMZ(500.0);
  cf.setIntensity(1000.0f);
  cf.setCharge(2);
  cf.setUniqueId(1000);

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
  cf.getPeptideIdentifications().push_back(pep_id1);

  cmap_in.push_back(cf);

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
  cmap_in.getUnassignedPeptideIdentifications().push_back(pep_id2);

  // Export features and PSMs
  auto features_table = ConsensusMapArrowIO::exportFeaturesToArrow(cmap_in);
  TEST_NOT_EQUAL(features_table, nullptr)

  auto psm_table = ConsensusMapArrowIO::exportPSMsToArrow(cmap_in);
  TEST_NOT_EQUAL(psm_table, nullptr)
  TEST_EQUAL(psm_table->num_rows(), 2)

  // Import: first features, then PSMs
  ConsensusMap cmap_out;
  cmap_out.setProteinIdentifications(cmap_in.getProteinIdentifications());

  bool ok = ConsensusMapArrowIO::importFeaturesFromArrow(features_table, cmap_out);
  TEST_EQUAL(ok, true)
  TEST_EQUAL(cmap_out.size(), 1)

  ok = ConsensusMapArrowIO::importPSMsFromArrow(psm_table, cmap_out);
  TEST_EQUAL(ok, true)

  // Verify feature PSMs
  TEST_EQUAL(cmap_out[0].getPeptideIdentifications().size(), 1)
  const auto& out_pid1 = cmap_out[0].getPeptideIdentifications()[0];
  TEST_EQUAL(out_pid1.getHits().size(), 1)
  TEST_EQUAL(out_pid1.getHits()[0].getSequence().toString(), "PEPTIDER")
  TEST_REAL_SIMILAR(out_pid1.getHits()[0].getScore(), 0.95)
  TEST_EQUAL(out_pid1.getHits()[0].getCharge(), 2)

  // Verify unassigned PSMs
  TEST_EQUAL(cmap_out.getUnassignedPeptideIdentifications().size(), 1)
  const auto& out_pid2 = cmap_out.getUnassignedPeptideIdentifications()[0];
  TEST_EQUAL(out_pid2.getHits().size(), 1)
  TEST_EQUAL(out_pid2.getHits()[0].getSequence().toString(), "ACDEFGHIK")
  TEST_REAL_SIMILAR(out_pid2.getHits()[0].getScore(), 0.80)
  TEST_EQUAL(out_pid2.getHits()[0].getCharge(), 3)
}
END_SECTION

/////////////////////////////////////////////////////////////
// PSM per-PSM higher_score_better + scan/reference_file_name round-trip
/////////////////////////////////////////////////////////////

START_SECTION(exportToParquet / importFromParquet - per-PSM higher_score_better independent from ProteinIdentification)
{
  ConsensusMap cmap;

  // ProteinIdentification with higher_score_better=false (run-level)
  ProteinIdentification prot_id;
  prot_id.setIdentifier("run_hsb_test");
  prot_id.setSearchEngine("Comet");
  prot_id.setScoreType("expect");
  prot_id.setHigherScoreBetter(false);
  cmap.setProteinIdentifications({prot_id});

  // ConsensusFeature with PeptideIdentification where higher_score_better=true (differs from run-level)
  ConsensusFeature cf;
  cf.setRT(100.0);
  cf.setMZ(500.0);
  cf.setIntensity(1000.0f);
  cf.setCharge(2);
  cf.setUniqueId(7001);

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
  cf.setPeptideIdentifications({pep_id1});
  cmap.push_back(cf);

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
  cmap.getUnassignedPeptideIdentifications().push_back(pep_id2);

  // --- Export and import ---
  String tmp_dir;
  NEW_TMP_FILE(tmp_dir)
  tmp_dir += ".cmd";

  TEST_EQUAL(ConsensusMapArrowIO::exportToParquet(cmap, tmp_dir), true)

  ConsensusMap imported;
  TEST_EQUAL(ConsensusMapArrowIO::importFromParquet(tmp_dir, imported), true)

  // Verify per-PSM higher_score_better: pep_id1 should be true (not run-level false)
  TEST_EQUAL(imported[0].getPeptideIdentifications().size(), 1)
  const PeptideIdentification& out_pid1 = imported[0].getPeptideIdentifications()[0];
  TEST_EQUAL(out_pid1.isHigherScoreBetter(), true)  // per-PSM value, NOT run-level false
  TEST_EQUAL(out_pid1.getScoreType(), "xcorr")

  // Verify unassigned PSM: higher_score_better should be false
  TEST_EQUAL(imported.getUnassignedPeptideIdentifications().size(), 1)
  const PeptideIdentification& out_pid2 = imported.getUnassignedPeptideIdentifications()[0];
  TEST_EQUAL(out_pid2.isHigherScoreBetter(), false)
  TEST_EQUAL(out_pid2.getScoreType(), "expect")

  // Verify scan survives round-trip via metavalue
  const PeptideHit& out_hit1 = out_pid1.getHits()[0];
  TEST_EQUAL(out_hit1.metaValueExists("scan"), true)
  TEST_EQUAL(static_cast<int>(out_hit1.getMetaValue("scan")), 99)
}
END_SECTION

/////////////////////////////////////////////////////////////
// Full Parquet directory round-trip test
/////////////////////////////////////////////////////////////

START_SECTION(exportToParquet / importFromParquet - full round-trip)
{
  ConsensusMap cmap;

  // --- Column headers ---
  ConsensusMap::ColumnHeader ch0;
  ch0.filename = "sample1.mzML";
  ch0.label = "light";
  ch0.size = 100;
  ch0.unique_id = 111;
  cmap.getColumnHeaders()[0] = ch0;

  ConsensusMap::ColumnHeader ch1;
  ch1.filename = "sample2.mzML";
  ch1.label = "heavy";
  ch1.size = 200;
  ch1.unique_id = 222;
  cmap.getColumnHeaders()[1] = ch1;

  cmap.setExperimentType("label-free");

  // --- Set ConsensusMap-level UniqueId ---
  cmap.setUniqueId(55555);

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

  cmap.setProteinIdentifications({prot_id});

  // --- ConsensusFeature 1: with handles, metavalues, and PSM ---
  ConsensusFeature cf1;
  cf1.setRT(100.0);
  cf1.setMZ(500.0);
  cf1.setIntensity(1000.0f);
  cf1.setCharge(2);
  cf1.setQuality(0.9f);
  cf1.setWidth(1.5f);
  cf1.setUniqueId(1001);
  cf1.setMetaValue("custom_val", 42);

  FeatureHandle fh0;
  fh0.setMapIndex(0);
  fh0.setUniqueId(100);
  fh0.setRT(99.5);
  fh0.setMZ(500.0);
  fh0.setIntensity(800.0f);
  fh0.setCharge(2);
  cf1.insert(fh0);

  FeatureHandle fh1;
  fh1.setMapIndex(1);
  fh1.setUniqueId(200);
  fh1.setRT(100.5);
  fh1.setMZ(500.1);
  fh1.setIntensity(1200.0f);
  fh1.setCharge(2);
  cf1.insert(fh1);

  PeptideIdentification pep_id1;
  pep_id1.setScoreType("expect");
  pep_id1.setHigherScoreBetter(false);
  pep_id1.setIdentifier("run_full_1");
  PeptideHit pep_hit1;
  pep_hit1.setSequence(AASequence::fromString("PEPTIDER"));
  pep_hit1.setScore(0.001);
  pep_hit1.setCharge(2);
  pep_id1.insertHit(pep_hit1);
  cf1.setPeptideIdentifications({pep_id1});

  cmap.push_back(cf1);

  // --- ConsensusFeature 2: simple, no PSMs ---
  ConsensusFeature cf2;
  cf2.setRT(200.0);
  cf2.setMZ(600.0);
  cf2.setIntensity(2000.0f);
  cf2.setCharge(3);
  cf2.setUniqueId(2001);
  cmap.push_back(cf2);

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
  cmap.setUnassignedPeptideIdentifications({unassigned});

  // --- Export to temp directory ---
  String tmp_dir;
  NEW_TMP_FILE(tmp_dir)
  tmp_dir += ".cmd";

  TEST_EQUAL(ConsensusMapArrowIO::exportToParquet(cmap, tmp_dir), true)

  // --- Import back ---
  ConsensusMap imported;
  TEST_EQUAL(ConsensusMapArrowIO::importFromParquet(tmp_dir, imported), true)

  // --- Verify column headers ---
  TEST_EQUAL(imported.getColumnHeaders().size(), 2)
  TEST_EQUAL(imported.getColumnHeaders().at(0).filename, "sample1.mzML")
  TEST_EQUAL(imported.getColumnHeaders().at(0).label, "light")
  TEST_EQUAL(imported.getColumnHeaders().at(0).size, 100)
  TEST_EQUAL(imported.getColumnHeaders().at(0).unique_id, 111)
  TEST_EQUAL(imported.getColumnHeaders().at(1).filename, "sample2.mzML")
  TEST_EQUAL(imported.getColumnHeaders().at(1).label, "heavy")

  // --- Verify experiment type ---
  TEST_EQUAL(imported.getExperimentType(), "label-free")

  // --- Verify ConsensusMap-level UniqueId ---
  TEST_EQUAL(imported.getUniqueId(), 55555)

  // --- Verify features ---
  TEST_EQUAL(imported.size(), 2)

  // Feature 1
  TEST_REAL_SIMILAR(imported[0].getRT(), 100.0)
  TEST_REAL_SIMILAR(imported[0].getMZ(), 500.0)
  TEST_REAL_SIMILAR(imported[0].getIntensity(), 1000.0)
  TEST_EQUAL(imported[0].getCharge(), 2)
  TEST_REAL_SIMILAR(imported[0].getQuality(), 0.9)
  TEST_REAL_SIMILAR(imported[0].getWidth(), 1.5)
  TEST_EQUAL(imported[0].getUniqueId(), 1001)
  TEST_EQUAL(int(imported[0].getMetaValue("custom_val")), 42)

  // Feature 1 handles
  TEST_EQUAL(imported[0].getFeatures().size(), 2)
  auto hit = imported[0].getFeatures().begin();
  TEST_EQUAL(hit->getMapIndex(), 0)
  TEST_EQUAL(hit->getUniqueId(), 100)
  TEST_REAL_SIMILAR(hit->getIntensity(), 800.0)
  ++hit;
  TEST_EQUAL(hit->getMapIndex(), 1)
  TEST_EQUAL(hit->getUniqueId(), 200)
  TEST_REAL_SIMILAR(hit->getIntensity(), 1200.0)

  // Feature 2
  TEST_REAL_SIMILAR(imported[1].getRT(), 200.0)
  TEST_EQUAL(imported[1].getCharge(), 3)
  TEST_EQUAL(imported[1].getUniqueId(), 2001)
  TEST_EQUAL(imported[1].getFeatures().size(), 0)

  // --- Verify PSMs ---
  TEST_EQUAL(imported[0].getPeptideIdentifications().size(), 1)
  TEST_EQUAL(imported[0].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "PEPTIDER")
  TEST_REAL_SIMILAR(imported[0].getPeptideIdentifications()[0].getHits()[0].getScore(), 0.001)

  TEST_EQUAL(imported.getUnassignedPeptideIdentifications().size(), 1)
  TEST_EQUAL(imported.getUnassignedPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "ACDEFGHIK")

  // --- Verify protein identifications ---
  TEST_EQUAL(imported.getProteinIdentifications().size(), 1)
  TEST_EQUAL(imported.getProteinIdentifications()[0].getIdentifier(), "run_full_1")
  TEST_EQUAL(imported.getProteinIdentifications()[0].getSearchEngine(), "Comet")
  TEST_EQUAL(imported.getProteinIdentifications()[0].getHits().size(), 1)
  TEST_EQUAL(imported.getProteinIdentifications()[0].getHits()[0].getAccession(), "P12345")
  TEST_EQUAL(imported.getProteinIdentifications()[0].getProteinGroups().size(), 1)
}
END_SECTION

/////////////////////////////////////////////////////////////
// ConsensusMap metadata round-trip tests
/////////////////////////////////////////////////////////////

START_SECTION(exportToParquet / importFromParquet - metadata round-trip (DocumentIdentifier + DataProcessing + ColumnHeaders + UniqueId + MetaValues))
{
  ConsensusMap cmap;

  // --- Set DocumentIdentifier fields ---
  cmap.setIdentifier("test_cmap_lsid_456");
  cmap.setLoadedFilePath("/data/experiments/test_run.consensusXML");

  // --- Set ConsensusMap-level UniqueId ---
  cmap.setUniqueId(77777);

  // --- Set ConsensusMap-level MetaValues ---
  cmap.setMetaValue("analysis_type", String("differential"));
  cmap.setMetaValue("num_runs", 5);
  cmap.setMetaValue("threshold", 0.01);

  // --- Set column headers with metavalues ---
  ConsensusMap::ColumnHeader ch0;
  ch0.filename = "run_A.mzML";
  ch0.label = "sample_A";
  ch0.size = 500;
  ch0.unique_id = 1111;
  ch0.setMetaValue("injection_order", 1);
  ch0.setMetaValue("sample_group", String("control"));
  cmap.getColumnHeaders()[0] = ch0;

  ConsensusMap::ColumnHeader ch1;
  ch1.filename = "run_B.mzML";
  ch1.label = "sample_B";
  ch1.size = 600;
  ch1.unique_id = 2222;
  ch1.setMetaValue("injection_order", 2);
  ch1.setMetaValue("sample_group", String("treatment"));
  cmap.getColumnHeaders()[1] = ch1;

  ConsensusMap::ColumnHeader ch2;
  ch2.filename = "run_C.mzML";
  ch2.label = "sample_C";
  ch2.size = 700;
  ch2.unique_id = 3333;
  // ch2 intentionally has no metavalues
  cmap.getColumnHeaders()[2] = ch2;

  cmap.setExperimentType("labeled_MS1");

  // --- Set DataProcessing ---
  DataProcessing dp1;
  dp1.getSoftware().setName("FeatureFinderCentroided");
  dp1.getSoftware().setVersion("3.2.0");
  dp1.setCompletionTime(DateTime::fromString("2025-06-15T14:30:00", "yyyy-MM-ddThh:mm:ss"));
  dp1.getProcessingActions().insert(DataProcessing::PEAK_PICKING);
  dp1.getProcessingActions().insert(DataProcessing::FILTERING);
  dp1.setMetaValue("parameter_file", String("params.ini"));
  dp1.setMetaValue("num_threads", 8);

  DataProcessing dp2;
  dp2.getSoftware().setName("FeatureLinkerUnlabeledQT");
  dp2.getSoftware().setVersion("3.2.0");
  dp2.setCompletionTime(DateTime::fromString("2025-06-15T15:00:00", "yyyy-MM-ddThh:mm:ss"));
  dp2.getProcessingActions().insert(DataProcessing::FEATURE_GROUPING);
  dp2.setMetaValue("max_rt_shift", 300.5);

  cmap.setDataProcessing({dp1, dp2});

  // --- Add a simple consensus feature so the map is non-empty ---
  ConsensusFeature cf;
  cf.setRT(100.0);
  cf.setMZ(500.0);
  cf.setIntensity(1000.0f);
  cf.setCharge(2);
  cf.setUniqueId(9999);
  cmap.push_back(cf);

  // Need ProteinIdentification for export
  ProteinIdentification prot_id;
  prot_id.setIdentifier("run_meta_test");
  cmap.setProteinIdentifications({prot_id});

  // --- Export ---
  String tmp_dir;
  NEW_TMP_FILE(tmp_dir)
  tmp_dir += ".cmd";

  TEST_EQUAL(ConsensusMapArrowIO::exportToParquet(cmap, tmp_dir), true)

  // --- Import ---
  ConsensusMap imported;
  TEST_EQUAL(ConsensusMapArrowIO::importFromParquet(tmp_dir, imported), true)

  // --- Verify DocumentIdentifier ---
  TEST_EQUAL(imported.getIdentifier(), "test_cmap_lsid_456")
  TEST_EQUAL(imported.getLoadedFilePath(), "/data/experiments/test_run.consensusXML")

  // --- Verify ConsensusMap-level UniqueId ---
  TEST_EQUAL(imported.getUniqueId(), 77777)

  // --- Verify ConsensusMap-level MetaValues ---
  TEST_EQUAL(String(imported.getMetaValue("analysis_type")), "differential")
  TEST_EQUAL(int(imported.getMetaValue("num_runs")), 5)
  TEST_REAL_SIMILAR(double(imported.getMetaValue("threshold")), 0.01)

  // --- Verify experiment type ---
  TEST_EQUAL(imported.getExperimentType(), "labeled_MS1")

  // --- Verify column headers ---
  TEST_EQUAL(imported.getColumnHeaders().size(), 3)
  TEST_EQUAL(imported.getColumnHeaders().at(0).filename, "run_A.mzML")
  TEST_EQUAL(imported.getColumnHeaders().at(0).label, "sample_A")
  TEST_EQUAL(imported.getColumnHeaders().at(0).size, 500)
  TEST_EQUAL(imported.getColumnHeaders().at(0).unique_id, 1111)
  TEST_EQUAL(imported.getColumnHeaders().at(1).filename, "run_B.mzML")
  TEST_EQUAL(imported.getColumnHeaders().at(1).label, "sample_B")
  TEST_EQUAL(imported.getColumnHeaders().at(2).filename, "run_C.mzML")
  TEST_EQUAL(imported.getColumnHeaders().at(2).label, "sample_C")

  // --- Verify column header metavalues ---
  TEST_EQUAL(int(imported.getColumnHeaders().at(0).getMetaValue("injection_order")), 1)
  TEST_EQUAL(String(imported.getColumnHeaders().at(0).getMetaValue("sample_group")), "control")
  TEST_EQUAL(int(imported.getColumnHeaders().at(1).getMetaValue("injection_order")), 2)
  TEST_EQUAL(String(imported.getColumnHeaders().at(1).getMetaValue("sample_group")), "treatment")

  // --- Verify DataProcessing ---
  TEST_EQUAL(imported.getDataProcessing().size(), 2)

  const auto& out_dp1 = imported.getDataProcessing()[0];
  TEST_EQUAL(out_dp1.getSoftware().getName(), "FeatureFinderCentroided")
  TEST_EQUAL(out_dp1.getSoftware().getVersion(), "3.2.0")
  TEST_EQUAL(out_dp1.getCompletionTime().toString("yyyy-MM-ddThh:mm:ss"), "2025-06-15T14:30:00")
  TEST_EQUAL(out_dp1.getProcessingActions().count(DataProcessing::PEAK_PICKING), 1)
  TEST_EQUAL(out_dp1.getProcessingActions().count(DataProcessing::FILTERING), 1)
  TEST_EQUAL(String(out_dp1.getMetaValue("parameter_file")), "params.ini")
  TEST_EQUAL(int(out_dp1.getMetaValue("num_threads")), 8)

  const auto& out_dp2 = imported.getDataProcessing()[1];
  TEST_EQUAL(out_dp2.getSoftware().getName(), "FeatureLinkerUnlabeledQT")
  TEST_EQUAL(out_dp2.getProcessingActions().count(DataProcessing::FEATURE_GROUPING), 1)
  TEST_REAL_SIMILAR(double(out_dp2.getMetaValue("max_rt_shift")), 300.5)
}
END_SECTION

END_TEST

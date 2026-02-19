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
#include <OpenMS/FORMAT/ProteinIdentificationArrowExport.h>
///////////////////////////

#include <OpenMS/config.h>

#ifdef WITH_PARQUET

#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/ProteinHit.h>

#include <arrow/api.h>

using namespace OpenMS;
using namespace std;

START_TEST(ProteinIdentificationArrowExport, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(static std::shared_ptr<arrow::Table> exportProteinsToArrow(...))
{
  // Create test data: one ProteinIdentification with 2 ProteinHits
  vector<ProteinIdentification> protein_ids;

  ProteinIdentification prot_id;
  prot_id.setIdentifier("run_1");
  prot_id.setSearchEngine("Mascot");
  prot_id.setSearchEngineVersion("2.7");
  prot_id.setScoreType("Mascot");
  prot_id.setHigherScoreBetter(true);
  prot_id.setSignificanceThreshold(0.05);
  prot_id.setPrimaryMSRunPath(StringList{"sample_1.mzML"});

  DateTime dt;
  dt.set("2024-01-15 10:30:00");
  prot_id.setDateTime(dt);

  // Hit 1: target protein with coverage and description
  ProteinHit hit1;
  hit1.setAccession("sp|P12345|PROT1_HUMAN");
  hit1.setScore(150.5);
  hit1.setRank(1);
  hit1.setCoverage(45.2);
  hit1.setSequence("MKWVTFISLLLLFSSAYS");
  hit1.setDescription("Serum albumin");
  hit1.setMetaValue("target_decoy", "target");
  hit1.setMetaValue("custom_score", 0.99);

  // Hit 2: decoy protein with unknown coverage (no description)
  ProteinHit hit2;
  hit2.setAccession("DECOY_sp|Q67890|PROT2_HUMAN");
  hit2.setScore(25.3);
  hit2.setRank(2);
  // coverage defaults to COVERAGE_UNKNOWN (-1)
  hit2.setMetaValue("target_decoy", "decoy");

  prot_id.insertHit(hit1);
  prot_id.insertHit(hit2);
  protein_ids.push_back(prot_id);

  // Export to Arrow table
  auto table = ProteinIdentificationArrowExport::exportProteinsToArrow(protein_ids);
  TEST_NOT_EQUAL(table, nullptr)

  // Verify number of rows (2 hits)
  TEST_EQUAL(table->num_rows(), 2)

  // Verify 19 columns
  TEST_EQUAL(table->num_columns(), 19)

  // Verify column names
  auto schema = table->schema();
  TEST_EQUAL(schema->field(0)->name(), "accession")
  TEST_EQUAL(schema->field(1)->name(), "score")
  TEST_EQUAL(schema->field(2)->name(), "score_type")
  TEST_EQUAL(schema->field(3)->name(), "higher_score_better")
  TEST_EQUAL(schema->field(4)->name(), "rank")
  TEST_EQUAL(schema->field(5)->name(), "coverage")
  TEST_EQUAL(schema->field(6)->name(), "sequence")
  TEST_EQUAL(schema->field(7)->name(), "description")
  TEST_EQUAL(schema->field(8)->name(), "is_decoy")
  TEST_EQUAL(schema->field(9)->name(), "run_identifier")
  TEST_EQUAL(schema->field(10)->name(), "reference_file_name")
  TEST_EQUAL(schema->field(11)->name(), "search_engine")
  TEST_EQUAL(schema->field(12)->name(), "search_engine_version")
  TEST_EQUAL(schema->field(13)->name(), "inference_engine")
  TEST_EQUAL(schema->field(14)->name(), "inference_engine_version")
  TEST_EQUAL(schema->field(15)->name(), "significance_threshold")
  TEST_EQUAL(schema->field(16)->name(), "date")
  TEST_EQUAL(schema->field(17)->name(), "modifications")
  TEST_EQUAL(schema->field(18)->name(), "metavalues")

  // Verify accession values
  auto acc_col = table->GetColumnByName("accession");
  auto acc_arr = std::static_pointer_cast<arrow::StringArray>(acc_col->chunk(0));
  TEST_EQUAL(acc_arr->GetString(0), "sp|P12345|PROT1_HUMAN")
  TEST_EQUAL(acc_arr->GetString(1), "DECOY_sp|Q67890|PROT2_HUMAN")

  // Verify score values
  auto score_col = table->GetColumnByName("score");
  auto score_arr = std::static_pointer_cast<arrow::DoubleArray>(score_col->chunk(0));
  TEST_REAL_SIMILAR(score_arr->Value(0), 150.5)
  TEST_REAL_SIMILAR(score_arr->Value(1), 25.3)

  // Verify score_type values
  auto st_col = table->GetColumnByName("score_type");
  auto st_arr = std::static_pointer_cast<arrow::StringArray>(st_col->chunk(0));
  TEST_EQUAL(st_arr->GetString(0), "Mascot")
  TEST_EQUAL(st_arr->GetString(1), "Mascot")

  // Verify higher_score_better values
  auto hsb_col = table->GetColumnByName("higher_score_better");
  auto hsb_arr = std::static_pointer_cast<arrow::BooleanArray>(hsb_col->chunk(0));
  TEST_EQUAL(hsb_arr->Value(0), true)
  TEST_EQUAL(hsb_arr->Value(1), true)

  // Verify rank values
  auto rank_col = table->GetColumnByName("rank");
  auto rank_arr = std::static_pointer_cast<arrow::Int32Array>(rank_col->chunk(0));
  TEST_EQUAL(rank_arr->Value(0), 1)
  TEST_EQUAL(rank_arr->Value(1), 2)

  // Verify coverage: first hit has valid coverage, second is null (COVERAGE_UNKNOWN)
  auto cov_col = table->GetColumnByName("coverage");
  auto cov_arr = std::static_pointer_cast<arrow::DoubleArray>(cov_col->chunk(0));
  TEST_EQUAL(cov_arr->IsNull(0), false)
  TEST_REAL_SIMILAR(cov_arr->Value(0), 45.2)
  TEST_EQUAL(cov_arr->IsNull(1), true)

  // Verify run_identifier
  auto rid_col = table->GetColumnByName("run_identifier");
  auto rid_arr = std::static_pointer_cast<arrow::StringArray>(rid_col->chunk(0));
  TEST_EQUAL(rid_arr->GetString(0), "run_1")
  TEST_EQUAL(rid_arr->GetString(1), "run_1")

  // Verify is_decoy: first is target (0), second is decoy (1)
  auto decoy_col = table->GetColumnByName("is_decoy");
  auto decoy_arr = std::static_pointer_cast<arrow::Int32Array>(decoy_col->chunk(0));
  TEST_EQUAL(decoy_arr->Value(0), 0)
  TEST_EQUAL(decoy_arr->Value(1), 1)

  // Verify description: first has value, second is null
  auto desc_col = table->GetColumnByName("description");
  auto desc_arr = std::static_pointer_cast<arrow::StringArray>(desc_col->chunk(0));
  TEST_EQUAL(desc_arr->IsNull(0), false)
  TEST_EQUAL(desc_arr->GetString(0), "Serum albumin")
  TEST_EQUAL(desc_arr->IsNull(1), true)

  // Verify sequence: first has value, second is null (empty)
  auto seq_col = table->GetColumnByName("sequence");
  auto seq_arr = std::static_pointer_cast<arrow::StringArray>(seq_col->chunk(0));
  TEST_EQUAL(seq_arr->IsNull(0), false)
  TEST_EQUAL(seq_arr->GetString(0), "MKWVTFISLLLLFSSAYS")
  TEST_EQUAL(seq_arr->IsNull(1), true)

  // Verify reference_file_name
  auto ref_col = table->GetColumnByName("reference_file_name");
  auto ref_arr = std::static_pointer_cast<arrow::StringArray>(ref_col->chunk(0));
  TEST_EQUAL(ref_arr->GetString(0), "sample_1.mzML")
  TEST_EQUAL(ref_arr->GetString(1), "sample_1.mzML")

  // Verify search_engine
  auto se_col = table->GetColumnByName("search_engine");
  auto se_arr = std::static_pointer_cast<arrow::StringArray>(se_col->chunk(0));
  TEST_EQUAL(se_arr->GetString(0), "Mascot")

  // Verify search_engine_version
  auto sev_col = table->GetColumnByName("search_engine_version");
  auto sev_arr = std::static_pointer_cast<arrow::StringArray>(sev_col->chunk(0));
  TEST_EQUAL(sev_arr->GetString(0), "2.7")

  // Verify inference_engine is null (not set)
  auto ie_col = table->GetColumnByName("inference_engine");
  auto ie_arr = std::static_pointer_cast<arrow::StringArray>(ie_col->chunk(0));
  TEST_EQUAL(ie_arr->IsNull(0), true)

  // Verify modifications: both should be null (no protein modifications set)
  auto mod_col = table->GetColumnByName("modifications");
  auto mod_arr = std::static_pointer_cast<arrow::ListArray>(mod_col->chunk(0));
  TEST_EQUAL(mod_arr->IsNull(0), true)
  TEST_EQUAL(mod_arr->IsNull(1), true)

  // Verify metavalues: first hit should have "custom_score" but not "Description" or "target_decoy"
  auto mv_col = table->GetColumnByName("metavalues");
  auto mv_arr = std::static_pointer_cast<arrow::ListArray>(mv_col->chunk(0));
  // First hit has custom_score metavalue (Description and target_decoy are excluded)
  TEST_EQUAL(mv_arr->value_length(0) >= 1, true)
  // Second hit has no extra metavalues (only target_decoy which is excluded)
  TEST_EQUAL(mv_arr->value_length(1), 0)

  // Verify data types
  TEST_EQUAL(schema->field(0)->type()->id(), arrow::Type::STRING)    // accession
  TEST_EQUAL(schema->field(1)->type()->id(), arrow::Type::DOUBLE)    // score
  TEST_EQUAL(schema->field(3)->type()->id(), arrow::Type::BOOL)      // higher_score_better
  TEST_EQUAL(schema->field(4)->type()->id(), arrow::Type::INT32)     // rank
  TEST_EQUAL(schema->field(5)->type()->id(), arrow::Type::DOUBLE)    // coverage
  TEST_EQUAL(schema->field(8)->type()->id(), arrow::Type::INT32)     // is_decoy
  TEST_EQUAL(schema->field(17)->type()->id(), arrow::Type::LIST)     // modifications
  TEST_EQUAL(schema->field(18)->type()->id(), arrow::Type::LIST)     // metavalues
}
END_SECTION

START_SECTION(Test empty protein identifications)
{
  // Test with empty input
  vector<ProteinIdentification> empty_ids;
  auto table = ProteinIdentificationArrowExport::exportProteinsToArrow(empty_ids);
  TEST_NOT_EQUAL(table, nullptr)
  TEST_EQUAL(table->num_rows(), 0)
  TEST_EQUAL(table->num_columns(), 19)
}
END_SECTION

START_SECTION(Test is_decoy null when target_decoy not set)
{
  // Test that is_decoy is null when target_decoy metavalue is absent
  vector<ProteinIdentification> protein_ids;
  ProteinIdentification prot_id;
  prot_id.setIdentifier("run_2");
  prot_id.setSearchEngine("TestEngine");
  prot_id.setScoreType("TestScore");

  ProteinHit hit;
  hit.setAccession("PROT_NO_TD");
  hit.setScore(10.0);
  // no target_decoy metavalue set

  prot_id.insertHit(hit);
  protein_ids.push_back(prot_id);

  auto table = ProteinIdentificationArrowExport::exportProteinsToArrow(protein_ids);
  TEST_NOT_EQUAL(table, nullptr)

  auto decoy_col = table->GetColumnByName("is_decoy");
  auto decoy_arr = std::static_pointer_cast<arrow::Int32Array>(decoy_col->chunk(0));
  TEST_EQUAL(decoy_arr->IsNull(0), true)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST

#else // WITH_PARQUET

START_TEST(ProteinIdentificationArrowExport, "$Id$")
END_TEST

#endif // WITH_PARQUET

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
#include <OpenMS/FORMAT/ProteinIdentificationArrowIO.h>
///////////////////////////

#include <OpenMS/config.h>

#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/ProteinHit.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>

#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/reader.h>

using namespace OpenMS;
using namespace std;

START_TEST(ProteinIdentificationArrowIO, "$Id$")

/////////////////////////////////////////////////////////////
// Export tests
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
  auto table = ProteinIdentificationArrowIO::exportProteinsToArrow(protein_ids);
  TEST_NOT_EQUAL(table, nullptr)

  // Verify number of rows (2 hits)
  TEST_EQUAL(table->num_rows(), 2)

  // Verify 10 columns
  TEST_EQUAL(table->num_columns(), 10)

  // Verify column names
  auto schema = table->schema();
  TEST_EQUAL(schema->field(0)->name(), "accession")
  TEST_EQUAL(schema->field(1)->name(), "score")
  TEST_EQUAL(schema->field(2)->name(), "rank")
  TEST_EQUAL(schema->field(3)->name(), "coverage")
  TEST_EQUAL(schema->field(4)->name(), "sequence")
  TEST_EQUAL(schema->field(5)->name(), "description")
  TEST_EQUAL(schema->field(6)->name(), "is_decoy")
  TEST_EQUAL(schema->field(7)->name(), "run_identifier")
  TEST_EQUAL(schema->field(8)->name(), "modifications")
  TEST_EQUAL(schema->field(9)->name(), "metavalues")

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

  // Verify is_decoy: first is target (false), second is decoy (true)
  auto decoy_col = table->GetColumnByName("is_decoy");
  auto decoy_arr = std::static_pointer_cast<arrow::BooleanArray>(decoy_col->chunk(0));
  TEST_EQUAL(decoy_arr->Value(0), false)
  TEST_EQUAL(decoy_arr->Value(1), true)

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
  TEST_EQUAL(schema->field(2)->type()->id(), arrow::Type::INT32)     // rank
  TEST_EQUAL(schema->field(3)->type()->id(), arrow::Type::DOUBLE)    // coverage
  TEST_EQUAL(schema->field(6)->type()->id(), arrow::Type::BOOL)      // is_decoy
  TEST_EQUAL(schema->field(8)->type()->id(), arrow::Type::LIST)      // modifications
  TEST_EQUAL(schema->field(9)->type()->id(), arrow::Type::LIST)      // metavalues
}
END_SECTION

START_SECTION(Test empty protein identifications)
{
  // Test with empty input
  vector<ProteinIdentification> empty_ids;
  auto table = ProteinIdentificationArrowIO::exportProteinsToArrow(empty_ids);
  TEST_NOT_EQUAL(table, nullptr)
  TEST_EQUAL(table->num_rows(), 0)
  TEST_EQUAL(table->num_columns(), 10)
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

  auto table = ProteinIdentificationArrowIO::exportProteinsToArrow(protein_ids);
  TEST_NOT_EQUAL(table, nullptr)

  auto decoy_col = table->GetColumnByName("is_decoy");
  auto decoy_arr = std::static_pointer_cast<arrow::BooleanArray>(decoy_col->chunk(0));
  TEST_EQUAL(decoy_arr->IsNull(0), true)
}
END_SECTION

START_SECTION(exportProteinGroupsToArrow())
{
  // Create test data: one ProteinIdentification with one protein_group and one indistinguishable group
  vector<ProteinIdentification> protein_ids;

  ProteinIdentification prot_id;
  prot_id.setIdentifier("run_groups_1");

  // Protein group: 2 accessions, probability 0.95
  ProteinIdentification::ProteinGroup pg;
  pg.probability = 0.95;
  pg.accessions = {"sp|P12345|PROT1_HUMAN", "sp|P12346|PROT2_HUMAN"};
  prot_id.insertProteinGroup(pg);

  // Indistinguishable group: 2 accessions, probability 0.0
  ProteinIdentification::ProteinGroup ig;
  ig.probability = 0.0;
  ig.accessions = {"sp|Q11111|PROTA_HUMAN", "sp|Q22222|PROTB_HUMAN"};
  prot_id.insertIndistinguishableProteins(ig);

  protein_ids.push_back(prot_id);

  // Export to Arrow table
  auto table = ProteinIdentificationArrowIO::exportProteinGroupsToArrow(protein_ids);
  TEST_NOT_EQUAL(table, nullptr)

  // Verify number of rows (1 protein_group + 1 indistinguishable = 2)
  TEST_EQUAL(table->num_rows(), 2)

  // Verify 8 columns
  TEST_EQUAL(table->num_columns(), 8)

  // Verify column names
  auto schema = table->schema();
  TEST_EQUAL(schema->field(0)->name(), "group_type")
  TEST_EQUAL(schema->field(1)->name(), "probability")
  TEST_EQUAL(schema->field(2)->name(), "accessions")
  TEST_EQUAL(schema->field(3)->name(), "run_identifier")
  TEST_EQUAL(schema->field(4)->name(), "group_index")
  TEST_EQUAL(schema->field(5)->name(), "float_data")
  TEST_EQUAL(schema->field(6)->name(), "string_data")
  TEST_EQUAL(schema->field(7)->name(), "integer_data")

  // Verify group_type values
  auto gt_col = table->GetColumnByName("group_type");
  auto gt_arr = std::static_pointer_cast<arrow::StringArray>(gt_col->chunk(0));
  TEST_EQUAL(gt_arr->GetString(0), "protein_group")
  TEST_EQUAL(gt_arr->GetString(1), "indistinguishable")

  // Verify probability values
  auto prob_col = table->GetColumnByName("probability");
  auto prob_arr = std::static_pointer_cast<arrow::DoubleArray>(prob_col->chunk(0));
  TEST_REAL_SIMILAR(prob_arr->Value(0), 0.95)
  TEST_REAL_SIMILAR(prob_arr->Value(1), 0.0)

  // Verify run_identifier
  auto rid_col = table->GetColumnByName("run_identifier");
  auto rid_arr = std::static_pointer_cast<arrow::StringArray>(rid_col->chunk(0));
  TEST_EQUAL(rid_arr->GetString(0), "run_groups_1")
  TEST_EQUAL(rid_arr->GetString(1), "run_groups_1")

  // Verify group_index (0 for both since they are independent counters)
  auto gi_col = table->GetColumnByName("group_index");
  auto gi_arr = std::static_pointer_cast<arrow::Int32Array>(gi_col->chunk(0));
  TEST_EQUAL(gi_arr->Value(0), 0)
  TEST_EQUAL(gi_arr->Value(1), 0)

  // Verify accessions list has correct values
  auto acc_col = table->GetColumnByName("accessions");
  auto acc_arr = std::static_pointer_cast<arrow::ListArray>(acc_col->chunk(0));
  // First group: 2 accessions
  TEST_EQUAL(acc_arr->value_length(0), 2)
  auto acc_values_0 = std::static_pointer_cast<arrow::StringArray>(acc_arr->value_slice(0));
  TEST_EQUAL(acc_values_0->GetString(0), "sp|P12345|PROT1_HUMAN")
  TEST_EQUAL(acc_values_0->GetString(1), "sp|P12346|PROT2_HUMAN")
  // Second group: 2 accessions
  TEST_EQUAL(acc_arr->value_length(1), 2)
  auto acc_values_1 = std::static_pointer_cast<arrow::StringArray>(acc_arr->value_slice(1));
  TEST_EQUAL(acc_values_1->GetString(0), "sp|Q11111|PROTA_HUMAN")
  TEST_EQUAL(acc_values_1->GetString(1), "sp|Q22222|PROTB_HUMAN")

  // Verify data array columns are null (no data arrays set)
  auto fd_col = table->GetColumnByName("float_data");
  auto fd_arr = std::static_pointer_cast<arrow::ListArray>(fd_col->chunk(0));
  TEST_EQUAL(fd_arr->IsNull(0), true)
  TEST_EQUAL(fd_arr->IsNull(1), true)

  auto sd_col = table->GetColumnByName("string_data");
  auto sd_arr = std::static_pointer_cast<arrow::ListArray>(sd_col->chunk(0));
  TEST_EQUAL(sd_arr->IsNull(0), true)
  TEST_EQUAL(sd_arr->IsNull(1), true)

  auto id_col = table->GetColumnByName("integer_data");
  auto id_arr = std::static_pointer_cast<arrow::ListArray>(id_col->chunk(0));
  TEST_EQUAL(id_arr->IsNull(0), true)
  TEST_EQUAL(id_arr->IsNull(1), true)

  // Verify data types
  TEST_EQUAL(schema->field(0)->type()->id(), arrow::Type::STRING)   // group_type
  TEST_EQUAL(schema->field(1)->type()->id(), arrow::Type::DOUBLE)   // probability
  TEST_EQUAL(schema->field(2)->type()->id(), arrow::Type::LIST)     // accessions
  TEST_EQUAL(schema->field(3)->type()->id(), arrow::Type::STRING)   // run_identifier
  TEST_EQUAL(schema->field(4)->type()->id(), arrow::Type::INT32)    // group_index
  TEST_EQUAL(schema->field(5)->type()->id(), arrow::Type::LIST)     // float_data
  TEST_EQUAL(schema->field(6)->type()->id(), arrow::Type::LIST)     // string_data
  TEST_EQUAL(schema->field(7)->type()->id(), arrow::Type::LIST)     // integer_data
}
END_SECTION

START_SECTION(exportSearchParamsToArrow())
{
  // Create test data: one ProteinIdentification with SearchParameters populated
  vector<ProteinIdentification> protein_ids;

  ProteinIdentification prot_id;
  prot_id.setIdentifier("run_search_1");
  prot_id.setSearchEngine("Comet");
  prot_id.setSearchEngineVersion("2023.01");
  prot_id.setScoreType("expect");
  prot_id.setHigherScoreBetter(false);
  prot_id.setSignificanceThreshold(0.01);
  prot_id.setPrimaryMSRunPath(StringList{"sample_1.mzML", "sample_2.mzML"});

  DateTime dt;
  dt.set("2024-06-15 14:30:00");
  prot_id.setDateTime(dt);

  // Configure SearchParameters
  ProteinIdentification::SearchParameters sp;
  sp.db = "uniprot_human";
  sp.db_version = "2024_01";
  sp.taxonomy = "Homo sapiens";
  sp.charges = "2:4";
  sp.mass_type = ProteinIdentification::PeakMassType::MONOISOTOPIC;
  sp.precursor_mass_tolerance = 10.0;
  sp.precursor_mass_tolerance_ppm = true;
  sp.fragment_mass_tolerance = 0.02;
  sp.fragment_mass_tolerance_ppm = false;
  sp.digestion_enzyme = *ProteaseDB::getInstance()->getEnzyme("Trypsin");
  sp.enzyme_term_specificity = EnzymaticDigestion::SPEC_FULL;
  sp.missed_cleavages = 2;
  sp.fixed_modifications = {"Carbamidomethyl (C)"};
  sp.variable_modifications = {"Oxidation (M)", "Acetyl (Protein N-term)"};

  prot_id.setSearchParameters(sp);
  protein_ids.push_back(prot_id);

  // Export to Arrow table
  auto table = ProteinIdentificationArrowIO::exportSearchParamsToArrow(protein_ids);
  TEST_NOT_EQUAL(table, nullptr)

  // Verify 1 row (one ProteinIdentification)
  TEST_EQUAL(table->num_rows(), 1)

  // Verify 26 columns
  TEST_EQUAL(table->num_columns(), 26)

  // Verify column names
  auto schema = table->schema();
  TEST_EQUAL(schema->field(0)->name(), "run_identifier")
  TEST_EQUAL(schema->field(1)->name(), "search_engine")
  TEST_EQUAL(schema->field(2)->name(), "search_engine_version")
  TEST_EQUAL(schema->field(3)->name(), "inference_engine")
  TEST_EQUAL(schema->field(4)->name(), "inference_engine_version")
  TEST_EQUAL(schema->field(5)->name(), "date")
  TEST_EQUAL(schema->field(6)->name(), "score_type")
  TEST_EQUAL(schema->field(7)->name(), "higher_score_better")
  TEST_EQUAL(schema->field(8)->name(), "significance_threshold")
  TEST_EQUAL(schema->field(9)->name(), "db")
  TEST_EQUAL(schema->field(10)->name(), "db_version")
  TEST_EQUAL(schema->field(11)->name(), "taxonomy")
  TEST_EQUAL(schema->field(12)->name(), "charges")
  TEST_EQUAL(schema->field(13)->name(), "mass_type")
  TEST_EQUAL(schema->field(14)->name(), "precursor_mass_tolerance")
  TEST_EQUAL(schema->field(15)->name(), "precursor_mass_tolerance_ppm")
  TEST_EQUAL(schema->field(16)->name(), "fragment_mass_tolerance")
  TEST_EQUAL(schema->field(17)->name(), "fragment_mass_tolerance_ppm")
  TEST_EQUAL(schema->field(18)->name(), "digestion_enzyme")
  TEST_EQUAL(schema->field(19)->name(), "enzyme_term_specificity")
  TEST_EQUAL(schema->field(20)->name(), "missed_cleavages")
  TEST_EQUAL(schema->field(21)->name(), "fixed_modifications")
  TEST_EQUAL(schema->field(22)->name(), "variable_modifications")
  TEST_EQUAL(schema->field(23)->name(), "primary_ms_run_paths")
  TEST_EQUAL(schema->field(24)->name(), "metavalues")
  TEST_EQUAL(schema->field(25)->name(), "sp_metavalues")

  // Verify run_identifier
  auto rid_col = table->GetColumnByName("run_identifier");
  auto rid_arr = std::static_pointer_cast<arrow::StringArray>(rid_col->chunk(0));
  TEST_EQUAL(rid_arr->GetString(0), "run_search_1")

  // Verify search_engine
  auto se_col = table->GetColumnByName("search_engine");
  auto se_arr = std::static_pointer_cast<arrow::StringArray>(se_col->chunk(0));
  TEST_EQUAL(se_arr->GetString(0), "Comet")

  // Verify search_engine_version
  auto sev_col = table->GetColumnByName("search_engine_version");
  auto sev_arr = std::static_pointer_cast<arrow::StringArray>(sev_col->chunk(0));
  TEST_EQUAL(sev_arr->GetString(0), "2023.01")

  // Verify inference_engine is null (not set)
  auto ie_col = table->GetColumnByName("inference_engine");
  auto ie_arr = std::static_pointer_cast<arrow::StringArray>(ie_col->chunk(0));
  TEST_EQUAL(ie_arr->IsNull(0), true)

  // Verify db
  auto db_col = table->GetColumnByName("db");
  auto db_arr = std::static_pointer_cast<arrow::StringArray>(db_col->chunk(0));
  TEST_EQUAL(db_arr->GetString(0), "uniprot_human")

  // Verify db_version
  auto dbv_col = table->GetColumnByName("db_version");
  auto dbv_arr = std::static_pointer_cast<arrow::StringArray>(dbv_col->chunk(0));
  TEST_EQUAL(dbv_arr->GetString(0), "2024_01")

  // Verify mass_type
  auto mt_col = table->GetColumnByName("mass_type");
  auto mt_arr = std::static_pointer_cast<arrow::StringArray>(mt_col->chunk(0));
  TEST_EQUAL(mt_arr->GetString(0), "MONOISOTOPIC")

  // Verify precursor_mass_tolerance
  auto pmt_col = table->GetColumnByName("precursor_mass_tolerance");
  auto pmt_arr = std::static_pointer_cast<arrow::DoubleArray>(pmt_col->chunk(0));
  TEST_REAL_SIMILAR(pmt_arr->Value(0), 10.0)

  // Verify precursor_mass_tolerance_ppm
  auto ppm_col = table->GetColumnByName("precursor_mass_tolerance_ppm");
  auto ppm_arr = std::static_pointer_cast<arrow::BooleanArray>(ppm_col->chunk(0));
  TEST_EQUAL(ppm_arr->Value(0), true)

  // Verify fragment_mass_tolerance
  auto fmt_col = table->GetColumnByName("fragment_mass_tolerance");
  auto fmt_arr = std::static_pointer_cast<arrow::DoubleArray>(fmt_col->chunk(0));
  TEST_REAL_SIMILAR(fmt_arr->Value(0), 0.02)

  // Verify fragment_mass_tolerance_ppm
  auto fppm_col = table->GetColumnByName("fragment_mass_tolerance_ppm");
  auto fppm_arr = std::static_pointer_cast<arrow::BooleanArray>(fppm_col->chunk(0));
  TEST_EQUAL(fppm_arr->Value(0), false)

  // Verify digestion_enzyme
  auto enz_col = table->GetColumnByName("digestion_enzyme");
  auto enz_arr = std::static_pointer_cast<arrow::StringArray>(enz_col->chunk(0));
  TEST_EQUAL(enz_arr->GetString(0), "Trypsin")

  // Verify enzyme_term_specificity
  auto spec_col = table->GetColumnByName("enzyme_term_specificity");
  auto spec_arr = std::static_pointer_cast<arrow::StringArray>(spec_col->chunk(0));
  TEST_EQUAL(spec_arr->GetString(0), "FULL")

  // Verify missed_cleavages
  auto mc_col = table->GetColumnByName("missed_cleavages");
  auto mc_arr = std::static_pointer_cast<arrow::Int32Array>(mc_col->chunk(0));
  TEST_EQUAL(mc_arr->Value(0), 2)

  // Verify higher_score_better
  auto hsb_col = table->GetColumnByName("higher_score_better");
  auto hsb_arr = std::static_pointer_cast<arrow::BooleanArray>(hsb_col->chunk(0));
  TEST_EQUAL(hsb_arr->Value(0), false)

  // Verify score_type
  auto st_col = table->GetColumnByName("score_type");
  auto st_arr = std::static_pointer_cast<arrow::StringArray>(st_col->chunk(0));
  TEST_EQUAL(st_arr->GetString(0), "expect")

  // Verify fixed_modifications list
  auto fm_col = table->GetColumnByName("fixed_modifications");
  auto fm_arr = std::static_pointer_cast<arrow::ListArray>(fm_col->chunk(0));
  TEST_EQUAL(fm_arr->value_length(0), 1)
  auto fm_values = std::static_pointer_cast<arrow::StringArray>(fm_arr->value_slice(0));
  TEST_EQUAL(fm_values->GetString(0), "Carbamidomethyl (C)")

  // Verify variable_modifications list
  auto vm_col = table->GetColumnByName("variable_modifications");
  auto vm_arr = std::static_pointer_cast<arrow::ListArray>(vm_col->chunk(0));
  TEST_EQUAL(vm_arr->value_length(0), 2)
  auto vm_values = std::static_pointer_cast<arrow::StringArray>(vm_arr->value_slice(0));
  TEST_EQUAL(vm_values->GetString(0), "Oxidation (M)")
  TEST_EQUAL(vm_values->GetString(1), "Acetyl (Protein N-term)")

  // Verify primary_ms_run_paths list
  auto msrp_col = table->GetColumnByName("primary_ms_run_paths");
  auto msrp_arr = std::static_pointer_cast<arrow::ListArray>(msrp_col->chunk(0));
  TEST_EQUAL(msrp_arr->value_length(0), 2)
  auto msrp_values = std::static_pointer_cast<arrow::StringArray>(msrp_arr->value_slice(0));
  TEST_EQUAL(msrp_values->GetString(0), "sample_1.mzML")
  TEST_EQUAL(msrp_values->GetString(1), "sample_2.mzML")

  // Verify data types
  TEST_EQUAL(schema->field(5)->type()->id(), arrow::Type::TIMESTAMP) // date
  TEST_EQUAL(schema->field(0)->type()->id(), arrow::Type::STRING)    // run_identifier
  TEST_EQUAL(schema->field(1)->type()->id(), arrow::Type::STRING)    // search_engine
  TEST_EQUAL(schema->field(7)->type()->id(), arrow::Type::BOOL)      // higher_score_better
  TEST_EQUAL(schema->field(8)->type()->id(), arrow::Type::DOUBLE)    // significance_threshold
  TEST_EQUAL(schema->field(13)->type()->id(), arrow::Type::STRING)   // mass_type
  TEST_EQUAL(schema->field(14)->type()->id(), arrow::Type::DOUBLE)   // precursor_mass_tolerance
  TEST_EQUAL(schema->field(15)->type()->id(), arrow::Type::BOOL)     // precursor_mass_tolerance_ppm
  TEST_EQUAL(schema->field(16)->type()->id(), arrow::Type::DOUBLE)   // fragment_mass_tolerance
  TEST_EQUAL(schema->field(17)->type()->id(), arrow::Type::BOOL)     // fragment_mass_tolerance_ppm
  TEST_EQUAL(schema->field(20)->type()->id(), arrow::Type::INT32)    // missed_cleavages
  TEST_EQUAL(schema->field(21)->type()->id(), arrow::Type::LIST)     // fixed_modifications
  TEST_EQUAL(schema->field(22)->type()->id(), arrow::Type::LIST)     // variable_modifications
  TEST_EQUAL(schema->field(23)->type()->id(), arrow::Type::LIST)     // primary_ms_run_paths
  TEST_EQUAL(schema->field(24)->type()->id(), arrow::Type::LIST)     // metavalues
  TEST_EQUAL(schema->field(25)->type()->id(), arrow::Type::LIST)     // sp_metavalues
}
END_SECTION

START_SECTION(exportProteinsToParquet())
{
  ProteinIdentification prot_id;
  prot_id.setIdentifier("run_parquet_1");
  prot_id.setSearchEngine("MS-GF+");
  prot_id.setScoreType("q-value");
  prot_id.setHigherScoreBetter(false);

  ProteinHit hit1;
  hit1.setAccession("P12345");
  hit1.setScore(0.001);
  prot_id.setHits({hit1});

  std::vector<ProteinIdentification> prot_ids = {prot_id};

  String tmp_file;
  NEW_TMP_FILE(tmp_file)
  tmp_file += ".parquet";

  TEST_EQUAL(ProteinIdentificationArrowIO::exportProteinsToParquet(prot_ids, tmp_file), true)

  // Read back and verify
  auto infile_result = arrow::io::ReadableFile::Open(std::string(tmp_file));
  TEST_EQUAL(infile_result.ok(), true)
  std::shared_ptr<arrow::io::ReadableFile> infile = *infile_result;

  auto reader_result = parquet::arrow::OpenFile(infile, arrow::default_memory_pool());
  TEST_EQUAL(reader_result.ok(), true)
  std::unique_ptr<parquet::arrow::FileReader> reader = std::move(reader_result.ValueOrDie());

  std::shared_ptr<arrow::Table> table;
  auto read_status = reader->ReadTable(&table);
  TEST_EQUAL(read_status.ok(), true)
  TEST_EQUAL(table->num_rows(), 1)
  TEST_EQUAL(table->num_columns(), 10)

  // Check file metadata
  auto metadata = table->schema()->metadata();
  TEST_NOT_EQUAL(metadata, nullptr)
  int idx = metadata->FindKey("file_type");
  TEST_EQUAL(idx >= 0, true)
  TEST_EQUAL(metadata->value(idx), "proteins")
}
END_SECTION

START_SECTION(exportProteinGroupsToParquet())
{
  ProteinIdentification prot_id;
  prot_id.setIdentifier("run_parquet_groups_1");

  ProteinIdentification::ProteinGroup pg;
  pg.probability = 0.95;
  pg.accessions = {"sp|P12345|PROT1_HUMAN", "sp|P12346|PROT2_HUMAN"};
  prot_id.insertProteinGroup(pg);

  ProteinIdentification::ProteinGroup ig;
  ig.probability = 0.0;
  ig.accessions = {"sp|Q11111|PROTA_HUMAN"};
  prot_id.insertIndistinguishableProteins(ig);

  std::vector<ProteinIdentification> prot_ids = {prot_id};

  String tmp_file;
  NEW_TMP_FILE(tmp_file)
  tmp_file += ".parquet";

  TEST_EQUAL(ProteinIdentificationArrowIO::exportProteinGroupsToParquet(prot_ids, tmp_file), true)

  // Read back and verify
  auto infile_result = arrow::io::ReadableFile::Open(std::string(tmp_file));
  TEST_EQUAL(infile_result.ok(), true)
  std::shared_ptr<arrow::io::ReadableFile> infile = *infile_result;

  auto reader_result = parquet::arrow::OpenFile(infile, arrow::default_memory_pool());
  TEST_EQUAL(reader_result.ok(), true)
  std::unique_ptr<parquet::arrow::FileReader> reader = std::move(reader_result.ValueOrDie());

  std::shared_ptr<arrow::Table> table;
  auto read_status = reader->ReadTable(&table);
  TEST_EQUAL(read_status.ok(), true)
  TEST_EQUAL(table->num_rows(), 2)
  TEST_EQUAL(table->num_columns(), 8)

  // Check file metadata
  auto metadata = table->schema()->metadata();
  TEST_NOT_EQUAL(metadata, nullptr)
  int idx = metadata->FindKey("file_type");
  TEST_EQUAL(idx >= 0, true)
  TEST_EQUAL(metadata->value(idx), "protein_groups")
}
END_SECTION

START_SECTION(exportSearchParamsToParquet())
{
  ProteinIdentification prot_id;
  prot_id.setIdentifier("run_parquet_search_1");
  prot_id.setSearchEngine("Comet");
  prot_id.setScoreType("expect");
  prot_id.setHigherScoreBetter(false);

  ProteinIdentification::SearchParameters sp;
  sp.db = "uniprot_human";
  sp.precursor_mass_tolerance = 10.0;
  sp.precursor_mass_tolerance_ppm = true;
  sp.fragment_mass_tolerance = 0.02;
  sp.fragment_mass_tolerance_ppm = false;
  sp.digestion_enzyme = *ProteaseDB::getInstance()->getEnzyme("Trypsin");
  sp.missed_cleavages = 2;
  sp.fixed_modifications = {"Carbamidomethyl (C)"};
  sp.variable_modifications = {"Oxidation (M)"};
  prot_id.setSearchParameters(sp);

  std::vector<ProteinIdentification> prot_ids = {prot_id};

  String tmp_file;
  NEW_TMP_FILE(tmp_file)
  tmp_file += ".parquet";

  TEST_EQUAL(ProteinIdentificationArrowIO::exportSearchParamsToParquet(prot_ids, tmp_file), true)

  // Read back and verify
  auto infile_result = arrow::io::ReadableFile::Open(std::string(tmp_file));
  TEST_EQUAL(infile_result.ok(), true)
  std::shared_ptr<arrow::io::ReadableFile> infile = *infile_result;

  auto reader_result = parquet::arrow::OpenFile(infile, arrow::default_memory_pool());
  TEST_EQUAL(reader_result.ok(), true)
  std::unique_ptr<parquet::arrow::FileReader> reader = std::move(reader_result.ValueOrDie());

  std::shared_ptr<arrow::Table> table;
  auto read_status = reader->ReadTable(&table);
  TEST_EQUAL(read_status.ok(), true)
  TEST_EQUAL(table->num_rows(), 1)
  TEST_EQUAL(table->num_columns(), 26)

  // Check file metadata
  auto metadata = table->schema()->metadata();
  TEST_NOT_EQUAL(metadata, nullptr)
  int idx = metadata->FindKey("file_type");
  TEST_EQUAL(idx >= 0, true)
  TEST_EQUAL(metadata->value(idx), "search_params")
}
END_SECTION

/////////////////////////////////////////////////////////////
// Import tests (round trips)
/////////////////////////////////////////////////////////////

START_SECTION(importSearchParamsFromArrow - round trip)
{
  // Create test data
  vector<ProteinIdentification> orig_ids;

  ProteinIdentification prot_id;
  prot_id.setIdentifier("run_search_1");
  prot_id.setSearchEngine("Comet");
  prot_id.setSearchEngineVersion("2023.01");
  prot_id.setScoreType("expect");
  prot_id.setHigherScoreBetter(false);
  prot_id.setSignificanceThreshold(0.01);
  prot_id.setPrimaryMSRunPath(StringList{"sample_1.mzML", "sample_2.mzML"});

  DateTime dt;
  dt.set("2024-06-15 14:30:00");
  prot_id.setDateTime(dt);

  ProteinIdentification::SearchParameters sp;
  sp.db = "uniprot_human";
  sp.db_version = "2024_01";
  sp.taxonomy = "Homo sapiens";
  sp.charges = "2:4";
  sp.mass_type = ProteinIdentification::PeakMassType::MONOISOTOPIC;
  sp.precursor_mass_tolerance = 10.0;
  sp.precursor_mass_tolerance_ppm = true;
  sp.fragment_mass_tolerance = 0.02;
  sp.fragment_mass_tolerance_ppm = false;
  sp.digestion_enzyme = *ProteaseDB::getInstance()->getEnzyme("Trypsin");
  sp.enzyme_term_specificity = EnzymaticDigestion::SPEC_FULL;
  sp.missed_cleavages = 2;
  sp.fixed_modifications = {"Carbamidomethyl (C)"};
  sp.variable_modifications = {"Oxidation (M)", "Acetyl (Protein N-term)"};
  sp.setMetaValue("sp_custom_param", "my_value");
  sp.setMetaValue("sp_custom_int", 123);
  prot_id.setSearchParameters(sp);

  // Set inference engine
  prot_id.setInferenceEngine("Percolator");
  prot_id.setInferenceEngineVersion("3.06");

  // Set a ProteinIdentification-level metavalue (should NOT end up on SearchParameters)
  prot_id.setMetaValue("prot_id_custom", "prot_level_value");

  orig_ids.push_back(prot_id);

  // Export to Arrow
  auto table = ProteinIdentificationArrowIO::exportSearchParamsToArrow(orig_ids);
  TEST_NOT_EQUAL(table, nullptr)

  // Import back
  vector<ProteinIdentification> imported_ids;
  TEST_EQUAL(ProteinIdentificationArrowIO::importSearchParamsFromArrow(table, imported_ids), true)
  TEST_EQUAL(imported_ids.size(), 1)

  const auto& imp = imported_ids[0];
  TEST_EQUAL(imp.getIdentifier(), "run_search_1")
  TEST_EQUAL(imp.getSearchEngine(), "Comet")
  TEST_EQUAL(imp.getSearchEngineVersion(), "2023.01")
  TEST_EQUAL(imp.getScoreType(), "expect")
  TEST_EQUAL(imp.isHigherScoreBetter(), false)
  TEST_REAL_SIMILAR(imp.getSignificanceThreshold(), 0.01)

  // Date
  TEST_EQUAL(imp.getDateTime().toString("yyyy-MM-ddThh:mm:ss"), "2024-06-15T14:30:00")

  // Primary MS run paths
  StringList ms_runs;
  imp.getPrimaryMSRunPath(ms_runs);
  TEST_EQUAL(ms_runs.size(), 2)
  TEST_EQUAL(ms_runs[0], "sample_1.mzML")
  TEST_EQUAL(ms_runs[1], "sample_2.mzML")

  // SearchParameters
  const auto& imp_sp = imp.getSearchParameters();
  TEST_EQUAL(imp_sp.db, "uniprot_human")
  TEST_EQUAL(imp_sp.db_version, "2024_01")
  TEST_EQUAL(imp_sp.taxonomy, "Homo sapiens")
  TEST_EQUAL(imp_sp.charges, "2:4")
  TEST_EQUAL(imp_sp.mass_type == ProteinIdentification::PeakMassType::MONOISOTOPIC, true)
  TEST_REAL_SIMILAR(imp_sp.precursor_mass_tolerance, 10.0)
  TEST_EQUAL(imp_sp.precursor_mass_tolerance_ppm, true)
  TEST_REAL_SIMILAR(imp_sp.fragment_mass_tolerance, 0.02)
  TEST_EQUAL(imp_sp.fragment_mass_tolerance_ppm, false)
  TEST_EQUAL(imp_sp.digestion_enzyme.getName(), "Trypsin")
  TEST_EQUAL(imp_sp.enzyme_term_specificity == EnzymaticDigestion::SPEC_FULL, true)
  TEST_EQUAL(imp_sp.missed_cleavages, 2)
  TEST_EQUAL(imp_sp.fixed_modifications.size(), 1)
  TEST_EQUAL(imp_sp.fixed_modifications[0], "Carbamidomethyl (C)")
  TEST_EQUAL(imp_sp.variable_modifications.size(), 2)
  TEST_EQUAL(imp_sp.variable_modifications[0], "Oxidation (M)")
  TEST_EQUAL(imp_sp.variable_modifications[1], "Acetyl (Protein N-term)")

  // Inference engine
  TEST_EQUAL(imp.getInferenceEngine(), "Percolator")
  TEST_EQUAL(imp.getInferenceEngineVersion(), "3.06")

  // SearchParameters metavalues should be on SearchParameters, not on ProteinIdentification
  TEST_EQUAL(String(imp_sp.getMetaValue("sp_custom_param")), "my_value")
  TEST_EQUAL(int(imp_sp.getMetaValue("sp_custom_int")), 123)

  // ProteinIdentification metavalue should be on ProteinIdentification
  TEST_EQUAL(String(imp.getMetaValue("prot_id_custom")), "prot_level_value")

  // Cross-check: sp_custom_param should NOT be on the ProteinIdentification object
  TEST_EQUAL(imp.metaValueExists("sp_custom_param"), false)
}
END_SECTION

START_SECTION(importProteinsFromArrow - round trip)
{
  // Create test data
  vector<ProteinIdentification> orig_ids;

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

  // Hit 2: decoy protein with unknown coverage
  ProteinHit hit2;
  hit2.setAccession("DECOY_sp|Q67890|PROT2_HUMAN");
  hit2.setScore(25.3);
  hit2.setRank(2);
  hit2.setMetaValue("target_decoy", "decoy");

  prot_id.insertHit(hit1);
  prot_id.insertHit(hit2);
  orig_ids.push_back(prot_id);

  // Export to Arrow
  auto table = ProteinIdentificationArrowIO::exportProteinsToArrow(orig_ids);
  TEST_NOT_EQUAL(table, nullptr)

  // Import back - first create shell via search params export/import (proper workflow)
  auto sp_table = ProteinIdentificationArrowIO::exportSearchParamsToArrow(orig_ids);
  vector<ProteinIdentification> imported_ids;
  TEST_EQUAL(ProteinIdentificationArrowIO::importSearchParamsFromArrow(sp_table, imported_ids), true)
  TEST_EQUAL(imported_ids.size(), 1)

  // Then import proteins into the shells
  TEST_EQUAL(ProteinIdentificationArrowIO::importProteinsFromArrow(table, imported_ids), true)
  TEST_EQUAL(imported_ids.size(), 1)

  const auto& imp = imported_ids[0];
  TEST_EQUAL(imp.getIdentifier(), "run_1")
  TEST_EQUAL(imp.getSearchEngine(), "Mascot")
  TEST_EQUAL(imp.getHits().size(), 2)

  // Hit 1
  const auto& imp_hit1 = imp.getHits()[0];
  TEST_EQUAL(imp_hit1.getAccession(), "sp|P12345|PROT1_HUMAN")
  TEST_REAL_SIMILAR(imp_hit1.getScore(), 150.5)
  TEST_EQUAL(imp_hit1.getRank(), 1)
  TEST_REAL_SIMILAR(imp_hit1.getCoverage(), 45.2)
  TEST_EQUAL(imp_hit1.getSequence(), "MKWVTFISLLLLFSSAYS")
  TEST_EQUAL(imp_hit1.getDescription(), "Serum albumin")
  TEST_EQUAL(String(imp_hit1.getMetaValue("target_decoy")), "target")
  TEST_REAL_SIMILAR(double(imp_hit1.getMetaValue("custom_score")), 0.99)

  // Hit 2
  const auto& imp_hit2 = imp.getHits()[1];
  TEST_EQUAL(imp_hit2.getAccession(), "DECOY_sp|Q67890|PROT2_HUMAN")
  TEST_REAL_SIMILAR(imp_hit2.getScore(), 25.3)
  TEST_EQUAL(imp_hit2.getRank(), 2)
  TEST_REAL_SIMILAR(imp_hit2.getCoverage(), ProteinHit::COVERAGE_UNKNOWN) // was null -> default
  TEST_EQUAL(String(imp_hit2.getMetaValue("target_decoy")), "decoy")
}
END_SECTION

START_SECTION(importProteinGroupsFromArrow - round trip)
{
  // Create test data
  vector<ProteinIdentification> orig_ids;

  ProteinIdentification prot_id;
  prot_id.setIdentifier("run_groups_1");

  // Protein group
  ProteinIdentification::ProteinGroup pg;
  pg.probability = 0.95;
  pg.accessions = {"sp|P12345|PROT1_HUMAN", "sp|P12346|PROT2_HUMAN"};
  prot_id.insertProteinGroup(pg);

  // Indistinguishable group
  ProteinIdentification::ProteinGroup ig;
  ig.probability = 0.0;
  ig.accessions = {"sp|Q11111|PROTA_HUMAN", "sp|Q22222|PROTB_HUMAN"};
  prot_id.insertIndistinguishableProteins(ig);

  orig_ids.push_back(prot_id);

  // Export to Arrow
  auto table = ProteinIdentificationArrowIO::exportProteinGroupsToArrow(orig_ids);
  TEST_NOT_EQUAL(table, nullptr)

  // Import back (create shell first to test run_id matching)
  vector<ProteinIdentification> imported_ids;
  ProteinIdentification shell;
  shell.setIdentifier("run_groups_1");
  imported_ids.push_back(shell);

  TEST_EQUAL(ProteinIdentificationArrowIO::importProteinGroupsFromArrow(table, imported_ids), true)
  TEST_EQUAL(imported_ids.size(), 1)

  // Verify protein groups
  const auto& imp = imported_ids[0];
  TEST_EQUAL(imp.getProteinGroups().size(), 1)
  TEST_REAL_SIMILAR(imp.getProteinGroups()[0].probability, 0.95)
  TEST_EQUAL(imp.getProteinGroups()[0].accessions.size(), 2)
  TEST_EQUAL(imp.getProteinGroups()[0].accessions[0], "sp|P12345|PROT1_HUMAN")
  TEST_EQUAL(imp.getProteinGroups()[0].accessions[1], "sp|P12346|PROT2_HUMAN")

  // Verify indistinguishable proteins
  TEST_EQUAL(imp.getIndistinguishableProteins().size(), 1)
  TEST_REAL_SIMILAR(imp.getIndistinguishableProteins()[0].probability, 0.0)
  TEST_EQUAL(imp.getIndistinguishableProteins()[0].accessions.size(), 2)
  TEST_EQUAL(imp.getIndistinguishableProteins()[0].accessions[0], "sp|Q11111|PROTA_HUMAN")
  TEST_EQUAL(imp.getIndistinguishableProteins()[0].accessions[1], "sp|Q22222|PROTB_HUMAN")
}
END_SECTION

START_SECTION(importSearchParamsFromParquet - file round trip)
{
  // Create test data
  ProteinIdentification prot_id;
  prot_id.setIdentifier("run_parquet_sp");
  prot_id.setSearchEngine("Comet");
  prot_id.setScoreType("expect");
  prot_id.setHigherScoreBetter(false);

  ProteinIdentification::SearchParameters sp;
  sp.db = "uniprot_human";
  sp.precursor_mass_tolerance = 10.0;
  sp.precursor_mass_tolerance_ppm = true;
  sp.fragment_mass_tolerance = 0.02;
  sp.fragment_mass_tolerance_ppm = false;
  sp.digestion_enzyme = *ProteaseDB::getInstance()->getEnzyme("Trypsin");
  sp.missed_cleavages = 2;
  sp.fixed_modifications = {"Carbamidomethyl (C)"};
  sp.variable_modifications = {"Oxidation (M)"};
  prot_id.setSearchParameters(sp);

  vector<ProteinIdentification> orig_ids = {prot_id};

  String tmp_file;
  NEW_TMP_FILE(tmp_file)
  tmp_file += ".parquet";

  TEST_EQUAL(ProteinIdentificationArrowIO::exportSearchParamsToParquet(orig_ids, tmp_file), true)

  vector<ProteinIdentification> imported_ids;
  TEST_EQUAL(ProteinIdentificationArrowIO::importSearchParamsFromParquet(tmp_file, imported_ids), true)
  TEST_EQUAL(imported_ids.size(), 1)
  TEST_EQUAL(imported_ids[0].getIdentifier(), "run_parquet_sp")
  TEST_EQUAL(imported_ids[0].getSearchParameters().db, "uniprot_human")
  TEST_EQUAL(imported_ids[0].getSearchParameters().digestion_enzyme.getName(), "Trypsin")
}
END_SECTION

START_SECTION(importProteinsFromParquet - file round trip)
{
  ProteinIdentification prot_id;
  prot_id.setIdentifier("run_parquet_prot");
  prot_id.setSearchEngine("MS-GF+");
  prot_id.setScoreType("q-value");
  prot_id.setHigherScoreBetter(false);

  ProteinHit hit1;
  hit1.setAccession("P12345");
  hit1.setScore(0.001);
  prot_id.setHits({hit1});

  vector<ProteinIdentification> orig_ids = {prot_id};

  String tmp_file;
  NEW_TMP_FILE(tmp_file)
  tmp_file += ".parquet";

  TEST_EQUAL(ProteinIdentificationArrowIO::exportProteinsToParquet(orig_ids, tmp_file), true)

  vector<ProteinIdentification> imported_ids;
  TEST_EQUAL(ProteinIdentificationArrowIO::importProteinsFromParquet(tmp_file, imported_ids), true)
  TEST_EQUAL(imported_ids.size(), 1)
  TEST_EQUAL(imported_ids[0].getHits().size(), 1)
  TEST_EQUAL(imported_ids[0].getHits()[0].getAccession(), "P12345")
  TEST_REAL_SIMILAR(imported_ids[0].getHits()[0].getScore(), 0.001)
}
END_SECTION

START_SECTION(importProteinGroupsFromParquet - file round trip)
{
  ProteinIdentification prot_id;
  prot_id.setIdentifier("run_parquet_pg");

  ProteinIdentification::ProteinGroup pg;
  pg.probability = 0.95;
  pg.accessions = {"P12345", "P12346"};
  prot_id.insertProteinGroup(pg);

  vector<ProteinIdentification> orig_ids = {prot_id};

  String tmp_file;
  NEW_TMP_FILE(tmp_file)
  tmp_file += ".parquet";

  TEST_EQUAL(ProteinIdentificationArrowIO::exportProteinGroupsToParquet(orig_ids, tmp_file), true)

  vector<ProteinIdentification> imported_ids;
  ProteinIdentification shell;
  shell.setIdentifier("run_parquet_pg");
  imported_ids.push_back(shell);

  TEST_EQUAL(ProteinIdentificationArrowIO::importProteinGroupsFromParquet(tmp_file, imported_ids), true)
  TEST_EQUAL(imported_ids.size(), 1)
  TEST_EQUAL(imported_ids[0].getProteinGroups().size(), 1)
  TEST_REAL_SIMILAR(imported_ids[0].getProteinGroups()[0].probability, 0.95)
}
END_SECTION

START_SECTION(importFromParquet - full combined round trip)
{
  // Create rich test data with 2 runs
  vector<ProteinIdentification> orig_ids;

  // Run 1
  {
    ProteinIdentification prot_id;
    prot_id.setIdentifier("run_combined_1");
    prot_id.setSearchEngine("Comet");
    prot_id.setSearchEngineVersion("2023.01");
    prot_id.setScoreType("expect");
    prot_id.setHigherScoreBetter(false);
    prot_id.setSignificanceThreshold(0.01);
    prot_id.setPrimaryMSRunPath(StringList{"sample_1.mzML"});

    DateTime dt;
    dt.set("2024-06-15 14:30:00");
    prot_id.setDateTime(dt);

    ProteinIdentification::SearchParameters sp;
    sp.db = "uniprot_human";
    sp.db_version = "2024_01";
    sp.mass_type = ProteinIdentification::PeakMassType::MONOISOTOPIC;
    sp.precursor_mass_tolerance = 10.0;
    sp.precursor_mass_tolerance_ppm = true;
    sp.fragment_mass_tolerance = 0.02;
    sp.fragment_mass_tolerance_ppm = false;
    sp.digestion_enzyme = *ProteaseDB::getInstance()->getEnzyme("Trypsin");
    sp.enzyme_term_specificity = EnzymaticDigestion::SPEC_FULL;
    sp.missed_cleavages = 2;
    sp.fixed_modifications = {"Carbamidomethyl (C)"};
    sp.variable_modifications = {"Oxidation (M)"};
    prot_id.setSearchParameters(sp);

    ProteinHit hit1;
    hit1.setAccession("sp|P12345|PROT1_HUMAN");
    hit1.setScore(0.001);
    hit1.setRank(1);
    hit1.setCoverage(45.2);
    hit1.setDescription("Serum albumin");
    hit1.setMetaValue("target_decoy", "target");
    prot_id.insertHit(hit1);

    ProteinHit hit2;
    hit2.setAccession("sp|Q67890|PROT2_HUMAN");
    hit2.setScore(0.005);
    hit2.setRank(2);
    hit2.setMetaValue("target_decoy", "target");
    prot_id.insertHit(hit2);

    ProteinIdentification::ProteinGroup pg;
    pg.probability = 0.99;
    pg.accessions = {"sp|P12345|PROT1_HUMAN", "sp|Q67890|PROT2_HUMAN"};
    prot_id.insertProteinGroup(pg);

    ProteinIdentification::ProteinGroup ig;
    ig.probability = 0.0;
    ig.accessions = {"sp|P12345|PROT1_HUMAN"};
    prot_id.insertIndistinguishableProteins(ig);

    orig_ids.push_back(prot_id);
  }

  // Run 2
  {
    ProteinIdentification prot_id;
    prot_id.setIdentifier("run_combined_2");
    prot_id.setSearchEngine("MS-GF+");
    prot_id.setScoreType("SpecEValue");
    prot_id.setHigherScoreBetter(false);

    ProteinIdentification::SearchParameters sp;
    sp.db = "uniprot_mouse";
    sp.mass_type = ProteinIdentification::PeakMassType::MONOISOTOPIC;
    sp.precursor_mass_tolerance = 20.0;
    sp.precursor_mass_tolerance_ppm = true;
    sp.fragment_mass_tolerance = 0.5;
    sp.fragment_mass_tolerance_ppm = false;
    sp.digestion_enzyme = *ProteaseDB::getInstance()->getEnzyme("Trypsin/P");
    sp.enzyme_term_specificity = EnzymaticDigestion::SPEC_SEMI;
    sp.missed_cleavages = 1;
    prot_id.setSearchParameters(sp);

    ProteinHit hit;
    hit.setAccession("MOUSE_PROT1");
    hit.setScore(1e-10);
    prot_id.insertHit(hit);

    orig_ids.push_back(prot_id);
  }

  // Export all 3 tables
  String proteins_file, groups_file, params_file;
  NEW_TMP_FILE(proteins_file)
  proteins_file += ".parquet";
  NEW_TMP_FILE(groups_file)
  groups_file += ".parquet";
  NEW_TMP_FILE(params_file)
  params_file += ".parquet";

  TEST_EQUAL(ProteinIdentificationArrowIO::exportProteinsToParquet(orig_ids, proteins_file), true)
  TEST_EQUAL(ProteinIdentificationArrowIO::exportProteinGroupsToParquet(orig_ids, groups_file), true)
  TEST_EQUAL(ProteinIdentificationArrowIO::exportSearchParamsToParquet(orig_ids, params_file), true)

  // Import combined
  vector<ProteinIdentification> imported_ids;
  TEST_EQUAL(ProteinIdentificationArrowIO::importFromParquet(
    proteins_file, groups_file, params_file, imported_ids), true)

  // Should have 2 runs
  TEST_EQUAL(imported_ids.size(), 2)

  // Find runs by identifier (order may differ)
  const ProteinIdentification* imp_run1 = nullptr;
  const ProteinIdentification* imp_run2 = nullptr;
  for (const auto& pid : imported_ids)
  {
    if (pid.getIdentifier() == "run_combined_1") imp_run1 = &pid;
    else if (pid.getIdentifier() == "run_combined_2") imp_run2 = &pid;
  }
  TEST_NOT_EQUAL(imp_run1, nullptr)
  TEST_NOT_EQUAL(imp_run2, nullptr)

  // Verify run 1
  TEST_EQUAL(imp_run1->getSearchEngine(), "Comet")
  TEST_EQUAL(imp_run1->getSearchEngineVersion(), "2023.01")
  TEST_EQUAL(imp_run1->getScoreType(), "expect")
  TEST_EQUAL(imp_run1->isHigherScoreBetter(), false)
  TEST_EQUAL(imp_run1->getHits().size(), 2)
  TEST_EQUAL(imp_run1->getProteinGroups().size(), 1)
  TEST_EQUAL(imp_run1->getIndistinguishableProteins().size(), 1)

  // Verify run 1 search params
  const auto& sp1 = imp_run1->getSearchParameters();
  TEST_EQUAL(sp1.db, "uniprot_human")
  TEST_REAL_SIMILAR(sp1.precursor_mass_tolerance, 10.0)
  TEST_EQUAL(sp1.digestion_enzyme.getName(), "Trypsin")

  // Verify run 1 hits
  bool found_p12345 = false;
  for (const auto& hit : imp_run1->getHits())
  {
    if (hit.getAccession() == "sp|P12345|PROT1_HUMAN")
    {
      found_p12345 = true;
      TEST_REAL_SIMILAR(hit.getScore(), 0.001)
      TEST_REAL_SIMILAR(hit.getCoverage(), 45.2)
    }
  }
  TEST_EQUAL(found_p12345, true)

  // Verify run 1 protein groups
  TEST_REAL_SIMILAR(imp_run1->getProteinGroups()[0].probability, 0.99)
  TEST_EQUAL(imp_run1->getProteinGroups()[0].accessions.size(), 2)

  // Verify run 2
  TEST_EQUAL(imp_run2->getSearchEngine(), "MS-GF+")
  TEST_EQUAL(imp_run2->getHits().size(), 1)
  TEST_EQUAL(imp_run2->getHits()[0].getAccession(), "MOUSE_PROT1")

  const auto& sp2 = imp_run2->getSearchParameters();
  TEST_EQUAL(sp2.db, "uniprot_mouse")
  TEST_REAL_SIMILAR(sp2.precursor_mass_tolerance, 20.0)
  TEST_EQUAL(sp2.digestion_enzyme.getName(), "Trypsin/P")
}
END_SECTION

START_SECTION(empty input round trip)
{
  // Test with empty input
  vector<ProteinIdentification> empty_ids;

  auto sp_table = ProteinIdentificationArrowIO::exportSearchParamsToArrow(empty_ids);
  auto prot_table = ProteinIdentificationArrowIO::exportProteinsToArrow(empty_ids);
  auto pg_table = ProteinIdentificationArrowIO::exportProteinGroupsToArrow(empty_ids);

  vector<ProteinIdentification> imported_ids;
  TEST_EQUAL(ProteinIdentificationArrowIO::importSearchParamsFromArrow(sp_table, imported_ids), true)
  TEST_EQUAL(imported_ids.size(), 0)

  TEST_EQUAL(ProteinIdentificationArrowIO::importProteinsFromArrow(prot_table, imported_ids), true)
  TEST_EQUAL(imported_ids.size(), 0)

  TEST_EQUAL(ProteinIdentificationArrowIO::importProteinGroupsFromArrow(pg_table, imported_ids), true)
  TEST_EQUAL(imported_ids.size(), 0)
}
END_SECTION

START_SECTION(nullable fields round trip)
{
  // Test that nullable fields survive round-trip correctly
  vector<ProteinIdentification> orig_ids;

  ProteinIdentification prot_id;
  prot_id.setIdentifier("run_nullable");
  prot_id.setSearchEngine("TestEngine");
  prot_id.setScoreType("TestScore");
  prot_id.setHigherScoreBetter(true);

  // Hit with minimal data (many nullable fields will be null)
  ProteinHit hit;
  hit.setAccession("MINIMAL_PROT");
  hit.setScore(1.0);
  // No coverage set (COVERAGE_UNKNOWN), no sequence, no description, no target_decoy
  prot_id.insertHit(hit);
  orig_ids.push_back(prot_id);

  auto table = ProteinIdentificationArrowIO::exportProteinsToArrow(orig_ids);
  TEST_NOT_EQUAL(table, nullptr)

  vector<ProteinIdentification> imported_ids;
  TEST_EQUAL(ProteinIdentificationArrowIO::importProteinsFromArrow(table, imported_ids), true)
  TEST_EQUAL(imported_ids.size(), 1)

  const auto& imp_hit = imported_ids[0].getHits()[0];
  TEST_EQUAL(imp_hit.getAccession(), "MINIMAL_PROT")
  TEST_REAL_SIMILAR(imp_hit.getScore(), 1.0)
  // Coverage should be default (COVERAGE_UNKNOWN since null was read)
  TEST_REAL_SIMILAR(imp_hit.getCoverage(), ProteinHit::COVERAGE_UNKNOWN)
  // Sequence should be empty
  TEST_EQUAL(imp_hit.getSequence(), "")
}
END_SECTION

START_SECTION(metavalue type preservation round trip)
{
  // Test that int, float, and string metavalues preserve their types
  vector<ProteinIdentification> orig_ids;

  ProteinIdentification prot_id;
  prot_id.setIdentifier("run_mv_types");
  prot_id.setSearchEngine("TestEngine");
  prot_id.setScoreType("TestScore");
  prot_id.setHigherScoreBetter(true);

  ProteinHit hit;
  hit.setAccession("PROT_MV");
  hit.setScore(1.0);
  hit.setMetaValue("my_int", 42);
  hit.setMetaValue("my_float", 3.14);
  hit.setMetaValue("my_string", "hello_world");
  hit.setMetaValue("test_int_list", DataValue(IntList{1, 2, 3}));
  hit.setMetaValue("test_double_list", DataValue(DoubleList{1.5, 2.5}));
  hit.setMetaValue("test_string_list", DataValue(StringList{"a", "b", "c"}));
  prot_id.insertHit(hit);
  orig_ids.push_back(prot_id);

  auto table = ProteinIdentificationArrowIO::exportProteinsToArrow(orig_ids);
  TEST_NOT_EQUAL(table, nullptr)

  vector<ProteinIdentification> imported_ids;
  TEST_EQUAL(ProteinIdentificationArrowIO::importProteinsFromArrow(table, imported_ids), true)
  TEST_EQUAL(imported_ids.size(), 1)

  const auto& imp_hit = imported_ids[0].getHits()[0];
  TEST_EQUAL(int(imp_hit.getMetaValue("my_int")), 42)
  TEST_REAL_SIMILAR(double(imp_hit.getMetaValue("my_float")), 3.14)
  TEST_EQUAL(String(imp_hit.getMetaValue("my_string")), "hello_world")

  // Check list metavalue types are preserved
  TEST_EQUAL(imp_hit.getMetaValue("test_int_list").valueType(), DataValue::INT_LIST)
  TEST_EQUAL(imp_hit.getMetaValue("test_int_list") == DataValue(IntList{1, 2, 3}), true)
  TEST_EQUAL(imp_hit.getMetaValue("test_double_list").valueType(), DataValue::DOUBLE_LIST)
  TEST_EQUAL(imp_hit.getMetaValue("test_double_list") == DataValue(DoubleList{1.5, 2.5}), true)
  TEST_EQUAL(imp_hit.getMetaValue("test_string_list").valueType(), DataValue::STRING_LIST)
  TEST_EQUAL(imp_hit.getMetaValue("test_string_list") == DataValue(StringList{"a", "b", "c"}), true)
}
END_SECTION

START_SECTION(protein groups with data arrays round trip)
{
  // Test that float_data, string_data, and integer_data survive round-trip
  vector<ProteinIdentification> orig_ids;

  ProteinIdentification prot_id;
  prot_id.setIdentifier("run_data_arrays");

  ProteinIdentification::ProteinGroup pg;
  pg.probability = 0.85;
  pg.accessions = {"PROT_A", "PROT_B"};

  // Add float data array
  ProteinIdentification::ProteinGroup::FloatDataArray fda;
  fda.setName("q_values");
  fda.push_back(0.01f);
  fda.push_back(0.02f);
  pg.getFloatDataArrays().push_back(fda);

  // Add string data array
  ProteinIdentification::ProteinGroup::StringDataArray sda;
  sda.setName("annotations");
  sda.push_back("group_leader");
  sda.push_back("member");
  pg.getStringDataArrays().push_back(sda);

  // Add integer data array
  ProteinIdentification::ProteinGroup::IntegerDataArray ida;
  ida.setName("counts");
  ida.push_back(10);
  ida.push_back(20);
  pg.getIntegerDataArrays().push_back(ida);

  prot_id.insertProteinGroup(pg);
  orig_ids.push_back(prot_id);

  auto table = ProteinIdentificationArrowIO::exportProteinGroupsToArrow(orig_ids);
  TEST_NOT_EQUAL(table, nullptr)

  vector<ProteinIdentification> imported_ids;
  ProteinIdentification shell;
  shell.setIdentifier("run_data_arrays");
  imported_ids.push_back(shell);

  TEST_EQUAL(ProteinIdentificationArrowIO::importProteinGroupsFromArrow(table, imported_ids), true)
  TEST_EQUAL(imported_ids[0].getProteinGroups().size(), 1)

  const auto& imp_pg = imported_ids[0].getProteinGroups()[0];
  TEST_REAL_SIMILAR(imp_pg.probability, 0.85)

  // Float data arrays
  TEST_EQUAL(imp_pg.getFloatDataArrays().size(), 1)
  TEST_EQUAL(imp_pg.getFloatDataArrays()[0].getName(), "q_values")
  TEST_EQUAL(imp_pg.getFloatDataArrays()[0].size(), 2)
  TEST_REAL_SIMILAR(imp_pg.getFloatDataArrays()[0][0], 0.01)
  TEST_REAL_SIMILAR(imp_pg.getFloatDataArrays()[0][1], 0.02)

  // String data arrays
  TEST_EQUAL(imp_pg.getStringDataArrays().size(), 1)
  TEST_EQUAL(imp_pg.getStringDataArrays()[0].getName(), "annotations")
  TEST_EQUAL(imp_pg.getStringDataArrays()[0].size(), 2)
  TEST_EQUAL(imp_pg.getStringDataArrays()[0][0], "group_leader")
  TEST_EQUAL(imp_pg.getStringDataArrays()[0][1], "member")

  // Integer data arrays
  TEST_EQUAL(imp_pg.getIntegerDataArrays().size(), 1)
  TEST_EQUAL(imp_pg.getIntegerDataArrays()[0].getName(), "counts")
  TEST_EQUAL(imp_pg.getIntegerDataArrays()[0].size(), 2)
  TEST_EQUAL(imp_pg.getIntegerDataArrays()[0][0], 10)
  TEST_EQUAL(imp_pg.getIntegerDataArrays()[0][1], 20)
}
END_SECTION

/////////////////////////////////////////////////////////////
// Identifier handling parity helpers (synthesize / rename / check unique)
/////////////////////////////////////////////////////////////

START_SECTION(([EXTRA] synthesizeRunIdentifiers + applyRunIdentifierRename - happy path))
{
  vector<ProteinIdentification> prot_ids(2);
  prot_ids[0].setIdentifier("orig_run_A");
  prot_ids[0].setSearchEngine("Comet");
  prot_ids[0].setDateTime(DateTime::fromString("2025-08-01T12:00:00"));
  prot_ids[1].setIdentifier("orig_run_B");
  prot_ids[1].setSearchEngine("Mascot");
  prot_ids[1].setDateTime(DateTime::fromString("2025-08-02T13:00:00"));

  PeptideIdentificationList pep_ids;
  PeptideIdentification p1; p1.setIdentifier("orig_run_A"); pep_ids.push_back(p1);
  PeptideIdentification p2; p2.setIdentifier("orig_run_B"); pep_ids.push_back(p2);
  PeptideIdentification p3; p3.setIdentifier("orphan");     pep_ids.push_back(p3);

  auto rename = ProteinIdentificationArrowIO::synthesizeRunIdentifiers(prot_ids);

  // Both ProtIDs get fresh identifiers; in-place mutation visible.
  TEST_EQUAL(rename.size(), 2)
  TEST_NOT_EQUAL(prot_ids[0].getIdentifier(), "orig_run_A")
  TEST_NOT_EQUAL(prot_ids[1].getIdentifier(), "orig_run_B")
  TEST_NOT_EQUAL(prot_ids[0].getIdentifier(), prot_ids[1].getIdentifier())

  // Synthesized identifier shape: <search_engine>_<date>_<UniqueIdGenerator>.
  TEST_TRUE(prot_ids[0].getIdentifier().hasPrefix("Comet_2025-08-01T12:00:00_"))
  TEST_TRUE(prot_ids[1].getIdentifier().hasPrefix("Mascot_2025-08-02T13:00:00_"))

  // Rename map agrees with the synthesized values.
  TEST_STRING_EQUAL(rename["orig_run_A"], prot_ids[0].getIdentifier())
  TEST_STRING_EQUAL(rename["orig_run_B"], prot_ids[1].getIdentifier())

  // applyRunIdentifierRename re-stamps matching pep_ids; orphans untouched.
  ProteinIdentificationArrowIO::applyRunIdentifierRename(rename, pep_ids);
  TEST_STRING_EQUAL(pep_ids[0].getIdentifier(), prot_ids[0].getIdentifier())
  TEST_STRING_EQUAL(pep_ids[1].getIdentifier(), prot_ids[1].getIdentifier())
  TEST_STRING_EQUAL(pep_ids[2].getIdentifier(), "orphan")
}
END_SECTION

START_SECTION(([EXTRA] synthesizeRunIdentifiers - edge cases (empty search engine, invalid date)))
{
  vector<ProteinIdentification> prot_ids(2);
  // Empty search engine -> prefix becomes "unknown_".
  prot_ids[0].setIdentifier("emptySE");
  prot_ids[0].setSearchEngine("");
  prot_ids[0].setDateTime(DateTime::fromString("2025-09-09T09:09:09"));
  // Invalid (default) DateTime -> placeholder "1900-01-01T00:00:00".
  prot_ids[1].setIdentifier("noDate");
  prot_ids[1].setSearchEngine("Sage");
  prot_ids[1].setDateTime(DateTime{});

  auto rename = ProteinIdentificationArrowIO::synthesizeRunIdentifiers(prot_ids);

  TEST_TRUE(prot_ids[0].getIdentifier().hasPrefix("unknown_2025-09-09T09:09:09_"))
  TEST_TRUE(prot_ids[1].getIdentifier().hasPrefix("Sage_1900-01-01T00:00:00_"))
  TEST_EQUAL(rename.size(), 2)
}
END_SECTION

START_SECTION(([EXTRA] synthesizeRunIdentifiers - duplicate stored identifiers collapse rename map with WARN))
{
  vector<ProteinIdentification> prot_ids(2);
  // Both share stored identifier "dup" (pre-fix-#2b corruption scenario).
  prot_ids[0].setIdentifier("dup");
  prot_ids[0].setSearchEngine("Comet");
  prot_ids[0].setDateTime(DateTime::fromString("2025-09-09T09:09:09"));
  prot_ids[1].setIdentifier("dup");
  prot_ids[1].setSearchEngine("Mascot");
  prot_ids[1].setDateTime(DateTime::fromString("2025-09-09T09:09:09"));

  auto rename = ProteinIdentificationArrowIO::synthesizeRunIdentifiers(prot_ids);

  // Each ProtID gets its own distinct synthesized identifier ...
  TEST_NOT_EQUAL(prot_ids[0].getIdentifier(), prot_ids[1].getIdentifier())
  TEST_TRUE(prot_ids[0].getIdentifier().hasPrefix("Comet_"))
  TEST_TRUE(prot_ids[1].getIdentifier().hasPrefix("Mascot_"))
  // ... but rename map has 1 entry, last-seen wins. Pep_ids with stored
  // identifier "dup" will all re-stamp to the LAST ProtID (the Mascot one).
  TEST_EQUAL(rename.size(), 1)
  TEST_STRING_EQUAL(rename["dup"], prot_ids[1].getIdentifier())
}
END_SECTION

START_SECTION(static void checkUniqueIdentifiers(const std::vector<ProteinIdentification>&))
{
  // No duplicates -> no throw. (If the call throws, the test framework propagates it as a failure.)
  vector<ProteinIdentification> ok(2);
  ok[0].setIdentifier("run_X");
  ok[1].setIdentifier("run_Y");
  ProteinIdentificationArrowIO::checkUniqueIdentifiers(ok);
  TEST_TRUE(true)  // reached only if the call above didn't throw

  // Duplicate identifiers -> throws Exception::InvalidValue with the canonical
  // message text used by XMLHandler::checkUniqueIdentifiers_ — log-grepping
  // "ProteinIdentification run identifiers are not unique" finds both lanes.
  vector<ProteinIdentification> dup(2);
  dup[0].setIdentifier("dup");
  dup[1].setIdentifier("dup");
  TEST_EXCEPTION_WITH_MESSAGE(
    Exception::InvalidValue,
    ProteinIdentificationArrowIO::checkUniqueIdentifiers(dup),
    "the value 'dup' was used but is not valid; ProteinIdentification run identifiers are not unique. This can lead to loss of unique PeptideIdentification assignment. Duplicated Protein-ID is:")
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST

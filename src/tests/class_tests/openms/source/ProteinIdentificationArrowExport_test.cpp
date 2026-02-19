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
#include <OpenMS/CHEMISTRY/ProteaseDB.h>

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
  auto table = ProteinIdentificationArrowExport::exportProteinGroupsToArrow(protein_ids);
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
  auto table = ProteinIdentificationArrowExport::exportSearchParamsToArrow(protein_ids);
  TEST_NOT_EQUAL(table, nullptr)

  // Verify 1 row (one ProteinIdentification)
  TEST_EQUAL(table->num_rows(), 1)

  // Verify 25 columns
  TEST_EQUAL(table->num_columns(), 25)

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
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST

#else // WITH_PARQUET

START_TEST(ProteinIdentificationArrowExport, "$Id$")
END_TEST

#endif // WITH_PARQUET

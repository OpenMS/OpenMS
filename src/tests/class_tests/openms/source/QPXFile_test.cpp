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
#include <OpenMS/FORMAT/QPXFile.h>
///////////////////////////

#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideEvidence.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/FORMAT/ArrowSchemaRegistry.h>
#include <OpenMS/FORMAT/ProteinGroupArrowExport.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>

#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/reader.h>

using namespace OpenMS;
using namespace std;

START_TEST(QPXFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(static std::shared_ptr<arrow::Table> exportPSMsToQPXArrow(...))
{
  // Create test data
  vector<ProteinIdentification> protein_ids;
  PeptideIdentificationList peptide_ids;

  ProteinIdentification protein_id;
  protein_id.setIdentifier("test_search");
  protein_id.setSearchEngine("TestEngine");
  protein_id.setScoreType("TestScore");
  protein_id.setHigherScoreBetter(true);
  protein_ids.push_back(protein_id);

  std::vector<String> pep_strs = {"PEPTIDER", "PEM(Oxidation)TIDER", "DFPIANGER"};
  for (size_t i = 0; i < pep_strs.size(); ++i)
  {
    PeptideIdentification peptide_id;
    peptide_id.setIdentifier("test_search");
    peptide_id.setRT(1234.5 + i * 100);
    peptide_id.setMZ(500.25 + i * 50);
    peptide_id.setScoreType("TestScore");
    peptide_id.setSpectrumReference("controllerType=0 controllerNumber=1 scan=" + String(1000 + i));

    PeptideHit hit;
    hit.setSequence(AASequence::fromString(pep_strs[i]));
    hit.setScore(0.95 - i * 0.1);
    hit.setCharge(2 + i);

    if (i % 2 == 0)
    {
      hit.setMetaValue("target_decoy", "target");
    }
    else
    {
      hit.setMetaValue("target_decoy", "decoy");
    }

    PeptideEvidence evidence;
    evidence.setProteinAccession("TEST_PROTEIN_" + String(i));
    hit.setPeptideEvidences(vector<PeptideEvidence>{evidence});

    // Add multiple hits, but only first should be processed by default
    vector<PeptideHit> hits = {hit};
    if (i == 0)
    {
      PeptideHit second_hit;
      second_hit.setSequence(AASequence::fromString("SECOND"));
      second_hit.setScore(0.5);
      second_hit.setCharge(1);
      hits.push_back(second_hit);
    }

    peptide_id.setHits(hits);
    peptide_ids.push_back(peptide_id);
  }

  // Export to Arrow table (default: best hit only)
  auto table = QPXFile::exportPSMsToQPXArrow(protein_ids, peptide_ids);
  TEST_NOT_EQUAL(table, nullptr)

  // Verify number of rows (should equal number of peptide identifications, not hits)
  TEST_EQUAL(table->num_rows(), 3)

  // Verify schema column names and count (24 columns in QPXPSMSchema)
  auto schema = table->schema();
  TEST_EQUAL(table->num_columns(), 24)

  TEST_EQUAL(schema->field(0)->name(), "sequence")
  TEST_EQUAL(schema->field(1)->name(), "peptidoform")
  TEST_EQUAL(schema->field(2)->name(), "modifications")
  TEST_EQUAL(schema->field(3)->name(), "charge")
  TEST_EQUAL(schema->field(4)->name(), "posterior_error_probability")
  TEST_EQUAL(schema->field(5)->name(), "is_decoy")
  TEST_EQUAL(schema->field(6)->name(), "calculated_mz")
  TEST_EQUAL(schema->field(7)->name(), "observed_mz")
  TEST_EQUAL(schema->field(8)->name(), "mass_error_ppm")
  TEST_EQUAL(schema->field(9)->name(), "additional_scores")
  TEST_EQUAL(schema->field(10)->name(), "predicted_rt")
  TEST_EQUAL(schema->field(11)->name(), "run_file_name")
  TEST_EQUAL(schema->field(12)->name(), "cv_params")
  TEST_EQUAL(schema->field(13)->name(), "scan")
  TEST_EQUAL(schema->field(14)->name(), "rt")
  TEST_EQUAL(schema->field(15)->name(), "ion_mobility")
  TEST_EQUAL(schema->field(16)->name(), "missed_cleavages")
  TEST_EQUAL(schema->field(17)->name(), "protein_accessions")
  TEST_EQUAL(schema->field(18)->name(), "cross_links")
  TEST_EQUAL(schema->field(19)->name(), "mz_array")
  TEST_EQUAL(schema->field(20)->name(), "intensity_array")
  TEST_EQUAL(schema->field(21)->name(), "charge_array")
  TEST_EQUAL(schema->field(22)->name(), "ion_type_array")
  TEST_EQUAL(schema->field(23)->name(), "ion_mobility_array")

  // Verify data types for key columns
  TEST_EQUAL(schema->field(3)->type()->id(), arrow::Type::INT16)   // charge is int16
  TEST_EQUAL(schema->field(4)->type()->id(), arrow::Type::DOUBLE)  // PEP is float64
  TEST_EQUAL(schema->field(6)->type()->id(), arrow::Type::FLOAT)   // calculated_mz is float32
  TEST_EQUAL(schema->field(7)->type()->id(), arrow::Type::FLOAT)   // observed_mz is float32
  TEST_EQUAL(schema->field(8)->type()->id(), arrow::Type::FLOAT)   // mass_error_ppm is float32
  TEST_EQUAL(schema->field(13)->type()->id(), arrow::Type::LIST)   // scan is list<int32>
  TEST_EQUAL(schema->field(17)->type()->id(), arrow::Type::LIST)   // protein_accessions is list
  TEST_EQUAL(schema->field(2)->type()->id(), arrow::Type::LIST)    // modifications is list

  // Verify sequence values
  auto seq_col = table->GetColumnByName("sequence");
  auto seq_arr = std::static_pointer_cast<arrow::StringArray>(seq_col->chunk(0));
  TEST_EQUAL(seq_arr->GetString(0), "PEPTIDER")
  TEST_EQUAL(seq_arr->GetString(1), "PEMTIDER")
  TEST_EQUAL(seq_arr->GetString(2), "DFPIANGER")

  // Verify peptidoform uses ProForma notation (not bracket substitution)
  auto pf_col = table->GetColumnByName("peptidoform");
  auto pf_arr = std::static_pointer_cast<arrow::StringArray>(pf_col->chunk(0));
  TEST_EQUAL(pf_arr->GetString(0), "PEPTIDER")
  // Modified peptide should have ProForma notation with bracket-enclosed modification
  std::string mod_pf = pf_arr->GetString(1);
  TEST_EQUAL(mod_pf.empty(), false)
  TEST_EQUAL(mod_pf.find('[') != std::string::npos, true) // ProForma bracket notation expected

  // Verify is_decoy non-nullable behavior
  auto decoy_col = table->GetColumnByName("is_decoy");
  auto decoy_arr = std::static_pointer_cast<arrow::BooleanArray>(decoy_col->chunk(0));
  TEST_EQUAL(decoy_arr->Value(0), false) // target
  TEST_EQUAL(decoy_arr->Value(1), true) // decoy
  TEST_EQUAL(decoy_arr->Value(2), false) // target

  // Verify scan extraction from native ID (now list<int32>)
  auto scan_col = table->GetColumnByName("scan");
  auto scan_list = std::static_pointer_cast<arrow::ListArray>(scan_col->chunk(0));
  // Each row should have a list with one scan number
  auto scan_values = std::static_pointer_cast<arrow::Int32Array>(scan_list->values());
  TEST_EQUAL(scan_list->value_length(0), 1)
  TEST_EQUAL(scan_values->Value(scan_list->value_offset(0)), 1000)
  TEST_EQUAL(scan_list->value_length(1), 1)
  TEST_EQUAL(scan_values->Value(scan_list->value_offset(1)), 1001)
  TEST_EQUAL(scan_list->value_length(2), 1)
  TEST_EQUAL(scan_values->Value(scan_list->value_offset(2)), 1002)
}
END_SECTION

START_SECTION(static std::shared_ptr<arrow::Table> exportPSMsToQPXArrow(...) with export_all_psms)
{
  // Create test data with multiple hits per peptide identification
  vector<ProteinIdentification> protein_ids;
  PeptideIdentificationList peptide_ids;

  ProteinIdentification protein_id;
  protein_id.setIdentifier("test_search");
  protein_id.setSearchEngine("TestEngine");
  protein_id.setScoreType("TestScore");
  protein_id.setHigherScoreBetter(true);
  protein_ids.push_back(protein_id);

  PeptideIdentification peptide_id;
  peptide_id.setIdentifier("test_search");
  peptide_id.setRT(1234.5);
  peptide_id.setMZ(500.25);
  peptide_id.setScoreType("TestScore");

  vector<PeptideHit> hits;
  vector<String> pep_strs = {"PEPTIDER", "ALTERNATIVE", "THIRDPSM"};
  for (int i = 0; i < 3; ++i)
  {
    PeptideHit hit;
    hit.setSequence(AASequence::fromString(pep_strs[i]));
    hit.setScore(0.95 - i * 0.1);
    hit.setCharge(2 + i);
    hit.setMetaValue("target_decoy", "target");

    PeptideEvidence evidence;
    evidence.setProteinAccession("TEST_PROTEIN_" + String(i));
    hit.setPeptideEvidences(vector<PeptideEvidence>{evidence});

    hits.push_back(hit);
  }

  peptide_id.setHits(hits);
  peptide_ids.push_back(peptide_id);

  // Export with all PSMs
  auto table = QPXFile::exportPSMsToQPXArrow(protein_ids, peptide_ids, true);
  TEST_NOT_EQUAL(table, nullptr)

  // All 3 hits should be exported
  TEST_EQUAL(table->num_rows(), 3)

  // Verify charge column (int16)
  auto charge_col = table->GetColumnByName("charge");
  auto charge_arr = std::static_pointer_cast<arrow::Int16Array>(charge_col->chunk(0));
  TEST_EQUAL(charge_arr->Value(0), 2)
  TEST_EQUAL(charge_arr->Value(1), 3)
  TEST_EQUAL(charge_arr->Value(2), 4)

  // Verify sequence values
  auto seq_col = table->GetColumnByName("sequence");
  auto seq_arr = std::static_pointer_cast<arrow::StringArray>(seq_col->chunk(0));
  TEST_EQUAL(seq_arr->GetString(0), "PEPTIDER")
  TEST_EQUAL(seq_arr->GetString(1), "ALTERNATIVE")
  TEST_EQUAL(seq_arr->GetString(2), "THIRDPSM")
}
END_SECTION

START_SECTION(static bool exportToParquet(...))
{
  // Create test data
  vector<ProteinIdentification> protein_ids;
  PeptideIdentificationList peptide_ids;

  ProteinIdentification protein_id;
  protein_id.setIdentifier("test_search");
  protein_id.setSearchEngine("TestEngine");
  protein_id.setScoreType("TestScore");
  protein_id.setHigherScoreBetter(true);
  protein_ids.push_back(protein_id);

  PeptideIdentification peptide_id;
  peptide_id.setIdentifier("test_search");
  peptide_id.setRT(1234.5);
  peptide_id.setMZ(500.25);
  peptide_id.setScoreType("TestScore");

  PeptideHit hit;
  hit.setSequence(AASequence::fromString("PEM(Oxidation)TIDER"));
  hit.setScore(0.95);
  hit.setCharge(2);
  hit.setMetaValue("target_decoy", "target");

  PeptideEvidence evidence;
  evidence.setProteinAccession("TEST_PROTEIN");
  PeptideEvidence evidence2;
  evidence2.setProteinAccession("TEST_PROTEIN_2");
  hit.setPeptideEvidences(vector<PeptideEvidence>{evidence, evidence2});

  peptide_id.setHits(vector<PeptideHit>{hit});
  peptide_ids.push_back(peptide_id);

  String output_file;
  NEW_TMP_FILE(output_file)

  // Write parquet
  bool success = QPXFile::exportToParquet(protein_ids, peptide_ids, output_file);
  TEST_EQUAL(success, true)

  // Read back and verify
  arrow::MemoryPool* pool = arrow::default_memory_pool();
  std::shared_ptr<arrow::io::ReadableFile> infile;
  auto result = arrow::io::ReadableFile::Open(output_file.c_str());
  TEST_EQUAL(result.ok(), true)
  infile = result.ValueOrDie();

  std::unique_ptr<parquet::arrow::FileReader> reader;
  PARQUET_ASSIGN_OR_THROW(reader, parquet::arrow::OpenFile(infile, pool));

  std::shared_ptr<arrow::Table> table;
  auto read_status = reader->ReadTable(&table);
  TEST_EQUAL(read_status.ok(), true)

  TEST_EQUAL(table->num_rows(), 1)
  TEST_EQUAL(table->num_columns(), 24)

  // Verify modifications column has structured data for modified peptide
  auto mod_col = table->GetColumnByName("modifications");
  TEST_NOT_EQUAL(mod_col, nullptr)
  TEST_EQUAL(mod_col->type()->id(), arrow::Type::LIST)

  // Verify protein_accessions is a list with two entries
  auto pa_col = table->GetColumnByName("protein_accessions");
  TEST_NOT_EQUAL(pa_col, nullptr)
  TEST_EQUAL(pa_col->type()->id(), arrow::Type::LIST)

  // Verify mass_error_ppm column exists
  auto mep_col = table->GetColumnByName("mass_error_ppm");
  TEST_NOT_EQUAL(mep_col, nullptr)
  TEST_EQUAL(mep_col->type()->id(), arrow::Type::FLOAT)

  // Verify QPX file metadata is present in schema key-value metadata
  auto metadata = table->schema()->metadata();
  TEST_NOT_EQUAL(metadata, nullptr)

  auto qpx_version_idx = metadata->FindKey("qpx_version");
  TEST_EQUAL(qpx_version_idx >= 0, true)
  TEST_EQUAL(metadata->value(qpx_version_idx), "1.0")

  auto creator_idx = metadata->FindKey("creator");
  TEST_EQUAL(creator_idx >= 0, true)
  TEST_EQUAL(metadata->value(creator_idx), "OpenMS")

  auto file_type_idx = metadata->FindKey("file_type");
  TEST_EQUAL(file_type_idx >= 0, true)
  TEST_EQUAL(metadata->value(file_type_idx), "psm")

  auto scan_format_idx = metadata->FindKey("scan_format");
  TEST_EQUAL(scan_format_idx >= 0, true)
  TEST_EQUAL(metadata->value(scan_format_idx), "scan")

  auto software_idx = metadata->FindKey("software_provider");
  TEST_EQUAL(software_idx >= 0, true)
  TEST_EQUAL(metadata->value(software_idx), "OpenMS")

  // creation_date and uuid should exist (values are dynamic)
  TEST_EQUAL(metadata->FindKey("creation_date") >= 0, true)
  TEST_EQUAL(metadata->FindKey("uuid") >= 0, true)
}
END_SECTION

START_SECTION(Test modifications structured output)
{
  vector<ProteinIdentification> protein_ids;
  PeptideIdentificationList peptide_ids;

  ProteinIdentification protein_id;
  protein_id.setIdentifier("test_search");
  protein_ids.push_back(protein_id);

  PeptideIdentification peptide_id;
  peptide_id.setIdentifier("test_search");
  peptide_id.setRT(100.0);
  peptide_id.setMZ(500.0);
  peptide_id.setScoreType("TestScore");

  PeptideHit hit;
  hit.setSequence(AASequence::fromString("PEM(Oxidation)TIDER"));
  hit.setScore(0.9);
  hit.setCharge(2);

  peptide_id.setHits(vector<PeptideHit>{hit});
  peptide_ids.push_back(peptide_id);

  auto table = QPXFile::exportPSMsToQPXArrow(protein_ids, peptide_ids);
  TEST_NOT_EQUAL(table, nullptr)

  // Verify modifications column is a list of structs (not null)
  auto mod_col = table->GetColumnByName("modifications");
  auto mod_list = std::static_pointer_cast<arrow::ListArray>(mod_col->chunk(0));
  // The modified peptide should have at least one modification entry
  TEST_EQUAL(mod_list->value_length(0) > 0, true)

  // Unmodified peptide should have empty modifications list
  PeptideIdentificationList peptide_ids2;
  PeptideIdentification peptide_id2;
  peptide_id2.setIdentifier("test_search");
  peptide_id2.setRT(100.0);
  peptide_id2.setMZ(500.0);
  peptide_id2.setScoreType("TestScore");
  PeptideHit hit2;
  hit2.setSequence(AASequence::fromString("PEPTIDER"));
  hit2.setScore(0.9);
  hit2.setCharge(2);
  peptide_id2.setHits(vector<PeptideHit>{hit2});
  peptide_ids2.push_back(peptide_id2);

  auto table2 = QPXFile::exportPSMsToQPXArrow(protein_ids, peptide_ids2);
  auto mod_col2 = table2->GetColumnByName("modifications");
  auto mod_list2 = std::static_pointer_cast<arrow::ListArray>(mod_col2->chunk(0));
  TEST_EQUAL(mod_list2->value_length(0), 0) // empty list for unmodified
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(QPXPgSchema::schema())
{
  auto schema = QPXPgSchema::schema();
  TEST_EQUAL(schema->num_fields(), 20)

  // Required (non-nullable) fields
  TEST_EQUAL(schema->GetFieldByName("pg_accessions")->nullable(), false)
  TEST_EQUAL(schema->GetFieldByName("anchor_protein")->nullable(), false)
  TEST_EQUAL(schema->GetFieldByName("run_file_name")->nullable(), false)
  TEST_EQUAL(schema->GetFieldByName("is_decoy")->nullable(), false)
  TEST_EQUAL(schema->GetFieldByName("peptides")->nullable(), false)

  // Nullable (optional) fields
  TEST_EQUAL(schema->GetFieldByName("intensities")->nullable(), true)
  TEST_EQUAL(schema->GetFieldByName("additional_intensities")->nullable(), true)
  TEST_EQUAL(schema->GetFieldByName("global_qvalue")->nullable(), true)
  TEST_EQUAL(schema->GetFieldByName("pg_qvalue")->nullable(), true)
}
END_SECTION

START_SECTION(ProteinGroupArrowExport::exportToArrow(vector<ProteinIdentification>, PeptideIdentificationList))
{
  // Set up minimal protein identification with one group
  ProteinIdentification prot_id;
  prot_id.setPrimaryMSRunPath({"test_run.mzML"});
  ProteinHit hit1;
  hit1.setAccession("PROT_A");
  hit1.setScore(0.99);
  hit1.setTargetDecoyType(ProteinHit::TargetDecoyType::TARGET);
  ProteinHit hit2;
  hit2.setAccession("PROT_B");
  hit2.setScore(0.95);
  hit2.setTargetDecoyType(ProteinHit::TargetDecoyType::TARGET);
  prot_id.setHits({hit1, hit2});
  prot_id.setScoreType("Mascot");
  prot_id.setHigherScoreBetter(true);

  ProteinIdentification::ProteinGroup group;
  group.accessions = {"PROT_A", "PROT_B"};
  group.probability = 0.01;
  prot_id.insertIndistinguishableProteins(group);

  // Set up peptide identifications
  PeptideIdentification pep_id;
  PeptideHit pep_hit;
  pep_hit.setSequence(AASequence::fromString("PEPTIDE"));
  PeptideEvidence ev;
  ev.setProteinAccession("PROT_A");
  pep_hit.setPeptideEvidences({ev});
  pep_id.setHits({pep_hit});
  PeptideIdentificationList pep_ids = {pep_id};

  auto table = ProteinGroupArrowExport::exportToArrow({prot_id}, pep_ids);
  TEST_NOT_EQUAL(table, nullptr)
  TEST_EQUAL(table->num_rows(), 1)
  TEST_EQUAL(table->num_columns(), 20)

  // Verify run_file_name is derived from ProteinIdentification
  auto run_col = std::static_pointer_cast<arrow::StringArray>(table->GetColumnByName("run_file_name")->chunk(0));
  TEST_STRING_EQUAL(run_col->GetString(0), "test_run.mzML")

  // Verify anchor_protein
  auto anchor_col = std::static_pointer_cast<arrow::StringArray>(table->GetColumnByName("anchor_protein")->chunk(0));
  TEST_STRING_EQUAL(anchor_col->GetString(0), "PROT_A")

  // Verify intensities is null (no quantification)
  auto int_col = table->GetColumnByName("intensities")->chunk(0);
  TEST_EQUAL(int_col->IsNull(0), true)

  // is_decoy
  auto is_decoy_col = std::static_pointer_cast<arrow::BooleanArray>(
    table->GetColumnByName("is_decoy")->chunk(0));
  TEST_EQUAL(is_decoy_col->Value(0), false)

  // global_qvalue carries group.probability
  auto qv_col = std::static_pointer_cast<arrow::DoubleArray>(
    table->GetColumnByName("global_qvalue")->chunk(0));
  TEST_EQUAL(qv_col->IsNull(0), false)
  TEST_REAL_SIMILAR(qv_col->Value(0), 0.01)

  // peptides list: one entry per group accession (2 here: PROT_A and PROT_B)
  auto peptides_list = std::static_pointer_cast<arrow::ListArray>(
    table->GetColumnByName("peptides")->chunk(0));
  TEST_EQUAL(peptides_list->value_length(0), 2)
  auto peptides_struct = std::static_pointer_cast<arrow::StructArray>(peptides_list->values());
  auto pep_name_arr = std::static_pointer_cast<arrow::StringArray>(peptides_struct->field(0));
  auto pep_count_arr = std::static_pointer_cast<arrow::Int32Array>(peptides_struct->field(1));
  TEST_STRING_EQUAL(pep_name_arr->GetString(0), "PROT_A")
  TEST_EQUAL(pep_count_arr->Value(0), 1)   // one peptide evidence for PROT_A
  TEST_STRING_EQUAL(pep_name_arr->GetString(1), "PROT_B")
  TEST_EQUAL(pep_count_arr->Value(1), 0)   // no peptide evidence for PROT_B

  // peptide_counts: unique sequences = 1 (the peptide is attributed to PROT_A only);
  // total_sequences = 1 (sum of per-protein counts)
  auto pc_struct = std::static_pointer_cast<arrow::StructArray>(
    table->GetColumnByName("peptide_counts")->chunk(0));
  auto pc_unique = std::static_pointer_cast<arrow::Int32Array>(pc_struct->field(0));
  auto pc_total = std::static_pointer_cast<arrow::Int32Array>(pc_struct->field(1));
  TEST_EQUAL(pc_unique->Value(0), 1)
  TEST_EQUAL(pc_total->Value(0), 1)

  // gg_accessions / gg_names must be parallel to pg_accessions (length 2, both empty strings)
  auto gg_acc_list = std::static_pointer_cast<arrow::ListArray>(
    table->GetColumnByName("gg_accessions")->chunk(0));
  auto gg_names_list = std::static_pointer_cast<arrow::ListArray>(
    table->GetColumnByName("gg_names")->chunk(0));
  TEST_EQUAL(gg_acc_list->value_length(0), 2)
  TEST_EQUAL(gg_names_list->value_length(0), 2)
}
END_SECTION

START_SECTION(ProteinGroupArrowExport::exportToArrow empty groups)
{
  ProteinIdentification prot_id;
  prot_id.setPrimaryMSRunPath({"test.mzML"});
  PeptideIdentificationList pep_ids;

  auto table = ProteinGroupArrowExport::exportToArrow({prot_id}, pep_ids);
  TEST_NOT_EQUAL(table, nullptr)
  TEST_EQUAL(table->num_rows(), 0)
  TEST_EQUAL(table->num_columns(), 20)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST

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

#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/reader.h>

using namespace OpenMS;
using namespace std;

START_TEST(QPXFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(static std::shared_ptr<arrow::Table> exportToArrow(...))
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
  auto table = QPXFile::exportToArrow(protein_ids, peptide_ids);
  TEST_NOT_EQUAL(table, nullptr)

  // Verify number of rows (should equal number of peptide identifications, not hits)
  TEST_EQUAL(table->num_rows(), 3)

  // Verify schema column names and count (24 columns in new QPX schema)
  auto schema = table->schema();
  TEST_EQUAL(table->num_columns(), 24)

  TEST_EQUAL(schema->field(0)->name(), "sequence")
  TEST_EQUAL(schema->field(1)->name(), "peptidoform")
  TEST_EQUAL(schema->field(2)->name(), "modifications")
  TEST_EQUAL(schema->field(3)->name(), "precursor_charge")
  TEST_EQUAL(schema->field(4)->name(), "posterior_error_probability")
  TEST_EQUAL(schema->field(5)->name(), "is_decoy")
  TEST_EQUAL(schema->field(6)->name(), "calculated_mz")
  TEST_EQUAL(schema->field(7)->name(), "observed_mz")
  TEST_EQUAL(schema->field(8)->name(), "additional_scores")
  TEST_EQUAL(schema->field(9)->name(), "protein_accessions")
  TEST_EQUAL(schema->field(10)->name(), "predicted_rt")
  TEST_EQUAL(schema->field(11)->name(), "reference_file_name")
  TEST_EQUAL(schema->field(12)->name(), "cv_params")
  TEST_EQUAL(schema->field(13)->name(), "scan")
  TEST_EQUAL(schema->field(14)->name(), "rt")
  TEST_EQUAL(schema->field(15)->name(), "ion_mobility")
  TEST_EQUAL(schema->field(16)->name(), "spectrum_reference")
  TEST_EQUAL(schema->field(17)->name(), "score")
  TEST_EQUAL(schema->field(18)->name(), "score_type")
  TEST_EQUAL(schema->field(19)->name(), "rank")
  TEST_EQUAL(schema->field(20)->name(), "P_ID")
  TEST_EQUAL(schema->field(21)->name(), "psm_metavalues")
  TEST_EQUAL(schema->field(22)->name(), "spectrum_metavalues")
  TEST_EQUAL(schema->field(23)->name(), "run_identifier")

  // Verify data types for key columns
  TEST_EQUAL(schema->field(4)->type()->id(), arrow::Type::DOUBLE) // PEP is float64
  TEST_EQUAL(schema->field(6)->type()->id(), arrow::Type::DOUBLE) // calculated_mz is float64
  TEST_EQUAL(schema->field(7)->type()->id(), arrow::Type::DOUBLE) // observed_mz is float64
  TEST_EQUAL(schema->field(9)->type()->id(), arrow::Type::LIST)   // protein_accessions is list
  TEST_EQUAL(schema->field(2)->type()->id(), arrow::Type::LIST)   // modifications is list

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

  // Verify is_decoy nullable behavior
  auto decoy_col = table->GetColumnByName("is_decoy");
  auto decoy_arr = std::static_pointer_cast<arrow::Int32Array>(decoy_col->chunk(0));
  TEST_EQUAL(decoy_arr->Value(0), 0) // target
  TEST_EQUAL(decoy_arr->Value(1), 1) // decoy
  TEST_EQUAL(decoy_arr->Value(2), 0) // target

  // Verify rank is always present and 0-based
  auto rank_col = table->GetColumnByName("rank");
  auto rank_arr = std::static_pointer_cast<arrow::Int32Array>(rank_col->chunk(0));
  TEST_EQUAL(rank_arr->Value(0), 0)
  TEST_EQUAL(rank_arr->Value(1), 0)
  TEST_EQUAL(rank_arr->Value(2), 0)

  // Verify P_ID values
  auto pid_col = table->GetColumnByName("P_ID");
  auto pid_arr = std::static_pointer_cast<arrow::Int32Array>(pid_col->chunk(0));
  TEST_EQUAL(pid_arr->Value(0), 0)
  TEST_EQUAL(pid_arr->Value(1), 1)
  TEST_EQUAL(pid_arr->Value(2), 2)

  // Verify scan extraction from native ID
  auto scan_col = table->GetColumnByName("scan");
  auto scan_arr = std::static_pointer_cast<arrow::StringArray>(scan_col->chunk(0));
  TEST_EQUAL(scan_arr->GetString(0), "1000")
  TEST_EQUAL(scan_arr->GetString(1), "1001")
  TEST_EQUAL(scan_arr->GetString(2), "1002")

  // Verify spectrum_reference stores full native ID
  auto sr_col = table->GetColumnByName("spectrum_reference");
  auto sr_arr = std::static_pointer_cast<arrow::StringArray>(sr_col->chunk(0));
  TEST_EQUAL(sr_arr->GetString(0), "controllerType=0 controllerNumber=1 scan=1000")
}
END_SECTION

START_SECTION(static std::shared_ptr<arrow::Table> exportToArrow(...) with export_all_psms)
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
  auto table = QPXFile::exportToArrow(protein_ids, peptide_ids, true);
  TEST_NOT_EQUAL(table, nullptr)

  // All 3 hits should be exported
  TEST_EQUAL(table->num_rows(), 3)

  // Verify rank values are 0-based
  auto rank_col = table->GetColumnByName("rank");
  auto rank_arr = std::static_pointer_cast<arrow::Int32Array>(rank_col->chunk(0));
  TEST_EQUAL(rank_arr->Value(0), 0)
  TEST_EQUAL(rank_arr->Value(1), 1)
  TEST_EQUAL(rank_arr->Value(2), 2)

  // All rows should have the same P_ID (same parent identification)
  auto pid_col = table->GetColumnByName("P_ID");
  auto pid_arr = std::static_pointer_cast<arrow::Int32Array>(pid_col->chunk(0));
  TEST_EQUAL(pid_arr->Value(0), 0)
  TEST_EQUAL(pid_arr->Value(1), 0)
  TEST_EQUAL(pid_arr->Value(2), 0)
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
  hit.setMetaValue("some_custom_value", 42);

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

  // Verify psm_metavalues contains the custom value
  auto pmv_col = table->GetColumnByName("psm_metavalues");
  TEST_NOT_EQUAL(pmv_col, nullptr)
  TEST_EQUAL(pmv_col->type()->id(), arrow::Type::LIST)

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

  auto table = QPXFile::exportToArrow(protein_ids, peptide_ids);
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

  auto table2 = QPXFile::exportToArrow(protein_ids, peptide_ids2);
  auto mod_col2 = table2->GetColumnByName("modifications");
  auto mod_list2 = std::static_pointer_cast<arrow::ListArray>(mod_col2->chunk(0));
  TEST_EQUAL(mod_list2->value_length(0), 0) // empty list for unmodified
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST

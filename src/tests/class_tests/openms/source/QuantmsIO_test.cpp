// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FORMAT/QuantmsIO.h>
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

START_TEST(QuantmsIO, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

QuantmsIO* ptr = nullptr;
QuantmsIO* null_ptr = nullptr;

START_SECTION(QuantmsIO())
{
  ptr = new QuantmsIO();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~QuantmsIO())
{
  delete ptr;
}
END_SECTION

START_SECTION((void store(const String& filename, const std::vector<ProteinIdentification>& protein_identifications, const PeptideIdentificationList& peptide_identifications)))
{
  QuantmsIO file;
  
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
  // Create multiple peptide identifications to test row count
  for (size_t i = 0; i < pep_strs.size(); ++i)
  {
    PeptideIdentification peptide_id;
    peptide_id.setIdentifier("test_search");
    peptide_id.setRT(1234.5 + i * 100);
    peptide_id.setMZ(500.25 + i * 50);
    peptide_id.setScoreType("TestScore");
    
    PeptideHit hit;
    hit.setSequence(AASequence::fromString(pep_strs[i]));
    hit.setScore(0.95 - i * 0.1);
    hit.setCharge(2 + i);
    
    // Add target_decoy metavalue for testing decoy detection
    if (i % 2 == 0)
    {
      hit.setMetaValue("target_decoy", "target");
    }
    else
    {
      hit.setMetaValue("target_decoy", "decoy");
    }
    
    // Add PEP score as metavalue for testing
    hit.setMetaValue("pep", 0.01 + i * 0.005);
    
    PeptideEvidence evidence;
    evidence.setProteinAccession("TEST_PROTEIN_" + String(i));
    hit.setPeptideEvidences(vector<PeptideEvidence>{evidence});
    
    // Add multiple hits, but only first should be processed
    vector<PeptideHit> hits = {hit};
    if (i == 0) {
      PeptideHit second_hit;
      second_hit.setSequence(AASequence::fromString("SECOND"));
      second_hit.setScore(0.5);
      second_hit.setCharge(1);
      hits.push_back(second_hit);
    }
    
    peptide_id.setHits(hits);
    peptide_ids.push_back(peptide_id);
  }
  
  String output_file;
  NEW_TMP_FILE(output_file)
  
  // Store the data
  TEST_NOT_EQUAL(peptide_ids.size(), 0)
  file.store(output_file, protein_ids, peptide_ids);
  
  // Read back the parquet file and verify rows and columns
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
  
  // Verify number of rows (should equal number of peptide identifications)
  TEST_EQUAL(table->num_rows(), 3)
  
  // Verify number of columns (should match the quantms.io PSM schema without file_metadata)
  TEST_EQUAL(table->num_columns(), 19)
  
  // Verify column names match the schema
  auto schema = table->schema();
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
  TEST_EQUAL(schema->field(16)->name(), "number_peaks")
  TEST_EQUAL(schema->field(17)->name(), "mz_array")
  TEST_EQUAL(schema->field(18)->name(), "intensity_array")
  
  // TODO: Add verification of file-level metadata when parquet metadata API is available
  // The file metadata should be stored at the Parquet file level, not as schema fields
}
END_SECTION

START_SECTION((void store(const String& filename, const std::vector<ProteinIdentification>& protein_identifications, const PeptideIdentificationList& peptide_identifications, bool export_all_psms)))
{
  QuantmsIO file;
  
  // Create test data with multiple hits per peptide identification
  vector<ProteinIdentification> protein_ids;
  PeptideIdentificationList peptide_ids;
  
  ProteinIdentification protein_id;
  protein_id.setIdentifier("test_search");
  protein_id.setSearchEngine("TestEngine");
  protein_id.setScoreType("TestScore");
  protein_id.setHigherScoreBetter(true);
  protein_ids.push_back(protein_id);

  // Create one peptide identification with multiple hits
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
    hit.setRank(i); // Set rank explicitly (0-based)
    hit.setMetaValue("target_decoy", "target");
    
    PeptideEvidence evidence;
    evidence.setProteinAccession("TEST_PROTEIN_" + String(i));
    hit.setPeptideEvidences(vector<PeptideEvidence>{evidence});
    
    hits.push_back(hit);
  }
  
  peptide_id.setHits(hits);
  peptide_ids.push_back(peptide_id);
  
  String output_file;
  NEW_TMP_FILE(output_file)
  
  // Store the data with export_all_psms = true
  file.store(output_file, protein_ids, peptide_ids, true);
  
  // Read back the parquet file and verify rows and columns
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
  
  // Verify number of rows (should equal number of peptide hits, not peptide identifications)
  TEST_EQUAL(table->num_rows(), 3)
  
  // Verify number of columns (should include rank column)
  TEST_EQUAL(table->num_columns(), 20) // 19 + 1 rank column
  
  // Verify rank column exists and is in the correct position (after precursor_charge)
  auto schema = table->schema();
  TEST_EQUAL(schema->field(4)->name(), "rank")
  
  // Verify rank values
  auto rank_column = table->GetColumnByName("rank");
  TEST_NOT_EQUAL(rank_column, nullptr)
  auto rank_array = std::static_pointer_cast<arrow::Int32Array>(rank_column->chunk(0));
  TEST_EQUAL(rank_array->Value(0), 1) // First hit should have rank 1
  TEST_EQUAL(rank_array->Value(1), 2) // Second hit should have rank 2
  TEST_EQUAL(rank_array->Value(2), 3) // Third hit should have rank 3
}
END_SECTION

START_SECTION((void store(const String& filename, const std::vector<ProteinIdentification>& protein_identifications, const PeptideIdentificationList& peptide_identifications, bool export_all_psms, const std::set<String>& meta_value_keys)))
{
  QuantmsIO file;
  
  // Create test data 
  vector<ProteinIdentification> protein_ids;
  PeptideIdentificationList peptide_ids;
  
  ProteinIdentification protein_id;
  protein_id.setIdentifier("test_search");
  protein_id.setSearchEngine("TestEngine");
  protein_id.setScoreType("TestScore");
  protein_id.setHigherScoreBetter(true);
  protein_ids.push_back(protein_id);

  // Create peptide identification with meta values
  PeptideIdentification peptide_id;
  peptide_id.setIdentifier("test_search");
  peptide_id.setRT(1234.5);
  peptide_id.setMZ(500.25);
  peptide_id.setScoreType("TestScore");
  
  PeptideHit hit;
  hit.setSequence(AASequence::fromString("PEPTIDER"));
  hit.setScore(0.95);
  hit.setCharge(2);
  hit.setMetaValue("target_decoy", "target");
  
  // Add meta values for testing
  hit.setMetaValue("confidence", 0.85);
  hit.setMetaValue("mass_error", 2.5);
  hit.setMetaValue("score_type", String("E-value"));
  
  PeptideEvidence evidence;
  evidence.setProteinAccession("TEST_PROTEIN");
  hit.setPeptideEvidences(vector<PeptideEvidence>{evidence});
  
  peptide_id.setHits(vector<PeptideHit>{hit});
  peptide_ids.push_back(peptide_id);
  
  // Define meta value keys to export
  std::set<String> meta_keys = {"confidence", "mass_error", "nonexistent_key"};
  
  String output_file;
  NEW_TMP_FILE(output_file)
  
  // Store the data with meta values
  file.store(output_file, protein_ids, peptide_ids, false, meta_keys);
  
  // Read back the parquet file
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
  
  // Verify number of columns (19 standard + 3 meta value columns)
  TEST_EQUAL(table->num_columns(), 22)
  
  // Verify meta value columns exist
  auto schema = table->schema();
  TEST_EQUAL(schema->field(19)->name(), "confidence")
  TEST_EQUAL(schema->field(20)->name(), "mass_error")
  TEST_EQUAL(schema->field(21)->name(), "nonexistent_key")
  
  // Verify meta value data and types
  auto confidence_column = table->GetColumnByName("confidence");
  TEST_NOT_EQUAL(confidence_column, nullptr)
  // confidence should be a double array since 0.85 is a double value
  auto confidence_array = std::static_pointer_cast<arrow::DoubleArray>(confidence_column->chunk(0));
  TEST_EQUAL(confidence_array->Value(0), 0.85)
  
  auto mass_error_column = table->GetColumnByName("mass_error");
  TEST_NOT_EQUAL(mass_error_column, nullptr)
  // mass_error should be a double array since 2.5 is a double value  
  auto mass_error_array = std::static_pointer_cast<arrow::DoubleArray>(mass_error_column->chunk(0));
  TEST_EQUAL(mass_error_array->Value(0), 2.5)
  
  auto nonexistent_column = table->GetColumnByName("nonexistent_key");
  TEST_NOT_EQUAL(nonexistent_column, nullptr)
  // nonexistent_key should be string type (default) with null value
  auto nonexistent_array = std::static_pointer_cast<arrow::StringArray>(nonexistent_column->chunk(0));
  TEST_EQUAL(nonexistent_array->IsNull(0), true)
}
END_SECTION

START_SECTION((Test meta value type detection and proper Arrow column types))
{
  QuantmsIO file;
  
  // Create test data with different meta value types
  vector<ProteinIdentification> protein_ids;
  PeptideIdentificationList peptide_ids;
  
  ProteinIdentification protein_id;
  protein_id.setIdentifier("test_search");
  protein_id.setSearchEngine("TestEngine");
  protein_id.setScoreType("TestScore");
  protein_id.setHigherScoreBetter(true);
  protein_ids.push_back(protein_id);

  // Create peptide identification with different types of meta values
  PeptideIdentification peptide_id;
  peptide_id.setIdentifier("test_search");
  peptide_id.setRT(1234.5);
  peptide_id.setMZ(500.25);
  peptide_id.setScoreType("TestScore");
  
  PeptideHit hit;
  hit.setSequence(AASequence::fromString("PEPTIDER"));
  hit.setScore(0.95);
  hit.setCharge(2);
  hit.setMetaValue("target_decoy", "target");
  
  // Add different types of meta values for testing
  hit.setMetaValue("string_value", String("test_string"));      // STRING_VALUE
  hit.setMetaValue("int_value", 42);                           // INT_VALUE
  hit.setMetaValue("double_value", 3.14159);                   // DOUBLE_VALUE
  hit.setMetaValue("string_list", StringList{"a", "b", "c"}); // STRING_LIST
  hit.setMetaValue("int_list", IntList{1, 2, 3});            // INT_LIST
  hit.setMetaValue("double_list", DoubleList{1.1, 2.2, 3.3}); // DOUBLE_LIST
  
  PeptideEvidence evidence;
  evidence.setProteinAccession("TEST_PROTEIN");
  hit.setPeptideEvidences(vector<PeptideEvidence>{evidence});
  
  peptide_id.setHits(vector<PeptideHit>{hit});
  peptide_ids.push_back(peptide_id);
  
  // Define meta value keys to export - test all types
  std::set<String> meta_keys = {"string_value", "int_value", "double_value", "string_list", "int_list", "double_list"};
  
  String output_file;
  NEW_TMP_FILE(output_file)
  
  // Store the data with meta values
  file.store(output_file, protein_ids, peptide_ids, false, meta_keys);
  
  // Read back the parquet file
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
  
  // Verify number of columns (19 standard + 6 meta value columns)
  TEST_EQUAL(table->num_columns(), 25)
  
  // Verify schema types for meta value columns
  auto schema = table->schema();
  
  // Check that string_value is utf8 type
  auto string_field = schema->GetFieldByName("string_value");
  TEST_NOT_EQUAL(string_field, nullptr)
  TEST_EQUAL(string_field->type()->id(), arrow::Type::STRING)
  
  // Check that int_value is int64 type
  auto int_field = schema->GetFieldByName("int_value");
  TEST_NOT_EQUAL(int_field, nullptr)  
  TEST_EQUAL(int_field->type()->id(), arrow::Type::INT64)
  
  // Check that double_value is float64 type
  auto double_field = schema->GetFieldByName("double_value");
  TEST_NOT_EQUAL(double_field, nullptr)
  TEST_EQUAL(double_field->type()->id(), arrow::Type::DOUBLE)
  
  // Check that string_list is list of string type
  auto string_list_field = schema->GetFieldByName("string_list");
  TEST_NOT_EQUAL(string_list_field, nullptr)
  TEST_EQUAL(string_list_field->type()->id(), arrow::Type::LIST)
  
  // Check that int_list is list of int64 type  
  auto int_list_field = schema->GetFieldByName("int_list");
  TEST_NOT_EQUAL(int_list_field, nullptr)
  TEST_EQUAL(int_list_field->type()->id(), arrow::Type::LIST)
  
  // Check that double_list is list of double type
  auto double_list_field = schema->GetFieldByName("double_list");
  TEST_NOT_EQUAL(double_list_field, nullptr)
  TEST_EQUAL(double_list_field->type()->id(), arrow::Type::LIST)
  
  // Verify actual data values
  auto string_column = table->GetColumnByName("string_value");
  auto string_array = std::static_pointer_cast<arrow::StringArray>(string_column->chunk(0));
  TEST_EQUAL(string_array->GetString(0), "test_string")
  
  auto int_column = table->GetColumnByName("int_value");
  auto int_array = std::static_pointer_cast<arrow::Int64Array>(int_column->chunk(0));
  TEST_EQUAL(int_array->Value(0), 42)
  
  auto double_column = table->GetColumnByName("double_value");
  auto double_array = std::static_pointer_cast<arrow::DoubleArray>(double_column->chunk(0));
  TEST_REAL_SIMILAR(double_array->Value(0), 3.14159)
}
END_SECTION

START_SECTION((Test meta value type conflict detection throws exception))
{
  QuantmsIO file;
  
  // Create test data with type conflicts for the same meta value key
  vector<ProteinIdentification> protein_ids;
  PeptideIdentificationList peptide_ids;
  
  ProteinIdentification protein_id;
  protein_id.setIdentifier("test_search");
  protein_id.setSearchEngine("TestEngine");
  protein_id.setScoreType("TestScore");
  protein_id.setHigherScoreBetter(true);
  protein_ids.push_back(protein_id);

  // Create two peptide identifications with different types for the same meta value key
  PeptideIdentification peptide_id1;
  peptide_id1.setIdentifier("test_search");
  peptide_id1.setRT(1234.5);
  peptide_id1.setMZ(500.25);
  peptide_id1.setScoreType("TestScore");
  
  PeptideHit hit1;
  hit1.setSequence(AASequence::fromString("PEPTIDER"));
  hit1.setScore(0.95);
  hit1.setCharge(2);
  hit1.setMetaValue("target_decoy", "target");
  hit1.setMetaValue("conflicting_value", 42);  // INT_VALUE
  
  PeptideEvidence evidence1;
  evidence1.setProteinAccession("TEST_PROTEIN_1");
  hit1.setPeptideEvidences(vector<PeptideEvidence>{evidence1});
  
  peptide_id1.setHits(vector<PeptideHit>{hit1});
  peptide_ids.push_back(peptide_id1);
  
  // Second peptide identification with different type for same key
  PeptideIdentification peptide_id2;
  peptide_id2.setIdentifier("test_search");
  peptide_id2.setRT(1334.5);
  peptide_id2.setMZ(600.25);
  peptide_id2.setScoreType("TestScore");
  
  PeptideHit hit2;
  hit2.setSequence(AASequence::fromString("ALTERNATIVE"));
  hit2.setScore(0.85);
  hit2.setCharge(2);
  hit2.setMetaValue("target_decoy", "target");
  hit2.setMetaValue("conflicting_value", String("text")); // STRING_VALUE - conflicts with INT_VALUE above
  
  PeptideEvidence evidence2;
  evidence2.setProteinAccession("TEST_PROTEIN_2");
  hit2.setPeptideEvidences(vector<PeptideEvidence>{evidence2});
  
  peptide_id2.setHits(vector<PeptideHit>{hit2});
  peptide_ids.push_back(peptide_id2);
  
  // Define meta value keys that include the conflicting key
  std::set<String> meta_keys = {"conflicting_value"};
  
  String output_file;
  NEW_TMP_FILE(output_file)
  
  // This should throw an InvalidParameter exception due to type conflict
  TEST_EXCEPTION(Exception::InvalidParameter, file.store(output_file, protein_ids, peptide_ids, false, meta_keys))
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST

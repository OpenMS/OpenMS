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

  std:vector<String> pep_strs = {"PEPTIDER", "PEM(Oxidation)TIDER", "DFPIANGER"};
  // Create multiple peptide identifications to test row count
  for (int i = 0; i < pep_strs.size(); ++i)
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
  
  // Verify number of columns (should match the quantms.io PSM schema)
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
  TEST_EQUAL(schema->field(9)->name(), "mp_accessions")
  TEST_EQUAL(schema->field(10)->name(), "predicted_rt")
  TEST_EQUAL(schema->field(11)->name(), "reference_file_name")
  TEST_EQUAL(schema->field(12)->name(), "cv_params")
  TEST_EQUAL(schema->field(13)->name(), "scan")
  TEST_EQUAL(schema->field(14)->name(), "rt")
  TEST_EQUAL(schema->field(15)->name(), "ion_mobility")
  TEST_EQUAL(schema->field(16)->name(), "num_peaks")
  TEST_EQUAL(schema->field(17)->name(), "mz_array")
  TEST_EQUAL(schema->field(18)->name(), "intensity_array")
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST

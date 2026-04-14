// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: OpenMS Team $
// $Authors: OpenMS Team $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FORMAT/BedRModFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/METADATA/ID/IdentificationData.h>
#include <OpenMS/CHEMISTRY/RibonucleotideDB.h>
#include <OpenMS/CHEMISTRY/NASequence.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(BedRModFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

BedRModFile* ptr = nullptr;
BedRModFile* nullPointer = nullptr;

START_SECTION((BedRModFile()))
{
  ptr = new BedRModFile();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((~BedRModFile()))
{
  delete ptr;
}
END_SECTION

START_SECTION((void store(const String& out_file, const IdentificationData& id_data, const String& chebi_mapping_file = "")))
{
  // Create test identification data
  IdentificationData id_data;
  
  // Add a parent sequence (RNA transcript)
  IdentificationData::ParentSequence parent;
  parent.accession = "test_rna";
  parent.sequence_type = IdentificationData::MoleculeType::RNA;
  parent.sequence = "AUGCAUGC"; // Simple test sequence
  auto parent_ref = id_data.registerParentSequence(parent);
  
  // Create a sequence with both unmodified and modified bases
  // Sequence: A-U-G-mC (where mC is methylcytidine)
  NASequence na_seq;
  auto ribo_db = RibonucleotideDB::getInstance();
  na_seq.append(ribo_db->getRibonucleotide("A"));
  na_seq.append(ribo_db->getRibonucleotide("U"));
  na_seq.append(ribo_db->getRibonucleotide("G"));
  na_seq.append(ribo_db->getRibonucleotide("mC")); // Modified cytidine
  
  // Create identified oligo
  IdentificationData::IdentifiedOligo oligo;
  oligo.sequence = na_seq;
  
  // Add parent match (positions 0-3 on parent)
  IdentificationData::ParentMatch parent_match;
  parent_match.start_pos = 0;
  parent_match.end_pos = 3;
  oligo.parent_matches[parent_ref].insert(parent_match);
  
  auto oligo_ref = id_data.registerIdentifiedOligo(oligo);
  
  // Create observation (spectrum)
  IdentificationData::InputFile input_file;
  input_file.name = "test.mzML";
  auto input_ref = id_data.registerInputFile(input_file);
  
  IdentificationData::Observation obs;
  obs.input_file = input_ref;
  auto obs_ref = id_data.registerObservation(obs);
  
  // Create observation match
  IdentificationData::ObservationMatch match;
  match.observation_ref = obs_ref;
  match.identified_molecule_var = oligo_ref;
  
  // Add a score
  auto score_ref = id_data.registerScoreType(IdentificationData::ScoreType("hyperscore"));
  match.addScore(score_ref, 100.0);
  
  id_data.registerObservationMatch(match);
  
  // Store to file
  String test_file = NEW_TMP_FILE();
  BedRModFile file;
  file.store(test_file, id_data);
  
  // Read and verify the output
  TextFile output;
  output.load(test_file);
  
  // Check header lines
  TEST_EQUAL(output[0], "#fileformat=bedRModv2")
  TEST_EQUAL(output[2], "#modification_type=RNA")
  
  // Count data lines (non-header)
  int data_lines = 0;
  bool found_unmodified = false;
  bool found_modified = false;
  
  for (Size i = 0; i < output.size(); ++i)
  {
    String line = output[i];
    if (line.empty() || line[0] == '#')
    {
      continue;
    }
    data_lines++;
    
    // Parse the line to check positions
    std::vector<String> fields;
    line.split('\t', fields);
    
    // Check if we have both modified and unmodified bases
    if (fields.size() > 3)
    {
      Int pos = fields[1].toInt();
      String chebi = fields[3];
      
      // Position 3 should be the modified base (mC)
      if (pos == 3 && chebi != "0")
      {
        found_modified = true;
      }
      // Positions 0, 1, 2 should be unmodified bases
      if (pos < 3)
      {
        found_unmodified = true;
      }
    }
  }
  
  // With the new behavior, we should have 4 data lines (A, U, G, mC)
  TEST_EQUAL(data_lines, 4)
  TEST_TRUE(found_unmodified) // Should find unmodified bases
  TEST_TRUE(found_modified)   // Should find modified base
}
END_SECTION

START_SECTION((void store - terminal modifications excluded))
{
  // Test that terminal modifications (5' and 3') are excluded
  IdentificationData id_data;
  
  // Add a parent sequence
  IdentificationData::ParentSequence parent;
  parent.accession = "test_rna_terminal";
  parent.sequence_type = IdentificationData::MoleculeType::RNA;
  parent.sequence = "AUGCAUGC";
  auto parent_ref = id_data.registerParentSequence(parent);
  
  // Create a sequence with a terminal modification
  NASequence na_seq;
  auto ribo_db = RibonucleotideDB::getInstance();
  
  // Add a 5' terminal phosphate
  na_seq.append(ribo_db->getRibonucleotide("5'-p"));
  
  // Add regular bases
  na_seq.append(ribo_db->getRibonucleotide("A"));
  na_seq.append(ribo_db->getRibonucleotide("U"));
  
  // Add a 3' terminal phosphate
  na_seq.append(ribo_db->getRibonucleotide("3'-p"));
  
  // Create identified oligo
  IdentificationData::IdentifiedOligo oligo;
  oligo.sequence = na_seq;
  
  // Add parent match
  IdentificationData::ParentMatch parent_match;
  parent_match.start_pos = 0;
  parent_match.end_pos = 3;
  oligo.parent_matches[parent_ref].insert(parent_match);
  
  auto oligo_ref = id_data.registerIdentifiedOligo(oligo);
  
  // Create observation
  IdentificationData::InputFile input_file;
  input_file.name = "test.mzML";
  auto input_ref = id_data.registerInputFile(input_file);
  
  IdentificationData::Observation obs;
  obs.input_file = input_ref;
  auto obs_ref = id_data.registerObservation(obs);
  
  // Create observation match
  IdentificationData::ObservationMatch match;
  match.observation_ref = obs_ref;
  match.identified_molecule_var = oligo_ref;
  
  auto score_ref = id_data.registerScoreType(IdentificationData::ScoreType("hyperscore"));
  match.addScore(score_ref, 100.0);
  
  id_data.registerObservationMatch(match);
  
  // Store to file
  String test_file = NEW_TMP_FILE();
  BedRModFile file;
  file.store(test_file, id_data);
  
  // Read and verify the output
  TextFile output;
  output.load(test_file);
  
  // Count data lines - should only have the 2 non-terminal bases (A and U)
  int data_lines = 0;
  for (Size i = 0; i < output.size(); ++i)
  {
    String line = output[i];
    if (line.empty() || line[0] == '#')
    {
      continue;
    }
    data_lines++;
  }
  
  // Should have 2 data lines (A and U), terminal modifications excluded
  TEST_EQUAL(data_lines, 2)
}
END_SECTION

END_TEST

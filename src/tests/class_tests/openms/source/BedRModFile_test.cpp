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
#include <OpenMS/DATASTRUCTURES/StringUtils.h>

using namespace OpenMS;
using namespace std;
namespace ID = IdentificationDataInternal;

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

START_SECTION((void store(const std::string& out_file, const IdentificationData& id_data, const std::string& chebi_mapping_file = "")))
{
  // Create test identification data
  IdentificationData id_data;

  // Add a parent sequence (RNA transcript)
  ID::ParentSequence parent("test_rna", ID::MoleculeType::RNA, "AUG[m5C]AUGC");
  auto parent_ref = id_data.registerParentSequence(parent);

  // Create a sequence with both unmodified and modified bases: AUG[m5C]
  NASequence na_seq = NASequence::fromString("AUG[m5C]");

  // Create identified oligo with parent match at positions 0-3
  ID::IdentifiedOligo oligo(na_seq);
  oligo.parent_matches[parent_ref].insert(ID::ParentMatch(0, 3));
  auto oligo_ref = id_data.registerIdentifiedOligo(oligo);

  // Create observation (spectrum)
  auto input_ref = id_data.registerInputFile(ID::InputFile("test.mzML"));
  ID::Observation obs("spectrum_1", input_ref, 100.0, 500.0);
  auto obs_ref = id_data.registerObservation(obs);

  // Create observation match with a hyperscore and a PSM-level q-value
  ID::ScoreType score_type("hyperscore", true);
  auto score_ref = id_data.registerScoreType(score_type);
  ID::ScoreType qvalue_type("PSM-level q-value", false);
  auto qvalue_ref = id_data.registerScoreType(qvalue_type);
  ID::ObservationMatch match(oligo_ref, obs_ref, 2);
  match.addScore(score_ref, 100.0);
  match.addScore(qvalue_ref, 0.01);
  id_data.registerObservationMatch(match);

  // Store to file
  std::string test_file;
  NEW_TMP_FILE(test_file);
  BedRModFile file;
  file.store(test_file, id_data);

  // Read and verify the output
  TextFile output;
  output.load(test_file);
  vector<std::string> lines(output.begin(), output.end());

  // Count data lines (non-header, non-empty)
  int data_lines = 0;
  bool found_mod_names_header = false;

  for (const auto& line : lines)
  {
    if (line.empty())
    {
      continue;
    }
    if (line.starts_with("#modification_names="))
    {
      // Should list all 4 bases including the modified one (m5C)
      found_mod_names_header = (line.find("m5C") != std::string::npos);
    }
    if (line[0] == '#')
    {
      continue;
    }
    data_lines++;
  }

  // With the new behavior, we should have 4 data lines (A, U, G, m5C)
  TEST_EQUAL(data_lines, 4)
  TEST_TRUE(found_mod_names_header) // Should list m5C in modification_names header
}
END_SECTION

START_SECTION((void store - terminal modifications excluded))
{
  // Test that terminal modifications (5' and 3') are excluded
  IdentificationData id_data;

  // Add a parent sequence
  ID::ParentSequence parent("test_rna_terminal", ID::MoleculeType::RNA, "AUGCAUGC");
  auto parent_ref = id_data.registerParentSequence(parent);

  // Create a sequence with terminal modifications: pAUp (5'-p + AU + 3'-p)
  NASequence na_seq = NASequence::fromString("pAUp");

  // Create identified oligo with parent match
  ID::IdentifiedOligo oligo(na_seq);
  oligo.parent_matches[parent_ref].insert(ID::ParentMatch(0, 1));
  auto oligo_ref = id_data.registerIdentifiedOligo(oligo);

  // Create observation
  auto input_ref = id_data.registerInputFile(ID::InputFile("test.mzML"));
  ID::Observation obs("spectrum_1", input_ref, 100.0, 500.0);
  auto obs_ref = id_data.registerObservation(obs);

  // Create observation match with a hyperscore and a PSM-level q-value
  ID::ScoreType score_type("hyperscore", true);
  auto score_ref = id_data.registerScoreType(score_type);
  ID::ScoreType qvalue_type("PSM-level q-value", false);
  auto qvalue_ref = id_data.registerScoreType(qvalue_type);
  ID::ObservationMatch match(oligo_ref, obs_ref, 2);
  match.addScore(score_ref, 100.0);
  match.addScore(qvalue_ref, 0.01);
  id_data.registerObservationMatch(match);

  // Store to file
  std::string test_file;
  NEW_TMP_FILE(test_file);
  BedRModFile file;
  file.store(test_file, id_data);

  // Read and verify the output
  TextFile output;
  output.load(test_file);
  vector<std::string> lines(output.begin(), output.end());

  // Count data lines - should only have the 2 non-terminal bases (A and U)
  int data_lines = 0;
  for (const auto& line : lines)
  {
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

START_SECTION((void store - matches without q-value are skipped))
{
  // Test that matches without a valid q-value score are excluded from output
  IdentificationData id_data;

  // Add a parent sequence
  ID::ParentSequence parent("test_rna_no_fdr", ID::MoleculeType::RNA, "AUGC");
  auto parent_ref = id_data.registerParentSequence(parent);

  // Create sequence
  NASequence na_seq = NASequence::fromString("AUG");

  // Create identified oligo
  ID::IdentifiedOligo oligo(na_seq);
  oligo.parent_matches[parent_ref].insert(ID::ParentMatch(0, 2));
  auto oligo_ref = id_data.registerIdentifiedOligo(oligo);

  // Create observation
  auto input_ref = id_data.registerInputFile(ID::InputFile("test.mzML"));
  ID::Observation obs("spectrum_1", input_ref, 100.0, 500.0);
  auto obs_ref = id_data.registerObservation(obs);

  // Create observation match with only a hyperscore (no q-value)
  ID::ScoreType score_type("hyperscore", true);
  auto score_ref = id_data.registerScoreType(score_type);
  ID::ObservationMatch match(oligo_ref, obs_ref, 2);
  match.addScore(score_ref, 100.0);
  // Deliberately NOT adding q-value score to simulate "FDR not run"
  id_data.registerObservationMatch(match);

  // Store to file
  std::string test_file;
  NEW_TMP_FILE(test_file);
  BedRModFile file;
  file.store(test_file, id_data);

  // Read and verify the output
  TextFile output;
  output.load(test_file);
  vector<std::string> lines(output.begin(), output.end());

  // Count data lines
  int data_lines = 0;
  for (const auto& line : lines)
  {
    if (line.empty() || line[0] == '#')
    {
      continue;
    }
    data_lines++;
  }

  // Should have 0 data lines since no q-value was provided
  TEST_EQUAL(data_lines, 0)
}
END_SECTION

START_SECTION((void store - sentinel q-values do not produce rows or affect frequency))
{
  IdentificationData id_data;
  auto parent_ref = id_data.registerParentSequence(
    ID::ParentSequence("test_rna_sentinel", ID::MoleculeType::RNA, "ACG"));

  ID::IdentifiedOligo valid_oligo(NASequence::fromString("ACG"));
  valid_oligo.parent_matches[parent_ref].insert(ID::ParentMatch(0, 2));
  auto valid_ref = id_data.registerIdentifiedOligo(valid_oligo);

  ID::IdentifiedOligo sentinel_oligo(NASequence::fromString("A[m5C]G"));
  sentinel_oligo.parent_matches[parent_ref].insert(ID::ParentMatch(0, 2));
  auto sentinel_ref = id_data.registerIdentifiedOligo(sentinel_oligo);

  auto input_ref = id_data.registerInputFile(ID::InputFile("test.mzML"));
  auto valid_obs = id_data.registerObservation(ID::Observation("valid", input_ref, 100.0, 500.0));
  auto sentinel_obs = id_data.registerObservation(ID::Observation("sentinel", input_ref, 101.0, 500.0));
  auto score_ref = id_data.registerScoreType(ID::ScoreType("hyperscore", true));
  auto qvalue_ref = id_data.registerScoreType(ID::ScoreType("PSM-level q-value", false));

  ID::ObservationMatch valid_match(valid_ref, valid_obs, 2);
  valid_match.addScore(score_ref, 100.0);
  valid_match.addScore(qvalue_ref, 0.01);
  id_data.registerObservationMatch(valid_match);

  ID::ObservationMatch sentinel_match(sentinel_ref, sentinel_obs, 2);
  sentinel_match.addScore(score_ref, 100.0);
  sentinel_match.addScore(qvalue_ref, -1.0);
  id_data.registerObservationMatch(sentinel_match);

  std::string test_file;
  NEW_TMP_FILE(test_file);
  BedRModFile().store(test_file, id_data);

  TextFile output;
  output.load(test_file);
  int data_lines = 0;
  for (const auto& line : output)
  {
    if (line.empty() || line[0] == '#')
    {
      continue;
    }
    vector<std::string> fields;
    StringUtils::split(line, '\t', fields);
    TEST_TRUE(fields.size() >= 11)
    if (fields.size() >= 11)
    {
      TEST_EQUAL(fields[4], "0.01")
      TEST_EQUAL(fields[10], "100")
    }
    ++data_lines;
  }
  TEST_EQUAL(data_lines, 3)
}
END_SECTION

START_SECTION((void store - target_mapping_count reflects target occurrences))
{
  // Test that target_mapping_count column shows count of target sequence occurrences
  IdentificationData id_data;

  // Add three parent sequences - 2 targets and 1 decoy
  ID::ParentSequence parent1("target_rna_1", ID::MoleculeType::RNA, "AUGCAUGC", "", 0.0, false);
  auto parent1_ref = id_data.registerParentSequence(parent1);

  ID::ParentSequence parent2("target_rna_2", ID::MoleculeType::RNA, "CAUGCAUG", "", 0.0, false);
  auto parent2_ref = id_data.registerParentSequence(parent2);

  ID::ParentSequence parent_decoy("DECOY_rna", ID::MoleculeType::RNA, "GCAUGCAU", "", 0.0, true);
  auto parent_decoy_ref = id_data.registerParentSequence(parent_decoy);

  // Create sequence that maps to multiple locations
  NASequence na_seq = NASequence::fromString("AUG");

  // Create identified oligo with matches to all three parents
  ID::IdentifiedOligo oligo(na_seq);
  oligo.parent_matches[parent1_ref].insert(ID::ParentMatch(0, 2));  // target 1
  oligo.parent_matches[parent1_ref].insert(ID::ParentMatch(4, 6));  // target 1 (second occurrence)
  oligo.parent_matches[parent2_ref].insert(ID::ParentMatch(1, 3));  // target 2
  oligo.parent_matches[parent_decoy_ref].insert(ID::ParentMatch(2, 4));  // decoy (should not count)
  auto oligo_ref = id_data.registerIdentifiedOligo(oligo);

  // Create observation
  auto input_ref = id_data.registerInputFile(ID::InputFile("test.mzML"));
  ID::Observation obs("spectrum_1", input_ref, 100.0, 500.0);
  auto obs_ref = id_data.registerObservation(obs);

  // Create observation match with scores
  ID::ScoreType score_type("hyperscore", true);
  auto score_ref = id_data.registerScoreType(score_type);
  ID::ScoreType qvalue_type("PSM-level q-value", false);
  auto qvalue_ref = id_data.registerScoreType(qvalue_type);
  ID::ObservationMatch match(oligo_ref, obs_ref, 2);
  match.addScore(score_ref, 100.0);
  match.addScore(qvalue_ref, 0.01);
  id_data.registerObservationMatch(match);

  // Store to file
  std::string test_file;
  NEW_TMP_FILE(test_file);
  BedRModFile file;
  file.store(test_file, id_data);

  // Read and verify the output
  TextFile output;
  output.load(test_file);
  vector<std::string> lines(output.begin(), output.end());

  // Check that target_mapping_count (column 12) equals 3 for all data lines
  // (2 occurrences in target 1 + 1 occurrence in target 2 = 3 total target occurrences)
  // The decoy occurrence should NOT be counted
  int data_lines_checked = 0;
  for (const auto& line : lines)
  {
    if (line.empty() || line[0] == '#')
    {
      continue;
    }
    
    // Parse the line (tab-separated)
    vector<std::string> fields;
    StringUtils::split(line, '\t', fields);
    
    // Check that we have enough fields
    TEST_TRUE(fields.size() >= 12)
    
    if (fields.size() >= 12)
    {
      // Column 12 (0-indexed: 11) is unique_mapping/target_mapping_count
      Int mapping_count = StringUtils::toInt32(fields[11]);
      
      // Should be 3 (only target occurrences)
      TEST_EQUAL(mapping_count, 3)
      data_lines_checked++;
    }
  }

  // Should have checked some data lines
  TEST_TRUE(data_lines_checked > 0)
}
END_SECTION

START_SECTION((void store - target_mapping_count equals 1 for unique mapping))
{
  // Test that target_mapping_count is 1 for a sequence with single target occurrence
  IdentificationData id_data;

  // Add one target parent sequence
  ID::ParentSequence parent("unique_target", ID::MoleculeType::RNA, "AUGC", "", 0.0, false);
  auto parent_ref = id_data.registerParentSequence(parent);

  // Create sequence
  NASequence na_seq = NASequence::fromString("AUG");

  // Create identified oligo with single match to one parent
  ID::IdentifiedOligo oligo(na_seq);
  oligo.parent_matches[parent_ref].insert(ID::ParentMatch(0, 2));
  auto oligo_ref = id_data.registerIdentifiedOligo(oligo);

  // Create observation
  auto input_ref = id_data.registerInputFile(ID::InputFile("test.mzML"));
  ID::Observation obs("spectrum_1", input_ref, 100.0, 500.0);
  auto obs_ref = id_data.registerObservation(obs);

  // Create observation match with scores
  ID::ScoreType score_type("hyperscore", true);
  auto score_ref = id_data.registerScoreType(score_type);
  ID::ScoreType qvalue_type("PSM-level q-value", false);
  auto qvalue_ref = id_data.registerScoreType(qvalue_type);
  ID::ObservationMatch match(oligo_ref, obs_ref, 2);
  match.addScore(score_ref, 100.0);
  match.addScore(qvalue_ref, 0.01);
  id_data.registerObservationMatch(match);

  // Store to file
  std::string test_file;
  NEW_TMP_FILE(test_file);
  BedRModFile file;
  file.store(test_file, id_data);

  // Read and verify the output
  TextFile output;
  output.load(test_file);
  vector<std::string> lines(output.begin(), output.end());

  // Check that target_mapping_count equals 1
  int data_lines_checked = 0;
  for (const auto& line : lines)
  {
    if (line.empty() || line[0] == '#')
    {
      continue;
    }
    
    vector<std::string> fields;
    StringUtils::split(line, '\t', fields);
    
    TEST_TRUE(fields.size() >= 12)
    
    if (fields.size() >= 12)
    {
      Int mapping_count = StringUtils::toInt32(fields[11]);
      TEST_EQUAL(mapping_count, 1)
      data_lines_checked++;
    }
  }

  TEST_TRUE(data_lines_checked > 0)
}
END_SECTION

START_SECTION((void store - decoy-only mappings have target_mapping_count of 0))
{
  // Test that sequences mapping only to decoys have target_mapping_count = 0
  IdentificationData id_data;

  // Add only decoy parent sequences
  ID::ParentSequence parent_decoy1("DECOY_rna_1", ID::MoleculeType::RNA, "AUGC", "", 0.0, true);
  auto parent_decoy1_ref = id_data.registerParentSequence(parent_decoy1);

  ID::ParentSequence parent_decoy2("DECOY_rna_2", ID::MoleculeType::RNA, "CAUGC", "", 0.0, true);
  auto parent_decoy2_ref = id_data.registerParentSequence(parent_decoy2);

  // Create sequence
  NASequence na_seq = NASequence::fromString("AUG");

  // Create identified oligo with matches only to decoys
  ID::IdentifiedOligo oligo(na_seq);
  oligo.parent_matches[parent_decoy1_ref].insert(ID::ParentMatch(0, 2));
  oligo.parent_matches[parent_decoy2_ref].insert(ID::ParentMatch(1, 3));
  auto oligo_ref = id_data.registerIdentifiedOligo(oligo);

  // Create observation
  auto input_ref = id_data.registerInputFile(ID::InputFile("test.mzML"));
  ID::Observation obs("spectrum_1", input_ref, 100.0, 500.0);
  auto obs_ref = id_data.registerObservation(obs);

  // Create observation match with scores
  ID::ScoreType score_type("hyperscore", true);
  auto score_ref = id_data.registerScoreType(score_type);
  ID::ScoreType qvalue_type("PSM-level q-value", false);
  auto qvalue_ref = id_data.registerScoreType(qvalue_type);
  ID::ObservationMatch match(oligo_ref, obs_ref, 2);
  match.addScore(score_ref, 100.0);
  match.addScore(qvalue_ref, 0.01);
  id_data.registerObservationMatch(match);

  // Store to file
  std::string test_file;
  NEW_TMP_FILE(test_file);
  BedRModFile file;
  file.store(test_file, id_data);

  // Read and verify the output
  TextFile output;
  output.load(test_file);
  vector<std::string> lines(output.begin(), output.end());

  // Check that target_mapping_count equals 0 for decoy-only mappings
  int data_lines_checked = 0;
  for (const auto& line : lines)
  {
    if (line.empty() || line[0] == '#')
    {
      continue;
    }
    
    vector<std::string> fields;
    StringUtils::split(line, '\t', fields);
    
    TEST_TRUE(fields.size() >= 12)
    
    if (fields.size() >= 12)
    {
      Int mapping_count = StringUtils::toInt32(fields[11]);
      TEST_EQUAL(mapping_count, 0)
      data_lines_checked++;
    }
  }

  TEST_TRUE(data_lines_checked > 0)
}
END_SECTION

END_TEST

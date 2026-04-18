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

START_SECTION((void store(const String& out_file, const IdentificationData& id_data, const String& chebi_mapping_file = "")))
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
  String test_file;
  NEW_TMP_FILE(test_file);
  BedRModFile file;
  file.store(test_file, id_data);

  // Read and verify the output
  TextFile output;
  output.load(test_file);
  vector<String> lines(output.begin(), output.end());

  // Count data lines (non-header, non-empty)
  int data_lines = 0;
  bool found_mod_names_header = false;

  for (const auto& line : lines)
  {
    if (line.empty())
    {
      continue;
    }
    if (line.hasPrefix("#modification_names="))
    {
      // Should list all 4 bases including the modified one (m5C)
      found_mod_names_header = (line.find("m5C") != String::npos);
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
  String test_file;
  NEW_TMP_FILE(test_file);
  BedRModFile file;
  file.store(test_file, id_data);

  // Read and verify the output
  TextFile output;
  output.load(test_file);
  vector<String> lines(output.begin(), output.end());

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

END_TEST

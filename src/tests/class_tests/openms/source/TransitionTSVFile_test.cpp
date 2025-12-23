// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>

#include <boost/assign/std/vector.hpp>

///////////////////////////
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionTSVFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>
///////////////////////////

#include <fstream>
#include <sstream>
#include <algorithm>

using namespace OpenMS;
using namespace std;

// Helper function to read TSV file into a vector of lines (sorted for comparison)
std::vector<std::string> readAndSortTSV(const std::string& filename)
{
  std::vector<std::string> lines;
  std::ifstream file(filename);
  std::string line;

  // Skip header
  std::getline(file, line);

  // Read all data lines
  while (std::getline(file, line))
  {
    if (!line.empty())
    {
      lines.push_back(line);
    }
  }

  // Sort for order-independent comparison
  std::sort(lines.begin(), lines.end());
  return lines;
}

// Helper function to parse a TSV line into fields
std::vector<std::string> parseTSVLine(const std::string& line)
{
  std::vector<std::string> fields;
  std::stringstream ss(line);
  std::string field;
  while (std::getline(ss, field, '\t'))
  {
    fields.push_back(field);
  }
  return fields;
}

// Helper to normalize NA/empty for comparison
std::string normalizeField(const std::string& field)
{
  if (field.empty() || field == "NA")
  {
    return "NA";
  }
  return field;
}

START_TEST(TransitionTSVFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

TransitionTSVFile* ptr = nullptr;
TransitionTSVFile* nullPointer = nullptr;

START_SECTION(TransitionTSVFile())
{
  ptr = new TransitionTSVFile();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~TransitionTSVFile())
{
  delete ptr;
}
END_SECTION

START_SECTION( void convertTargetedExperimentToTSV(const char * filename, OpenMS::TargetedExperiment & targeted_exp))
{
  // Test roundtrip: TraML -> TargetedExperiment -> TSV
  std::string input_file = OPENMS_GET_TEST_DATA_PATH("MRMAssay_detectingTransistionCompound_input.TraML");
  std::string output_file;
  NEW_TMP_FILE(output_file);

  TraMLFile traml;
  TargetedExperiment targeted_exp;
  traml.load(input_file, targeted_exp);

  TransitionTSVFile tsv_writer;
  tsv_writer.convertTargetedExperimentToTSV(output_file.c_str(), targeted_exp);

  // Verify output file exists and has content
  std::ifstream ifs(output_file);
  TEST_EQUAL(ifs.good(), true)

  std::string line;
  int line_count = 0;
  while (std::getline(ifs, line)) { line_count++; }
  TEST_EQUAL(line_count > 1, true) // Header + at least one data line
}
END_SECTION

START_SECTION( void convertTSVToTargetedExperiment(const char * filename, OpenMS::TargetedExperiment & targeted_exp))
{
  // Load from TraML and write to TSV, then read TSV back
  std::string input_file = OPENMS_GET_TEST_DATA_PATH("MRMAssay_detectingTransistionCompound_input.TraML");
  std::string tsv_file;
  NEW_TMP_FILE(tsv_file);

  TraMLFile traml;
  TargetedExperiment targeted_exp_in;
  traml.load(input_file, targeted_exp_in);

  TransitionTSVFile tsv_handler;
  tsv_handler.convertTargetedExperimentToTSV(tsv_file.c_str(), targeted_exp_in);

  // Now read the TSV back
  TargetedExperiment targeted_exp_out;
  tsv_handler.convertTSVToTargetedExperiment(tsv_file.c_str(), FileTypes::TSV, targeted_exp_out);

  // Verify we loaded some transitions
  TEST_EQUAL(targeted_exp_out.getTransitions().size() > 0, true)
  // Note: Peptides may not be loaded if the TraML only contains compounds
  TEST_EQUAL(targeted_exp_out.getCompounds().size() > 0 || targeted_exp_out.getPeptides().size() > 0, true)
}
END_SECTION

START_SECTION( void validateTargetedExperiment(OpenMS::TargetedExperiment & targeted_exp))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION([EXTRA] Light path TSV roundtrip)
{
  // Test roundtrip: TraML -> Heavy -> TSV -> Light -> TSV
  std::string input_file = OPENMS_GET_TEST_DATA_PATH("MRMAssay_detectingTransistionCompound_input.TraML");
  std::string tsv_file, output_file;
  NEW_TMP_FILE(tsv_file);
  NEW_TMP_FILE(output_file);

  // Load TraML to heavy, convert to TSV
  TraMLFile traml;
  TargetedExperiment targeted_exp;
  traml.load(input_file, targeted_exp);

  TransitionTSVFile tsv_handler;
  tsv_handler.convertTargetedExperimentToTSV(tsv_file.c_str(), targeted_exp);

  // Read TSV into LightTargetedExperiment
  OpenSwath::LightTargetedExperiment light_exp;
  tsv_handler.convertTSVToTargetedExperiment(tsv_file.c_str(), FileTypes::TSV, light_exp);

  // Verify we loaded some data
  TEST_EQUAL(light_exp.transitions.size() > 0, true)
  TEST_EQUAL(light_exp.compounds.size() > 0, true)

  // Write back to TSV
  tsv_handler.convertLightTargetedExperimentToTSV(output_file.c_str(), light_exp);

  // Verify output file exists and has content
  std::ifstream ifs(output_file);
  TEST_EQUAL(ifs.good(), true)

  std::string line;
  int line_count = 0;
  while (std::getline(ifs, line)) { line_count++; }
  TEST_EQUAL(line_count > 1, true)
}
END_SECTION

START_SECTION([EXTRA] Light vs Heavy path equivalence via PQP roundtrip)
{
  // Test that light and heavy paths produce equivalent output when both go through PQP
  // This tests: TraML -> Heavy -> PQP -> Heavy -> TSV vs TraML -> Light -> PQP -> Light -> TSV

  std::string input_file = OPENMS_GET_TEST_DATA_PATH("MRMAssay_detectingTransistionCompound_input.TraML");
  std::string heavy_pqp, light_pqp, heavy_tsv, light_tsv;
  NEW_TMP_FILE(heavy_pqp);
  NEW_TMP_FILE(light_pqp);
  NEW_TMP_FILE(heavy_tsv);
  NEW_TMP_FILE(light_tsv);

  TransitionTSVFile tsv_file;
  TransitionPQPFile pqp_file;
  TraMLFile traml;

  // Load TraML
  TargetedExperiment targeted_exp;
  traml.load(input_file, targeted_exp);

  // Heavy path: TargetedExperiment -> PQP -> TargetedExperiment -> TSV
  {
    TargetedExperiment heavy_exp_out;
    pqp_file.convertTargetedExperimentToPQP(heavy_pqp.c_str(), targeted_exp);
    pqp_file.convertPQPToTargetedExperiment(heavy_pqp.c_str(), heavy_exp_out);
    tsv_file.convertTargetedExperimentToTSV(heavy_tsv.c_str(), heavy_exp_out);
  }

  // Light path: TargetedExperiment -> LightTargetedExperiment -> PQP -> LightTargetedExperiment -> TSV
  {
    OpenSwath::LightTargetedExperiment light_exp_in, light_exp_out;
    OpenSwathDataAccessHelper::convertTargetedExp(targeted_exp, light_exp_in);
    pqp_file.convertLightTargetedExperimentToPQP(light_pqp.c_str(), light_exp_in);
    pqp_file.convertPQPToTargetedExperiment(light_pqp.c_str(), light_exp_out);
    tsv_file.convertLightTargetedExperimentToTSV(light_tsv.c_str(), light_exp_out);
  }

  // Compare outputs (sorted for order-independent comparison)
  std::vector<std::string> heavy_lines = readAndSortTSV(heavy_tsv);
  std::vector<std::string> light_lines = readAndSortTSV(light_tsv);

  // Should have same number of transitions
  TEST_EQUAL(heavy_lines.size(), light_lines.size())

  // Compare field by field for each line
  // Fields that should match after PQP roundtrip (indices based on TSV header):
  // 0: PrecursorMz, 1: ProductMz, 2: PrecursorCharge, 3: ProductCharge
  // 4: LibraryIntensity, 5: NormalizedRetentionTime, 6: PeptideSequence
  // 7: ModifiedPeptideSequence, 8: PeptideGroupLabel
  // 22: TransitionGroupId, 23: TransitionId
  // 24: Decoy, 25: DetectingTransition, 26: IdentifyingTransition, 27: QuantifyingTransition
  //
  // Note: Fields like FragmentType (17) may differ between empty and "NA"
  // due to different handling - we normalize these for comparison

  std::vector<size_t> numeric_fields = {0, 1, 4, 5};
  std::vector<size_t> string_fields = {2, 3, 6, 7, 8, 22, 23, 24, 25, 26, 27};

  for (size_t i = 0; i < std::min(heavy_lines.size(), light_lines.size()); ++i)
  {
    std::vector<std::string> heavy_fields = parseTSVLine(heavy_lines[i]);
    std::vector<std::string> light_fields = parseTSVLine(light_lines[i]);

    // Compare numeric fields
    for (size_t field_idx : numeric_fields)
    {
      if (field_idx < heavy_fields.size() && field_idx < light_fields.size())
      {
        double heavy_val = std::stod(heavy_fields[field_idx]);
        double light_val = std::stod(light_fields[field_idx]);
        TEST_REAL_SIMILAR(heavy_val, light_val)
      }
    }

    // Compare string fields (normalize empty/"NA" values)
    for (size_t field_idx : string_fields)
    {
      if (field_idx < heavy_fields.size() && field_idx < light_fields.size())
      {
        std::string heavy_norm = normalizeField(heavy_fields[field_idx]);
        std::string light_norm = normalizeField(light_fields[field_idx]);
        TEST_STRING_EQUAL(heavy_norm, light_norm)
      }
    }
  }
}
END_SECTION

START_SECTION([EXTRA] Light path preserves transition flags)
{
  // Test that detecting/quantifying/identifying flags are preserved in light path
  std::string input_file = OPENMS_GET_TEST_DATA_PATH("MRMAssay_detectingTransistionCompound_input.TraML");
  std::string pqp_file_path;
  NEW_TMP_FILE(pqp_file_path);

  TraMLFile traml;
  TransitionPQPFile pqp_file;
  TargetedExperiment targeted_exp;
  OpenSwath::LightTargetedExperiment light_exp_in, light_exp_out;

  // Load TraML and convert to light
  traml.load(input_file, targeted_exp);
  OpenSwathDataAccessHelper::convertTargetedExp(targeted_exp, light_exp_in);

  // Store original flag values
  std::vector<bool> original_detecting, original_quantifying, original_identifying, original_decoy;
  for (const auto& tr : light_exp_in.transitions)
  {
    original_detecting.push_back(tr.isDetectingTransition());
    original_quantifying.push_back(tr.isQuantifyingTransition());
    original_identifying.push_back(tr.isIdentifyingTransition());
    original_decoy.push_back(tr.getDecoy());
  }

  // Write to PQP and read back
  pqp_file.convertLightTargetedExperimentToPQP(pqp_file_path.c_str(), light_exp_in);
  pqp_file.convertPQPToTargetedExperiment(pqp_file_path.c_str(), light_exp_out);

  // Verify flags are preserved
  TEST_EQUAL(light_exp_out.transitions.size(), original_detecting.size())

  for (size_t i = 0; i < light_exp_out.transitions.size(); ++i)
  {
    TEST_EQUAL(light_exp_out.transitions[i].isDetectingTransition(), original_detecting[i])
    TEST_EQUAL(light_exp_out.transitions[i].isQuantifyingTransition(), original_quantifying[i])
    TEST_EQUAL(light_exp_out.transitions[i].isIdentifyingTransition(), original_identifying[i])
    TEST_EQUAL(light_exp_out.transitions[i].getDecoy(), original_decoy[i])
  }
}
END_SECTION

START_SECTION([EXTRA] Light path preserves fragment annotation)
{
  // Test that fragment type/number/charge are preserved in light path
  std::string input_file = OPENMS_GET_TEST_DATA_PATH("MRMAssay_detectingTransistionCompound_input.TraML");
  std::string pqp_file_path;
  NEW_TMP_FILE(pqp_file_path);

  TraMLFile traml;
  TransitionPQPFile pqp_file;
  TargetedExperiment targeted_exp;
  OpenSwath::LightTargetedExperiment light_exp_in, light_exp_out;

  // Load TraML and convert to light
  traml.load(input_file, targeted_exp);
  OpenSwathDataAccessHelper::convertTargetedExp(targeted_exp, light_exp_in);

  // Store original fragment info
  std::vector<std::string> original_annotations;
  for (const auto& tr : light_exp_in.transitions)
  {
    original_annotations.push_back(tr.getAnnotation());
  }

  // Write to PQP and read back
  pqp_file.convertLightTargetedExperimentToPQP(pqp_file_path.c_str(), light_exp_in);
  pqp_file.convertPQPToTargetedExperiment(pqp_file_path.c_str(), light_exp_out);

  // Verify annotations are preserved
  TEST_EQUAL(light_exp_out.transitions.size(), original_annotations.size())

  for (size_t i = 0; i < light_exp_out.transitions.size(); ++i)
  {
    TEST_STRING_EQUAL(light_exp_out.transitions[i].getAnnotation(), original_annotations[i])
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

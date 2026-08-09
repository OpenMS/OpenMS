// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathLibraryPreparation.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionTSVFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>
#include <OpenMS/SYSTEM/File.h>

#include <algorithm>
#include <string>
#include <vector>

using namespace OpenMS;
using namespace std;

namespace
{
  std::string toppDataPath_(const std::string& filename)
  {
    return File::absolutePath(std::string(OPENMS_GET_TEST_DATA_PATH("")) + "../../../topp/" + filename);
  }

  std::vector<std::string> sortedLines_(const std::string& filename)
  {
    TextFile text_file;
    text_file.load(filename);
    std::vector<std::string> lines(text_file.begin(), text_file.end());
    std::sort(lines.begin(), lines.end());
    return lines;
  }

  void testSortedFilesEqual_(const std::string& actual_file, const std::string& expected_file)
  {
    const auto actual_lines = sortedLines_(actual_file);
    const auto expected_lines = sortedLines_(expected_file);
    TEST_EQUAL(actual_lines.size(), expected_lines.size())
    for (Size i = 0; i < actual_lines.size(); ++i)
    {
      TEST_EQUAL(actual_lines[i], expected_lines[i])
    }
  }

  OpenSwathLibraryPreparation::AssayGeneratorParameters makeIPFTestParameters_()
  {
    OpenSwathLibraryPreparation::AssayGeneratorParameters params;
    params.allowed_fragment_charges = {1, 2, 3, 4};
    params.enable_ipf = true;
    params.test_mode = true;
    params.unimod_file = toppDataPath_("OpenSwathAssayGenerator_input_2_unimod.xml");
    return params;
  }

  OpenSwathLibraryPreparation::DecoyGeneratorParameters makeDeterministicDecoyParameters_()
  {
    OpenSwathLibraryPreparation::DecoyGeneratorParameters params;
    params.method = "pseudo-reverse";
    params.min_decoy_fraction = 0.4;
    params.enable_detection_specific_losses = true;
    params.enable_detection_unspecific_losses = true;
    return params;
  }
}

START_TEST(OpenSwathLibraryPreparation, "$Id$")

OpenSwathLibraryPreparation* ptr = nullptr;
OpenSwathLibraryPreparation* null_ptr = nullptr;

START_SECTION(OpenSwathLibraryPreparation())
{
  ptr = new OpenSwathLibraryPreparation();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~OpenSwathLibraryPreparation())
{
  delete ptr;
}
END_SECTION

START_SECTION([EXTRA] normalizeLibraryToPQP normalizes a prepared decoy-containing TSV library to PQP)
{
  OpenSwathLibraryPreparation prep;
  prep.setLogType(ProgressLogger::NONE);

  std::string output_pqp;
  NEW_TMP_FILE(output_pqp);
  output_pqp += ".pqp";

  const auto stats = prep.normalizeLibraryToPQP(
    toppDataPath_("OpenSwathDecoyGenerator_output_6_light.tsv"),
    FileTypes::TSV,
    output_pqp);

  TEST_TRUE(File::exists(output_pqp))
  TEST_TRUE(stats.transition_count > 0)
  TEST_TRUE(stats.compound_count > 0)
  TEST_TRUE(stats.hasDecoys())

  TransitionPQPFile pqp_reader;
  OpenSwath::LightTargetedExperiment light_exp;
  pqp_reader.convertPQPToTargetedExperiment(output_pqp.c_str(), light_exp, true);

  const Size decoy_count = static_cast<Size>(std::count_if(
    light_exp.transitions.begin(), light_exp.transitions.end(),
    [](const auto& transition)
    {
      return transition.getDecoy();
    }));

  TEST_EQUAL(light_exp.transitions.size(), stats.transition_count)
  TEST_EQUAL(light_exp.compounds.size(), stats.compound_count)
  TEST_EQUAL(decoy_count, stats.decoy_transition_count)
}
END_SECTION

START_SECTION([EXTRA] prepareAssays remains deterministic for IPF test-mode output through the shared helper)
{
  OpenSwathLibraryPreparation prep;
  prep.setLogType(ProgressLogger::NONE);

  std::string output_pqp_1;
  NEW_TMP_FILE(output_pqp_1);
  output_pqp_1 += ".pqp";

  std::string output_pqp_2;
  NEW_TMP_FILE(output_pqp_2);
  output_pqp_2 += ".pqp";

  const auto stats_1 = prep.prepareAssays(
    toppDataPath_("OpenSwathAssayGenerator_input_4.pqp"),
    FileTypes::PQP,
    output_pqp_1,
    FileTypes::PQP,
    makeIPFTestParameters_());

  const auto stats_2 = prep.prepareAssays(
    toppDataPath_("OpenSwathAssayGenerator_input_4.pqp"),
    FileTypes::PQP,
    output_pqp_2,
    FileTypes::PQP,
    makeIPFTestParameters_());

  TEST_TRUE(File::exists(output_pqp_1))
  TEST_TRUE(File::exists(output_pqp_2))
  TEST_TRUE(stats_1.transition_count > 0)
  TEST_TRUE(stats_1.identifying_transition_count > 0)
  TEST_EQUAL(stats_1.transition_count, stats_2.transition_count)
  TEST_EQUAL(stats_1.identifying_transition_count, stats_2.identifying_transition_count)
  TEST_EQUAL(stats_1.compound_count, stats_2.compound_count)
  TEST_EQUAL(stats_1.protein_count, stats_2.protein_count)

  TransitionPQPFile pqp_reader;
  TargetedExperiment targeted_exp_1;
  TargetedExperiment targeted_exp_2;
  pqp_reader.convertPQPToTargetedExperiment(output_pqp_1.c_str(), targeted_exp_1, true);
  pqp_reader.convertPQPToTargetedExperiment(output_pqp_2.c_str(), targeted_exp_2, true);

  std::string output_tsv_1;
  NEW_TMP_FILE(output_tsv_1);
  std::string output_tsv_2;
  NEW_TMP_FILE(output_tsv_2);
  TransitionTSVFile tsv_writer;
  tsv_writer.convertTargetedExperimentToTSV(output_tsv_1.c_str(), targeted_exp_1);
  tsv_writer.convertTargetedExperimentToTSV(output_tsv_2.c_str(), targeted_exp_2);

  testSortedFilesEqual_(output_tsv_1, output_tsv_2);
}
END_SECTION

START_SECTION([EXTRA] prepareEmpiricalLibraryToPQP runs assay preparation plus decoy generation and remains deterministic with deterministic decoys)
{
  OpenSwathLibraryPreparation prep;
  prep.setLogType(ProgressLogger::NONE);

  auto assay_params = makeIPFTestParameters_();
  assay_params.enable_ipf = false;
  assay_params.unimod_file.clear();
  const auto decoy_params = makeDeterministicDecoyParameters_();

  std::string output_pqp_1;
  NEW_TMP_FILE(output_pqp_1);
  output_pqp_1 += ".pqp";

  std::string output_pqp_2;
  NEW_TMP_FILE(output_pqp_2);
  output_pqp_2 += ".pqp";

  const auto stats_1 = prep.prepareEmpiricalLibraryToPQP(
    toppDataPath_("OpenSwathDecoyGenerator_input_4.tsv"),
    FileTypes::TSV,
    output_pqp_1,
    assay_params,
    decoy_params,
    Param(),
    File::path(output_pqp_1));

  const auto stats_2 = prep.prepareEmpiricalLibraryToPQP(
    toppDataPath_("OpenSwathDecoyGenerator_input_4.tsv"),
    FileTypes::TSV,
    output_pqp_2,
    assay_params,
    decoy_params,
    Param(),
    File::path(output_pqp_2));

  TEST_TRUE(File::exists(output_pqp_1))
  TEST_TRUE(File::exists(output_pqp_2))
  TEST_TRUE(stats_1.transition_count > 0)
  TEST_TRUE(stats_1.compound_count > 0)
  TEST_TRUE(stats_1.hasDecoys())
  TEST_EQUAL(stats_1.transition_count, stats_2.transition_count)
  TEST_EQUAL(stats_1.decoy_transition_count, stats_2.decoy_transition_count)

  TransitionPQPFile pqp_reader;
  TargetedExperiment targeted_exp_1;
  TargetedExperiment targeted_exp_2;
  pqp_reader.convertPQPToTargetedExperiment(output_pqp_1.c_str(), targeted_exp_1, true);
  pqp_reader.convertPQPToTargetedExperiment(output_pqp_2.c_str(), targeted_exp_2, true);

  std::string output_tsv_1;
  NEW_TMP_FILE(output_tsv_1);
  std::string output_tsv_2;
  NEW_TMP_FILE(output_tsv_2);
  TransitionTSVFile tsv_writer;
  tsv_writer.convertTargetedExperimentToTSV(output_tsv_1.c_str(), targeted_exp_1);
  tsv_writer.convertTargetedExperimentToTSV(output_tsv_2.c_str(), targeted_exp_2);

  testSortedFilesEqual_(output_tsv_1, output_tsv_2);
}
END_SECTION

START_SECTION([EXTRA] prepareEmpiricalLibraryToPQP preserves decoy flags when heavy TraML fallback is normalized to PQP)
{
  OpenSwathLibraryPreparation prep;
  prep.setLogType(ProgressLogger::NONE);

  auto assay_params = makeIPFTestParameters_();
  assay_params.enable_ipf = false;
  assay_params.unimod_file.clear();

  auto decoy_params = makeDeterministicDecoyParameters_();
  decoy_params.min_decoy_fraction = 0.0;

  std::string output_pqp;
  NEW_TMP_FILE(output_pqp);
  output_pqp += ".pqp";

  const auto stats = prep.prepareEmpiricalLibraryToPQP(
    toppDataPath_("OpenSwathWorkflow_1_input.TraML"),
    FileTypes::TRAML,
    output_pqp,
    assay_params,
    decoy_params);

  TEST_TRUE(File::exists(output_pqp))
  TEST_TRUE(stats.transition_count > 0)
  TEST_TRUE(stats.hasDecoys())

  TransitionPQPFile pqp_reader;
  OpenSwath::LightTargetedExperiment light_exp;
  pqp_reader.convertPQPToTargetedExperiment(output_pqp.c_str(), light_exp, true);

  const Size decoy_count = static_cast<Size>(std::count_if(
    light_exp.transitions.begin(), light_exp.transitions.end(),
    [](const auto& transition)
    {
      return transition.getDecoy();
    }));

  TEST_TRUE(decoy_count > 0)
  TEST_EQUAL(decoy_count, stats.decoy_transition_count)
}
END_SECTION

END_TEST

// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/config.h>

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathLibraryPreparation.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionTSVFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/TempFiles.h>

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

  std::string classTestDataPath_(const std::string& filename)
  {
    return File::absolutePath(std::string(OPENMS_GET_TEST_DATA_PATH("")) + filename);
  }

  std::vector<std::string> sortedLines_(const std::string& filename)
  {
    TextFile text_file;
    text_file.load(filename);
    std::vector<std::string> lines(text_file.begin(), text_file.end());
    std::sort(lines.begin(), lines.end());
    return lines;
  }

  // The scratch directory must hold nothing but what the caller put there: no
  // leftover intermediate and no nested directory (which would lengthen paths).
  void testScratchDirectoryUntouched_(const std::string& scratch_dir, const std::string& expected_file)
  {
    StringList files;
    File::fileList(scratch_dir, "*", files);
    TEST_EQUAL(files.size(), 1)
    if (!files.empty())
    {
      TEST_EQUAL(files[0], expected_file)
    }
    TEST_EQUAL(File::listDirectories(scratch_dir).size(), 0)
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

START_SECTION([EXTRA] ensureUnimodLoaded initializes the requested IPF modification source before library parsing)
{
  OpenSwathLibraryPreparation prep;
  prep.setLogType(ProgressLogger::NONE);
  prep.ensureUnimodLoaded(makeIPFTestParameters_());
}
END_SECTION

START_SECTION([EXTRA] normalizeLibraryToPQP normalizes a prepared decoy-containing TSV library to PQP)
{
  OpenSwathLibraryPreparation prep;
  prep.setLogType(ProgressLogger::NONE);

  std::string output_pqp;
  NEW_TMP_FILE_EXT(output_pqp, ".pqp");

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
  NEW_TMP_FILE_EXT(output_pqp_1, ".pqp");

  std::string output_pqp_2;
  NEW_TMP_FILE_EXT(output_pqp_2, ".pqp");

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

START_SECTION([EXTRA] prepareAssays rejects fallback UIS windows that cannot be constructed)
{
  OpenSwathLibraryPreparation prep;
  prep.setLogType(ProgressLogger::NONE);

  std::string output_pqp;
  NEW_TMP_FILE_EXT(output_pqp, ".pqp");

  auto params = makeIPFTestParameters_();
  params.enable_swath_specifity = false;
  params.precursor_mz_threshold = 0.0;
  TEST_EXCEPTION_WITH_MESSAGE(Exception::InvalidParameter,
    prep.prepareAssays(toppDataPath_("OpenSwathAssayGenerator_input_4.pqp"), FileTypes::PQP, output_pqp, FileTypes::PQP, params),
    "AssayGeneratorParameters::precursor_mz_threshold must be greater than zero when constructing fallback UIS SWATH windows.")

  params = makeIPFTestParameters_();
  params.enable_swath_specifity = false;
  params.precursor_lower_mz_limit = 1200.0;
  params.precursor_upper_mz_limit = 400.0;
  TEST_EXCEPTION_WITH_MESSAGE(Exception::InvalidParameter,
    prep.prepareAssays(toppDataPath_("OpenSwathAssayGenerator_input_4.pqp"), FileTypes::PQP, output_pqp, FileTypes::PQP, params),
    "AssayGeneratorParameters::precursor_upper_mz_limit must not be below precursor_lower_mz_limit when constructing fallback UIS SWATH windows.")
}
END_SECTION

#ifndef WITH_ONNX
START_SECTION([EXTRA] preparePredictedLibraryToPQP requires an ONNX-enabled build)
{
  OpenSwathLibraryPreparation prep;
  prep.setLogType(ProgressLogger::NONE);
  TEST_EXCEPTION(Exception::Precondition,
    prep.preparePredictedLibraryToPQP("unused.fasta", "unused.pqp",
                                      OpenSwathLibraryPreparation::AssayGeneratorParameters(),
                                      OpenSwathLibraryPreparation::DecoyGeneratorParameters(),
                                      OpenSwathLibraryPreparation::PredictedLibraryParameters()))
}
END_SECTION
#endif

START_SECTION([EXTRA] prepareEmpiricalLibraryToPQP runs assay preparation plus decoy generation and remains deterministic with deterministic decoys)
{
  OpenSwathLibraryPreparation prep;
  prep.setLogType(ProgressLogger::NONE);

  auto assay_params = makeIPFTestParameters_();
  assay_params.enable_ipf = false;
  assay_params.unimod_file.clear();
  const auto decoy_params = makeDeterministicDecoyParameters_();

  // A caller-provided scratch directory is a shared parent, not an invocation-owned
  // workspace. A pre-existing prepared_assays.pqp must never be overwritten or removed.
  TempDir shared_scratch_parent;
  const std::string scratch_sentinel = shared_scratch_parent.getPath() + "/prepared_assays.pqp";
  TextFile sentinel;
  sentinel.addLine("must survive");
  sentinel.store(scratch_sentinel);

  std::string output_pqp_1;
  NEW_TMP_FILE_EXT(output_pqp_1, ".pqp");

  std::string output_pqp_2;
  NEW_TMP_FILE_EXT(output_pqp_2, ".pqp");

  // Use a fixture that remains decoyable after assay preparation at the normal
  // 40% coverage gate. The former TSV fixture relied on the now-removed raw-library
  // fallback after assay filtering, so it is intentionally covered by the separate
  // fail-closed regression below instead of serving as a positive-path fixture.
  const auto stats_1 = prep.prepareEmpiricalLibraryToPQP(
    classTestDataPath_("MRMDecoyGenerator_input.TraML"),
    FileTypes::TRAML,
    output_pqp_1,
    assay_params,
    decoy_params,
    Param(),
    shared_scratch_parent.getPath());

  const auto sentinel_after_first = sortedLines_(scratch_sentinel);
  TEST_EQUAL(sentinel_after_first.size(), 1)
  TEST_EQUAL(sentinel_after_first[0], "must survive")

  const auto stats_2 = prep.prepareEmpiricalLibraryToPQP(
    classTestDataPath_("MRMDecoyGenerator_input.TraML"),
    FileTypes::TRAML,
    output_pqp_2,
    assay_params,
    decoy_params,
    Param(),
    shared_scratch_parent.getPath());

  const auto sentinel_after_second = sortedLines_(scratch_sentinel);
  TEST_EQUAL(sentinel_after_second.size(), 1)
  TEST_EQUAL(sentinel_after_second[0], "must survive")
  testScratchDirectoryUntouched_(shared_scratch_parent.getPath(), "prepared_assays.pqp");

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

START_SECTION([EXTRA] prepareEmpiricalLibraryToPQP fails closed when assay preparation yields zero transitions)
{
  OpenSwathLibraryPreparation prep;
  prep.setLogType(ProgressLogger::NONE);

  auto assay_params = makeIPFTestParameters_();
  assay_params.enable_ipf = false;
  assay_params.unimod_file.clear();
  assay_params.min_transitions = 1000;

  const auto decoy_params = makeDeterministicDecoyParameters_();

  TempDir shared_scratch_parent;
  const std::string scratch_sentinel = shared_scratch_parent.getPath() + "/prepared_assays.pqp";
  TextFile sentinel;
  sentinel.addLine("must survive failure");
  sentinel.store(scratch_sentinel);

  std::string output_pqp;
  NEW_TMP_FILE_EXT(output_pqp, ".pqp");

  TEST_EXCEPTION(Exception::Precondition,
    prep.prepareEmpiricalLibraryToPQP(
      classTestDataPath_("MRMDecoyGenerator_input.TraML"),
      FileTypes::TRAML,
      output_pqp,
      assay_params,
      decoy_params,
      Param(),
      shared_scratch_parent.getPath()))
  TEST_FALSE(File::exists(output_pqp))

  const auto sentinel_after_failure = sortedLines_(scratch_sentinel);
  TEST_EQUAL(sentinel_after_failure.size(), 1)
  TEST_EQUAL(sentinel_after_failure[0], "must survive failure")
  testScratchDirectoryUntouched_(shared_scratch_parent.getPath(), "prepared_assays.pqp");
}
END_SECTION

START_SECTION([EXTRA] prepareEmpiricalLibraryToPQP fails closed when decoy generation yields zero decoys)
{
  OpenSwathLibraryPreparation prep;
  prep.setLogType(ProgressLogger::NONE);

  auto assay_params = makeIPFTestParameters_();
  assay_params.enable_ipf = false;
  assay_params.unimod_file.clear();

  // Disable the coverage gate and request no decoys, so decoy generation succeeds and
  // writes a target-only library before the zero-decoy check rejects it.
  auto decoy_params = makeDeterministicDecoyParameters_();
  decoy_params.min_decoy_fraction = 0.0;
  decoy_params.aim_decoy_fraction = 0.0;

  TempDir scratch_dir;
  std::string output_pqp;
  NEW_TMP_FILE_EXT(output_pqp, ".pqp");

  TEST_EXCEPTION(Exception::Precondition,
    prep.prepareEmpiricalLibraryToPQP(
      classTestDataPath_("MRMDecoyGenerator_input.TraML"),
      FileTypes::TRAML,
      output_pqp,
      assay_params,
      decoy_params,
      Param(),
      scratch_dir.getPath()))
  TEST_FALSE(File::exists(output_pqp))

  StringList leftover_files;
  File::fileList(scratch_dir.getPath(), "*", leftover_files);
  TEST_EQUAL(leftover_files.size(), 0)
}
END_SECTION


START_SECTION([EXTRA] generateDecoys preserves the historical light-path requirement for protein annotations)
{
  OpenSwathLibraryPreparation prep;
  prep.setLogType(ProgressLogger::NONE);

  TransitionTSVFile tsv;
  OpenSwath::LightTargetedExperiment light_exp;
  tsv.convertTSVToTargetedExperiment(
    toppDataPath_("OpenSwathWorkflow_23_input.tsv").c_str(),
    FileTypes::TSV,
    light_exp);
  light_exp.proteins.clear();
  for (auto& compound : light_exp.compounds)
  {
    compound.protein_refs.clear();
  }

  std::string proteinless_tsv;
  NEW_TMP_FILE_EXT(proteinless_tsv, ".tsv");
  tsv.convertLightTargetedExperimentToTSV(proteinless_tsv.c_str(), light_exp);

  std::string output_tsv;
  NEW_TMP_FILE_EXT(output_tsv, ".tsv");

  auto decoy_params = makeDeterministicDecoyParameters_();
  decoy_params.min_decoy_fraction = 0.0;

  TEST_EXCEPTION(Exception::IllegalArgument,
    prep.generateDecoys(
      proteinless_tsv,
      FileTypes::TSV,
      output_tsv,
      FileTypes::TSV,
      decoy_params))
}
END_SECTION

START_SECTION([EXTRA] generateDecoys preserves decoy flags when heavy TraML is normalized to PQP)
{
  OpenSwathLibraryPreparation prep;
  prep.setLogType(ProgressLogger::NONE);

  const auto decoy_params = makeDeterministicDecoyParameters_();

  std::string output_pqp;
  NEW_TMP_FILE_EXT(output_pqp, ".pqp");

  const auto stats = prep.generateDecoys(
    classTestDataPath_("MRMDecoyGenerator_input.TraML"),
    FileTypes::TRAML,
    output_pqp,
    FileTypes::PQP,
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

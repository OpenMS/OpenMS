// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/test_config.h>

#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/FORMAT/DATAACCESS/MSChromatogramParquetConsumer.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/FORMAT/XICParquetFile.h>
#include <OpenMS/SYSTEM/File.h>

using namespace OpenMS;
using namespace std;

START_TEST(MSChromatogramParquetConsumer, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(MSChromatogramParquetConsumer_basic)
{
  // Load transitions and chromatograms from TOPP test data
  TargetedExperiment targeted_exp;
  TraMLFile().load(OPENMS_GET_TEST_DATA_PATH("MSChromatogramParquetConsumer_1_input.TraML"), targeted_exp);

  OpenSwath::LightTargetedExperiment light_exp;
  OpenSwathDataAccessHelper::convertTargetedExp(targeted_exp, light_exp);

  MSExperiment exp;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MSChromatogramParquetConsumer_1_output.chrom.mzML"), exp);
  TEST_EQUAL(exp.getChromatograms().empty(), false)

  std::string tmp;
  NEW_TMP_FILE(tmp);
  std::string out = tmp + ".xic";
  {
    MSChromatogramParquetConsumer consumer(out, 1, "test_source", light_exp);
    for (auto& chrom : exp.getChromatograms())
    {
      consumer.consumeChromatogram(chrom);
    }
  }

  TEST_EQUAL(File::exists(out), true)

  XICParquetFile xic(out);
  std::vector<XICChromatogram> chroms;
  xic.load(chroms);
  TEST_EQUAL(chroms.size(), exp.getChromatograms().size())
  TEST_EQUAL(chroms[0].rt.empty(), false)
  TEST_EQUAL(chroms[0].intensity.empty(), false)
}
END_SECTION

START_SECTION(MSChromatogramParquetConsumer_empty_chromatograms)
{
  TargetedExperiment targeted_exp;
  TraMLFile().load(OPENMS_GET_TEST_DATA_PATH("MSChromatogramParquetConsumer_1_input.TraML"), targeted_exp);

  OpenSwath::LightTargetedExperiment light_exp;
  OpenSwathDataAccessHelper::convertTargetedExp(targeted_exp, light_exp);

  std::string tmp;
  NEW_TMP_FILE(tmp);
  std::string out = tmp + ".xic";
  {
    MSChromatogramParquetConsumer consumer(out, 1, "test_source", light_exp);
    // Do not consume any chromatograms
  }

  XICParquetFile xic(out);
  std::vector<XICChromatogram> chroms;
  xic.load(chroms);
  TEST_EQUAL(chroms.size(), 0)
}
END_SECTION

START_SECTION(MSChromatogramParquetConsumer_trailing_empty_chromatograms_are_written)
{
  OpenSwath::LightTargetedExperiment light_exp;

  OpenSwath::LightCompound compound;
  compound.id = "pep1";
  compound.sequence = "PEPTIDE";
  compound.charge = 2;
  light_exp.compounds.push_back(compound);

  OpenSwath::LightTransition transition;
  transition.transition_name = "tr1";
  transition.peptide_ref = "pep1";
  transition.precursor_mz = 600.2;
  transition.product_mz = 500.2;
  transition.fragment_charge = 1;
  transition.fragment_nr = 7;
  transition.setFragmentType("y");
  transition.setDetectingTransition(true);
  light_exp.transitions.push_back(transition);

  MSChromatogram nonempty_chrom;
  nonempty_chrom.setNativeID("tr1");
  ChromatogramPeak peak;
  peak.setRT(100.0);
  peak.setIntensity(1000.0);
  nonempty_chrom.push_back(peak);

  MSChromatogram empty_chrom;
  empty_chrom.setNativeID("tr1");

  std::string tmp;
  NEW_TMP_FILE(tmp);
  std::string out = tmp + ".xic";
  {
    MSChromatogramParquetConsumer consumer(out, 1, "test_source", light_exp);
    consumer.consumeChromatogram(nonempty_chrom);
    consumer.consumeChromatogram(empty_chrom);
  }

  XICParquetFile xic(out);
  std::vector<XICChromatogram> chroms;
  xic.load(chroms);
  TEST_EQUAL(chroms.size(), 2)
  TEST_EQUAL(chroms[0].rt.empty(), false)
  TEST_EQUAL(chroms[1].rt.empty(), true)
  TEST_EQUAL(chroms[1].intensity.empty(), true)
}
END_SECTION

START_SECTION(MSChromatogramParquetConsumer_destructor_no_throw)
{
  TargetedExperiment targeted_exp;
  TraMLFile().load(OPENMS_GET_TEST_DATA_PATH("MSChromatogramParquetConsumer_1_input.TraML"), targeted_exp);

  OpenSwath::LightTargetedExperiment light_exp;
  OpenSwathDataAccessHelper::convertTargetedExp(targeted_exp, light_exp);

  std::string out = File::getTempDirectory() + "/openms_missing_dir/xic_out.xic";
  bool caught = false;
  try
  {
    MSChromatogramParquetConsumer consumer(out, 1, "test_source", light_exp);
  }
  catch (...)
  {
    caught = true;
  }
  TEST_EQUAL(caught, false)
}
END_SECTION

START_SECTION(MSChromatogramParquetConsumer_finalize_surfaces_errors)
{
  // issue #9488 (FORM-129): unlike the destructor (which must not throw and only logs),
  // calling finalize() explicitly surfaces write errors as exceptions so callers can react
  // during normal control flow. OpenSwathWorkflow now calls finalize() before delete.
  TargetedExperiment targeted_exp;
  TraMLFile().load(OPENMS_GET_TEST_DATA_PATH("MSChromatogramParquetConsumer_1_input.TraML"), targeted_exp);

  OpenSwath::LightTargetedExperiment light_exp;
  OpenSwathDataAccessHelper::convertTargetedExp(targeted_exp, light_exp);

  // Unwritable path (non-existent directory) -> finalize()/write_() must throw, NOT swallow.
  std::string out = File::getTempDirectory() + "/openms_missing_dir/xic_finalize.xic";
  MSChromatogramParquetConsumer consumer(out, 1, "test_source", light_exp);
  TEST_EXCEPTION(Exception::FileNotWritable, consumer.finalize())
}
END_SECTION

START_SECTION(MSChromatogramParquetConsumer_mixed_target_and_decoy_identifying_transitions_keep_single_target_analyte)
{
  OpenSwath::LightTargetedExperiment light_exp;

  OpenSwath::LightCompound compound;
  compound.id = "2407";
  compound.sequence = "TSHVENDYIDNPSLALT(UniMod:21)TGPK";
  compound.charge = 3;
  light_exp.compounds.push_back(compound);

  OpenSwath::LightTransition target_transition;
  target_transition.transition_name = "tr_target";
  target_transition.peptide_ref = "2407";
  target_transition.precursor_mz = 784.6967;
  target_transition.product_mz = 326.1459;
  target_transition.fragment_charge = 1;
  target_transition.fragment_nr = 3;
  target_transition.setFragmentType("b");
  target_transition.setDetectingTransition(true);
  target_transition.setIdentifyingTransition(false);
  target_transition.setDecoy(false);
  light_exp.transitions.push_back(target_transition);

  OpenSwath::LightTransition identifying_decoy = target_transition;
  identifying_decoy.transition_name = "tr_ident_decoy";
  identifying_decoy.product_mz = 302.1710;
  identifying_decoy.setDetectingTransition(false);
  identifying_decoy.setIdentifyingTransition(true);
  identifying_decoy.setDecoy(true);
  light_exp.transitions.push_back(identifying_decoy);

  MSChromatogram target_chrom;
  target_chrom.setNativeID("tr_target");
  ChromatogramPeak peak;
  peak.setRT(100.0);
  peak.setIntensity(1000.0);
  target_chrom.push_back(peak);

  MSChromatogram decoy_chrom;
  decoy_chrom.setNativeID("tr_ident_decoy");
  decoy_chrom.push_back(peak);

  std::string tmp;
  NEW_TMP_FILE(tmp);
  std::string out = tmp + ".xic";
  {
    MSChromatogramParquetConsumer consumer(out, 1, "test_source", light_exp);
    consumer.consumeChromatogram(target_chrom);
    consumer.consumeChromatogram(decoy_chrom);
  }

  XICParquetFile xic(out);
  std::vector<XICChromatogram> chroms;
  xic.load(chroms);
  TEST_EQUAL(chroms.size(), 2)
  for (const auto& chrom : chroms)
  {
    TEST_EQUAL(chrom.has_precursor_decoy, true)
    TEST_EQUAL(chrom.precursor_decoy, 0)
  }

  std::vector<XICAnalyte> analytes;
  xic.getAnalytes(analytes);
  TEST_EQUAL(analytes.size(), 1)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

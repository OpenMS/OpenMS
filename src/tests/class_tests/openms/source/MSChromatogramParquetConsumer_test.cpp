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

  String tmp;
  NEW_TMP_FILE(tmp);
  String out = tmp + ".xic";
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

  String tmp;
  NEW_TMP_FILE(tmp);
  String out = tmp + ".xic";
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

START_SECTION(MSChromatogramParquetConsumer_destructor_no_throw)
{
  TargetedExperiment targeted_exp;
  TraMLFile().load(OPENMS_GET_TEST_DATA_PATH("MSChromatogramParquetConsumer_1_input.TraML"), targeted_exp);

  OpenSwath::LightTargetedExperiment light_exp;
  OpenSwathDataAccessHelper::convertTargetedExp(targeted_exp, light_exp);

  String out = File::getTempDirectory() + "/openms_missing_dir/xic_out.xic";
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

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

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

#include <OpenMS/FORMAT/DATAACCESS/MobilogramParquetConsumer.h>
#include <OpenMS/KERNEL/Mobilogram.h>
#include <OpenMS/KERNEL/MobilityPeak1D.h>
#include <OpenMS/SYSTEM/File.h>

using namespace OpenMS;
using namespace std;

START_TEST(MobilogramParquetConsumer, "$Id$")

START_SECTION(MobilogramParquetConsumer_basic)
{
#ifdef WITH_PARQUET
  // create a small mobilogram
  Mobilogram m;
  m.setRT(123.45);
  m.push_back(MobilityPeak1D(1.0, 10.0));
  m.push_back(MobilityPeak1D(2.0, 20.0));

  OpenSwath::LightTargetedExperiment light_exp; // empty

  String tmp;
  NEW_TMP_FILE(tmp);
  String out = tmp + ".xim";
  {
    MobilogramParquetConsumer consumer(out, 1, "test_source", light_exp);
    consumer.consumeMobilogram(m);
  }

  TEST_EQUAL(File::exists(out), true)
#else
  OpenSwath::LightTargetedExperiment light_exp;
  TEST_EXCEPTION(Exception::NotImplemented, MobilogramParquetConsumer("dummy.xim", 1, "x", light_exp))
#endif
}
END_SECTION

START_SECTION(MobilogramParquetConsumer_empty_mobilograms)
{
#ifdef WITH_PARQUET
  OpenSwath::LightTargetedExperiment light_exp;

  String tmp;
  NEW_TMP_FILE(tmp);
  String out = tmp + ".xim";
  {
    MobilogramParquetConsumer consumer(out, 1, "test_source", light_exp);
    // Do not consume any mobilograms
  }

  TEST_EQUAL(File::exists(out), true)
#endif
}
END_SECTION

START_SECTION(MobilogramParquetConsumer_destructor_no_throw)
{
#ifdef WITH_PARQUET
  OpenSwath::LightTargetedExperiment light_exp;

  String out = File::getTempDirectory() + "/openms_missing_dir/xim_out.xim";
  bool caught = false;
  try
  {
    MobilogramParquetConsumer consumer(out, 1, "test_source", light_exp);
  }
  catch (...)
  {
    caught = true;
  }
  TEST_EQUAL(caught, false)
#endif
}
END_SECTION

END_TEST

// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

#include <OpenMS/config.h>

#ifdef WITH_THERMO_RAW

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/ThermoRawFile.h>
#include <OpenMS/SYSTEM/File.h>

using namespace OpenMS;
using namespace std;

START_TEST(ThermoRawFile, "$Id$")

START_SECTION(void load(const String& path, MSExperiment& exp))
{
  ThermoRawFile file;
  MSExperiment exp;

  TEST_EXCEPTION(Exception::FileNotFound, file.load("/definitely/not/present.raw", exp))
}
END_SECTION

#ifdef THERMO_RAW_TEST_DATA
START_SECTION(round-trip load raw -> mzML -> reload MSExperiment)
{
  ThermoRawFile file;
  MSExperiment original;
  file.load(THERMO_RAW_TEST_DATA, original);

  TEST_NOT_EQUAL(original.size(), 0)
  TEST_EQUAL(original.getSourceFiles().empty(), false)

  String tmp_mzml = File::getTempDirectory() + "/" + File::getUniqueName() + "_thermo_roundtrip.mzML";
  MzMLFile().store(tmp_mzml, original);

  MSExperiment reloaded;
  MzMLFile().load(tmp_mzml, reloaded);
  File::remove(tmp_mzml);

  TEST_EQUAL(original.size(), reloaded.size())
  TEST_EQUAL(original.getSourceFiles().size(), reloaded.getSourceFiles().size())

  if (!original.empty())
  {
    TEST_EQUAL(original[0].getMSLevel(), reloaded[0].getMSLevel())
    TEST_REAL_SIMILAR(original[0].getRT(), reloaded[0].getRT())
    TEST_EQUAL(original[0].size(), reloaded[0].size())
  }

  bool found_msn_precursor = false;
  for (Size i = 0; i < reloaded.size(); ++i)
  {
    TEST_EQUAL(original[i].getMSLevel(), reloaded[i].getMSLevel())
    if (reloaded[i].getMSLevel() > 1 && !reloaded[i].getPrecursors().empty())
    {
      TEST_NOT_EQUAL(reloaded[i].getPrecursors()[0].getMZ(), 0.0)
      found_msn_precursor = true;
      break;
    }
  }
  TEST_EQUAL(found_msn_precursor, true)
}
END_SECTION
#endif

END_TEST

#endif // WITH_THERMO_RAW
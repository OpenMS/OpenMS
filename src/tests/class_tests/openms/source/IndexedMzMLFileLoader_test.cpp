// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FORMAT/IndexedMzMLFileLoader.h>
#include <OpenMS/KERNEL/OnDiscMSExperiment.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(IndexedMzMLFileLoader, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

IndexedMzMLFileLoader* ptr = nullptr;
IndexedMzMLFileLoader* nullPointer = nullptr;
START_SECTION((IndexedMzMLFileLoader()))
  ptr = new IndexedMzMLFileLoader;
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~IndexedMzMLFileLoader()))
  delete ptr;
END_SECTION

START_SECTION(const PeakFileOptions& getOptions() const)
  IndexedMzMLFileLoader file;
  const PeakFileOptions& options = file.getOptions();
  TEST_EQUAL(options.hasMSLevels(),false)
END_SECTION

START_SECTION(PeakFileOptions& getOptions())
  IndexedMzMLFileLoader file;
  file.getOptions().addMSLevel(1);
  TEST_EQUAL(file.getOptions().hasMSLevels(),true);
END_SECTION

START_SECTION(void setOptions(const PeakFileOptions &))
  IndexedMzMLFileLoader file;
  const PeakFileOptions& options = file.getOptions();
  TEST_EQUAL(options.hasMSLevels(),false)
  TEST_EQUAL(file.getOptions().hasMSLevels(),false);
  PeakFileOptions new_options(options);
  new_options.addMSLevel(1);
  file.setOptions(new_options);
  TEST_EQUAL(file.getOptions().hasMSLevels(),true);
END_SECTION

TOLERANCE_ABSOLUTE(0.01)

START_SECTION(bool load(const String& filename, OnDiscPeakMap& exp))
{
  IndexedMzMLFileLoader file;
  OnDiscPeakMap exp;
  file.load(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"),exp);

  PeakMap exp2;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"),exp2);

  TEST_EQUAL(exp.getNrSpectra(), exp2.getSpectra().size())
  TEST_EQUAL(exp.getNrChromatograms(), exp2.getChromatograms().size())
  TEST_EQUAL(exp.getNrSpectra(), 2)
  TEST_EQUAL(exp.getNrChromatograms(), 1)
  TEST_EQUAL(exp.getSpectrum(0) == exp2.getSpectra()[0], true)

  for (Size i = 0; i < exp.getNrSpectra(); i++)
  {
    TEST_EQUAL(exp.getSpectrum(i) == exp2.getSpectra()[i], true)
  }
  for (Size i = 0; i < exp.getNrChromatograms(); i++)
  {
    TEST_EQUAL(exp.getChromatogram(i) == exp2.getChromatograms()[i], true)
  }

  TEST_EQUAL(*exp.getExperimentalSettings() == (OpenMS::ExperimentalSettings)exp2, true)
}
END_SECTION

START_SECTION([EXTRA]CheckParsing)
{
  // Check return value of load
  IndexedMzMLFileLoader file;
  OnDiscPeakMap exp;
  bool success;
  success = file.load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"), exp);
  TEST_EQUAL(success, false)
  success = file.load(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"),exp);
  TEST_EQUAL(success, true)
}
END_SECTION

START_SECTION(void store(const String& filename, OnDiscPeakMap& exp))
{
  IndexedMzMLFileLoader file;
  OnDiscPeakMap exp, exp_;
  file.load(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"),exp_);
  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  std::cout << "Storing in file " << tmp_filename << std::endl;
  file.store(tmp_filename,exp_);

  bool success = file.load(tmp_filename,exp);
  TEST_EQUAL(success, true)

  PeakMap exp2;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"),exp2);

  TEST_EQUAL(exp.getNrSpectra(), exp2.getSpectra().size())
  TEST_EQUAL(exp.getNrChromatograms(), exp2.getChromatograms().size())
  TEST_EQUAL(exp.getNrSpectra(), 2)
  TEST_EQUAL(exp.getNrChromatograms(), 1)
  TEST_EQUAL(exp.getSpectrum(0) == exp2.getSpectra()[0], true)

  for (Size i = 0; i < exp.getNrSpectra(); i++)
  {
    TEST_EQUAL(exp.getSpectrum(i) == exp2.getSpectra()[i], true)
  }
  for (Size i = 0; i < exp.getNrChromatograms(); i++)
  {
    TEST_EQUAL(exp.getChromatogram(i) == exp2.getChromatograms()[i], true)
  }

  TEST_EQUAL(*exp.getExperimentalSettings() == (OpenMS::ExperimentalSettings)exp2, true)
}
END_SECTION

START_SECTION(void store(const String& filename, PeakMap& exp))
{
  IndexedMzMLFileLoader file;
  OnDiscPeakMap exp;
  PeakMap exp2;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"),exp2);

  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  std::cout << "Storing in file " << tmp_filename << std::endl;
  file.store(tmp_filename,exp2);

  bool success = file.load(tmp_filename,exp);
  TEST_EQUAL(success, true)

  TEST_EQUAL(exp.getNrSpectra(), exp2.getSpectra().size())
  TEST_EQUAL(exp.getNrChromatograms(), exp2.getChromatograms().size())
  TEST_EQUAL(exp.getNrSpectra(), 2)
  TEST_EQUAL(exp.getNrChromatograms(), 1)
  TEST_EQUAL(exp.getSpectrum(0) == exp2.getSpectra()[0], true)

  for (Size i = 0; i < exp.getNrSpectra(); i++)
  {
    TEST_EQUAL(exp.getSpectrum(i) == exp2.getSpectra()[i], true)
  }
  for (Size i = 0; i < exp.getNrChromatograms(); i++)
  {
    TEST_EQUAL(exp.getChromatogram(i) == exp2.getChromatograms()[i], true)
  }

  TEST_EQUAL(*exp.getExperimentalSettings() == (OpenMS::ExperimentalSettings)exp2, true)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


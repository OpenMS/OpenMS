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
#include <OpenMS/FORMAT/HANDLERS/MzMLSqliteSwathHandler.h>
///////////////////////////

#include <OpenMS/FORMAT/SqMassFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <QFile>

using namespace OpenMS;
using namespace OpenMS::Internal;
using namespace std;

///////////////////////////

START_TEST(SqMassFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MzMLSqliteSwathHandler* ptr = nullptr;
MzMLSqliteSwathHandler* nullPointer = nullptr;

std::string tmp_filename;
NEW_TMP_FILE(tmp_filename);
MSExperiment exp_orig;
MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("SwathFile.mzML"), exp_orig);
SqMassFile file;
file.store(tmp_filename, exp_orig);

START_SECTION((MzMLSqliteSwathHandler()))
  ptr = new MzMLSqliteSwathHandler(tmp_filename);
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~MzMLSqliteSwathHandler()))
  delete ptr;
END_SECTION

TOLERANCE_RELATIVE(1.0005)

START_SECTION(std::vector<OpenSwath::SwathMap> readSwathWindows())
{
  MzMLSqliteSwathHandler handler(tmp_filename);

  std::vector<OpenSwath::SwathMap> maps = handler.readSwathWindows();
  TEST_EQUAL(maps.size(), 5)

  TEST_EQUAL(maps[0].ms1, false)
  TEST_REAL_SIMILAR(maps[0].lower, 400.0)
  TEST_REAL_SIMILAR(maps[0].center, 412.5)
  TEST_REAL_SIMILAR(maps[0].upper, 425.0)
  TEST_REAL_SIMILAR(maps[1].lower, 425.0)
  TEST_REAL_SIMILAR(maps[1].upper, 450.0)
  TEST_REAL_SIMILAR(maps[4].lower, 500.0)
  TEST_REAL_SIMILAR(maps[4].upper, 525.0)
}
END_SECTION

START_SECTION(std::vector<int> readMS1Spectra())
{
  MzMLSqliteSwathHandler handler(tmp_filename);

  TEST_EQUAL(handler.readMS1Spectra().size(), 19)
  TEST_EQUAL(handler.readMS1Spectra()[0], 0)
  TEST_EQUAL(handler.readMS1Spectra()[18], 108)
}
END_SECTION

START_SECTION(std::vector<int> readSpectraForWindow(const OpenSwath::SwathMap& swath_map))
{
  MzMLSqliteSwathHandler handler(tmp_filename);

  std::vector<OpenSwath::SwathMap> maps = handler.readSwathWindows();
  TEST_EQUAL(maps.size(), 5)

  TEST_EQUAL(handler.readSpectraForWindow(maps[0]).size(), 19)
  TEST_EQUAL(handler.readSpectraForWindow(maps[0])[0], 1)
  TEST_EQUAL(handler.readSpectraForWindow(maps[0])[18], 109)
  TEST_EQUAL(handler.readSpectraForWindow(maps[1]).size(), 19)
  TEST_EQUAL(handler.readSpectraForWindow(maps[1])[0], 2)
  TEST_EQUAL(handler.readSpectraForWindow(maps[1])[18], 110)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


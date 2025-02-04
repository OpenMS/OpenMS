// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FORMAT/MS2File.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/MSExperiment.h>

using namespace OpenMS;
using namespace std;

///////////////////////////
///////////////////////////

START_TEST(MS2File, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


MS2File * ptr = nullptr;
MS2File* nullPointer = nullptr;
START_SECTION((MS2File()))
ptr = new MS2File;
TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~MS2File()))
delete ptr;
END_SECTION

TOLERANCE_ABSOLUTE(0.01)

START_SECTION((template <typename MapType> void load(const String &filename, MapType & exp)))
MS2File file;
PeakMap exp;
file.load(OPENMS_GET_TEST_DATA_PATH("MS2File_test_spectra.ms2"), exp);

//test DocumentIdentifier addition
TEST_STRING_EQUAL(exp.getLoadedFilePath(), OPENMS_GET_TEST_DATA_PATH("MS2File_test_spectra.ms2"));
TEST_STRING_EQUAL(FileTypes::typeToName(exp.getLoadedFileType()), "ms2");

TEST_EQUAL(exp.size(), 2)

TEST_EQUAL(exp[0].size(), 4)
TEST_EQUAL(exp[1].size(), 4)

TEST_STRING_EQUAL(exp[0].getNativeID(), "index=0")
TEST_STRING_EQUAL(exp[1].getNativeID(), "index=1")

TEST_REAL_SIMILAR(exp[0].getPrecursors()[0].getMZ(), 444.44)
TEST_REAL_SIMILAR(exp[1].getPrecursors()[0].getMZ(), 555.555)

END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

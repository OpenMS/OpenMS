// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

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


MS2File* ptr = 0;
MS2File* nullPointer = 0;
START_SECTION((MS2File()))
	ptr = new MS2File;
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~MS2File()))
	delete ptr;
END_SECTION

TOLERANCE_ABSOLUTE(0.01)

START_SECTION((template <typename MapType> void load(const String& filename, MapType& exp)))
	MS2File file;
	PeakMap exp;
	file.load(OPENMS_GET_TEST_DATA_PATH("MS2File_test_spectra.ms2"), exp);

	//test DocumentIdentifier addition
	TEST_STRING_EQUAL(exp.getLoadedFilePath(), OPENMS_GET_TEST_DATA_PATH("MS2File_test_spectra.ms2"));
	TEST_STRING_EQUAL(FileHandler::typeToName(exp.getLoadedFileType()),"ms2");

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


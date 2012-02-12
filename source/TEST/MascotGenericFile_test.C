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
#include <OpenMS/FORMAT/MascotGenericFile.h>
#include <sstream>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MascotGenericFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MascotGenericFile* ptr = 0;
MascotGenericFile* nullPointer = 0;
START_SECTION(MascotGenericFile())
{
	ptr = new MascotGenericFile();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(virtual ~MascotGenericFile())
{
	delete ptr;
}
END_SECTION

ptr = new MascotGenericFile();

START_SECTION((template < typename MapType > void load(const String &filename, MapType &exp)))
{
  PeakMap exp;
  ptr->load(OPENMS_GET_TEST_DATA_PATH("MascotInfile_test.mascot_in"), exp);
  TEST_EQUAL(exp.size(), 1)

  TEST_EQUAL(exp.begin()->size(), 9)
}
END_SECTION

START_SECTION((void store(std::ostream &os, const String &filename, const PeakMap &experiment)))
{
	PeakMap exp;
	ptr->load(OPENMS_GET_TEST_DATA_PATH("MascotInfile_test.mascot_in"), exp);
	
	stringstream ss;
	ptr->store(ss, "bla", exp);

	// BEGIN IONS
	// TITLE=Testtitle
	// PEPMASS=1998
	// RTINSECONDS=25.379
	// 1 1
	// 2 4
	// 3 9
	// 4 16
	// 5 25
	// 6 36
	// 7 49
	// 8 64
	// 9 81
	// END IONS
	vector<String> strings;
	strings.push_back("BEGIN IONS");
	strings.push_back("TITLE=Testtitle");
	strings.push_back("PEPMASS=1998");
	strings.push_back("RTINSECONDS=25.37");
	strings.push_back("1 1");
	strings.push_back("2 4");
	strings.push_back("3 9");
	strings.push_back("4 16");
	strings.push_back("5 25");
	strings.push_back("6 36");
	strings.push_back("7 49");
	strings.push_back("8 64");
	strings.push_back("9 81");
	strings.push_back("END IONS");

	String mgf_file(ss.str());
	for (Size i = 0; i != strings.size(); ++i)
	{
		TEST_EQUAL(mgf_file.hasSubstring(mgf_file[i]), true)
	}	
}
END_SECTION

START_SECTION((void store(const String &filename, const PeakMap &experiment)))
{
  String tmp_name("MascotGenericFile_1.tmp");
	NEW_TMP_FILE(tmp_name)
	PeakMap exp;
  ptr->load(OPENMS_GET_TEST_DATA_PATH("MascotInfile_test.mascot_in"), exp);


	ptr->store(tmp_name, exp);

	PeakMap exp2;
	ptr->load(tmp_name, exp2);
	TEST_EQUAL(exp.size() == exp2.size(), true)
	TEST_EQUAL(exp.begin()->size() == exp2.begin()->size(), true)
	TEST_REAL_SIMILAR(exp.begin()->getRT(), exp2.begin()->getRT())
	TEST_REAL_SIMILAR(exp.begin()->getPrecursors().begin()->getMZ(), exp2.begin()->getPrecursors().begin()->getMZ())
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST




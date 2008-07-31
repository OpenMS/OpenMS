// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

using namespace OpenMS;

START_TEST(ClassTest, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

CHECK()
RESULT

CHECK()
RESULT

CHECK()
	std::string tmp_filename;
	NEW_TMP_FILE(tmp_filename);
	TEST::this_test = (tmp_filename != "");
RESULT


CHECK()
	PRECISION(0.5)
	TEST_EQUAL(TEST::precision, 0.5)
RESULT

CHECK()
	PRECISION(0.011)
	TEST_REAL_EQUAL(1.0, 1.01)
	TEST_REAL_EQUAL(1.0, 1.0)
	TEST_REAL_EQUAL(-1.0, -1.01)
	TEST_REAL_EQUAL(-1.01, -1.0)
RESULT

CHECK()
	TEST_EQUAL(1.0, 1.0)
	TEST_EQUAL('A', 'A')
RESULT

CHECK()
	TEST_NOT_EQUAL(0, 1)
	TEST_NOT_EQUAL('A', 'B')
RESULT

CHECK()
	TEST_EXCEPTION(Exception::NullPointer, throw Exception::NullPointer(__FILE__, __LINE__, __PRETTY_FUNCTION__))
RESULT

CHECK()
	TEST_EXCEPTION_WITH_MESSAGE(Exception::NullPointer, throw Exception::NullPointer(__FILE__, __LINE__, __PRETTY_FUNCTION__), "a null pointer was specified")
RESULT

CHECK()
	struct Dummy
	{
		std::string f_dummy(double, float,int,unsigned,long,unsigned long,char) { return __PRETTY_FUNCTION__; }
	} dummy;

	STATUS("\n\n\tExample for usage of __PRETTY_FUNCTION__ inside a member function of a nested class in main():\n\t" << dummy.f_dummy(0,0,0,0,0,0,0) << '\n')
RESULT

CHECK()
	STATUS("status message")
RESULT

CHECK()
	TEST_FILE("data/class_test_infile.txt", "data/class_test_template.txt")
RESULT

CHECK()
	ABORT_IF(true)
	TEST_EQUAL(1, 0)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

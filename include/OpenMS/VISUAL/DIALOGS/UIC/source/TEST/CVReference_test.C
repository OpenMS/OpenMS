// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/DATASTRUCTURES/CVReference.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(CVReference, "$Id: CVReference_test.C 5908 2009-08-26 13:44:26Z marc_sturm $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

CVReference* ptr = 0;
START_SECTION(CVReference())
{
	ptr = new CVReference();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(virtual ~CVReference())
{
	delete ptr;
}
END_SECTION

ptr = new CVReference();

START_SECTION((CVReference(const CVReference &rhs)))
{
 	CVReference cvr;
	TEST_STRING_EQUAL(CVReference(cvr).getName(), cvr.getName())
	TEST_STRING_EQUAL(CVReference(cvr).getIdentifier(), cvr.getIdentifier())

	cvr.setName("my_test_name");
	TEST_STRING_EQUAL(CVReference(cvr).getName(), "my_test_name")

	cvr.setIdentifier("my_test_identifier");
	TEST_STRING_EQUAL(CVReference(cvr).getIdentifier(), "my_test_identifier")
}
END_SECTION

START_SECTION((CVReference& operator=(const CVReference &rhs)))
{
  CVReference cvr, cvr_copy;
	cvr_copy = cvr;
	TEST_STRING_EQUAL(cvr_copy.getName(), "")
	TEST_STRING_EQUAL(cvr_copy.getIdentifier(), "")

	cvr.setName("my_test_name");
	cvr_copy = cvr;
	TEST_STRING_EQUAL(cvr_copy.getName(), "my_test_name")

	cvr.setIdentifier("my_test_identifier");
	cvr_copy = cvr;
	TEST_STRING_EQUAL(cvr_copy.getIdentifier(), "my_test_identifier")
}
END_SECTION

START_SECTION((bool operator == (const CVReference& rhs) const))
{
  CVReference cvr, cvr_copy;
	TEST_EQUAL(cvr == cvr_copy, true)
  cvr_copy = cvr;
	TEST_EQUAL(cvr == cvr_copy, true)

  cvr.setName("my_test_name");
	TEST_EQUAL(cvr == cvr_copy, false)
  cvr_copy = cvr;
	TEST_EQUAL(cvr == cvr_copy, true)

  cvr.setIdentifier("my_test_identifier");
	TEST_EQUAL(cvr == cvr_copy, false)
  cvr_copy = cvr;
	TEST_EQUAL(cvr == cvr_copy, true)
}
END_SECTION

START_SECTION((bool operator != (const CVReference& rhs) const))
{
  CVReference cvr, cvr_copy;
  TEST_EQUAL(cvr != cvr_copy, false)
  cvr_copy = cvr;
  TEST_EQUAL(cvr != cvr_copy, false)

  cvr.setName("my_test_name");
  TEST_EQUAL(cvr != cvr_copy, true)
  cvr_copy = cvr;
  TEST_EQUAL(cvr != cvr_copy, false)

  cvr.setIdentifier("my_test_identifier");
  TEST_EQUAL(cvr != cvr_copy, true)
  cvr_copy = cvr;
  TEST_EQUAL(cvr != cvr_copy, false)
}
END_SECTION

START_SECTION((void setName(const String &name)))
{
  ptr->setName("my_test_name");
	TEST_STRING_EQUAL(ptr->getName(), "my_test_name")
}
END_SECTION

START_SECTION((const String& getName() const ))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void setIdentifier(const String &identifier)))
{
  ptr->setIdentifier("my_test_identifier");
	TEST_STRING_EQUAL(ptr->getIdentifier(), "my_test_identifier")
}
END_SECTION

START_SECTION((const String& getIdentifier() const ))
{
  NOT_TESTABLE
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST




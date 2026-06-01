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
#include <OpenMS/DATASTRUCTURES/CVReference.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(CVReference, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

CVReference* ptr = nullptr;
CVReference* nullPointer = nullptr;
START_SECTION(CVReference())
{
	ptr = new CVReference();
	TEST_NOT_EQUAL(ptr, nullPointer)
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
	TEST_TRUE(cvr == cvr_copy)
  cvr_copy = cvr;
	TEST_TRUE(cvr == cvr_copy)

  cvr.setName("my_test_name");
	TEST_EQUAL(cvr == cvr_copy, false)
  cvr_copy = cvr;
	TEST_TRUE(cvr == cvr_copy)

  cvr.setIdentifier("my_test_identifier");
	TEST_EQUAL(cvr == cvr_copy, false)
  cvr_copy = cvr;
	TEST_TRUE(cvr == cvr_copy)
}
END_SECTION

START_SECTION((bool operator != (const CVReference& rhs) const))
{
  CVReference cvr, cvr_copy;
  TEST_EQUAL(cvr != cvr_copy, false)
  cvr_copy = cvr;
  TEST_EQUAL(cvr != cvr_copy, false)

  cvr.setName("my_test_name");
  TEST_FALSE(cvr == cvr_copy)
  cvr_copy = cvr;
  TEST_EQUAL(cvr != cvr_copy, false)

  cvr.setIdentifier("my_test_identifier");
  TEST_FALSE(cvr == cvr_copy)
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

delete ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST




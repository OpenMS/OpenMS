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
#include <OpenMS/DATASTRUCTURES/CVMappingTerm.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(CVMappingTerm, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

CVMappingTerm* ptr = 0;
CVMappingTerm* nullPointer = 0;
START_SECTION(CVMappingTerm())
{
	ptr = new CVMappingTerm();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(virtual ~CVMappingTerm())
{
	delete ptr;
}
END_SECTION

ptr = new CVMappingTerm();

START_SECTION((CVMappingTerm(const CVMappingTerm &rhs)))
{
 	CVMappingTerm cvmt;
	cvmt.setAccession("my_test_accession");
	TEST_STRING_EQUAL(CVMappingTerm(cvmt).getAccession(), "my_test_accession")

	cvmt.setUseTermName(true);
	TEST_EQUAL(CVMappingTerm(cvmt).getUseTermName(), true)
	cvmt.setUseTermName(false);
	TEST_EQUAL(CVMappingTerm(cvmt).getUseTermName(), false)

	cvmt.setUseTerm(true);
	TEST_EQUAL(CVMappingTerm(cvmt).getUseTerm(), true)
	cvmt.setUseTerm(false);
	TEST_EQUAL(CVMappingTerm(cvmt).getUseTerm(), false)

	cvmt.setTermName("my_test_termname");
	TEST_STRING_EQUAL(CVMappingTerm(cvmt).getTermName(), "my_test_termname")

	cvmt.setIsRepeatable(true);
	TEST_EQUAL(CVMappingTerm(cvmt).getIsRepeatable(), true)
	cvmt.setIsRepeatable(false);
	TEST_EQUAL(CVMappingTerm(cvmt).getIsRepeatable(), false)

	cvmt.setAllowChildren(true);
	TEST_EQUAL(CVMappingTerm(cvmt).getAllowChildren(), true)
	cvmt.setAllowChildren(false);
	TEST_EQUAL(CVMappingTerm(cvmt).getAllowChildren(), false)

	cvmt.setCVIdentifierRef("my_test_cvidentifierref");
	TEST_STRING_EQUAL(CVMappingTerm(cvmt).getCVIdentifierRef(), "my_test_cvidentifierref")
}
END_SECTION

START_SECTION((CVMappingTerm& operator=(const CVMappingTerm &rhs)))
{
	CVMappingTerm cvmt, cvmt_copy;

  cvmt.setAccession("my_test_accession");
	cvmt_copy = cvmt;
  TEST_STRING_EQUAL(cvmt_copy.getAccession(), "my_test_accession")

  cvmt.setUseTermName(true);
	cvmt_copy = cvmt;
  TEST_EQUAL(cvmt_copy.getUseTermName(), true)
  cvmt.setUseTermName(false);
	cvmt_copy = cvmt;
  TEST_EQUAL(cvmt_copy.getUseTermName(), false)

  cvmt.setUseTerm(true);
	cvmt_copy = cvmt;
  TEST_EQUAL(cvmt_copy.getUseTerm(), true)
  cvmt.setUseTerm(false);
	cvmt_copy = cvmt;
  TEST_EQUAL(cvmt_copy.getUseTerm(), false)

  cvmt.setTermName("my_test_termname");
	cvmt_copy = cvmt;
  TEST_STRING_EQUAL(cvmt_copy.getTermName(), "my_test_termname")

  cvmt.setIsRepeatable(true);
	cvmt_copy = cvmt;
  TEST_EQUAL(cvmt_copy.getIsRepeatable(), true)
  cvmt.setIsRepeatable(false);
	cvmt_copy = cvmt;
  TEST_EQUAL(cvmt_copy.getIsRepeatable(), false)

  cvmt.setAllowChildren(true);
	cvmt_copy = cvmt;
  TEST_EQUAL(cvmt_copy.getAllowChildren(), true)
  cvmt.setAllowChildren(false);
	cvmt_copy = cvmt;
  TEST_EQUAL(cvmt_copy.getAllowChildren(), false)

  cvmt.setCVIdentifierRef("my_test_cvidentifierref");
	cvmt_copy = cvmt;
  TEST_STRING_EQUAL(cvmt_copy.getCVIdentifierRef(), "my_test_cvidentifierref")
}
END_SECTION


START_SECTION((bool operator == (const CVMappingTerm& rhs) const))
{
  CVMappingTerm cvmt, cvmt_copy;

  cvmt.setAccession("my_test_accession");
	TEST_EQUAL(cvmt_copy == cvmt, false)
  cvmt_copy = cvmt;
	TEST_EQUAL(cvmt_copy == cvmt, true)

  cvmt.setUseTermName(true);
	TEST_EQUAL(cvmt_copy == cvmt, false)
  cvmt_copy = cvmt;
	TEST_EQUAL(cvmt_copy == cvmt, true)
  cvmt.setUseTermName(false);
	TEST_EQUAL(cvmt_copy == cvmt, false)
  cvmt_copy = cvmt;
	TEST_EQUAL(cvmt_copy == cvmt, true)

  cvmt.setUseTerm(true);
	TEST_EQUAL(cvmt_copy == cvmt, false)
  cvmt_copy = cvmt;
	TEST_EQUAL(cvmt_copy == cvmt, true)
  cvmt.setUseTerm(false);
	TEST_EQUAL(cvmt_copy == cvmt, false)
  cvmt_copy = cvmt;
	TEST_EQUAL(cvmt_copy == cvmt, true)

  cvmt.setTermName("my_test_termname");
	TEST_EQUAL(cvmt_copy == cvmt, false)
  cvmt_copy = cvmt;
	TEST_EQUAL(cvmt_copy == cvmt, true)
  
	cvmt.setIsRepeatable(true);
	TEST_EQUAL(cvmt_copy == cvmt, false)
  cvmt_copy = cvmt;
	TEST_EQUAL(cvmt_copy == cvmt, true)
  cvmt.setIsRepeatable(false);
	TEST_EQUAL(cvmt_copy == cvmt, false)
  cvmt_copy = cvmt;
	TEST_EQUAL(cvmt_copy == cvmt, true)

  cvmt.setAllowChildren(true);
	TEST_EQUAL(cvmt_copy == cvmt, false)
  cvmt_copy = cvmt;
	TEST_EQUAL(cvmt_copy == cvmt, true)
  cvmt.setAllowChildren(false);
	TEST_EQUAL(cvmt_copy == cvmt, false)
  cvmt_copy = cvmt;
	TEST_EQUAL(cvmt_copy == cvmt, true)

  cvmt.setCVIdentifierRef("my_test_cvidentifierref");
	TEST_EQUAL(cvmt_copy == cvmt, false)
  cvmt_copy = cvmt;
	TEST_EQUAL(cvmt_copy == cvmt, true)
}
END_SECTION

START_SECTION((bool operator != (const CVMappingTerm& rhs) const))
{
  CVMappingTerm cvmt, cvmt_copy;

  cvmt.setAccession("my_test_accession");
  TEST_EQUAL(cvmt_copy == cvmt, false)
  cvmt_copy = cvmt;
  TEST_EQUAL(cvmt_copy == cvmt, true)

  cvmt.setUseTermName(true);
  TEST_EQUAL(cvmt_copy == cvmt, false)
  cvmt_copy = cvmt;
  TEST_EQUAL(cvmt_copy == cvmt, true)
  cvmt.setUseTermName(false);
  TEST_EQUAL(cvmt_copy == cvmt, false)
  cvmt_copy = cvmt;
  TEST_EQUAL(cvmt_copy == cvmt, true)

  cvmt.setUseTerm(true);
  TEST_EQUAL(cvmt_copy == cvmt, false)
  cvmt_copy = cvmt;
  TEST_EQUAL(cvmt_copy == cvmt, true)
  cvmt.setUseTerm(false);
  TEST_EQUAL(cvmt_copy == cvmt, false)
  cvmt_copy = cvmt;
  TEST_EQUAL(cvmt_copy == cvmt, true)

  cvmt.setTermName("my_test_termname");
  TEST_EQUAL(cvmt_copy == cvmt, false)
  cvmt_copy = cvmt;
  TEST_EQUAL(cvmt_copy == cvmt, true)

  cvmt.setIsRepeatable(true);
  TEST_EQUAL(cvmt_copy == cvmt, false)
  cvmt_copy = cvmt;
  TEST_EQUAL(cvmt_copy == cvmt, true)
  cvmt.setIsRepeatable(false);
  TEST_EQUAL(cvmt_copy == cvmt, false)
  cvmt_copy = cvmt;
  TEST_EQUAL(cvmt_copy == cvmt, true)

  cvmt.setAllowChildren(true);
  TEST_EQUAL(cvmt_copy == cvmt, false)
  cvmt_copy = cvmt;
  TEST_EQUAL(cvmt_copy == cvmt, true)
  cvmt.setAllowChildren(false);
  TEST_EQUAL(cvmt_copy == cvmt, false)
  cvmt_copy = cvmt;
  TEST_EQUAL(cvmt_copy == cvmt, true)

  cvmt.setCVIdentifierRef("my_test_cvidentifierref");
  TEST_EQUAL(cvmt_copy == cvmt, false)
  cvmt_copy = cvmt;
  TEST_EQUAL(cvmt_copy == cvmt, true)
}
END_SECTION

START_SECTION((void setAccession(const String &accession)))
{
  ptr->setAccession("my_test_accession");
	TEST_STRING_EQUAL(ptr->getAccession(), "my_test_accession")
}
END_SECTION

START_SECTION((const String& getAccession() const ))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void setUseTermName(bool use_term_name)))
{
  ptr->setUseTermName(true);
	TEST_EQUAL(ptr->getUseTermName(), true)
	ptr->setUseTermName(false);
	TEST_EQUAL(ptr->getUseTermName(), false)
}
END_SECTION

START_SECTION((bool getUseTermName() const ))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void setUseTerm(bool use_term)))
{
  ptr->setUseTerm(true);
	TEST_EQUAL(ptr->getUseTerm(), true)
	ptr->setUseTerm(false);
	TEST_EQUAL(ptr->getUseTerm(), false)
}
END_SECTION

START_SECTION((bool getUseTerm() const ))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void setTermName(const String &term_name)))
{
  ptr->setTermName("my_test_termname");
	TEST_EQUAL(ptr->getTermName(), "my_test_termname")
}
END_SECTION

START_SECTION((const String& getTermName() const ))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void setIsRepeatable(bool is_repeatable)))
{
  ptr->setIsRepeatable(true);
	TEST_EQUAL(ptr->getIsRepeatable(), true)
	ptr->setIsRepeatable(false);
	TEST_EQUAL(ptr->getIsRepeatable(), false)
}
END_SECTION

START_SECTION((bool getIsRepeatable() const ))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void setAllowChildren(bool allow_children)))
{
  ptr->setAllowChildren(true);
	TEST_EQUAL(ptr->getAllowChildren(), true)
	ptr->setAllowChildren(false);
	TEST_EQUAL(ptr->getAllowChildren(), false)
}
END_SECTION

START_SECTION((bool getAllowChildren() const ))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void setCVIdentifierRef(const String &cv_identifier_ref)))
{
  ptr->setCVIdentifierRef("my_test_cvidentifierref");
	TEST_EQUAL(ptr->getCVIdentifierRef(), "my_test_cvidentifierref")
}
END_SECTION

START_SECTION((const String& getCVIdentifierRef() const ))
{
  NOT_TESTABLE
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST




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
// $Maintainer: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/METADATA/DocumentIdentifier.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(DocumentIdentifier, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

DocumentIdentifier* ptr = 0;
START_SECTION(DocumentIdentifier())
{
	ptr = new DocumentIdentifier();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(~DocumentIdentifier())
{
	delete ptr;
}
END_SECTION

START_SECTION((DocumentIdentifier(const DocumentIdentifier &source)))
{
  DocumentIdentifier di1;
	di1.setIdentifier("this is a test");
	DocumentIdentifier di2(di1);
	TEST_EQUAL(di2.getIdentifier(), "this is a test")
}
END_SECTION

START_SECTION((DocumentIdentifier& operator=(const DocumentIdentifier &source)))
{
  DocumentIdentifier di1;
	di1.setIdentifier("this is a test");
	DocumentIdentifier di2 = di1;
	TEST_EQUAL(di2.getIdentifier(), "this is a test")
}
END_SECTION

START_SECTION((bool operator==(const DocumentIdentifier &rhs) const))
{
  DocumentIdentifier di1;
	di1.setIdentifier("this is a test");
	DocumentIdentifier di2(di1);
	TEST_EQUAL(di1==di2, true)
}
END_SECTION

START_SECTION((virtual ~DocumentIdentifier()))
{
	NOT_TESTABLE
}
END_SECTION

START_SECTION((void setIdentifier(const String &id)))
{
  DocumentIdentifier di1;
	di1.setIdentifier("this is a test");
	TEST_EQUAL(di1.getIdentifier(), "this is a test")
}
END_SECTION

START_SECTION((const String& getIdentifier() const))
{
	// tested above	
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void swap(DocumentIdentifier& from)))
{
  DocumentIdentifier di1;
	di1.setIdentifier("this is a test");
	DocumentIdentifier di2;
	di1.swap(di2);
	TEST_EQUAL(di1.getIdentifier() == "", true)					
	TEST_EQUAL(di2.getIdentifier() == "this is a test", true)
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST




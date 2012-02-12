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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/FORMAT/XMLFile.h>

///////////////////////////

START_TEST(XMLFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace OpenMS::Internal;
using namespace std;

XMLFile* ptr = 0;
XMLFile* nullPointer = 0;

START_SECTION(XMLFile())
	ptr = new XMLFile();
  TEST_NOT_EQUAL(ptr,nullPointer)
END_SECTION

START_SECTION(~XMLFile())
	delete ptr;
END_SECTION

START_SECTION(XMLFile(const String &schema_location, const String &version))
	NOT_TESTABLE
END_SECTION

START_SECTION(bool isValid(const String &filename,  std::ostream& os = std::cerr) )
	XMLFile f("","");
	TEST_EXCEPTION(Exception::NotImplemented, f.isValid(""))
END_SECTION

START_SECTION(const String& getVersion() const)
	XMLFile f("","1.567");
	TEST_EQUAL( f.getVersion(),"1.567")
END_SECTION

START_SECTION(([EXTRA] void writeXMLEscape(const String& to_escape, std::ostream& os)))
	stringstream ss1, ss2, ss3;
	String s1("nothing_to_escape. Just a regular string...");
	String s2("This string contains an ampersand, &, which must be escaped.");
	String s3("This string also contains characters which is not allowed, and must be escaped; the characters are '>' and '<'");

	writeXMLEscape(s1, ss1);
	TEST_STRING_EQUAL(ss1.str(), "nothing_to_escape. Just a regular string...")
	
	writeXMLEscape(s2, ss2);
	TEST_STRING_EQUAL(ss2.str(), "This string contains an ampersand, &amp;, which must be escaped.")

	writeXMLEscape(s3, ss3);
	TEST_STRING_EQUAL(ss3.str(), "This string also contains characters which is not allowed, and must be escaped; the characters are &apos;&gt;&apos; and &apos;&lt;&apos;");

END_SECTION

START_SECTION(([EXTRA] String writeXMLEscape(const String& to_escape)))
  String s1("nothing_to_escape. Just a regular string...");
  String s2("This string contains an ampersand, &, which must be escaped.");
  String s3("This string also contains characters which is not allowed, and must be escaped; the characters are '>' and '<'");

  TEST_STRING_EQUAL(writeXMLEscape(s1), "nothing_to_escape. Just a regular string...")
  TEST_STRING_EQUAL(writeXMLEscape(s2), "This string contains an ampersand, &amp;, which must be escaped.")
  TEST_STRING_EQUAL(writeXMLEscape(s3), "This string also contains characters which is not allowed, and must be escaped; the characters are &apos;&gt;&apos; and &apos;&lt;&apos;");
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST

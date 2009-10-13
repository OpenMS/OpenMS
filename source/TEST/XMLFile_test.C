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
// $Maintainer: $
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

XMLFile* ptr;
START_SECTION(XMLFile())
	ptr = new XMLFile();
	TEST_NOT_EQUAL(ptr,0)
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

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST

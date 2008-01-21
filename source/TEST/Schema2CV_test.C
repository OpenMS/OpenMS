// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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

#include<OpenMS/FORMAT/Schema2CV.h>

///////////////////////////

START_TEST(Schema2CV, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

Schema2CV* ptr = 0;
CHECK((Schema2CV()))
	ptr = new Schema2CV();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((~Schema2CV()))
	delete ptr;
RESULT

Schema2CV map;

CHECK(void loadFromFile(const String& filename) throw (Exception::FileNotFound, Exception::ParseError))
	map.loadFromFile("data/Schema2CV.xml");
	NOT_TESTABLE
RESULT

CHECK(const std::vector<CVDesc>& getCVs() const)
	TEST_EQUAL(map.getCVs().size(),2)
	
	TEST_EQUAL(map.getCVs()[0].name,"Auto")
	TEST_EQUAL(map.getCVs()[0].version,"1")
	TEST_EQUAL(map.getCVs()[0].uri,"3")
	TEST_EQUAL(map.getCVs()[0].id,"A")
	TEST_EQUAL(map.getCVs()[0].format,"OBO")
	
	TEST_EQUAL(map.getCVs()[1].name,"Klettern")
	TEST_EQUAL(map.getCVs()[1].version,"2")
	TEST_EQUAL(map.getCVs()[1].uri,"4")
	TEST_EQUAL(map.getCVs()[1].id,"K")
	TEST_EQUAL(map.getCVs()[1].format,"OBO")
RESULT

CHECK(const std::vector<LocDesc>& getLocations() const)
	TEST_EQUAL(map.getLocations().size(),2)
	
	TEST_EQUAL(map.getLocations()[0].location, "//Kletterurlaub/Auto")
	TEST_EQUAL(map.getLocations()[0].strict, true)
	TEST_EQUAL(map.getLocations()[0].terms.size(), 3)
	
	TEST_EQUAL(map.getLocations()[0].terms[0].accession, "Auto:4")
	TEST_EQUAL(map.getLocations()[0].terms[0].cv,"A")
	TEST_EQUAL(map.getLocations()[0].terms[0].allowSelf, true)
	TEST_EQUAL(map.getLocations()[0].terms[0].allowChildren, false)
	TEST_EQUAL(map.getLocations()[0].terms[0].repeatable, true)
	
	TEST_EQUAL(map.getLocations()[0].terms[1].accession, "Auto:5")
	TEST_EQUAL(map.getLocations()[0].terms[1].cv,"A")
	TEST_EQUAL(map.getLocations()[0].terms[1].allowSelf, true)
	TEST_EQUAL(map.getLocations()[0].terms[1].allowChildren, false)
	TEST_EQUAL(map.getLocations()[0].terms[1].repeatable, true)
	
	TEST_EQUAL(map.getLocations()[0].terms[2].accession, "Auto:6")
	TEST_EQUAL(map.getLocations()[0].terms[2].cv,"A")
	TEST_EQUAL(map.getLocations()[0].terms[2].allowSelf, true)
	TEST_EQUAL(map.getLocations()[0].terms[2].allowChildren, false)
	TEST_EQUAL(map.getLocations()[0].terms[2].repeatable, true)

	TEST_EQUAL(map.getLocations()[1].location, "//Kletterurlaub/Gebiet")
	TEST_EQUAL(map.getLocations()[1].strict, false)
	TEST_EQUAL(map.getLocations()[1].terms.size(), 1)
	
	TEST_EQUAL(map.getLocations()[1].terms[0].accession, "Klettern:1")
	TEST_EQUAL(map.getLocations()[1].terms[0].cv,"K")
	TEST_EQUAL(map.getLocations()[1].terms[0].allowSelf, false)
	TEST_EQUAL(map.getLocations()[1].terms[0].allowChildren, true)
	TEST_EQUAL(map.getLocations()[1].terms[0].repeatable, true)

RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

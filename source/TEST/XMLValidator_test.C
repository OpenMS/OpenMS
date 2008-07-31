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

///////////////////////////

#include <OpenMS/FORMAT/XMLValidator.h>

///////////////////////////

START_TEST(XMLValidator, "XMLValidator")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

XMLValidator* ptr = 0;
CHECK((XMLValidator()))
	ptr = new XMLValidator;
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(([EXTRA]~XMLValidator()))
	delete ptr;
RESULT

CHECK((bool isValid(const String &filename, const String &schema) ))
	XMLValidator v;
	
	TEST_EQUAL(v.isValid("data/XMLValidator_valid.xml","data/XMLValidator.xsd"), true);
	
	TEST_EQUAL(v.isValid("data/XMLValidator_missing_element.xml","data/XMLValidator.xsd"), false);
	
	TEST_EQUAL(v.isValid("data/XMLValidator_missing_attribute.xml","data/XMLValidator.xsd"), false);
	
	TEST_EQUAL(v.isValid("data/XMLValidator_syntax.xml","data/XMLValidator.xsd"), false);
	
	//check vaild fail again to make sure internal states are ok
	TEST_EQUAL(v.isValid("data/XMLValidator_valid.xml","data/XMLValidator.xsd"), true);
	
	//test exception
	TEST_EXCEPTION(Exception::FileNotFound, v.isValid("data/this_file_does_not_exist.for_sure","data/XMLValidator.xsd"));
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

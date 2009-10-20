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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/FORMAT/VALIDATORS/XMLValidator.h>

///////////////////////////

START_TEST(XMLValidator, "XMLValidator")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

XMLValidator* ptr = 0;
START_SECTION((XMLValidator()))
	ptr = new XMLValidator;
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION(([EXTRA]~XMLValidator()))
	delete ptr;
END_SECTION

START_SECTION((bool isValid(const String &filename, const String &schema,  std::ostream& os = std::cerr) ))
	XMLValidator v;
	
	TEST_EQUAL(v.isValid(OPENMS_GET_TEST_DATA_PATH("XMLValidator_valid.xml"),OPENMS_GET_TEST_DATA_PATH("XMLValidator.xsd")), true);
	
	TEST_EQUAL(v.isValid(OPENMS_GET_TEST_DATA_PATH("XMLValidator_missing_element.xml"),OPENMS_GET_TEST_DATA_PATH("XMLValidator.xsd")), false);
	
	TEST_EQUAL(v.isValid(OPENMS_GET_TEST_DATA_PATH("XMLValidator_missing_attribute.xml"),OPENMS_GET_TEST_DATA_PATH("XMLValidator.xsd")), false);
	
	TEST_EQUAL(v.isValid(OPENMS_GET_TEST_DATA_PATH("XMLValidator_syntax.xml"),OPENMS_GET_TEST_DATA_PATH("XMLValidator.xsd")), false);
	
	//check vaild fail again to make sure internal states are ok
	TEST_EQUAL(v.isValid(OPENMS_GET_TEST_DATA_PATH("XMLValidator_valid.xml"),OPENMS_GET_TEST_DATA_PATH("XMLValidator.xsd")), true);
	
	//test exception
	TEST_EXCEPTION(Exception::FileNotFound, v.isValid(OPENMS_GET_TEST_DATA_PATH("this_file_does_not_exist.for_sure"),OPENMS_GET_TEST_DATA_PATH("XMLValidator.xsd")));
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

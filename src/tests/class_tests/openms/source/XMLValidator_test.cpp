// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FORMAT/VALIDATORS/XMLValidator.h>

///////////////////////////

START_TEST(XMLValidator, "XMLValidator")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

XMLValidator* ptr = nullptr;
XMLValidator* nullPointer = nullptr;
START_SECTION((XMLValidator()))
	ptr = new XMLValidator;
	TEST_NOT_EQUAL(ptr, nullPointer)
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

// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/FORMAT/VALIDATORS/SemanticValidator.h>
#include <OpenMS/DATASTRUCTURES/CVMappings.h>
#include <OpenMS/FORMAT/CVMappingFile.h>
#include <OpenMS/FORMAT/ControlledVocabulary.h>
#include <OpenMS/SYSTEM/File.h>

///////////////////////////

START_TEST(SemanticValidator, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace OpenMS::Internal;
using namespace std;

CVMappings mapping;
CVMappingFile().load(OPENMS_GET_TEST_DATA_PATH("SemanticValidator_mapping.xml"),mapping);

ControlledVocabulary cv;
cv.loadFromOBO("PSI",OPENMS_GET_TEST_DATA_PATH("SemanticValidator_cv.obo"));
cv.loadFromOBO("PATO",File::find("/CV/quality.obo"));
cv.loadFromOBO("UO",File::find("/CV/unit.obo"));
cv.loadFromOBO("brenda",File::find("/CV/brenda.obo"));
cv.loadFromOBO("GO",File::find("/CV/goslim_goa.obo"));

SemanticValidator* ptr = nullptr;
SemanticValidator* nullPointer = nullptr;
START_SECTION((SemanticValidator(const CVMappings& mapping, const ControlledVocabulary& cv)))
	ptr = new SemanticValidator(mapping,cv);
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~SemanticValidator()))
	delete ptr;
END_SECTION

START_SECTION((void setTag(const String& tag)))
	NOT_TESTABLE
END_SECTION
	
START_SECTION((void setAccessionAttribute(const String& accession)))
	NOT_TESTABLE
END_SECTION
	
START_SECTION((void setNameAttribute(const String& name)))
	NOT_TESTABLE
END_SECTION
	
START_SECTION((void setValueAttribute(const String& value)))
	NOT_TESTABLE
END_SECTION

START_SECTION((bool validate(const String &filename, StringList &errors, StringList &warnings)))
	StringList errors, warnings;
	
	//----------------------------------------------------------------------------------------
	//test exceptions
	SemanticValidator sv(mapping, cv);
	TEST_EXCEPTION(Exception::FileNotFound, sv.validate("/does/not/exist", errors, warnings));

	//----------------------------------------------------------------------------------------
	//test of valid file
	TEST_EQUAL(sv.validate(OPENMS_GET_TEST_DATA_PATH("SemanticValidator_valid.xml"), errors, warnings),true);
	TEST_EQUAL(errors.size(),0)
	TEST_EQUAL(warnings.size(),0)

	//----------------------------------------------------------------------------------------
	//test of corrupt file
	TEST_EQUAL(sv.validate(OPENMS_GET_TEST_DATA_PATH("SemanticValidator_corrupt.xml"), errors, warnings),false);
	TEST_EQUAL(errors.size(),5)
	TEST_STRING_EQUAL(errors[0],"Violated mapping rule 'R3' at element '/mzML/fileDescription/sourceFileList/sourceFile', 2 term(s) should be present, 1 found!")
	TEST_STRING_EQUAL(errors[1],"Name of CV term not correct: 'MS:1000554 - LCQ Deca2 - invalid repeat' should be 'LCQ Deca'")
	TEST_STRING_EQUAL(errors[2],"CV term used in invalid element: 'MS:1000030 - vendor' at element '/mzML/instrumentConfigurationList/instrumentConfiguration'")
	TEST_STRING_EQUAL(errors[3],"Violated mapping rule 'R6a' number of term repeats at element '/mzML/instrumentConfigurationList/instrumentConfiguration'")
	TEST_STRING_EQUAL(errors[4],"Violated mapping rule 'R17a' at element '/mzML/run/spectrumList/spectrum/spectrumDescription', 1 term(s) should be present, 0 found!")
	TEST_EQUAL(warnings.size(),4)
	TEST_STRING_EQUAL(warnings[0],"Unknown CV term: 'MS:1111569 - SHA-1' at element '/mzML/fileDescription/sourceFileList/sourceFile'")
	TEST_STRING_EQUAL(warnings[1],"Obsolete CV term: 'MS:1000030 - vendor' at element '/mzML/instrumentConfigurationList/instrumentConfiguration'")
	TEST_STRING_EQUAL(warnings[2],"No mapping rule found for element '/mzML/acquisitionSettingsList/acquisitionSettings/targetList/target'")
	TEST_STRING_EQUAL(warnings[3],"No mapping rule found for element '/mzML/acquisitionSettingsList/acquisitionSettings/targetList/target'")

END_SECTION

START_SECTION((void setCheckTermValueTypes(bool check)))
	SemanticValidator sv(mapping, cv);
	sv.setCheckTermValueTypes(true);
	NOT_TESTABLE
END_SECTION

START_SECTION((void setCheckUnits(bool check)))
	SemanticValidator sv(mapping, cv);
	sv.setCheckUnits(true);
	NOT_TESTABLE
END_SECTION

START_SECTION((void setUnitAccessionAttribute(const String &accession)))
	SemanticValidator sv(mapping, cv);
	sv.setUnitAccessionAttribute("unitAccession");
	NOT_TESTABLE
END_SECTION

START_SECTION((void setUnitNameAttribute(const String &name)))
	SemanticValidator sv(mapping, cv);
	sv.setUnitNameAttribute("unit");
	NOT_TESTABLE
END_SECTION




/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

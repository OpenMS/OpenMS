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

#include <OpenMS/FORMAT/SemanticValidator.h>
#include <OpenMS/FORMAT/CVMappings.h>
#include <OpenMS/FORMAT/CVMappingFile.h>
#include <OpenMS/FORMAT/ControlledVocabulary.h>
#include <OpenMS/SYSTEM/File.h>

///////////////////////////

START_TEST(SemanticValidator, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

CVMappings mapping;
CVMappingFile().load("data/SemanticValidator_mapping.xml",mapping);

ControlledVocabulary cv;
cv.loadFromOBO("PSI","data/SemanticValidator_cv.obo");
cv.loadFromOBO("PATO",File::find("/CV/quality.obo"));
cv.loadFromOBO("UO",File::find("/CV/unit.obo"));
cv.loadFromOBO("brenda",File::find("/CV/brenda.obo"));
cv.loadFromOBO("GO",File::find("/CV/goslim_goa.obo"));

SemanticValidator* ptr = 0;
CHECK(SemanticValidator(const CVMappings& mapping, const ControlledVocabulary& cv))
	ptr = new SemanticValidator(mapping,cv);
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((virtual ~SemanticValidator()))
	delete ptr;
RESULT

CHECK(void setTag(const String& tag))
	NOT_TESTABLE
RESULT
	
CHECK(void setAccessionAttribute(const String& accession))
	NOT_TESTABLE
RESULT
	
CHECK(void setNameAttribute(const String& name))
	NOT_TESTABLE
RESULT
	
CHECK(void setValueAttribute(const String& value))
	NOT_TESTABLE
RESULT

CHECK(bool validate(const String& filename, ValidationOutput& output))
	SemanticValidator::ValidationOutput out;
	
	//----------------------------------------------------------------------------------------
	//test exceptions
	SemanticValidator sv(mapping, cv);
	TEST_EXCEPTION(Exception::FileNotFound, sv.validate("/does/not/exist", out));

	//----------------------------------------------------------------------------------------
	//test of valid file
	TEST_EQUAL(sv.validate("data/SemanticValidator_valid.mzML", out),true);
	TEST_EQUAL(out.unknown_terms.size(),0)
	TEST_EQUAL(out.obsolete_terms.size(),0)
	TEST_EQUAL(out.invalid_location.size(),0)
	TEST_EQUAL(out.no_mapping.size(),0)
	TEST_EQUAL(out.violated.size(),0)
	TEST_EQUAL(out.violated_repeats.size(),0)

	//----------------------------------------------------------------------------------------
	//test of corrupt file
	TEST_EQUAL(sv.validate("data/SemanticValidator_corrupt.mzML", out),false);
	
	//unkonwn
	TEST_EQUAL(out.unknown_terms.size(),1)
	TEST_STRING_EQUAL(out.unknown_terms[0].path,"/mzML/fileDescription/sourceFileList/sourceFile/cvParam/@accession")
	TEST_STRING_EQUAL(out.unknown_terms[0].accession,"MS:1111569")
	TEST_STRING_EQUAL(out.unknown_terms[0].name,"SHA-1")
	TEST_STRING_EQUAL(out.unknown_terms[0].value,"81be39fb2700ab2f3c8b2234b91274968b6899b1")
	
	//obsolete
	TEST_EQUAL(out.obsolete_terms.size(),1)
	TEST_STRING_EQUAL(out.obsolete_terms[0].path,"/mzML/instrumentConfigurationList/instrumentConfiguration/cvParam/@accession")
	TEST_STRING_EQUAL(out.obsolete_terms[0].accession,"MS:1000030")
	TEST_STRING_EQUAL(out.obsolete_terms[0].name,"vendor")
	TEST_STRING_EQUAL(out.obsolete_terms[0].value,"Thermo")

	//invalid location
	TEST_EQUAL(out.invalid_location.size(),2)
	TEST_STRING_EQUAL(out.invalid_location[0].path,"/mzML/fileDescription/sourceFileList/sourceFile/cvParam/@accession")
	TEST_STRING_EQUAL(out.invalid_location[0].accession,"MS:1111569")
	TEST_STRING_EQUAL(out.invalid_location[1].path,"/mzML/instrumentConfigurationList/instrumentConfiguration/cvParam/@accession")
	TEST_STRING_EQUAL(out.invalid_location[1].accession,"MS:1000030")

	//no mapping rule found
	TEST_EQUAL(out.no_mapping.size(),2)
	TEST_STRING_EQUAL(out.no_mapping[0].path,"/mzML/acquisitionSettingsList/acquisitionSettings/targetList/target/cvParam/@accession")
	TEST_STRING_EQUAL(out.no_mapping[0].accession,"MS:1000040")
	TEST_STRING_EQUAL(out.no_mapping[1].path,"/mzML/acquisitionSettingsList/acquisitionSettings/targetList/target/cvParam/@accession")
	TEST_STRING_EQUAL(out.no_mapping[1].accession,"MS:1000040")

	//violated rules
	TEST_EQUAL(out.violated.size(),2)
	TEST_EQUAL(out.violated[0],"R3")
	TEST_EQUAL(out.violated[1],"R17a")

	//Repeats
	TEST_EQUAL(out.violated_repeats.size(),1)
	TEST_EQUAL(out.violated_repeats[0],"R6a")

RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

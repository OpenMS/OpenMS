// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Chris Bauer $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/FORMAT/PTMXMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/PTMXMLHandler.h>

#include <vector>

///////////////////////////

START_TEST(PTMXMLFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

PTMXMLFile* ptr = nullptr;
PTMXMLFile* nullPointer = nullptr;
PTMXMLFile xml_file;

START_SECTION((PTMXMLFile()))
	ptr = new PTMXMLFile();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((void load(const std::string& filename, std::map<std::string, std::pair<std::string, std::string> >& ptm_informations)))

	map<std::string, pair<std::string, std::string> > ptm_informations;
	xml_file.load(OPENMS_GET_TEST_DATA_PATH("PTMs.xml"), ptm_informations);
	
	TEST_EQUAL(ptm_informations["TEST"].first, "N2O2-CH3")
	TEST_EQUAL(ptm_informations["TEST"].second, "KLR")
END_SECTION


START_SECTION((void store(std::string filename, std::map<std::string, std::pair<std::string, std::string> > &ptm_informations) const))

	map<std::string, pair<std::string, std::string> > ptm_informations;
	xml_file.load(OPENMS_GET_TEST_DATA_PATH("PTMs.xml"), ptm_informations);
	std::string temp_filename;
	NEW_TMP_FILE(temp_filename)
	xml_file.store(temp_filename, ptm_informations);
	ptm_informations.clear();
	xml_file.load(temp_filename, ptm_informations);
	
	TEST_EQUAL(ptm_informations["TEST"].first, "N2O2-CH3")
	TEST_EQUAL(ptm_informations["TEST"].second, "KLR")
END_SECTION

delete ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FORMAT/PepNovoInfile.h>

#include <iostream>
#include <vector>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(PepNovoInfile, "$Id$")

/////////////////////////////////////////////////////////////

PepNovoInfile* ptr = nullptr;
PepNovoInfile* nullPointer = nullptr;
START_SECTION(PepNovoInfile())
	ptr = new PepNovoInfile();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~PepNovoInfile())
	delete ptr;
END_SECTION

StringList fix_mods, var_mods;
map<String,String>keys_and_mods;
fix_mods.push_back("Phospho (C)");
var_mods.push_back("Phospho (D)");
var_mods.push_back("Ethanolamine (C-term)");
//var_mods.push_back("TMT6plex (N-term)");

START_SECTION((bool operator==(const PepNovoInfile &pepnovo_infile) const))
  PepNovoInfile pepnovo_infile1;
  pepnovo_infile1.setModifications(fix_mods, var_mods);
  PepNovoInfile pepnovo_infile2;
  pepnovo_infile2 = pepnovo_infile1;
  TEST_EQUAL(( pepnovo_infile1 == pepnovo_infile2 ), true)
END_SECTION

START_SECTION((PepNovoInfile& operator=(const PepNovoInfile& pepnovo_infile)))
  PepNovoInfile pepnovo_infile1;
  pepnovo_infile1.setModifications(fix_mods, var_mods);
  PepNovoInfile pepnovo_infile2;
  pepnovo_infile2 = pepnovo_infile1;
  TEST_EQUAL(( pepnovo_infile1 == pepnovo_infile2 ), true)
END_SECTION


START_SECTION((PepNovoInfile(const PepNovoInfile &pepnovo_infile)))
  PepNovoInfile pepnovo_infile1;
  pepnovo_infile1.setModifications(fix_mods, var_mods);
  PepNovoInfile pepnovo_infile2;
  pepnovo_infile2 = pepnovo_infile1;
  TEST_EQUAL(( pepnovo_infile1 == pepnovo_infile2 ), true)
END_SECTION

START_SECTION((void setModifications(const StringList &fixed_mods, const StringList &variable_mods)))
	NOT_TESTABLE // will be tested in next section
END_SECTION

START_SECTION(void getModifications(std::map<String,String>& modification_key_map) const)
	PepNovoInfile pepnovo_infile;
	pepnovo_infile.setModifications(fix_mods, var_mods);
	pepnovo_infile.getModifications(keys_and_mods);

	//TEST_EQUAL(keys_and_mods.size(), 4)
  TEST_EQUAL(keys_and_mods.size(), 3)


  if(keys_and_mods.size()==4)
  {
    map<String, String>::iterator mod_it=keys_and_mods.begin();
    TEST_EQUAL((mod_it++)->first, "$+43")
    TEST_EQUAL((mod_it++)->first, "C+80")
    TEST_EQUAL((mod_it++)->first, "D+80")
//    TEST_EQUAL((mod_it)->first, "^+229")
  }
END_SECTION

START_SECTION(void store(const String& filename))
  PepNovoInfile pepnovo_infile;
  pepnovo_infile.setModifications(fix_mods, var_mods);
	String filename;
	NEW_TMP_FILE(filename)

	// test actual program
	pepnovo_infile.store(filename);
//	pepnovo_infile.store("test_infile.txt");
  

	TEST_FILE_EQUAL(filename.c_str(), OPENMS_GET_TEST_DATA_PATH("PepNovoInfile_test_template.txt"));
	// if the comparison fails because the unimod.xml has been replaced, remove non-ascii characters
	// from the unimod.xml file. E.g. registrated trademark symbol
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

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
// $Maintainer: Martin Langwisch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/FORMAT/InspectInfile.h>

#include <iostream>
#include <vector>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(InspectInfile, "$Id$")

/////////////////////////////////////////////////////////////

InspectInfile* ptr = 0;
CHECK(InspectInfile())
	ptr = new InspectInfile();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~InspectInfile())
	delete ptr;
RESULT

InspectInfile file;

CHECK(void setSpectra(const string& spectra))
	file.setSpectra("dummy4712");
	TEST_EQUAL(file.getSpectra(), "dummy4712")
RESULT

CHECK(const string& getSpectra() const)
	TEST_EQUAL(file.getSpectra(), "dummy4712")
RESULT


CHECK(void setDb(const String& db))
	file.setDb("dummy4711");
	TEST_EQUAL(file.getDb(), "dummy4711");
RESULT

CHECK(const String& getDb() const)
	TEST_EQUAL(file.getDb(), "dummy4711");
RESULT


CHECK(void setProtease(const String& protease))
	file.setProtease("Trypsin");
	TEST_EQUAL(file.getProtease(), "Trypsin")
RESULT

CHECK(const String& getProtease() const)
	TEST_EQUAL(file.getProtease(), "Trypsin")
RESULT


CHECK(void setMod(const vector< vector< String > >& mod))
	vector< vector< String > > mod;
	mod.push_back(vector< String >());
	mod.back().push_back("+57");
	mod.back().push_back("C");
	mod.back().push_back("fix");
	mod.back().push_back("carbamidomethylation");
	file.setMod(mod);
	TEST_EQUAL((file.getMod() == mod), true)
RESULT

CHECK(void addMod(const vector< String >& mod))
	vector< String > mod;
	mod.push_back("80");
	mod.push_back("STY");
	mod.push_back("opt");
	mod.push_back("phosphorylation");
	file.addMod(mod);
	TEST_EQUAL((file.getMod().back() == mod), true)
RESULT

CHECK(vector< vector< String > >& getMod())
	vector< vector< String > > mod;
	mod.push_back(vector< String >());
	mod.back().push_back("+57");
	mod.back().push_back("C");
	mod.back().push_back("fix");
	mod.back().push_back("carbamidomethylation");
	mod.push_back(vector< String >());
	mod.back().push_back("80");
	mod.back().push_back("STY");
	mod.back().push_back("opt");
	mod.back().push_back("phosphorylation");
	TEST_EQUAL((file.getMod() == mod), true)
RESULT


CHECK(void setMods(int mods))
	file.setMods(2);
	TEST_EQUAL(file.getMods(), 2)
RESULT

CHECK(const int getMods() const)
	TEST_EQUAL(file.getMods(), 2)
RESULT


CHECK(void setBlind(unsigned int blind))
	file.setBlind(1);
	TEST_EQUAL(file.getBlind(), 1)
RESULT

CHECK(const unsigned int getBlind() const)
	TEST_EQUAL(file.getBlind(), 1)
RESULT


CHECK(void setMaxPTMsize(DoubleReal maxptmsize))
	file.setMaxPTMsize(250);
	TEST_EQUAL(file.getMaxPTMsize(), 250)
RESULT

CHECK(const DoubleReal getMaxPTMsize() const)
	TEST_EQUAL(file.getMaxPTMsize(), 250)
RESULT


CHECK(void setPMTolerance(DoubleReal PM_tolerance))
	file.setPMTolerance(1.3);
	TEST_EQUAL(file.getPMTolerance(), 1.3)
RESULT

CHECK(const DoubleReal getPMTolerance() const)
	TEST_EQUAL(file.getPMTolerance(), 1.3)
RESULT


CHECK(void setIonTolerance(DoubleReal ion_tolerance))
	file.setIonTolerance(0.3);
	TEST_EQUAL(file.getIonTolerance(), 0.3)
RESULT

CHECK(const DoubleReal getIonTolerance() const)
	TEST_EQUAL(file.getIonTolerance(), 0.3)
RESULT


CHECK(void setMulticharge(unsigned int multicharge))
	file.setMulticharge(1);
	TEST_EQUAL(file.getMulticharge(), 1)
RESULT

CHECK(const unsigned int getMulticharge() const)
	TEST_EQUAL(file.getMulticharge(), 1)
RESULT


CHECK(void setInstrument(const String& instrument))
	file.setInstrument("ESI-ION-TRAP");
	TEST_EQUAL(file.getInstrument(), "ESI-ION-TRAP")
RESULT

CHECK(const String& getInstrument() const)
	TEST_EQUAL(file.getInstrument(), "ESI-ION-TRAP")
RESULT


CHECK(void setTagCount(int TagCount))
	file.setTagCount(1);
	TEST_EQUAL(file.getTagCount(), 1)
RESULT

CHECK(const int getTagCount() const)
	TEST_EQUAL(file.getTagCount(), 1)
RESULT


CHECK(void store(const String& filename) throw (Exception::UnableToCreateFile))
	file.store("InspectInfile_test.txt");
	TEST_FILE("InspectInfile_test.txt", "data/InspectInfile_test_template1.txt");
	remove("InspectInfile_test.txt");
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

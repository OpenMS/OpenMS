// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Nico Pfeifer $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/FORMAT/MascotInfile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <iostream>
#include <vector>

using namespace OpenMS;
using namespace std;


///////////////////////////

START_TEST(MascotInfile, "$Id$")

/////////////////////////////////////////////////////////////

//DPeakArray (dummy for spectrum)
DPeakArrayNonPolymorphic<1> spec;
DPeak<1> tmp;
vector<SignedInt> charges;
charges.push_back(2);
for (UnsignedInt i=1;i<10;i+=1)
{
	tmp.setPosition(DPosition<1>(i));
	tmp.setIntensity(i*i);
	spec.push_back(tmp);	
}

MascotInfile* ptr = 0;
CHECK(MascotInfile())
	ptr = new MascotInfile(spec,1998.0f,charges,"TestTitle", 25.379);
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~MascotInfile())
	delete ptr;
RESULT

MascotInfile file(spec,1998.0f,charges,"TestTitle", 25.379);

CHECK(setBoundary() / getBoundary())
	file.setBoundary("ABCDEFGHIJKMNOPQRSTUVWXYZ");
	TEST_EQUAL(file.getBoundary() , "ABCDEFGHIJKMNOPQRSTUVWXYZ")
RESULT

CHECK( write( filename) (defaults) )
	// here a fixed name has to be used as it has to be in the tamplate
	file.write("MascotInfile_test.txt");
	TEST_FILE("MascotInfile_test.txt", "data/MascotInfile_test_template1.txt");
RESULT

CHECK(setDB(...) / getDB())
	file.setDB("DB_TEST");
	TEST_EQUAL(file.getDB() , "DB_TEST")
RESULT

CHECK(setSearchType(...) / getSearchType())
	file.setSearchType("SearchType_TEST");
	TEST_EQUAL(file.getSearchType() , "SearchType_TEST")
RESULT

CHECK(setHits(...) / getHits())
	file.setHits("Hits_TEST");
	TEST_EQUAL(file.getHits() , "Hits_TEST")
RESULT

CHECK(setCleavage(...) / getCleavage())
	file.setCleavage("Cleavage_TEST");
	TEST_EQUAL(file.getCleavage() , "Cleavage_TEST")
RESULT

CHECK(setMassType(...) / getMassType())
	file.setMassType("MassType_TEST");
	TEST_EQUAL(file.getMassType() , "MassType_TEST")
RESULT

CHECK(setInstrument(...) / getInstrument())
	file.setInstrument("Instrument_TEST");
	TEST_EQUAL(file.getInstrument() , "Instrument_TEST")
RESULT

CHECK(setMissedCleavages(...) / getMissedCleavages())
	file.setMissedCleavages(4711);
	TEST_EQUAL(file.getMissedCleavages() , 4711)
RESULT

CHECK(setPrecursorMassTolerance(...) / getPrecursorMassTolerance())
	file.setPrecursorMassTolerance(4711.1f);
	TEST_REAL_EQUAL(file.getPrecursorMassTolerance() , 4711.1f)
RESULT

CHECK(setPeakMassTolerance(...) / getPeakMassTolerance())
	file.setPeakMassTolerance(4711.2f);
	TEST_REAL_EQUAL(file.getPeakMassTolerance() , 4711.2f)
RESULT

CHECK(setTaxonomy(...) / getTaxonomy())
	file.setTaxonomy("Taxonomy_TEST");
	TEST_EQUAL(file.getTaxonomy() , "Taxonomy_TEST")
RESULT

CHECK(setFormVersion(...) / getFormVersion())
	file.setFormVersion("FormVersion_TEST");
	TEST_EQUAL(file.getFormVersion() , "FormVersion_TEST")
RESULT

CHECK(setModifications(...) / getModifications())
	vector<String> mods;
	mods.push_back("Modifiactions_TEST_1");
	mods.push_back("Modifiactions_TEST_2");
	file.setModifications(mods);
	TEST_EQUAL(file.getModifications() == mods, true)
RESULT

CHECK( write( filename) (settings) )
	// here a fixed name has to be used as it has to be in the tamplate
	file.write("MascotInfile_test.txt");
	TEST_FILE("MascotInfile_test.txt", "data/MascotInfile_test_template2.txt");
	remove("MascotInfile_test.txt");
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

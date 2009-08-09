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
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <limits>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ModificationsDB, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ModificationsDB* ptr = 0;
START_SECTION(ModificationsDB* getInstance())
{
	ptr = ModificationsDB::getInstance();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(Size getNumberOfModifications() const)
	// range because data may change over time
	TEST_EQUAL(ptr->getNumberOfModifications() > 100, true);
END_SECTION

START_SECTION(const ResidueModification& getModification(Size index) const)
	//TEST_EQUAL(ptr->getModification(0).getId().size() > 0, true)
END_SECTION

START_SECTION(std::set<String> searchModifications(const String &name) const)
	set<String> mods = ptr->searchModifications("Phosphorylation");
	//TEST_EQUAL(mods.find("MOD:00046") != mods.end(), true)
	//TEST_EQUAL(mods.find("MOD:00047") != mods.end(), true)
	//TEST_EQUAL(mods.find("MOD:00048") != mods.end(), true)
END_SECTION

START_SECTION(const ResidueModification& getModification(const String &name) const)
	TEST_EQUAL(ptr->getModification("Carboxymethyl Cystenyl").getId(), "MOD:01062")
END_SECTION

START_SECTION(const ResidueModification& getModification(const String &residue_name, const String &mod_name) const)
	//TEST_EQUAL(ptr->getModification("S", "Phosphorylation").getId(), "MOD:00046")
END_SECTION

START_SECTION(Size findModificationIndex(const String &mod_name) const)
	Size index = numeric_limits<Size>::max();
	//index = ptr->findModificationIndex("MOD:00046");
	//TEST_NOT_EQUAL(index, numeric_limits<Size>::max())
END_SECTION
    
START_SECTION(void getModificationsByDiffMonoMass(std::vector< String > &mods, double mass, double error=0.0))
	vector<String> mods;
	ptr->getModificationsByDiffMonoMass(mods, 80.0, 0.1);
	set<String> uniq_mods;
	for (vector<String>::const_iterator it = mods.begin(); it != mods.end(); ++it)
	{
		uniq_mods.insert(*it);
	}

	//TEST_EQUAL(uniq_mods.find("MOD:00046") != uniq_mods.end(), true)
	//TEST_EQUAL(uniq_mods.find("MOD:00047") != uniq_mods.end(), true)
	//TEST_EQUAL(uniq_mods.find("MOD:00048") != uniq_mods.end(), true)
	//TEST_EQUAL(uniq_mods.find("MOD:00180") != uniq_mods.end(), true)
END_SECTION

START_SECTION(void getModificationsByDiffMonoMass(std::vector< String > &mods, const String &residue, double mass, double error=0.0))
	vector<String> mods;
	ptr->getModificationsByDiffMonoMass(mods, "S", 80.0, 0.1);
	set<String> uniq_mods;
	for (vector<String>::const_iterator it = mods.begin(); it != mods.end(); ++it)
	{
		uniq_mods.insert(*it);
	}

	//TEST_EQUAL(uniq_mods.find("MOD:00046") != uniq_mods.end(), true)
	//TEST_EQUAL(uniq_mods.find("MOD:00366") != uniq_mods.end(), true)
END_SECTION

START_SECTION(void readFromOBOFile(const String &filename))
	// implicitely tested above
	NOT_TESTABLE
END_SECTION

START_SECTION(void readFromUnimodXMLFile(const String &filename))
	// just provided for convenience at the moment
	NOT_TESTABLE
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST




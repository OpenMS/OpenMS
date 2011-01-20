// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Stephan Aiche $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <limits>
#include <algorithm>
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
	TEST_EQUAL(ptr->getModification(0).getId().size() > 0, true)
END_SECTION

START_SECTION(void searchModifications(std::set<const ResidueModification*>& mods, const String& orgin, const String& mod_name, ResidueModification::Term_Specificity term_spec) const)
	set<const ResidueModification*> mods;
	ptr->searchModifications(mods, "T", "Phosphorylation", ResidueModification::ANYWHERE);
	TEST_EQUAL(mods.size(), 1)
	TEST_STRING_EQUAL((*mods.begin())->getFullId(), "Phospho (T)")
END_SECTION

START_SECTION((void searchTerminalModifications(std::set< const ResidueModification * > &mods, const String &name, ResidueModification::Term_Specificity term_spec) const))
	set<const ResidueModification*> mods;
	ptr->searchTerminalModifications(mods, "NIC", ResidueModification::N_TERM);
	TEST_EQUAL(mods.size(), 1)
END_SECTION

START_SECTION((void searchModifications(std::set< const ResidueModification * > &mods, const String &mod_name, ResidueModification::Term_Specificity term_spec) const ))
{
  set<const ResidueModification*> mods;
  ptr->searchTerminalModifications(mods, "Label:18O(1)", ResidueModification::ANYWHERE);

  TEST_EQUAL(mods.size(), 4)
  ABORT_IF(mods.size() != 4)

  set<const ResidueModification*>::const_iterator mod_it = mods.begin();

  TEST_STRING_EQUAL((*mod_it)->getOrigin(), "Y")
  TEST_STRING_EQUAL((*mod_it)->getId(), "Label:18O(1)")
  TEST_EQUAL((*mod_it)->getTermSpecificity(), ResidueModification::ANYWHERE)
  ++mod_it;

  TEST_STRING_EQUAL((*mod_it)->getOrigin(), "T")
  TEST_STRING_EQUAL((*mod_it)->getId(), "Label:18O(1)")
  TEST_EQUAL((*mod_it)->getTermSpecificity(), ResidueModification::ANYWHERE)
  ++mod_it;

  TEST_STRING_EQUAL((*mod_it)->getOrigin(), "S")
  TEST_STRING_EQUAL((*mod_it)->getId(), "Label:18O(1)")
  TEST_EQUAL((*mod_it)->getTermSpecificity(), ResidueModification::ANYWHERE)
  ++mod_it;

  TEST_STRING_EQUAL((*mod_it)->getOrigin(), "C-term")
  TEST_STRING_EQUAL((*mod_it)->getId(), "Label:18O(1)")
  TEST_EQUAL((*mod_it)->getTermSpecificity(), ResidueModification::C_TERM)

  mods.clear();
  ptr->searchTerminalModifications(mods, "Label:18O(1)", ResidueModification::C_TERM);

  TEST_EQUAL(mods.size(), 1)
  ABORT_IF(mods.size() != 1)

  mod_it = mods.begin();
  TEST_STRING_EQUAL((*mod_it)->getOrigin(), "C-term")
  TEST_STRING_EQUAL((*mod_it)->getId(), "Label:18O(1)")
  TEST_EQUAL((*mod_it)->getTermSpecificity(), ResidueModification::C_TERM)

  mods.clear();
  ptr->searchTerminalModifications(mods, "Label:18O(1)", ResidueModification::N_TERM);

  TEST_EQUAL(mods.size(), 0)
  ABORT_IF(mods.size() != 0)
}
END_SECTION


START_SECTION((void getTerminalModificationsByDiffMonoMass(std::vector< String > &mods, DoubleReal mass, DoubleReal error, ResidueModification::Term_Specificity term_spec)))
	vector<String> mods;
	ptr->getTerminalModificationsByDiffMonoMass(mods, 42, 0.1, ResidueModification::N_TERM);
	set<String> uniq_mods;
	for (vector<String>::const_iterator it = mods.begin(); it != mods.end(); ++it)
	{
		uniq_mods.insert(*it);
	}
	TEST_EQUAL(mods.size(), 16)
	TEST_EQUAL(uniq_mods.size(), 16)
	TEST_EQUAL(uniq_mods.find("Acetyl (N-term)") != uniq_mods.end(), true)

END_SECTION

START_SECTION(const ResidueModification& getModification(const String & modification) const)
	TEST_EQUAL(ptr->getModification("Carboxymethyl (C)").getFullId(), "Carboxymethyl (C)")
	TEST_EQUAL(ptr->getModification("Carboxymethyl (C)").getId(), "Carboxymethyl")
END_SECTION

START_SECTION((const ResidueModification& getModification(const String &residue_name, const String &mod_name, ResidueModification::Term_Specificity term_spec) const))
	TEST_EQUAL(ptr->getModification("S", "Phosphorylation", ResidueModification::ANYWHERE).getId(), "Phospho")
	TEST_EQUAL(ptr->getModification("S", "Phosphorylation", ResidueModification::ANYWHERE).getFullId(), "Phospho (S)")
END_SECTION

START_SECTION(Size findModificationIndex(const String &mod_name) const)
	Size index = numeric_limits<Size>::max();
	index = ptr->findModificationIndex("Phospho (T)");
	TEST_NOT_EQUAL(index, numeric_limits<Size>::max())
END_SECTION
    
START_SECTION(void getModificationsByDiffMonoMass(std::vector< String > &mods, DoubleReal mass, DoubleReal error=0.0))
	vector<String> mods;
	ptr->getModificationsByDiffMonoMass(mods, 80.0, 0.1);
	set<String> uniq_mods;
	for (vector<String>::const_iterator it = mods.begin(); it != mods.end(); ++it)
	{
		uniq_mods.insert(*it);
	}

	TEST_EQUAL(uniq_mods.find("Phospho (S)") != uniq_mods.end(), true)
	TEST_EQUAL(uniq_mods.find("Phospho (T)") != uniq_mods.end(), true)
	TEST_EQUAL(uniq_mods.find("Phospho (Y)") != uniq_mods.end(), true)
	TEST_EQUAL(uniq_mods.find("Sulfo (S)") != uniq_mods.end(), true)
END_SECTION

START_SECTION(void readFromOBOFile(const String &filename))
	// implicitely tested above
	NOT_TESTABLE
END_SECTION

START_SECTION(void readFromUnimodXMLFile(const String &filename))
	// just provided for convenience at the moment
	NOT_TESTABLE
END_SECTION

START_SECTION((const ResidueModification& getTerminalModification(const String &name, ResidueModification::Term_Specificity term_spec) const))
	TEST_EQUAL(ptr->getTerminalModification("NIC", ResidueModification::N_TERM).getId(), "NIC")
	TEST_EQUAL(ptr->getTerminalModification("NIC", ResidueModification::N_TERM).getFullId(), "NIC (N-term)")
	TEST_EQUAL(ptr->getTerminalModification("Acetyl", ResidueModification::N_TERM).getFullId(), "Acetyl (N-term)")	
END_SECTION

START_SECTION((void getModificationsByDiffMonoMass(std::vector< String > &mods, const String &residue, DoubleReal mass, DoubleReal error=0.0)))
	vector<String> mods;
	ptr->getModificationsByDiffMonoMass(mods, "S", 80.0, 0.1);
	TEST_EQUAL(find(mods.begin(), mods.end(), "Phospho (S)") != mods.end(), true)
	TEST_EQUAL(find(mods.begin(), mods.end(), "Sulfo (S)") != mods.end(), true)
END_SECTION

START_SECTION((void getAllSearchModifications(std::vector< String > &modifications)))
	vector<String> mods;
	ptr->getAllSearchModifications(mods);
	TEST_EQUAL(find(mods.begin(), mods.end(), "Phospho (S)") != mods.end(), true)
	TEST_EQUAL(find(mods.begin(), mods.end(), "Sulfo (S)") != mods.end(), true)
	TEST_EQUAL(find(mods.begin(), mods.end(), "NIC (N-term)") != mods.end(), true)
	TEST_EQUAL(find(mods.begin(), mods.end(), "Phospho") != mods.end(), false)
	TEST_EQUAL(find(mods.begin(), mods.end(), "Dehydrated (N-term C)") != mods.end(), true)
END_SECTION



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST




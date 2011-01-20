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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/CHEMISTRY/ModificationDefinitionsSet.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ModificationDefinitionsSet, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ModificationDefinitionsSet* ptr = 0;
START_SECTION(ModificationDefinitionsSet())
{
	ptr = new ModificationDefinitionsSet();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION((ModificationDefinitionsSet(const ModificationDefinitionsSet &rhs)))
{
  ModificationDefinitionsSet mod_set;
	mod_set.setMaxModifications(2);
	ModificationDefinition mod_def, mod_def2;
	mod_def.setModification("Phospho (S)");
	mod_def.setFixedModification(true);
	mod_def2.setModification("Phospho (T)");
	mod_def2.setFixedModification(false);
	mod_def2.setMaxOccurences(10);
	ModificationDefinitionsSet mod_set2(mod_set);

	TEST_EQUAL(mod_set == mod_set2, true)
}
END_SECTION

START_SECTION((ModificationDefinitionsSet(const String &fixed_modifications, const String &variable_modifications="")))
{
  ModificationDefinitionsSet mod_set("Phospho (S),Phospho (T),Phospho (Y)", "Carbamidomethyl (C)");
	set<String> fixed_mods;
	fixed_mods.insert("Phospho (S)");
	fixed_mods.insert("Phospho (T)");
	fixed_mods.insert("Phospho (Y)");

	set<String> var_mods;
	var_mods.insert("Carbamidomethyl (C)");

	TEST_EQUAL(mod_set.getFixedModificationNames() == fixed_mods, true)
	TEST_EQUAL(mod_set.getVariableModificationNames() == var_mods, true)
}
END_SECTION

START_SECTION((virtual ~ModificationDefinitionsSet()))
{
  delete ptr;
}
END_SECTION

START_SECTION((void setMaxModifications(Size max_mod)))
{
  ModificationDefinitionsSet mod_set;
	mod_set.setMaxModifications(1);
	TEST_EQUAL(mod_set.getMaxModifications(), 1)
	mod_set.setMaxModifications(2);
	TEST_EQUAL(mod_set.getMaxModifications(), 2)
}
END_SECTION

START_SECTION((Size getMaxModifications() const ))
{
  // tested above
	NOT_TESTABLE
}
END_SECTION

START_SECTION((Size getNumberOfModifications() const ))
{
  ModificationDefinitionsSet mod_set("Phospho (S),Phospho (T),Phospho (Y)", "Carbamidomethyl (C)");
	TEST_EQUAL(mod_set.getNumberOfModifications(), 4)
	ModificationDefinitionsSet mod_set2("", "Carbamidomethyl (C)");
	TEST_EQUAL(mod_set2.getNumberOfModifications(), 1)

	ModificationDefinitionsSet mod_set3("Phospho (S)");
	TEST_EQUAL(mod_set3.getNumberOfModifications(), 1)
}
END_SECTION

START_SECTION((Size getNumberOfFixedModifications() const ))
{
  ModificationDefinitionsSet mod_set("Phospho (S),Phospho (T),Phospho (Y)", "Carbamidomethyl (C)");
  TEST_EQUAL(mod_set.getNumberOfFixedModifications(), 3)
  ModificationDefinitionsSet mod_set2("", "Carbamidomethyl (C)");
  TEST_EQUAL(mod_set2.getNumberOfFixedModifications(), 0)

  ModificationDefinitionsSet mod_set3("Phospho (S)");
  TEST_EQUAL(mod_set3.getNumberOfFixedModifications(), 1)
}
END_SECTION

START_SECTION((Size getNumberOfVariableModifications() const ))
{
  ModificationDefinitionsSet mod_set("Phospho (S),Phospho (T)", "Carbamidomethyl (C),Phospho (Y)");
  TEST_EQUAL(mod_set.getNumberOfVariableModifications(), 2)
  ModificationDefinitionsSet mod_set2("", "Carbamidomethyl (C)");
  TEST_EQUAL(mod_set2.getNumberOfVariableModifications(), 1)

  ModificationDefinitionsSet mod_set3("Phospho (S)");
  TEST_EQUAL(mod_set3.getNumberOfVariableModifications(), 0)
}
END_SECTION

START_SECTION((void addModification(const ModificationDefinition& mod_def)))
{
  ModificationDefinition mod_def;
	mod_def.setModification("Phospho (Y)");
	mod_def.setFixedModification(true);


	ModificationDefinitionsSet mod_set;
	mod_set.addModification(mod_def);

	ModificationDefinitionsSet mod_set2;

	TEST_EQUAL(mod_set != mod_set2, true)

	TEST_EQUAL(mod_set.getNumberOfModifications(), 1)
	TEST_EQUAL(mod_set.getNumberOfFixedModifications(), 1)
	TEST_EQUAL(mod_set.getNumberOfVariableModifications(), 0)

	ModificationDefinition mod_def3;
	mod_def3.setModification("Phospho (T)");
	mod_def3.setFixedModification(false);

	ModificationDefinitionsSet mod_set3;
	mod_set3.addModification(mod_def3);
	mod_set3.addModification(mod_def);

	TEST_EQUAL(mod_set3.getNumberOfModifications(), 2)
	TEST_EQUAL(mod_set3.getNumberOfFixedModifications(), 1)
	TEST_EQUAL(mod_set3.getNumberOfVariableModifications(), 1)
}
END_SECTION

START_SECTION((void setModifications(const std::set<ModificationDefinition>& mod_defs)))
{
  ModificationDefinition mod_def1, mod_def2;
	mod_def1.setModification("Phospho (T)");
	mod_def1.setFixedModification(true);
	mod_def2.setModification("Phospho (S)");
	mod_def2.setFixedModification(false);
	set<ModificationDefinition> mod_defs;
	mod_defs.insert(mod_def1);
	mod_defs.insert(mod_def2);

	ModificationDefinitionsSet mod_set;
	mod_set.setModifications(mod_defs);
	TEST_EQUAL(mod_set.getNumberOfModifications(), 2)
	TEST_EQUAL(mod_set.getNumberOfFixedModifications(), 1)
	TEST_EQUAL(mod_set.getNumberOfVariableModifications(), 1)
}
END_SECTION

START_SECTION((void setModifications(const String& fixed_modifications, const String& variable_modifications)))
{
  ModificationDefinitionsSet mod_set1("Phospho (S),Phospho (T),Phospho (Y)", "Carbamidomethyl (C)");
	ModificationDefinitionsSet mod_set2;
	mod_set2.setModifications("Phospho (S),Phospho (T),Phospho (Y)", "Carbamidomethyl (C)");

	TEST_EQUAL(mod_set1.getFixedModificationNames() == mod_set2.getFixedModificationNames(), true)
	TEST_EQUAL(mod_set1.getVariableModificationNames() == mod_set2.getVariableModificationNames(), true)
	TEST_EQUAL(mod_set1.getModificationNames() == mod_set2.getModificationNames(), true)
	TEST_EQUAL(mod_set1 == mod_set2, true)

	mod_set1.setModifications("Phospho (S)", "Carbamidomethyl (C)");
	TEST_EQUAL(mod_set1.getNumberOfModifications(), 2)
	TEST_EQUAL(mod_set1.getNumberOfFixedModifications(), 1)
	TEST_EQUAL(mod_set1.getNumberOfVariableModifications(), 1)
}
END_SECTION

START_SECTION((std::set<ModificationDefinition> getModifications() const ))
{
  ModificationDefinitionsSet mod_set1("Phospho (S),Phospho (T),Phospho (Y)", "Carbamidomethyl (C)");
	set<String> fixed_mods, var_mods;
	fixed_mods.insert("Phospho (S)");
	fixed_mods.insert("Phospho (T)");
	fixed_mods.insert("Phospho (Y)");
	var_mods.insert("Carbamidomethyl (C)");

	set<ModificationDefinition> mod_defs = mod_set1.getModifications();
	for (set<ModificationDefinition>::const_iterator it = mod_defs.begin(); it != mod_defs.end(); ++it)
	{
		if (it->isFixedModification())
		{
			TEST_EQUAL(fixed_mods.find(it->getModification()) != fixed_mods.end(), true)
		}
		else
		{
			TEST_EQUAL(var_mods.find(it->getModification()) != var_mods.end(), true)
		}
	}

}
END_SECTION

START_SECTION(const std::set<ModificationDefinition>& getFixedModifications() const)
  ModificationDefinitionsSet mod_set1("Phospho (S),Phospho (T),Phospho (Y)", "Carbamidomethyl (C)");
  set<String> fixed_mods;
  fixed_mods.insert("Phospho (S)");
  fixed_mods.insert("Phospho (T)");
  fixed_mods.insert("Phospho (Y)");
	
	set<ModificationDefinition> mod_defs = mod_set1.getFixedModifications();
	TEST_EQUAL(mod_defs.size(), 3)
  for (set<ModificationDefinition>::const_iterator it = mod_defs.begin(); it != mod_defs.end(); ++it)
  {
		TEST_EQUAL(it->isFixedModification(), true)
		TEST_EQUAL(fixed_mods.find(it->getModification()) != fixed_mods.end(), true)
  }
END_SECTION

START_SECTION(const std::set<ModificationDefinition>& getVariableModifications() const)
	ModificationDefinitionsSet mod_set1("Phospho (S),Phospho (T),Phospho (Y)", "Carbamidomethyl (C),Phospho (S)");
  set<String> mods;
  mods.insert("Phospho (S)");
  mods.insert("Carbamidomethyl (C)");

  set<ModificationDefinition> mod_defs = mod_set1.getVariableModifications();
  TEST_EQUAL(mod_defs.size(), 2)
  for (set<ModificationDefinition>::const_iterator it = mod_defs.begin(); it != mod_defs.end(); ++it)
  {
    TEST_EQUAL(it->isFixedModification(), false)
    TEST_EQUAL(mods.find(it->getModification()) != mods.end(), true)
  }
END_SECTION


START_SECTION((std::set<String> getModificationNames() const ))
{
  ModificationDefinitionsSet mod_set1("Phospho (S),Phospho (T),Phospho (Y)", "Carbamidomethyl (C)");
  set<String> mods;
  mods.insert("Phospho (S)");
  mods.insert("Phospho (T)");
  mods.insert("Phospho (Y)");
  mods.insert("Carbamidomethyl (C)");

	TEST_EQUAL(mod_set1.getModificationNames() == mods, true)
}
END_SECTION

START_SECTION((std::set<String> getFixedModificationNames() const ))
{
  ModificationDefinitionsSet mod_set1("Phospho (S),Phospho (T),Phospho (Y)", "Carbamidomethyl (C)");
  set<String> mods;
  mods.insert("Phospho (S)");
  mods.insert("Phospho (T)");
  mods.insert("Phospho (Y)");
	TEST_EQUAL(mod_set1.getFixedModificationNames() == mods, true)
}
END_SECTION

START_SECTION((std::set<String> getVariableModificationNames() const ))
{
  ModificationDefinitionsSet mod_set1("Phospho (S),Phospho (T)", "Phospho (Y),Carbamidomethyl (C)");
  set<String> mods;
  mods.insert("Carbamidomethyl (C)");
  mods.insert("Phospho (Y)");

	TEST_EQUAL(mod_set1.getVariableModificationNames() == mods, true)
}
END_SECTION

START_SECTION((ModificationDefinitionsSet& operator=(const ModificationDefinitionsSet& element)))
{
  ModificationDefinitionsSet mod_set1, mod_set2;
	mod_set1.setModifications("Phospho (S),Phospho (T),Phospho (Y)", "");
	TEST_EQUAL(mod_set1 == mod_set2, false)
	mod_set2 = mod_set1;
	TEST_EQUAL(mod_set1 == mod_set2, true)

	mod_set1.setMaxModifications(3);
	TEST_EQUAL(mod_set1 == mod_set2, false)
	mod_set2 = mod_set1;
	TEST_EQUAL(mod_set1 == mod_set2, true)

	mod_set1.setModifications("Phospho (S),Phospho (T),Phospho (Y)", "Carbamidomethyl (C)");
	TEST_EQUAL(mod_set1 == mod_set2, false)
	mod_set2 = mod_set1;
	TEST_EQUAL(mod_set1 == mod_set2, true)
}
END_SECTION

START_SECTION((bool operator==(const ModificationDefinitionsSet& rhs) const))
{
  ModificationDefinitionsSet mod_set1, mod_set2;
  mod_set1.setModifications("Phospho (S),Phospho (T),Phospho (Y)", "");
  TEST_EQUAL(mod_set1 == mod_set2, false)
  mod_set2 = mod_set1;
  TEST_EQUAL(mod_set1 == mod_set2, true)

  mod_set1.setMaxModifications(3);
  TEST_EQUAL(mod_set1 == mod_set2, false)
  mod_set2 = mod_set1;
  TEST_EQUAL(mod_set1 == mod_set2, true)

  mod_set1.setModifications("Phospho (S),Phospho (T),Phospho (Y)", "Carbamidomethyl (C)");
  TEST_EQUAL(mod_set1 == mod_set2, false)
  mod_set2 = mod_set1;
  TEST_EQUAL(mod_set1 == mod_set2, true)
}
END_SECTION

START_SECTION((bool operator!=(const ModificationDefinitionsSet& rhs) const))
{
  ModificationDefinitionsSet mod_set1, mod_set2;
  mod_set1.setModifications("Phospho (S),Phospho (T),Phospho (Y)", "");
  TEST_EQUAL(mod_set1 != mod_set2, true)
  mod_set2 = mod_set1;
  TEST_EQUAL(mod_set1 != mod_set2, false)

  mod_set1.setMaxModifications(3);
  TEST_EQUAL(mod_set1 != mod_set2, true)
  mod_set2 = mod_set1;
  TEST_EQUAL(mod_set1 != mod_set2, false)

  mod_set1.setModifications("Phospho (S),Phospho (T),Phospho (Y)", "Carbamidomethyl (C)");
  TEST_EQUAL(mod_set1 != mod_set2, true)
  mod_set2 = mod_set1;
  TEST_EQUAL(mod_set1 != mod_set2, false)
}
END_SECTION


START_SECTION((ModificationDefinitionsSet(const StringList &fixed_modifications, const StringList &variable_modifications=StringList::create(""))))
  ModificationDefinitionsSet mod_set(StringList::create("Phospho (S),Phospho (T),Phospho (Y)"), StringList::create("Carbamidomethyl (C)"));
  set<String> fixed_mods;
  fixed_mods.insert("Phospho (S)");
  fixed_mods.insert("Phospho (T)");
  fixed_mods.insert("Phospho (Y)");

  set<String> var_mods;
  var_mods.insert("Carbamidomethyl (C)");

  TEST_EQUAL(mod_set.getFixedModificationNames() == fixed_mods, true)
  TEST_EQUAL(mod_set.getVariableModificationNames() == var_mods, true)
END_SECTION


START_SECTION((void setModifications(const StringList &fixed_modifications, const StringList &variable_modifications)))
  ModificationDefinitionsSet mod_set;
  mod_set.setModifications(StringList::create("Phospho (T)"), StringList::create("Phospho (S)"));
  TEST_EQUAL(mod_set.getNumberOfModifications(), 2)
  TEST_EQUAL(mod_set.getNumberOfFixedModifications(), 1)
  TEST_EQUAL(mod_set.getNumberOfVariableModifications(), 1)
END_SECTION

START_SECTION((bool isCompatible(const AASequence &peptide) const))
	ModificationDefinitionsSet mod_set(StringList::create("Carbamidomethyl (C)"), StringList::create("Phospho (S),Phospho (T),Phospho (Y)"));
	AASequence pep1("CCTKPESER");
	AASequence pep2("C(Carbamidomethyl)CTKPESER");
	AASequence pep3("C(Carbamidomethyl)C(Carbamidomethyl)TKPESER");
	AASequence pep4("C(Carbamidomethyl)C(Carbamidomethyl)T(Phospho)TKPESER");
	AASequence pep5("(Acetyl)CCTKPESER");
	AASequence pep6("(Acetyl)C(Carbamidomethyl)C(Carbamidomethyl)TKPES(Phospho)ER");
	AASequence pep7("(Acetyl)C(Carbamidomethyl)C(Carbamidomethyl)T(Phospho)KPES(Phospho)ER");

	TEST_EQUAL(mod_set.isCompatible(pep1), false);
	TEST_EQUAL(mod_set.isCompatible(pep2), false);
	TEST_EQUAL(mod_set.isCompatible(pep3), true);
	TEST_EQUAL(mod_set.isCompatible(pep4), true);
	TEST_EQUAL(mod_set.isCompatible(pep5), false);
	TEST_EQUAL(mod_set.isCompatible(pep6), false);
	TEST_EQUAL(mod_set.isCompatible(pep7), false);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST




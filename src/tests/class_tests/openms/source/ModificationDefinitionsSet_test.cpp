// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/CHEMISTRY/ModificationDefinitionsSet.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ModificationDefinitionsSet, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ModificationDefinitionsSet* ptr = nullptr;
ModificationDefinitionsSet* nullPointer = nullptr;
START_SECTION(ModificationDefinitionsSet())
{
  ptr = new ModificationDefinitionsSet();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((ModificationDefinitionsSet(const ModificationDefinitionsSet& rhs)))
{
  ModificationDefinitionsSet mod_set;
  mod_set.setMaxModifications(2);
  ModificationDefinition mod_def, mod_def2;
  mod_def.setModification("Phospho (S)");
  mod_def.setFixedModification(true);
  mod_def2.setModification("Phospho (T)");
  mod_def2.setFixedModification(false);
  mod_def2.setMaxOccurrences(10);
  ModificationDefinitionsSet mod_set2(mod_set);

  TEST_TRUE(mod_set == mod_set2)
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

START_SECTION((Size getMaxModifications() const))
{
  // tested above
  NOT_TESTABLE
}
END_SECTION

START_SECTION((Size getNumberOfModifications() const))
{
  ModificationDefinitionsSet mod_set(ListUtils::create<String>("Phospho (S),Phospho (T),Phospho (Y)"), ListUtils::create<String>("Carbamidomethyl (C)"));
  TEST_EQUAL(mod_set.getNumberOfModifications(), 4)
  ModificationDefinitionsSet mod_set2(ListUtils::create<String>(""), ListUtils::create<String>("Carbamidomethyl (C)"));
  TEST_EQUAL(mod_set2.getNumberOfModifications(), 1)

  ModificationDefinitionsSet mod_set3(ListUtils::create<String>("Phospho (S)"), ListUtils::create<String>(""));
  TEST_EQUAL(mod_set3.getNumberOfModifications(), 1)
}
END_SECTION

START_SECTION((Size getNumberOfFixedModifications() const))
{
  ModificationDefinitionsSet mod_set(ListUtils::create<String>("Phospho (S),Phospho (T),Phospho (Y)"), ListUtils::create<String>("Carbamidomethyl (C)"));
  TEST_EQUAL(mod_set.getNumberOfFixedModifications(), 3)
  ModificationDefinitionsSet mod_set2(ListUtils::create<String>(""), ListUtils::create<String>("Carbamidomethyl (C)"));
  TEST_EQUAL(mod_set2.getNumberOfFixedModifications(), 0)

  ModificationDefinitionsSet mod_set3(ListUtils::create<String>("Phospho (S)"), ListUtils::create<String>(""));
  TEST_EQUAL(mod_set3.getNumberOfFixedModifications(), 1)
}
END_SECTION

START_SECTION((Size getNumberOfVariableModifications() const))
{
  ModificationDefinitionsSet mod_set(ListUtils::create<String>("Phospho (S),Phospho (T)"), ListUtils::create<String>("Carbamidomethyl (C),Phospho (Y)"));
  TEST_EQUAL(mod_set.getNumberOfVariableModifications(), 2)
  ModificationDefinitionsSet mod_set2(ListUtils::create<String>(""), ListUtils::create<String>("Carbamidomethyl (C)"));
  TEST_EQUAL(mod_set2.getNumberOfVariableModifications(), 1)

  ModificationDefinitionsSet mod_set3(ListUtils::create<String>("Phospho (S)"), ListUtils::create<String>(""));
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

  TEST_FALSE(mod_set == mod_set2)

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
  ModificationDefinitionsSet mod_set1(ListUtils::create<String>("Phospho (S),Phospho (T),Phospho (Y)"), ListUtils::create<String>("Carbamidomethyl (C)"));
  ModificationDefinitionsSet mod_set2;
  mod_set2.setModifications("Phospho (S),Phospho (T),Phospho (Y)", "Carbamidomethyl (C)");

  TEST_EQUAL(mod_set1.getFixedModificationNames() == mod_set2.getFixedModificationNames(), true)
  TEST_EQUAL(mod_set1.getVariableModificationNames() == mod_set2.getVariableModificationNames(), true)
  TEST_EQUAL(mod_set1.getModificationNames() == mod_set2.getModificationNames(), true)
  TEST_TRUE(mod_set1 == mod_set2)

  mod_set1.setModifications("Phospho (S)", "Carbamidomethyl (C)");
  TEST_EQUAL(mod_set1.getNumberOfModifications(), 2)
  TEST_EQUAL(mod_set1.getNumberOfFixedModifications(), 1)
  TEST_EQUAL(mod_set1.getNumberOfVariableModifications(), 1)
}
END_SECTION

START_SECTION((std::set<ModificationDefinition> getModifications() const))
{
  ModificationDefinitionsSet mod_set1(ListUtils::create<String>("Phospho (S),Phospho (T),Phospho (Y)"), ListUtils::create<String>("Carbamidomethyl (C)"));
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
      TEST_EQUAL(fixed_mods.find(it->getModificationName()) != fixed_mods.end(), true)
    }
    else
    {
      TEST_EQUAL(var_mods.find(it->getModificationName()) != var_mods.end(), true)
    }
  }

}
END_SECTION

START_SECTION(const std::set<ModificationDefinition>& getFixedModifications() const)
{
  ModificationDefinitionsSet mod_set1(ListUtils::create<String>("Phospho (S),Phospho (T),Phospho (Y)"), ListUtils::create<String>("Carbamidomethyl (C)"));
  set<String> fixed_mods;
  fixed_mods.insert("Phospho (S)");
  fixed_mods.insert("Phospho (T)");
  fixed_mods.insert("Phospho (Y)");

  set<ModificationDefinition> mod_defs = mod_set1.getFixedModifications();
  TEST_EQUAL(mod_defs.size(), 3)
  for (set<ModificationDefinition>::const_iterator it = mod_defs.begin(); it != mod_defs.end(); ++it)
  {
    TEST_EQUAL(it->isFixedModification(), true)
    TEST_EQUAL(fixed_mods.find(it->getModificationName()) != fixed_mods.end(), true)
  }
}
END_SECTION

START_SECTION(const std::set<ModificationDefinition>& getVariableModifications() const)
{
  ModificationDefinitionsSet mod_set1(ListUtils::create<String>("Phospho (S),Phospho (T),Phospho (Y)"), ListUtils::create<String>("Carbamidomethyl (C),Phospho (S)"));
  set<String> mods;
  mods.insert("Phospho (S)");
  mods.insert("Carbamidomethyl (C)");

  set<ModificationDefinition> mod_defs = mod_set1.getVariableModifications();
  TEST_EQUAL(mod_defs.size(), 2)
  for (set<ModificationDefinition>::const_iterator it = mod_defs.begin(); it != mod_defs.end(); ++it)
  {
    TEST_EQUAL(it->isFixedModification(), false)
    TEST_EQUAL(mods.find(it->getModificationName()) != mods.end(), true)
  }
}
END_SECTION

START_SECTION((std::set<String> getModificationNames() const))
{
  ModificationDefinitionsSet mod_set1(ListUtils::create<String>("Phospho (S),Phospho (T),Phospho (Y)"), ListUtils::create<String>("Carbamidomethyl (C)"));
  set<String> mods;
  mods.insert("Phospho (S)");
  mods.insert("Phospho (T)");
  mods.insert("Phospho (Y)");
  mods.insert("Carbamidomethyl (C)");

  TEST_EQUAL(mod_set1.getModificationNames() == mods, true)
}
END_SECTION

START_SECTION((void getModificationNames(StringList& fixed_modifications, StringList& variable_modifications) const ))
{
  StringList fixed_mods = ListUtils::create<String>("Phospho (S),Phospho (T),Phospho (Y)");
  StringList var_mods = ListUtils::create<String>("Carbamidomethyl (C)");
  ModificationDefinitionsSet mod_set(fixed_mods, var_mods);

  StringList fixed_mods_out, var_mods_out;
  mod_set.getModificationNames(fixed_mods_out, var_mods_out);

  TEST_STRING_EQUAL(ListUtils::concatenate<String>(fixed_mods, ","), 
                    ListUtils::concatenate<String>(fixed_mods_out, ","));
  TEST_STRING_EQUAL(ListUtils::concatenate<String>(var_mods, ","), 
                    ListUtils::concatenate<String>(var_mods_out, ","));
}
END_SECTION

START_SECTION((std::set<String> getFixedModificationNames() const))
{
  ModificationDefinitionsSet mod_set1(ListUtils::create<String>("Phospho (S),Phospho (T),Phospho (Y)"), ListUtils::create<String>("Carbamidomethyl (C)"));
  set<String> mods;
  mods.insert("Phospho (S)");
  mods.insert("Phospho (T)");
  mods.insert("Phospho (Y)");
  TEST_EQUAL(mod_set1.getFixedModificationNames() == mods, true)
}
END_SECTION

START_SECTION((std::set<String> getVariableModificationNames() const))
{
  ModificationDefinitionsSet mod_set1(ListUtils::create<String>("Phospho (S),Phospho (T)"), ListUtils::create<String>("Phospho (Y),Carbamidomethyl (C)"));
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
  TEST_TRUE(mod_set1 == mod_set2)

  mod_set1.setMaxModifications(3);
  TEST_EQUAL(mod_set1 == mod_set2, false)
  mod_set2 = mod_set1;
  TEST_TRUE(mod_set1 == mod_set2)

  mod_set1.setModifications("Phospho (S),Phospho (T),Phospho (Y)", "Carbamidomethyl (C)");
  TEST_EQUAL(mod_set1 == mod_set2, false)
  mod_set2 = mod_set1;
  TEST_TRUE(mod_set1 == mod_set2)
}
END_SECTION

START_SECTION((bool operator==(const ModificationDefinitionsSet& rhs) const))
{
  ModificationDefinitionsSet mod_set1, mod_set2;
  mod_set1.setModifications("Phospho (S),Phospho (T),Phospho (Y)", "");
  TEST_EQUAL(mod_set1 == mod_set2, false)
  mod_set2 = mod_set1;
  TEST_TRUE(mod_set1 == mod_set2)

  mod_set1.setMaxModifications(3);
  TEST_EQUAL(mod_set1 == mod_set2, false)
  mod_set2 = mod_set1;
  TEST_TRUE(mod_set1 == mod_set2)

  mod_set1.setModifications("Phospho (S),Phospho (T),Phospho (Y)", "Carbamidomethyl (C)");
  TEST_EQUAL(mod_set1 == mod_set2, false)
  mod_set2 = mod_set1;
  TEST_TRUE(mod_set1 == mod_set2)
}
END_SECTION

START_SECTION((bool operator!=(const ModificationDefinitionsSet& rhs) const))
{
  ModificationDefinitionsSet mod_set1, mod_set2;
  mod_set1.setModifications("Phospho (S),Phospho (T),Phospho (Y)", "");
  TEST_FALSE(mod_set1 == mod_set2)
  mod_set2 = mod_set1;
  TEST_EQUAL(mod_set1 != mod_set2, false)

  mod_set1.setMaxModifications(3);
  TEST_FALSE(mod_set1 == mod_set2)
  mod_set2 = mod_set1;
  TEST_EQUAL(mod_set1 != mod_set2, false)

  mod_set1.setModifications("Phospho (S),Phospho (T),Phospho (Y)", "Carbamidomethyl (C)");
  TEST_FALSE(mod_set1 == mod_set2)
  mod_set2 = mod_set1;
  TEST_EQUAL(mod_set1 != mod_set2, false)
}
END_SECTION

START_SECTION((ModificationDefinitionsSet(const StringList &fixed_modifications, const StringList &variable_modifications)))
{
  ModificationDefinitionsSet mod_set(ListUtils::create<String>("Phospho (S),Phospho (T),Phospho (Y)"), ListUtils::create<String>("Carbamidomethyl (C)"));
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

START_SECTION((void setModifications(const StringList &fixed_modifications, const StringList &variable_modifications)))
{
  ModificationDefinitionsSet mod_set;
  mod_set.setModifications(ListUtils::create<String>("Phospho (T)"), ListUtils::create<String>("Phospho (S)"));
  TEST_EQUAL(mod_set.getNumberOfModifications(), 2)
  TEST_EQUAL(mod_set.getNumberOfFixedModifications(), 1)
  TEST_EQUAL(mod_set.getNumberOfVariableModifications(), 1)
}
END_SECTION

START_SECTION((bool isCompatible(const AASequence &peptide) const))
{
  ModificationDefinitionsSet mod_set(ListUtils::create<String>("Carbamidomethyl (C)"), ListUtils::create<String>("Phospho (S),Phospho (T),Phospho (Y)"));
  AASequence pep1 = AASequence::fromString("CCTKPESER");
  AASequence pep2 = AASequence::fromString("C(Carbamidomethyl)CTKPESER");
  AASequence pep3 = AASequence::fromString("C(Carbamidomethyl)C(Carbamidomethyl)TKPESER");
  AASequence pep4 = AASequence::fromString("C(Carbamidomethyl)C(Carbamidomethyl)T(Phospho)TKPESER");
  AASequence pep5 = AASequence::fromString("(Acetyl)CCTKPESER");
  AASequence pep6 = AASequence::fromString("(Acetyl)C(Carbamidomethyl)C(Carbamidomethyl)TKPES(Phospho)ER");
  AASequence pep7 = AASequence::fromString("(Acetyl)C(Carbamidomethyl)C(Carbamidomethyl)T(Phospho)KPES(Phospho)ER");

  TEST_EQUAL(mod_set.isCompatible(pep1), false);
  TEST_EQUAL(mod_set.isCompatible(pep2), false);
  TEST_EQUAL(mod_set.isCompatible(pep3), true);
  TEST_EQUAL(mod_set.isCompatible(pep4), true);
  TEST_EQUAL(mod_set.isCompatible(pep5), false);
  TEST_EQUAL(mod_set.isCompatible(pep6), false);
  TEST_EQUAL(mod_set.isCompatible(pep7), false);
}
END_SECTION

START_SECTION((void findMatches(multimap<double, ModificationDefinition>& matches, double mass, const String& residue, ResidueModification::TermSpecificity term_spec, bool consider_fixed, bool consider_variable, bool is_delta, double tolerance) const))
{
  ModificationDefinitionsSet mod_set;
  mod_set.setModifications("Gln->pyro-Glu (N-term Q)", "Glu->pyro-Glu (N-term E),Oxidation (M)");
  multimap<double, ModificationDefinition> matches;
  // nothing to consider:
  TEST_EXCEPTION(Exception::IllegalArgument, mod_set.findMatches(matches, -18, "E", ResidueModification::N_TERM, false, false, true, 0.1));
  // wrong term. spec.:
  mod_set.findMatches(matches, -18, "E", ResidueModification::ANYWHERE, true, true, true, 0.1);
  TEST_EQUAL(matches.empty(), true);
  // wrong residue:
  mod_set.findMatches(matches, -18, "Q", ResidueModification::N_TERM, true, true, true, 0.1);
  TEST_EQUAL(matches.empty(), true);
  // wrong fixed/variable:
  mod_set.findMatches(matches, -18, "E", ResidueModification::N_TERM, true, false, true, 0.1);
  TEST_EQUAL(matches.empty(), true);
  // residue, low tolerance:
  mod_set.findMatches(matches, -18, "E", ResidueModification::N_TERM, true, true, true, 0.1);
  TEST_EQUAL(matches.size(), 1);
  TEST_EQUAL(matches.begin()->second.getModificationName(), "Glu->pyro-Glu (N-term E)");
  // no residue, low tolerance:
  mod_set.findMatches(matches, -18, "", ResidueModification::N_TERM, true, true, true, 0.1);
  TEST_EQUAL(matches.size(), 1); 
  TEST_EQUAL(matches.begin()->second.getModificationName(), "Glu->pyro-Glu (N-term E)");
  // no residue, high tolerance:
  mod_set.findMatches(matches, -18, "", ResidueModification::N_TERM, true, true, true, 2);
  TEST_EQUAL(matches.size(), 2);
  TEST_EQUAL(matches.begin()->second.getModificationName(), "Glu->pyro-Glu (N-term E)");
  TEST_EQUAL((++matches.begin())->second.getModificationName(), "Gln->pyro-Glu (N-term Q)");
}
END_SECTION

START_SECTION(void inferFromPeptides(const vector<PeptideIdentification>& peptides))
{
  vector<PeptideIdentification> peptides(2);
  PeptideHit hit;
  hit.setSequence(AASequence::fromString("AC(Carbamidomethyl)M"));
  peptides[0].insertHit(hit);
  hit.setSequence(AASequence::fromString("(Acetyl)AEM"));
  peptides[0].insertHit(hit);
  hit.setSequence(AASequence::fromString("AC(Carbamidomethyl)M(Oxidation)"));
  peptides[1].insertHit(hit);

  ModificationDefinitionsSet mod_defs;
  mod_defs.inferFromPeptides(peptides);
  set<String> mods = mod_defs.getFixedModificationNames();
  TEST_EQUAL(mods.size(), 1);
  set<String>::const_iterator it = mods.begin();
  TEST_STRING_EQUAL(*it, "Carbamidomethyl (C)");
  mods = mod_defs.getVariableModificationNames();
  TEST_EQUAL(mods.size(), 2);
  it = mods.begin();
  TEST_STRING_EQUAL(*it, "Acetyl (N-term)");
  ++it;
  TEST_STRING_EQUAL(*it, "Oxidation (M)");
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

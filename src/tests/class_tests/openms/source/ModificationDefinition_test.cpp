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
#include <OpenMS/CHEMISTRY/ModificationDefinition.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <unordered_set>
#include <unordered_map>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ModificationDefinition, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ModificationDefinition* ptr = nullptr;
ModificationDefinition* nullPointer = nullptr;
START_SECTION(ModificationDefinition())
{
  ptr = new ModificationDefinition();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((virtual ~ModificationDefinition()))
{
  delete ptr;
}
END_SECTION

ptr = new ModificationDefinition();

START_SECTION((ModificationDefinition(const ModificationDefinition& rhs)))
{
  ModificationDefinition mod_def;
  mod_def.setFixedModification(true);
  ModificationDefinition copy(mod_def);
  TEST_EQUAL(mod_def.isFixedModification(), copy.isFixedModification())

  mod_def.setFixedModification(false);
  ModificationDefinition copy2(mod_def);
  TEST_EQUAL(mod_def.isFixedModification(), copy2.isFixedModification())
}
END_SECTION

START_SECTION((ModificationDefinition(const String& mod, bool fixed = true, UInt max_occur = 0)))
{
  ModificationDefinition mod1("Acetyl (N-term)");
  TEST_EQUAL(mod1.getModificationName(), "Acetyl (N-term)");
  ModificationDefinition mod2("Oxidation (M)");
  TEST_EQUAL(mod2.getModificationName(), "Oxidation (M)");
  TEST_EQUAL(mod2.isFixedModification(), true);
  TEST_EQUAL(mod2.getMaxOccurrences(), 0);
  ModificationDefinition mod3("Carboxymethyl (C)", false, 2);
  TEST_EQUAL(mod3.getModificationName(), "Carboxymethyl (C)");
  TEST_EQUAL(mod3.isFixedModification(), false);
  TEST_EQUAL(mod3.getMaxOccurrences(), 2);
}
END_SECTION

START_SECTION((ModificationDefinition(const ResidueModification& mod, bool fixed = true, UInt max_occur = 0)))
{
  const ResidueModification res_mod1 = *ModificationsDB::getInstance()->getModification("Acetyl (N-term)");
  ModificationDefinition mod1(res_mod1);
  TEST_EQUAL(mod1.getModificationName(), "Acetyl (N-term)");
  TEST_EQUAL(mod1.isFixedModification(), true);
  TEST_EQUAL(mod1.getMaxOccurrences(), 0);
  const ResidueModification res_mod2 = *ModificationsDB::getInstance()->getModification("Oxidation (M)");
  ModificationDefinition mod2(res_mod2, false, 2);
  TEST_EQUAL(mod2.isFixedModification(), false);
  TEST_EQUAL(mod2.getMaxOccurrences(), 2);
}
END_SECTION

START_SECTION((void setFixedModification(bool fixed)))
{
  ptr->setFixedModification(true);
  TEST_EQUAL(ptr->isFixedModification(), true)
  ptr->setFixedModification(false);
  TEST_EQUAL(ptr->isFixedModification(), false)
}
END_SECTION

START_SECTION((bool isFixedModification() const))
{
  // tested above
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void setMaxOccurrences(UInt num)))
{
  ptr->setMaxOccurrences(1);
  TEST_EQUAL(ptr->getMaxOccurrences(), 1)
  ptr->setMaxOccurrences(1000);
  TEST_EQUAL(ptr->getMaxOccurrences(), 1000)
}
END_SECTION

START_SECTION((UInt getMaxOccurrences() const))
{
  // tested above
  NOT_TESTABLE
}
END_SECTION

START_SECTION((String getModificationName() const))
{
  ModificationDefinition mod1;
  mod1.setModification("Acetyl (N-term)");
  TEST_EQUAL(mod1.getModificationName(), "Acetyl (N-term)")
  mod1.setModification("Oxidation (M)");
  TEST_EQUAL(mod1.getModificationName(), "Oxidation (M)")
}
END_SECTION

START_SECTION((String getModification() const))
{
  const ResidueModification* rm = ModificationsDB::getInstance()->getModification("Acetyl (N-term)");
  ModificationDefinition mod1;
  mod1.setModification(rm->getFullId());
  TEST_EQUAL(rm, &(mod1.getModification()));
}
END_SECTION

START_SECTION((void setModification(const String& modification)))
{
  // tested above
  NOT_TESTABLE
}
END_SECTION

START_SECTION((ModificationDefinition& operator=(const ModificationDefinition& element)))
{
  ModificationDefinition mod_def;
  mod_def.setFixedModification(true);
  *ptr = mod_def;
  TEST_EQUAL(mod_def.isFixedModification(), ptr->isFixedModification())

  mod_def.setFixedModification(false);
  *ptr = mod_def;
  TEST_EQUAL(mod_def.isFixedModification(), ptr->isFixedModification())
  
}
END_SECTION

START_SECTION((bool operator==(const ModificationDefinition& rhs) const))
{
  ModificationDefinition m1, m2;
  TEST_TRUE(m1 == m2)
  m1.setFixedModification(false);
  TEST_EQUAL(m1 == m2, false)
  m1.setFixedModification(true);
  m1.setMaxOccurrences(15);
  TEST_EQUAL(m1 == m2, false)
  m1.setMaxOccurrences(0);
  m1.setModification("Oxidation (M)");
  TEST_EQUAL(m1 == m2, false)
  m2.setModification("Oxidation (M)");
  TEST_TRUE(m1 == m2)
}
END_SECTION

START_SECTION((bool operator!=(const ModificationDefinition& rhs) const))
{
  ModificationDefinition m1, m2;
  TEST_EQUAL(m1 != m2, false)
  m1.setFixedModification(false);
  TEST_FALSE(m1 == m2)
  m1.setFixedModification(true);
  m1.setMaxOccurrences(15);
  TEST_FALSE(m1 == m2)
  m1.setMaxOccurrences(0);
  m1.setModification("Oxidation (M)");
  TEST_FALSE(m1 == m2)
  m2.setModification("Oxidation (M)");
  TEST_EQUAL(m1 != m2, false)
}
END_SECTION

START_SECTION((bool operator<(const OpenMS::ModificationDefinition& rhs) const))
{
  ModificationDefinition m1, m2;
  m1.setModification("Oxidation (M)");
  m2.setModification("Carboxymethyl (C)");
  TEST_EQUAL(m1 < m2, false)
  TEST_EQUAL(m1 < m1, false)
  TEST_EQUAL(m2 < m1, true)
}
END_SECTION

delete ptr;

/////////////////////////////////////////////////////////////
// Hash tests
/////////////////////////////////////////////////////////////

START_SECTION(([EXTRA] std::hash<ModificationDefinition>))
{
  // Test 1: Equal objects must have equal hashes
  ModificationDefinition m1("Oxidation (M)", true, 3);
  ModificationDefinition m2("Oxidation (M)", true, 3);
  TEST_EQUAL(m1 == m2, true)
  TEST_EQUAL(std::hash<ModificationDefinition>{}(m1), std::hash<ModificationDefinition>{}(m2))

  // Test 2: Different modifications should (likely) have different hashes
  ModificationDefinition m3("Acetyl (N-term)", true, 0);
  // Not guaranteed, but extremely likely to be different
  TEST_NOT_EQUAL(std::hash<ModificationDefinition>{}(m1), std::hash<ModificationDefinition>{}(m3))

  // Test 3: Different fixed_modification should produce different hashes
  ModificationDefinition m4("Oxidation (M)", false, 3);
  TEST_EQUAL(m1 == m4, false)
  TEST_NOT_EQUAL(std::hash<ModificationDefinition>{}(m1), std::hash<ModificationDefinition>{}(m4))

  // Test 4: Different max_occurrences should produce different hashes
  ModificationDefinition m5("Oxidation (M)", true, 5);
  TEST_EQUAL(m1 == m5, false)
  TEST_NOT_EQUAL(std::hash<ModificationDefinition>{}(m1), std::hash<ModificationDefinition>{}(m5))

  // Test 5: Hash consistency - same object hashed multiple times
  std::size_t hash1 = std::hash<ModificationDefinition>{}(m1);
  std::size_t hash2 = std::hash<ModificationDefinition>{}(m1);
  TEST_EQUAL(hash1, hash2)
}
END_SECTION

START_SECTION(([EXTRA] ModificationDefinition in std::unordered_set))
{
  std::unordered_set<ModificationDefinition> mod_set;

  ModificationDefinition m1("Oxidation (M)", true, 0);
  ModificationDefinition m2("Acetyl (N-term)", false, 2);
  ModificationDefinition m3("Carboxymethyl (C)", true, 1);

  // Insert modifications
  mod_set.insert(m1);
  mod_set.insert(m2);
  mod_set.insert(m3);
  TEST_EQUAL(mod_set.size(), 3)

  // Insert duplicate - should not increase size
  ModificationDefinition m1_copy("Oxidation (M)", true, 0);
  mod_set.insert(m1_copy);
  TEST_EQUAL(mod_set.size(), 3)

  // Find operations
  TEST_EQUAL(mod_set.find(m1) != mod_set.end(), true)
  TEST_EQUAL(mod_set.find(m2) != mod_set.end(), true)
  TEST_EQUAL(mod_set.count(m1_copy), 1)

  // Non-existent element
  ModificationDefinition m4("Phospho (S)", true, 0);
  TEST_EQUAL(mod_set.find(m4) == mod_set.end(), true)
  TEST_EQUAL(mod_set.count(m4), 0)
}
END_SECTION

START_SECTION(([EXTRA] ModificationDefinition in std::unordered_map))
{
  std::unordered_map<ModificationDefinition, std::string> mod_map;

  ModificationDefinition m1("Oxidation (M)", true, 0);
  ModificationDefinition m2("Acetyl (N-term)", false, 2);

  // Insert key-value pairs
  mod_map[m1] = "oxidation_value";
  mod_map[m2] = "acetyl_value";
  TEST_EQUAL(mod_map.size(), 2)

  // Access by key
  TEST_EQUAL(mod_map[m1], "oxidation_value")
  TEST_EQUAL(mod_map[m2], "acetyl_value")

  // Update value
  mod_map[m1] = "updated_value";
  TEST_EQUAL(mod_map[m1], "updated_value")
  TEST_EQUAL(mod_map.size(), 2)

  // Find using equal object
  ModificationDefinition m1_copy("Oxidation (M)", true, 0);
  auto it = mod_map.find(m1_copy);
  TEST_EQUAL(it != mod_map.end(), true)
  TEST_EQUAL(it->second, "updated_value")
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST




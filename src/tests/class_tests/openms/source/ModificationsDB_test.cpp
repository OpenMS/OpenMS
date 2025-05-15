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
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <limits>
#include <algorithm>
///////////////////////////

using namespace OpenMS;
using namespace std;

struct ResidueModificationOriginCmp
{
  bool operator() (const ResidueModification* a, const ResidueModification* b) const
  {
    return a->getOrigin() < b->getOrigin();
  }
};

START_TEST(ModificationsDB, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(bool ModificationsDB::isInstantiated())
{
  bool instantiated = ModificationsDB::isInstantiated();
  TEST_EQUAL(instantiated, false);
}
END_SECTION

ModificationsDB* ptr = nullptr;
ModificationsDB* nullPointer = nullptr;

START_SECTION(ModificationsDB* getInstance())
{
	ptr = ModificationsDB::getInstance();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(bool ModificationsDB::isInstantiated())
{
  bool instantiated = ModificationsDB::isInstantiated();
  TEST_EQUAL(instantiated, true);
}
END_SECTION

START_SECTION(Size getNumberOfModifications() const)
	// range because data may change over time
	TEST_EQUAL(ptr->getNumberOfModifications() > 100, true);
END_SECTION

START_SECTION(const ResidueModification& getModification(Size index) const)
	TEST_EQUAL(!ptr->getModification(0)->getId().empty(), true)
END_SECTION

  START_SECTION((void searchModifications(std::set<const ResidueModification*>& mods, const String& mod_name, const String& residue, ResidueModification::TermSpecificity term_spec) const))
{
  set<const ResidueModification*> mods;
  ptr->searchModifications(mods, "Phosphorylation", "T", ResidueModification::ANYWHERE);
  TEST_EQUAL(mods.size(), 1);
  TEST_STRING_EQUAL((*mods.begin())->getFullId(), "Phospho (T)");
  // terminal mod:
  ptr->searchModifications(mods, "NIC", "", ResidueModification::N_TERM);
  TEST_EQUAL(mods.size(), 1);

  // Phosphorylation Decoy custom mod
  ptr->searchModifications(mods, "Phosphorylation Decoy", "A", ResidueModification::ANYWHERE); 
  TEST_EQUAL(mods.size(), 1);

  // Phosphorylation Decoy custom search by title  
  ptr->searchModifications(mods, "PhosphoDecoy");
  TEST_EQUAL(mods.size(), 3);
  
  ptr->searchModifications(mods, "Label:18O(1)");
  TEST_EQUAL(mods.size(), 4);
  ABORT_IF(mods.size() != 4);

  // Create a sorted set (sorted by getOrigin) instead of sorted by pointer
  // value -> this is more robust on different platforms.
  set<const ResidueModification*, ResidueModificationOriginCmp> mods_sorted;

  // Copy the mods into our sorted container
  set<const ResidueModification*>::const_iterator mod_it_ = mods.begin();
  for (; mod_it_ != mods.end(); mod_it_++)
  {
    mods_sorted.insert((*mod_it_));
  }

  set<const ResidueModification*, ResidueModificationOriginCmp>::const_iterator mod_it = mods_sorted.begin();

  TEST_EQUAL((*mod_it)->getOrigin(), 'S')
  TEST_STRING_EQUAL((*mod_it)->getId(), "Label:18O(1)")
  TEST_EQUAL((*mod_it)->getTermSpecificity(), ResidueModification::ANYWHERE)
  ++mod_it;

  TEST_EQUAL((*mod_it)->getOrigin(), 'T')
  TEST_STRING_EQUAL((*mod_it)->getId(), "Label:18O(1)")
  TEST_EQUAL((*mod_it)->getTermSpecificity(), ResidueModification::ANYWHERE)
  ++mod_it;

  TEST_EQUAL((*mod_it)->getOrigin(), 'X')
  TEST_STRING_EQUAL((*mod_it)->getId(), "Label:18O(1)")
  TEST_EQUAL((*mod_it)->getTermSpecificity(), ResidueModification::C_TERM)
  ++mod_it;

  TEST_EQUAL((*mod_it)->getOrigin(), 'Y')
  TEST_STRING_EQUAL((*mod_it)->getId(), "Label:18O(1)")
  TEST_EQUAL((*mod_it)->getTermSpecificity(), ResidueModification::ANYWHERE)

  ptr->searchModifications(mods, "Label:18O(1)", "", ResidueModification::C_TERM);

  TEST_EQUAL(mods.size(), 1)
  ABORT_IF(mods.size() != 1)

  // Copy the mods into our sorted container
  mods_sorted.clear();
  for (mod_it_ = mods.begin(); mod_it_ != mods.end(); mod_it_++)
  {
    mods_sorted.insert((*mod_it_));
  }

  mod_it = mods_sorted.begin();
  TEST_EQUAL((*mod_it)->getOrigin(), 'X')
  TEST_STRING_EQUAL((*mod_it)->getId(), "Label:18O(1)")
  TEST_EQUAL((*mod_it)->getTermSpecificity(), ResidueModification::C_TERM)

  // no match, thus mods should be empty
  ptr->searchModifications(mods, "Label:18O(1)", "", ResidueModification::N_TERM);

  TEST_EQUAL(mods.size(), 0)
}
END_SECTION


START_SECTION((void searchModificationsByDiffMonoMass(std::vector<String>& mods, double mass, double max_error, const String& residue, ResidueModification::TermSpecificity term_spec)))
{
  vector<String> mods;
  ptr->searchModificationsByDiffMonoMass(mods, 80.0, 0.1, "S");
  TEST_EQUAL(find(mods.begin(), mods.end(), "Phospho (S)") != mods.end(), true);
  TEST_EQUAL(find(mods.begin(), mods.end(), "Sulfo (S)") != mods.end(), true);

  // something exotic.. mods should return empty (without clearing it before)
  ptr->searchModificationsByDiffMonoMass(mods, 800000000.0, 0.1, "S");
  TEST_EQUAL(mods.size(), 0);

  // terminal mod:
  ptr->searchModificationsByDiffMonoMass(mods, 42, 0.1, "", ResidueModification::N_TERM);
  set<String> uniq_mods;
  for (vector<String>::const_iterator it = mods.begin(); it != mods.end(); ++it)
  {
    uniq_mods.insert(*it);
  }
  TEST_EQUAL(mods.size(), 18);
  TEST_EQUAL(uniq_mods.size(), 18);
  TEST_EQUAL(uniq_mods.find("Acetyl (N-term)") != uniq_mods.end(), true);

  // something exotic.. mods should return empty (without clearing it before)
  ptr->searchModificationsByDiffMonoMass(mods, 4200000, 0.1, "", ResidueModification::N_TERM);
  TEST_EQUAL(mods.size(), 0);

  ptr->searchModificationsByDiffMonoMass(mods, 80.0, 0.1);
  uniq_mods.clear();
  for (vector<String>::const_iterator it = mods.begin(); it != mods.end(); ++it)
  {
    uniq_mods.insert(*it);
  }

  TEST_EQUAL(uniq_mods.find("Phospho (S)") != uniq_mods.end(), true);
  TEST_EQUAL(uniq_mods.find("Phospho (T)") != uniq_mods.end(), true);
  TEST_EQUAL(uniq_mods.find("Phospho (Y)") != uniq_mods.end(), true);
  TEST_EQUAL(uniq_mods.find("Sulfo (S)") != uniq_mods.end(), true);

  // something exotic.. mods should return empty (without clearing it before)
  ptr->searchModificationsByDiffMonoMass(mods, 800000000.0, 0.1);
  TEST_EQUAL(mods.size(), 0);

  // make sure the common ones are also found for integer masses (this is how
  // integer mass search is done)
  mods.clear();
  ptr->searchModificationsByDiffMonoMass(mods, 80.0, 1.0, "S");
  TEST_EQUAL(mods.empty(), false)
  TEST_EQUAL(mods[0], "Phospho (S)")
  mods.clear();
  ptr->searchModificationsByDiffMonoMass(mods, 80.0, 1.0, "T");
  TEST_EQUAL(mods.empty(), false)
  TEST_EQUAL(mods[0], "Phospho (T)")
  mods.clear();
  ptr->searchModificationsByDiffMonoMass(mods, 80.0, 1.0, "Y");
  TEST_EQUAL(mods.empty(), false)
  TEST_EQUAL(mods[0], "Phospho (Y)")
  mods.clear();
  ptr->searchModificationsByDiffMonoMass(mods, 16.0, 1.0, "M");
  TEST_EQUAL(mods.empty(), false)
  TEST_EQUAL(mods[0], "Oxidation (M)")
  ptr->searchModificationsByDiffMonoMass(mods, 0.98, 0.1, "N");
  TEST_EQUAL(mods.empty(), false)
  TEST_EQUAL(mods[0], "Deamidated (N)")
  ptr->searchModificationsByDiffMonoMass(mods, 0.98, 1.0, "Q");
  TEST_EQUAL(mods.empty(), false)
  TEST_EQUAL(mods[0], "Deamidated (Q)")
}
END_SECTION

START_SECTION((const ResidueModification& getModification(const String& mod_name, const String& residue, ResidueModification::TermSpecificity term_spec) const))
{
  TEST_EQUAL(ptr->getModification("Carboxymethyl (C)")->getFullId(), "Carboxymethyl (C)");
  TEST_EQUAL(ptr->getModification("Carboxymethyl (C)")->getId(), "Carboxymethyl");

  TEST_EQUAL(ptr->getModification("Phosphorylation", "S", ResidueModification::ANYWHERE)->getId(), "Phospho");
  TEST_EQUAL(ptr->getModification("Phosphorylation", "S", ResidueModification::ANYWHERE)->getFullId(), "Phospho (S)");

  // terminal mod:
  TEST_EQUAL(ptr->getModification("NIC", "", ResidueModification::N_TERM)->getId(), "NIC");
  TEST_EQUAL(ptr->getModification("NIC", "", ResidueModification::N_TERM)->getFullId(), "NIC (N-term)");
  TEST_EQUAL(ptr->getModification("Acetyl", "", ResidueModification::N_TERM)->getFullId(), "Acetyl (N-term)");

  // missing modification (exception)
  TEST_EXCEPTION(Exception::InvalidValue, ptr->getModification("MISSING"));
  TEST_EXCEPTION(Exception::InvalidValue, ptr->getModification("MISSING", "", ResidueModification::N_TERM));
  TEST_EXCEPTION(Exception::InvalidValue, ptr->getModification("MISSING", "", ResidueModification::C_TERM));
}
END_SECTION

START_SECTION((Size findModificationIndex(const String& mod_name) const))
{
  Size index = numeric_limits<Size>::max();
  index = ptr->findModificationIndex("Phospho (T)");
  TEST_NOT_EQUAL(index, numeric_limits<Size>::max());
}
END_SECTION

START_SECTION(void readFromOBOFile(const String& filename))
	// implicitely tested above
	NOT_TESTABLE
END_SECTION

START_SECTION(void readFromUnimodXMLFile(const String& filename))
	// just provided for convenience at the moment
	NOT_TESTABLE
END_SECTION

START_SECTION((void getAllSearchModifications(std::vector<String>& modifications)))
{
  vector<String> mods;
  ptr->getAllSearchModifications(mods);
  TEST_EQUAL(find(mods.begin(), mods.end(), "Phospho (S)") != mods.end(), true);
  TEST_EQUAL(find(mods.begin(), mods.end(), "Sulfo (S)") != mods.end(), true);
  TEST_EQUAL(find(mods.begin(), mods.end(), "NIC (N-term)") != mods.end(), true);
  TEST_EQUAL(find(mods.begin(), mods.end(), "Phospho") != mods.end(), false);
  TEST_EQUAL(find(mods.begin(), mods.end(), "Dehydrated (N-term C)") != mods.end(), true);

  // repeat search .. return size should be the same
  Size old_size=mods.size();
  ptr->getAllSearchModifications(mods);
  TEST_EQUAL(mods.size(), old_size);
}
END_SECTION

START_SECTION((bool addModification(ResidueModification* modification)))
{
  TEST_EQUAL(ptr->has("Phospho (A)"), false);
  std::unique_ptr<ResidueModification> modification(new ResidueModification());
  modification->setFullId("Phospho (A)");
  ptr->addModification(std::move(modification));
  TEST_EQUAL(ptr->has("Phospho (A)"), true);
}
END_SECTION

START_SECTION([EXTRA] multithreaded example)
{
  // All measurements are best of three (wall time, Linux, 8 threads)
  //
  // Serial execution of code:
  // 1e6 iterations -> 6.36 seconds
  // Parallel execution of code:
  // 1e6 iterations -> 9.79 seconds with boost::shared_mutex
  // 1e6 iterations -> 6.28 seconds with std::mutex
  // 1e6 iterations -> 4.64 seconds with pragma critical

   static ModificationsDB* mdb = ModificationsDB::getInstance();

   int nr_iterations (1e4), test (0);
#pragma omp parallel for reduction (+: test)
  for (int k = 1; k < nr_iterations + 1; k++)
  {
    int mod_id = k;
    String modname = "mod" + String(mod_id);
    std::unique_ptr<ResidueModification> new_mod(new ResidueModification());
    new_mod->setFullId(modname);
    new_mod->setMonoMass( 0.11 * mod_id);
    new_mod->setAverageMass(1.0);
    new_mod->setDiffMonoMass( 0.05 * mod_id);
    mdb->addModification(std::move(new_mod));
    int tmp = (int)mdb->getModification(modname)->getAverageMass();
    test += tmp;
  }
  TEST_EQUAL(test, nr_iterations*1.0)

   // Every modification is the same
  test = 0;
  #pragma omp parallel for reduction (+: test)
  for (int k = 1; k < nr_iterations + 1; k++)
  {
    int mod_id = 42;
    String modname = "mod" + String(mod_id);
    if (!mdb->has(modname))
    {
      std::unique_ptr<ResidueModification> new_mod(new ResidueModification());
      new_mod->setFullId(modname);
      new_mod->setMonoMass( 0.11 * mod_id);
      new_mod->setAverageMass(1.0);
      new_mod->setDiffMonoMass( 0.05 * mod_id);
      mdb->addModification(std::move(new_mod));
    }
    int tmp = (int)mdb->getModification(modname)->getAverageMass();
    test += tmp;
  }
  TEST_EQUAL(test, nr_iterations*1.0)

 }
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/CHEMISTRY/CrossLinksDB.h>
#include <limits>
///////////////////////////

using namespace OpenMS;
using namespace std;

struct ResidueModificationOriginCmp
{
  bool operator() (const ResidueModification* a, const ResidueModification* b) const
  {
    if (a->getOrigin() == b->getOrigin())
    {
      return a->getTermSpecificity() < b->getTermSpecificity();
    }
    return a->getOrigin() < b->getOrigin();
  }
};

START_TEST(CrossLinksDB, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

CrossLinksDB* ptr = nullptr;
CrossLinksDB* nullPointer = nullptr;
START_SECTION(CrossLinksDB* getInstance())
{
	ptr = CrossLinksDB::getInstance();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(Size getNumberOfModifications() const)
	// range because data may change over time
	TEST_EQUAL(ptr->getNumberOfModifications() > 10, true);
END_SECTION

START_SECTION(const ResidueModification& getModification(Size index) const)
        TEST_EQUAL(!ptr->getModification(0)->getId().empty(), true)
END_SECTION

START_SECTION((void searchModifications(std::set<const ResidueModification*>& mods, const String& mod_name, const String& residue, ResidueModification::TermSpecificity term_spec) const))
{
  set<const ResidueModification*> mods;
  ptr->searchModifications(mods, "DSS", "K", ResidueModification::ANYWHERE);
  TEST_EQUAL(mods.size(), 1);
  TEST_STRING_EQUAL((*mods.begin())->getFullId(), "DSS (K)");
  // terminal mod:
  ptr->searchModifications(mods, "DSS", "", ResidueModification::N_TERM);
  TEST_EQUAL(mods.size(), 1);

  ptr->searchModifications(mods, "EDC");

  TEST_EQUAL(mods.size(), 8);
  ABORT_IF(mods.size() != 8);

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

  // EDC is a heterobifunctional cross-linker: one reactive site binds to C-term, D and E, the other to N-term, K, S, T
  // no distinction between the two sites implemented in CrossLinksDB or ResidueModification, the search engine has to take care of that for now
  TEST_EQUAL((*mod_it)->getOrigin(), 'D')
  TEST_STRING_EQUAL((*mod_it)->getId(), "EDC")
  TEST_EQUAL((*mod_it)->getTermSpecificity(), ResidueModification::ANYWHERE)
  ++mod_it;

  TEST_EQUAL((*mod_it)->getOrigin(), 'E')
  TEST_STRING_EQUAL((*mod_it)->getId(), "EDC")
  TEST_EQUAL((*mod_it)->getTermSpecificity(), ResidueModification::ANYWHERE)
  ++mod_it;

  TEST_EQUAL((*mod_it)->getOrigin(), 'K')
  TEST_STRING_EQUAL((*mod_it)->getId(), "EDC")
  TEST_EQUAL((*mod_it)->getTermSpecificity(), ResidueModification::ANYWHERE)
  ++mod_it;

  TEST_EQUAL((*mod_it)->getOrigin(), 'S')
  TEST_STRING_EQUAL((*mod_it)->getId(), "EDC")
  TEST_EQUAL((*mod_it)->getTermSpecificity(), ResidueModification::ANYWHERE)
  ++mod_it;

  TEST_EQUAL((*mod_it)->getOrigin(), 'T')
  TEST_STRING_EQUAL((*mod_it)->getId(), "EDC")
  TEST_EQUAL((*mod_it)->getTermSpecificity(), ResidueModification::ANYWHERE)
  ++mod_it;

  TEST_EQUAL((*mod_it)->getOrigin(), 'X')
  TEST_STRING_EQUAL((*mod_it)->getId(), "EDC")
  TEST_EQUAL((*mod_it)->getTermSpecificity(), ResidueModification::C_TERM)
  ++mod_it;

  TEST_EQUAL((*mod_it)->getOrigin(), 'X')
  TEST_STRING_EQUAL((*mod_it)->getId(), "EDC")
  TEST_EQUAL((*mod_it)->getTermSpecificity(), ResidueModification::N_TERM)
  ++mod_it;

  TEST_EQUAL((*mod_it)->getOrigin(), 'Y')
  TEST_STRING_EQUAL((*mod_it)->getId(), "EDC")
  TEST_EQUAL((*mod_it)->getTermSpecificity(), ResidueModification::ANYWHERE)

  ptr->searchModifications(mods, "EDC", "", ResidueModification::C_TERM);

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
  TEST_STRING_EQUAL((*mod_it)->getId(), "EDC")
  TEST_EQUAL((*mod_it)->getTermSpecificity(), ResidueModification::C_TERM)

  // no match, thus mods should be empty
  ptr->searchModifications(mods, "EDC", "R", ResidueModification::ANYWHERE);

  TEST_EQUAL(mods.size(), 0)
}
END_SECTION


START_SECTION((void searchModificationsByDiffMonoMass(std::vector<String>& mods, double mass, double max_error, const String& residue, ResidueModification::TermSpecificity term_spec)))
{
  vector<String> mods;
  // these two cross-linkers have exactly the same mass / structure after the cross-linking reaction
  ptr->searchModificationsByDiffMonoMass(mods, 138.06807961, 0.00001, "K");
  TEST_EQUAL(find(mods.begin(), mods.end(), "DSS (K)") != mods.end(), true);
  TEST_EQUAL(find(mods.begin(), mods.end(), "BS3 (K)") != mods.end(), true);

  // something exotic.. mods should return empty (without clearing it before)
  ptr->searchModificationsByDiffMonoMass(mods, 800000000.0, 0.1, "S");
  TEST_EQUAL(mods.size(), 0);

  // terminal mod:
  ptr->searchModificationsByDiffMonoMass(mods, 138.068, 0.01, "", ResidueModification::N_TERM);
  set<String> uniq_mods;
  for (vector<String>::const_iterator it = mods.begin(); it != mods.end(); ++it)
  {
    uniq_mods.insert(*it);
  }
  TEST_EQUAL(mods.size(), 2);
  TEST_EQUAL(uniq_mods.size(), 2);
  TEST_EQUAL(uniq_mods.find("BS3 (N-term)") != uniq_mods.end(), true);

  // something exotic.. mods should return empty (without clearing it before)
  ptr->searchModificationsByDiffMonoMass(mods, 4200000, 0.1, "", ResidueModification::N_TERM);
  TEST_EQUAL(mods.size(), 0);

  ptr->searchModificationsByDiffMonoMass(mods, 138.068, 0.01);
  uniq_mods.clear();
  for (vector<String>::const_iterator it = mods.begin(); it != mods.end(); ++it)
  {
    uniq_mods.insert(*it);
  }

  TEST_EQUAL(uniq_mods.find("DSS (K)") != uniq_mods.end(), true);
  TEST_EQUAL(uniq_mods.find("BS3 (K)") != uniq_mods.end(), true);
  TEST_EQUAL(uniq_mods.find("BS3 (N-term)") != uniq_mods.end(), true);
  TEST_EQUAL(uniq_mods.find("DSS (S)") != uniq_mods.end(), true);

  // something exotic.. mods should return empty (without clearing it before)
  ptr->searchModificationsByDiffMonoMass(mods, 800000000.0, 0.1);
  TEST_EQUAL(mods.size(), 0);
}
END_SECTION

START_SECTION((const ResidueModification& getModification(const String& mod_name, const String& residue, ResidueModification::TermSpecificity term_spec) const))
{
  TEST_EQUAL(ptr->getModification("EDC (E)")->getFullId(), "EDC (E)");
  TEST_EQUAL(ptr->getModification("EDC (E)")->getId(), "EDC");

  TEST_EQUAL(ptr->getModification("DSS", "S", ResidueModification::ANYWHERE)->getId(), "DSS");
  TEST_EQUAL(ptr->getModification("DSS", "S", ResidueModification::ANYWHERE)->getFullId(), "DSS (S)");

  // terminal mod:
  TEST_EQUAL(ptr->getModification("DSS", "", ResidueModification::N_TERM)->getId(), "DSS");
  TEST_EQUAL(ptr->getModification("BS3", "", ResidueModification::N_TERM)->getFullId(), "BS3 (N-term)");
  TEST_EQUAL(ptr->getModification("EDC", "", ResidueModification::N_TERM)->getFullId(), "EDC (N-term)");
}
END_SECTION

START_SECTION((Size findModificationIndex(const String& mod_name) const))
{
  Size index = numeric_limits<Size>::max();
  index = ptr->findModificationIndex("EDC (T)");
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
  TEST_EQUAL(find(mods.begin(), mods.end(), "EDC (S)") != mods.end(), true);
  TEST_EQUAL(find(mods.begin(), mods.end(), "DSS (K)") != mods.end(), true);
  TEST_EQUAL(find(mods.begin(), mods.end(), "BS3 (N-term)") != mods.end(), true);
  TEST_EQUAL(find(mods.begin(), mods.end(), "DSS") != mods.end(), false);
  TEST_EQUAL(find(mods.begin(), mods.end(), "EDC (E)") != mods.end(), true);

  // repeat search .. return size should be the same
  Size old_size=mods.size();
  ptr->getAllSearchModifications(mods);
  TEST_EQUAL(mods.size(), old_size);
}
END_SECTION

START_SECTION((bool addModification(ResidueModification* modification)))
{
  TEST_EQUAL(ptr->has("DSS (C-term)"), false);
  std::unique_ptr<ResidueModification> modification(new ResidueModification());
  modification->setFullId("DSS (C-term)");
  ptr->addModification(std::move(modification));
  TEST_EQUAL(ptr->has("DSS (C-term)"), true);
}
END_SECTION
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST




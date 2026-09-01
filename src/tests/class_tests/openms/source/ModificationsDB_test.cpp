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
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <limits>
#include <sstream>
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

const ModificationsDB* ptr = nullptr;
const ModificationsDB* nullPointer = nullptr;

START_SECTION(const ModificationsDB* getInstance())
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

  START_SECTION((void searchModifications(std::set<const ResidueModification*>& mods, const std::string& mod_name, const std::string& residue, ResidueModification::TermSpecificity term_spec) const))
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


START_SECTION((void searchModificationsByDiffMonoMass(std::vector<std::string>& mods, double mass, double max_error, const std::string& residue, ResidueModification::TermSpecificity term_spec)))
{
  vector<std::string> mods;
  ptr->searchModificationsByDiffMonoMass(mods, 80.0, 0.1, "S");
  TEST_EQUAL(find(mods.begin(), mods.end(), "Phospho (S)") != mods.end(), true);
  TEST_EQUAL(find(mods.begin(), mods.end(), "Sulfo (S)") != mods.end(), true);

  // something exotic.. mods should return empty (without clearing it before)
  ptr->searchModificationsByDiffMonoMass(mods, 800000000.0, 0.1, "S");
  TEST_EQUAL(mods.size(), 0);

  // terminal mod:
  ptr->searchModificationsByDiffMonoMass(mods, 42, 0.1, "", ResidueModification::N_TERM);
  set<std::string> uniq_mods;
  for (vector<std::string>::const_iterator it = mods.begin(); it != mods.end(); ++it)
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
  for (vector<std::string>::const_iterator it = mods.begin(); it != mods.end(); ++it)
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

  // Modifications without a mass difference (PSI-MOD "residue" terms such as
  // 'MOD:00026 L-threonine residue') cannot explain a mass shift and must not be reported
  // for one, no matter how large the tolerance is (issue #10029).
  std::vector<const ResidueModification*> mod_ptrs;
  ptr->searchModificationsByDiffMonoMass(mod_ptrs, -0.5, 1.0, "T");
  bool found_zero_diff = false;
  for (const ResidueModification* m : mod_ptrs)
  {
    if (m->getDiffMonoMass() == 0.0) found_zero_diff = true;
  }
  TEST_EQUAL(found_zero_diff, false)
  // ... they are still found when a zero mass difference is what was asked for
  ptr->searchModificationsByDiffMonoMass(mods, 0.0, 0.001, "T");
  TEST_EQUAL(find(mods.begin(), mods.end(), "L-threonine residue (T)") != mods.end(), true)
}
END_SECTION

START_SECTION((const ResidueModification* getBestModificationByDiffMonoMass(double mass, double max_error, const std::string& residue, ResidueModification::TermSpecificity term_spec) const))
{
  TEST_EQUAL(ptr->getBestModificationByDiffMonoMass(15.994915, 0.001, "M", ResidueModification::ANYWHERE)->getFullId(), "Oxidation (M)")

  // A modification with a mass difference of zero minimizes |diff mass - mass| for every small
  // mass difference and used to be returned as the best match for it (issue #10029).
  const ResidueModification* best = ptr->getBestModificationByDiffMonoMass(-0.5, 1.0, "T", ResidueModification::ANYWHERE);
  TEST_EQUAL(best != nullptr && best->getDiffMonoMass() == 0.0, false)
  TEST_EQUAL(ptr->getBestModificationByDiffMonoMass(0.001, 0.005, "T", ResidueModification::ANYWHERE) == nullptr, true)

  // ... but a zero mass difference is still matched when explicitly searched for
  best = ptr->getBestModificationByDiffMonoMass(0.0, 0.001, "T", ResidueModification::ANYWHERE);
  TEST_EQUAL(best == nullptr, false)
  TEST_REAL_SIMILAR(best->getDiffMonoMass(), 0.0)
}
END_SECTION

START_SECTION((const ResidueModification& getModification(const std::string& mod_name, const std::string& residue, ResidueModification::TermSpecificity term_spec) const))
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

START_SECTION([EXTRA] diff formula of a UniMod entry matches its diff mono mass)
{
  // ICAT-D (UNIMOD:13) and ICAT-D:2H(8) (UNIMOD:12) are the only unimod.xml entries with <Ignore>
  // blocks. Their <element> children used to be folded into the modification's own formula, leaving
  // it ~2250 Da heavier than the delta mass that getDiffMonoMass() reports.
  // See https://github.com/OpenMS/OpenMS/issues/10030
  for (const std::string& acc : {std::string("UNIMOD:12"), std::string("UNIMOD:13")})
  {
    const ResidueModification* mod = ptr->getModification(acc);
    TEST_REAL_SIMILAR(mod->getDiffFormula().getMonoWeight(), mod->getDiffMonoMass());
  }

  TEST_REAL_SIMILAR(ptr->getModification("UNIMOD:12")->getDiffMonoMass(), 450.275205);
  TEST_REAL_SIMILAR(ptr->getModification("UNIMOD:13")->getDiffMonoMass(), 442.224991);
}
END_SECTION

START_SECTION((Size findModificationIndex(const std::string& mod_name) const))
{
  Size index = numeric_limits<Size>::max();
  index = ptr->findModificationIndex("Phospho (T)");
  TEST_NOT_EQUAL(index, numeric_limits<Size>::max());
}
END_SECTION


START_SECTION((void getAllSearchModifications(std::vector<std::string>& modifications)))
{
  vector<std::string> mods;
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

   static const ModificationsDB* mdb = ModificationsDB::getInstance();

   int nr_iterations (1e4), test (0);
#pragma omp parallel for reduction (+: test)
  for (int k = 1; k < nr_iterations + 1; k++)
  {
    int mod_id = k;
    std::string modname = "mod" + StringUtils::toStr(mod_id);
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
    std::string modname = "mod" + StringUtils::toStr(mod_id);
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

START_SECTION([EXTRA] registerDefinition - same name is the same modification and re-registration is idempotent)
{
  const ModificationsDB* db = ModificationsDB::getInstance();
  ResidueModification d;
  d.setId("TestDef:Idem");
  d.setOrigin('K');
  d.setTermSpecificity(ResidueModification::ANYWHERE);
  d.setFullId();
  d.setDiffFormula(EmpiricalFormula("C2H2O"));
  d.setDiffMonoMass(EmpiricalFormula("C2H2O").getMonoWeight());

  const Size before = db->getNumberOfModifications();
  const ResidueModification* first = db->registerDefinition(d);
  TEST_TRUE(first != nullptr)
  TEST_EQUAL(db->getNumberOfModifications(), before + 1)
  if (first != nullptr)
  {
    TEST_EQUAL(first->getProvenance(), ResidueModification::DEFINED)
  }

  // Re-registration must be SILENT: addModification() logs "already exists" on its duplicate path,
  // and registerDefinition's has() pre-check exists precisely so that re-reading a file does not
  // emit that once per definition.
  std::ostringstream captured;
  OPENMS_LOG_WARN.insert(captured);
  const ResidueModification* second = db->registerDefinition(d);
  OPENMS_LOG_WARN.remove(captured);
  TEST_TRUE(first == second)
  TEST_EQUAL(db->getNumberOfModifications(), before + 1)
  TEST_TRUE(captured.str().find("already exists") == std::string::npos)
  TEST_TRUE(captured.str().find("disagrees") == std::string::npos)

  TEST_TRUE(db->hasDefinedModification("TestDef:Idem"))
  TEST_TRUE(db->hasDefinedModification("TestDef:Idem (K)"))
  TEST_FALSE(db->hasDefinedModification("Oxidation"))
  TEST_FALSE(db->hasDefinedModification("Oxidation (M)"))
  TEST_FALSE(db->hasDefinedModification("TestDef:NeverRegistered"))

  // the same FullId with different chemistry keeps the first definition
  ResidueModification conflict(d);
  conflict.setDiffFormula(EmpiricalFormula("C3H2O"));
  conflict.setDiffMonoMass(EmpiricalFormula("C3H2O").getMonoWeight());
  std::ostringstream conflict_log;
  OPENMS_LOG_WARN.insert(conflict_log);
  const ResidueModification* kept = db->registerDefinition(conflict);
  OPENMS_LOG_WARN.remove(conflict_log);
  TEST_TRUE(kept == first)
  TEST_TRUE(conflict_log.str().find("disagrees") != std::string::npos) // ... but a real conflict is loud
  if (kept != nullptr)
  {
    TEST_EQUAL(kept->getDiffFormula().toString(), EmpiricalFormula("C2H2O").toString())
  }
  TEST_EQUAL(db->getNumberOfModifications(), before + 1)

  ResidueModification nameless;
  TEST_EXCEPTION(Exception::MissingInformation, db->registerDefinition(nameless))
}
END_SECTION

START_SECTION([EXTRA] registerDefinition - one name on two residues is two entries that resolve at their own site)
{
  const ModificationsDB* db = ModificationsDB::getInstance();
  const double delta = EmpiricalFormula("C2H2O").getMonoWeight();
  const Size before = db->getNumberOfModifications();
  for (const char origin : {'K', 'Y'})
  {
    ResidueModification d;
    d.setId("TestDef:TwoSites");
    d.setOrigin(origin);
    d.setTermSpecificity(ResidueModification::ANYWHERE);
    d.setFullId();
    d.setDiffFormula(EmpiricalFormula("C2H2O"));
    d.setDiffMonoMass(delta);
    TEST_TRUE(db->registerDefinition(d) != nullptr)
  }
  TEST_EQUAL(db->getNumberOfModifications(), before + 2)

  const ResidueModification* on_k = db->getModification("TestDef:TwoSites", "K", ResidueModification::ANYWHERE);
  const ResidueModification* on_y = db->getModification("TestDef:TwoSites", "Y", ResidueModification::ANYWHERE);
  TEST_TRUE(on_k != nullptr)
  TEST_TRUE(on_y != nullptr)
  TEST_TRUE(on_k != on_y)
  if (on_k != nullptr && on_y != nullptr)
  {
    TEST_EQUAL(on_k->getOrigin(), 'K')
    TEST_EQUAL(on_y->getOrigin(), 'Y')
    TEST_EQUAL(on_k->getFullId(), "TestDef:TwoSites (K)")
    TEST_EQUAL(on_y->getFullId(), "TestDef:TwoSites (Y)")
  }

  // The sequence round trip that catches the origin-'X' trap: the residue letter must survive.
  AASequence seq = AASequence::fromString("PEPTK(TestDef:TwoSites)IDY(TestDef:TwoSites)E");
  TEST_EQUAL(seq.toString(), "PEPTK(TestDef:TwoSites)IDY(TestDef:TwoSites)E")
  TEST_TRUE(AASequence::fromString(seq.toString()) == seq)
  TEST_REAL_SIMILAR(seq.getMonoWeight(), AASequence::fromString("PEPTKIDYE").getMonoWeight() + 2.0 * delta)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
